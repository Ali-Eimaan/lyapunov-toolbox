"""Bridge from Python to the Julia SOS routines in ``julia/src``.

Why a bridge at all
-------------------
SumOfSquares.jl is the best-maintained SOS modelling layer in any language, and
reimplementing it in Python would mean either depending on a commercial solver
interface or writing a worse version of it. The rest of the portfolio is Python,
so the boundary lives here: this module is the *only* place in the repository
that knows both languages.

Backends
--------
Two are supported, selected by :func:`get_session` or the
``LYAPUNOV_TOOLBOX_BRIDGE`` environment variable:

``"juliacall"`` (default)
    ``juliacall`` / ``PythonCall.jl``. Manages its own ``libjulia``, works with a
    statically linked CPython, and does not need a special launcher.

``"pyjulia"``
    The ``julia`` package named in the repository specification. Kept because it
    is what the original design called for and because some environments already
    have it configured. On a statically linked CPython it requires either
    ``python-jl`` or ``julia.install()``; the error message for that case must say
    so explicitly rather than surfacing the raw segfault-adjacent failure.

Both are exercised in CI. If they ever disagree on a certificate, that is a bug
in this file, not a solver artefact.

Environment
-----------
``LYAPUNOV_TOOLBOX_JULIA_PROJECT``
    Path to the ``julia/`` project directory. Defaults to the sibling ``julia``
    directory resolved from this file, not from the working directory, so that
    the bridge works when imported from another repository.
``LYAPUNOV_TOOLBOX_BRIDGE``
    ``"juliacall"`` or ``"pyjulia"``.

Data crossing the boundary
--------------------------
Polynomials do **not** cross as objects. They cross as
``(monomial_exponents, coefficients)`` pairs — an integer array of shape
``(n_terms, n_vars)`` and a float array of shape ``(n_terms,)``. That
representation is explicit, cheap to serialise, trivially testable, and survives
a version bump of either polynomial library. Passing opaque wrapped Julia
objects into Python code that then evaluates them through the bridge is much
slower and produces stack traces nobody can read.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]
Backend = Literal["juliacall", "pyjulia"]

#: Default location of the Julia project, resolved relative to this file.
DEFAULT_JULIA_PROJECT = Path(__file__).resolve().parents[2] / "julia"


class JuliaBackendError(RuntimeError):
    """Raised when the Julia backend is unavailable or a Julia call fails.

    Always chains the original exception (``raise ... from exc``) and includes,
    in the message: which backend was attempted, the resolved project path, and
    the concrete next command the user should run. "Julia not found" with no
    further detail costs an hour of somebody's afternoon.
    """


@dataclass(frozen=True)
class RoaCertificate:
    """A region-of-attraction certificate, marshalled from Julia.

    This is a *value object*: it holds no live Julia handle, is picklable, and
    can be round-tripped through JSON. Downstream repositories store these
    alongside their results, which is what makes the "computed with
    lyapunov-toolbox at commit abc123" citation checkable rather than decorative.

    Attributes
    ----------
    exponents:
        ``(n_terms, n_vars)`` non-negative integers; row ``i`` is the exponent
        vector of monomial ``i`` in ``V``.
    coefficients:
        ``(n_terms,)`` coefficients aligned with ``exponents``.
    rho:
        Certified level: the set is ``{x : V(x) <= rho}``.
    variables:
        Names of the state variables, in the same column order as ``exponents``.
    epsilon:
        The positive-definiteness margin used in the SOS program.
    system_name, solver_status, solve_time, sdp_size:
        Provenance, copied verbatim from the Julia ``LyapunovCertificate``.
    bounded, boundedness_evidence:
        Whether the sublevel set was shown to be bounded, and how. A certificate
        with ``bounded=False`` is not a certificate; construction rejects it.
    verified:
        Whether Julia's ``verify_certificate`` passed all three of its checks at
        marshalling time. Never set this by hand.
    """

    exponents: NDArray[np.int64]
    coefficients: Array
    rho: float
    variables: tuple[str, ...]
    epsilon: float
    system_name: str
    solver_status: str
    solve_time: float
    sdp_size: int
    bounded: bool
    boundedness_evidence: str
    verified: bool

    def __post_init__(self) -> None:
        """Validate shapes, non-negative exponents, ``rho > 0``, and ``bounded``."""
        if self.exponents.ndim != 2:
            raise ValueError(f"exponents must be 2-D, got shape {self.exponents.shape}")
        if self.coefficients.ndim != 1:
            raise ValueError(f"coefficients must be 1-D, got shape {self.coefficients.shape}")
        if self.exponents.shape[0] != self.coefficients.shape[0]:
            raise ValueError(
                f"exponents ({self.exponents.shape[0]} terms) and coefficients "
                f"({self.coefficients.shape[0]} terms) must have the same length"
            )
        if np.any(self.exponents < 0):
            raise ValueError("exponents must be non-negative")
        if self.rho <= 0:
            raise ValueError(f"rho must be positive, got {self.rho}")
        if not self.bounded:
            raise ValueError("certificate is not bounded — not a valid ROA certificate")

    # ------------------------------------------------------------------ shape
    @property
    def n_vars(self) -> int:
        """Number of state variables."""
        return self.exponents.shape[1]

    @property
    def degree(self) -> int:
        """Total degree of ``V``."""
        return int(np.max(np.sum(self.exponents, axis=1)))

    # ------------------------------------------------------------- evaluation
    def evaluate(self, x: Array) -> Array:
        """Evaluate ``V`` at ``x``.

        Accepts ``(n_vars,)`` or ``(k, n_vars)`` and returns a scalar or ``(k,)``.
        """
        x = np.atleast_2d(x)
        if x.shape[1] != self.n_vars:
            raise ValueError(
                f"expected x with {self.n_vars} columns, got shape {x.shape}"
            )
        # monomial_values[i, k] = ∏_j x[k, j]^exponents[i, j]
        monomial_vals = np.prod(
            x[:, None, :] ** self.exponents[None, :, :], axis=-1
        )  # shape: (k, n_terms)
        result = monomial_vals @ self.coefficients  # shape: (k,)
        if result.shape[0] == 1:
            return float(result[0])
        return result

    def contains(self, x: Array, tol: float = 0.0) -> NDArray[np.bool_] | bool:
        """Membership in ``{V <= rho + tol}``."""
        vals = self.evaluate(x)
        if np.ndim(vals) == 0:
            return bool(vals <= self.rho + tol)
        return vals <= self.rho + tol

    def boundary_2d(
        self,
        axes: tuple[int, int] = (0, 1),
        fixed: dict[int, float] | None = None,
        n: int = 400,
    ) -> Array:
        """Points on ``{V = rho}`` in a 2-D slice, shape ``(n, 2)``, angularly ordered."""
        if fixed is None:
            fixed = {}
        angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
        boundary = np.zeros((n, 2))
        for i, theta in enumerate(angles):
            d = np.array([np.cos(theta), np.sin(theta)])
            lo, hi = 0.0, 100.0
            # Bisect along ray.
            for _ in range(60):
                mid = (lo + hi) / 2
                pt = np.zeros(self.n_vars)
                pt[list(axes)] = mid * d
                for k, v in fixed.items():
                    pt[k] = v
                if self.evaluate(pt) <= self.rho:
                    lo = mid
                else:
                    hi = mid
            boundary[i] = lo * d
        return boundary

    # ---------------------------------------------------------- serialisation
    def to_dict(self) -> dict[str, Any]:
        """JSON-serialisable dict with every field, plus a ``"format_version"`` key."""
        return {
            "format_version": 1,
            "exponents": self.exponents.tolist(),
            "coefficients": self.coefficients.tolist(),
            "rho": self.rho,
            "variables": list(self.variables),
            "epsilon": self.epsilon,
            "system_name": self.system_name,
            "solver_status": self.solver_status,
            "solve_time": self.solve_time,
            "sdp_size": self.sdp_size,
            "bounded": self.bounded,
            "boundedness_evidence": self.boundedness_evidence,
            "verified": self.verified,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> RoaCertificate:
        """Inverse of :meth:`to_dict`; validates ``format_version``."""
        if data.get("format_version") != 1:
            raise ValueError(
                f"unsupported format_version: {data.get('format_version')}; expected 1"
            )
        return cls(
            exponents=np.array(data["exponents"], dtype=np.int64),
            coefficients=np.array(data["coefficients"], dtype=np.float64),
            rho=float(data["rho"]),
            variables=tuple(data["variables"]),
            epsilon=float(data["epsilon"]),
            system_name=str(data["system_name"]),
            solver_status=str(data["solver_status"]),
            solve_time=float(data["solve_time"]),
            sdp_size=int(data["sdp_size"]),
            bounded=bool(data["bounded"]),
            boundedness_evidence=str(data["boundedness_evidence"]),
            verified=bool(data["verified"]),
        )

    def save(self, path: str | Path) -> None:
        """Write :meth:`to_dict` as indented JSON."""
        import json

        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2))


class JuliaSession:
    """A started Julia runtime with ``LyapunovToolbox`` loaded."""

    _instance: JuliaSession | None = None

    def __init__(
        self,
        backend: Backend | None = None,
        project: str | Path | None = None,
        instantiate: bool = False,
    ) -> None:
        if JuliaSession._instance is not None:
            raise RuntimeError(
                "A JuliaSession is already active. Call close() on it first, "
                "or use get_session() to retrieve the existing one."
            )
        import os

        self._backend: Backend = backend or os.environ.get("LYAPUNOV_TOOLBOX_BRIDGE", "juliacall")
        if project is None:
            project = os.environ.get("LYAPUNOV_TOOLBOX_JULIA_PROJECT", str(DEFAULT_JULIA_PROJECT))
        self._project = Path(project)
        self._jl = None

        if self._backend == "juliacall":
            self._jl = self._start_juliacall()
        elif self._backend == "pyjulia":
            self._jl = self._start_pyjulia()
        else:
            raise ValueError(
                f"Unknown backend '{self._backend}'; expected 'juliacall' or 'pyjulia'"
            )

        if instantiate:
            self.eval("using Pkg; Pkg.instantiate()")
        else:
            try:
                self.eval("using LyapunovToolbox")
            except Exception as exc:
                raise JuliaBackendError(
                    f"Failed to load LyapunovToolbox. The Julia project at "
                    f"'{self._project}' may not be instantiated. Run:\n"
                    f"  julia --project={self._project} -e 'using Pkg; Pkg.instantiate()'"
                ) from exc

        JuliaSession._instance = self

    def _start_juliacall(self) -> Any:
        """Start via ``juliacall``."""
        import os
        os.environ["PYTHON_JULIACALL_PROJECT"] = str(self._project)
        try:
            import juliacall
            return juliacall
        except ImportError as exc:
            raise JuliaBackendError(
                "juliacall is not installed. Install it with:\n"
                "  pip install juliacall"
            ) from exc

    def _start_pyjulia(self) -> Any:
        """Start via ``julia`` (PyJulia)."""
        try:
            import julia
        except ImportError as exc:
            raise JuliaBackendError(
                "pyjulia (the 'julia' package) is not installed. Install it with:\n"
                "  pip install julia"
            ) from exc
        try:
            jl = julia.Julia(compiled_modules=False)
            from julia import Main as _main
            _main.eval(f'using Pkg; Pkg.activate("{self._project}")')
            return jl
        except Exception as exc:
            raise JuliaBackendError(
                f"Failed to start PyJulia. On a statically-linked CPython, "
                f"use 'python-jl' or call 'julia.install()' first. "
                f"Original error: {exc}"
            ) from exc

    def eval(self, code: str) -> Any:
        """Evaluate Julia source in the ``Main`` module."""
        if self._backend == "juliacall":
            from juliacall import Main as jl_main
            return jl_main.seval(code)
        else:
            from julia import Main as jl_main
            return jl_main.eval(code)

    def call(self, fname: str, *args: Any, **kwargs: Any) -> Any:
        """Call ``LyapunovToolbox.<fname>`` with marshalled arguments."""
        try:
            if self._backend == "juliacall":
                from juliacall import Main as jl_main
                func = jl_main.seval(f"LyapunovToolbox.{fname}")
            else:
                from julia import Main as jl_main
                func = jl_main.eval(f"LyapunovToolbox.{fname}")
            return func(*args, **kwargs)
        except Exception as exc:
            raise JuliaBackendError(
                f"Call to LyapunovToolbox.{fname} failed: {exc}"
            ) from exc

    @property
    def version(self) -> str:
        """Julia version string."""
        return str(self.eval("VERSION"))

    def close(self) -> None:
        """Release the session."""
        JuliaSession._instance = None
        self._jl = None


# ------------------------------------------------------------------- session
def get_session(
    backend: Backend | None = None,
    project: str | Path | None = None,
    instantiate: bool = False,
) -> JuliaSession:
    """Return the process-wide :class:`JuliaSession`, starting it on first use."""
    if JuliaSession._instance is not None:
        return JuliaSession._instance
    return JuliaSession(backend=backend, project=project, instantiate=instantiate)


def julia_available(backend: Backend | None = None) -> bool:
    """Can the backend be started and ``LyapunovToolbox`` loaded?"""
    try:
        session = get_session(backend=backend, instantiate=False)
        session.eval("using LyapunovToolbox")
        return True
    except Exception:
        return False


# --------------------------------------------------------------- top-level API
def certify_roa(
    dynamics: Sequence[str] | Callable[..., Sequence[Any]],
    variables: Sequence[str],
    *,
    V_degree: int = 2,
    multiplier_degree: int = 2,
    epsilon: float = 1e-6,
    constraints: Sequence[str] = (),
    invariants: Sequence[str] = (),
    solver: str | None = None,
    verbose: bool = False,
) -> RoaCertificate:
    """Certify an ROA for a polynomial vector field given as strings."""
    if isinstance(dynamics, (list, tuple)):
        for expr in dynamics:
            _check_expression_safety(expr, variables)
        for expr in constraints:
            _check_expression_safety(expr, variables)
        for expr in invariants:
            _check_expression_safety(expr, variables)

    session = get_session()
    var_str = "[" + ", ".join(variables) + "]"
    dyn_str = "[" + ", ".join(dynamics) + "]"
    constr_str = "[" + ", ".join(constraints) + "]" if constraints else "[]"
    inv_str = "[" + ", ".join(invariants) + "]" if invariants else "[]"

    code = f"""
    using DynamicPolynomials
    @polyvar {var_str}
    x = [{", ".join(variables)}]
    f = {dyn_str}
    g = {constr_str}
    h = {inv_str}
    sys = LyapunovToolbox.PolynomialSystem("custom", f, x, zeros(length(x));
                                            constraints = g, equalities = h)
    opts = LyapunovToolbox.SOSOptions(ε = {epsilon},
                                       V_degree = {V_degree},
                                       multiplier_degree = {multiplier_degree},
                                       verbose = {verbose})
    result = LyapunovToolbox.synthesize_lyapunov(sys; opts = opts)
    cert = result.certificate
    LyapunovToolbox.verify_certificate(cert, sys)
    cert
    """
    try:
        jl_cert = session.eval(code)
    except Exception as exc:
        raise JuliaBackendError(
            f"ROA certification failed. Dynamics: {dynamics}, "
            f"variables: {variables}. Error: {exc}"
        ) from exc

    # Marshal to Python.
    exponents, coefficients = _marshal_polynomial(jl_cert.V, variables)
    is_verified = bool(session.eval(
        f"LyapunovToolbox.verify_certificate({jl_cert}, sys)"
    ))
    return RoaCertificate(
        exponents=exponents,
        coefficients=coefficients,
        rho=float(jl_cert.ρ),
        variables=tuple(variables),
        epsilon=float(epsilon),
        system_name=str(jl_cert.system_name),
        solver_status=str(jl_cert.status),
        solve_time=float(jl_cert.solve_time),
        sdp_size=int(jl_cert.gram_dim),
        bounded=bool(jl_cert.bounded),
        boundedness_evidence=str(jl_cert.boundedness_evidence),
        verified=is_verified,
    )


def synthesize(
    system: str,
    *,
    V_degree: int = 4,
    multiplier_degree: int = 2,
    vs_iterations: int = 10,
    shape: str | None = None,
    **kwargs: Any,
) -> tuple[RoaCertificate, dict[str, Any]]:
    """Run the full V–s iteration on one of the built-in benchmark systems."""
    session = get_session()
    system_map = {
        "van_der_pol": "LyapunovToolbox.van_der_pol(reversed=true)",
        "inverted_pendulum": "LyapunovToolbox.inverted_pendulum()",
        "unicycle_with_obstacle": "LyapunovToolbox.unicycle_with_obstacle()",
        "consensus_network": "LyapunovToolbox.consensus_network()",
    }
    if system not in system_map:
        raise ValueError(
            f"Unknown system '{system}'. Available: {list(system_map.keys())}"
        )

    code = f"""
    sys = {system_map[system]}
    opts = LyapunovToolbox.SOSOptions(V_degree = {V_degree},
                                       multiplier_degree = {multiplier_degree},
                                       vs_iterations = {vs_iterations})
    result = LyapunovToolbox.synthesize_lyapunov(sys; opts = opts)
    (result.certificate, result)
    """
    try:
        jl_cert, jl_result = session.eval(code)
    except Exception as exc:
        raise JuliaBackendError(f"Synthesis for '{system}' failed: {exc}") from exc

    variables = [str(v) for v in session.eval("sys.x")]
    exponents, coefficients = _marshal_polynomial(jl_cert.V, variables)
    is_verified = bool(session.eval(
        f"LyapunovToolbox.verify_certificate({jl_cert}, sys)"
    ))
    cert = RoaCertificate(
        exponents=exponents,
        coefficients=coefficients,
        rho=float(jl_cert.ρ),
        variables=tuple(variables),
        epsilon=float(kwargs.get("epsilon", 1e-6)),
        system_name=str(jl_cert.system_name),
        solver_status=str(jl_cert.status),
        solve_time=float(jl_cert.solve_time),
        sdp_size=int(jl_cert.gram_dim),
        bounded=bool(jl_cert.bounded),
        boundedness_evidence=str(jl_cert.boundedness_evidence),
        verified=is_verified,
    )
    history = {
        "rho_history": [float(v) for v in jl_result.ρ_history],
        "beta_history": [float(v) for v in jl_result.β_history],
        "volume_history": [float(v) for v in jl_result.volume_history],
        "iterations": int(jl_result.iterations),
        "converged": bool(jl_result.converged),
        "n_sdps": int(jl_result.total_sdps),
        "total_solve_time": float(jl_result.total_time),
    }
    return cert, history


def load_certificate(path: str | Path) -> RoaCertificate:
    """Load a certificate saved by :meth:`RoaCertificate.save`. Needs no Julia."""
    import json

    path = Path(path)
    data = json.loads(path.read_text())
    return RoaCertificate.from_dict(data)


def _marshal_polynomial(julia_poly: Any, variables: Sequence[str]) -> tuple[NDArray[np.int64], Array]:
    """Convert a Julia polynomial into ``(exponents, coefficients)``."""
    import juliacall

    # Convert to Julia dict of monomial => coefficient.
    terms_dict = juliacall.Main.seval(
        f"DynamicPolynomials.monomials({julia_poly}) => "
        f"DynamicPolynomials.coefficients({julia_poly})"
    )
    n_vars = len(variables)
    exponents_list = []
    coeffs_list = []
    for monomial, coeff in terms_dict.items():
        exp_vec = np.zeros(n_vars, dtype=np.int64)
        for i, var_name in enumerate(variables):
            exp_vec[i] = int(juliacall.Main.seval(f"DynamicPolynomials.degree({monomial}, {var_name})"))
        exponents_list.append(exp_vec)
        coeffs_list.append(float(coeff))

    if not exponents_list:
        exponents = np.zeros((1, n_vars), dtype=np.int64)
        coefficients = np.zeros(1, dtype=np.float64)
    else:
        exponents = np.array(exponents_list, dtype=np.int64)
        coefficients = np.array(coeffs_list, dtype=np.float64)
    return exponents, coefficients


def _check_expression_safety(expr: str, variables: Sequence[str]) -> None:
    """Reject a dynamics/constraint string containing anything but arithmetic."""
    import re

    allowed_pattern = re.compile(
        r'^[\d\s\+\-\*\/\^\(\)\.' +
        ''.join(re.escape(v) for v in variables) +
        r']+$'
    )
    if not allowed_pattern.match(expr):
        # Find the offending token.
        forbidden = {"using", "import", "include", "run", "ccall"}
        for token in forbidden:
            if token in expr:
                raise ValueError(
                    f"Forbidden token '{token}' in expression '{expr}'. "
                    f"Only arithmetic on {variables} is allowed."
                )
        raise ValueError(
            f"Expression '{expr}' contains characters outside the allowed set. "
            f"Only arithmetic (+, -, *, /, ^) on {variables} is allowed."
        )
