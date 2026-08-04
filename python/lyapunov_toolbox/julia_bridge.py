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
import re
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
    def evaluate(self, x: Array) -> Array | float:
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
        level: float | None = None,
    ) -> Array:
        """Points on ``{V = level}`` in a 2-D slice, shape ``(n, 2)``, angularly ordered.

        ``level`` defaults to ``self.rho``; pass a smaller value to draw an
        interior contour. This is the only correct way to draw ``{V = c}`` for a
        non-quadratic ``V`` — radial rescaling of the outer boundary is wrong
        whenever ``V`` is not a homogeneous quadratic.

        Matches Julia's ``roa_boundary_2d``: marches along rays from the origin,
        bisecting on the radius for ``V(r·d) = level``. Ray-marching is valid for a
        star-shaped sublevel set; where it is not (a non-star-shaped set at high
        degree) a ray crosses ``{V = level}`` more than once, and the innermost
        crossing would silently draw the wrong boundary. So, like the Julia side,
        this first coarse-scans each ray for the number of sign changes of
        ``V - level`` and raises ``ValueError`` if any ray crosses more than once.
        """
        if level is None:
            level = self.rho
        if fixed is None:
            fixed = {}
        angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
        boundary = np.zeros((n, 2))
        for i, theta in enumerate(angles):
            d = np.array([np.cos(theta), np.sin(theta)])

            def ray_point(r: float) -> NDArray[np.float64]:
                pt = np.zeros(self.n_vars)
                pt[list(axes)] = r * d
                for k, v in fixed.items():
                    pt[k] = v
                return pt

            # Coarse scan for multiple crossings along this ray, as in Julia's
            # roa_boundary_2d. Sign changes of (V - level) on [0, 100].
            ns = 100
            crossings = 0
            prev: int | None = None
            for k in range(ns + 1):
                r = 0.0 + 100.0 * k / ns
                val = float(self.evaluate(ray_point(r))) - level
                s = 1 if val > 0 else (-1 if val < 0 else 0)
                if prev is not None and s != prev and s != 0:
                    crossings += 1
                if s != 0:
                    prev = s
            if crossings > 1:
                raise ValueError(
                    f"boundary_2d: ray at angle {theta} crosses the level set "
                    f"{crossings} times; the sublevel set may not be star-shaped. "
                    f"Ray-marching is not valid in this case."
                )

            lo, hi = 0.0, 100.0
            # Bisect along ray.
            for _ in range(60):
                mid = (lo + hi) / 2
                if float(self.evaluate(ray_point(mid))) <= level:
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

        self._backend: Backend = backend or os.environ.get("LYAPUNOV_TOOLBOX_BRIDGE", "juliacall")  # type: ignore[assignment]
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
        try:
            import juliacall
        except ImportError as exc:
            raise JuliaBackendError(
                "juliacall is not installed. Install it with:\n"
                "  pip install juliacall"
            ) from exc

        # Activate the toolbox project inside the running Julia process rather
        # than via PYTHON_JULIACALL_PROJECT: newer juliacall versions refuse to
        # start when that variable is set without PYTHON_JULIACALL_EXE, and
        # activating at runtime avoids leaking state into the process environ.
        from juliacall import Main as jl_main

        jl_main.seval(f'using Pkg; Pkg.activate("{self._project}")')
        return juliacall

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
    """Return the process-wide :class:`JuliaSession`, starting it on first use.

    If a session is already active, its configuration must match the requested
    ``backend``/``project``; otherwise this raises ``RuntimeError`` naming both
    values. Returning the live session while silently ignoring a conflicting
    request makes the caller believe it is using ``backend="pyjulia"`` with
    project A when the process is actually running ``"juliacall"`` against
    project B — the hardest kind of bug to find, since everything except the
    provenance fields still looks right.
    """
    import os

    if JuliaSession._instance is not None:
        live = JuliaSession._instance
        want_backend = backend or os.environ.get("LYAPUNOV_TOOLBOX_BRIDGE", "juliacall")
        want_project = Path(
            project if project is not None
            else os.environ.get("LYAPUNOV_TOOLBOX_JULIA_PROJECT",
                                str(DEFAULT_JULIA_PROJECT))
        ).resolve()
        live_backend = live._backend
        live_project = live._project.resolve()
        if want_backend != live_backend or want_project != live_project:
            raise RuntimeError(
                f"An active JuliaSession uses backend='{live_backend}' with "
                f"project='{live_project}', but backend='{want_backend}' with "
                f"project='{want_project}' was requested. Call close() on the "
                f"existing session first, or use get_session() with no arguments "
                f"to keep the running configuration."
            )
        return live
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
    # Validate variable names are valid Julia identifiers before
    # interpolating into a Julia f-string (R1-17).
    _validate_variable_names(variables)
    if callable(dynamics):
        # A Python callable mapping the state variables to the RHS components.
        dyn_str = _vector_literal(dynamics(variables), variables)
    else:
        dyn_str = _vector_literal(dynamics, variables)
    constr_str = _vector_literal(constraints, variables) if constraints else "[]"
    inv_str = _vector_literal(invariants, variables) if invariants else "[]"

    session = get_session()
    var_str = " ".join(variables)

    # `solver` is a JuMP optimizer factory in Julia, not a name string, and the
    # default (`nothing` → default_solver() → Clarabel) is the only one the
    # bridge can construct without importing additional packages into the
    # generated Julia source. Accept only that; reject anything else loudly
    # rather than silently ignoring the argument.
    if solver is not None:
        raise ValueError(
            f"certify_roa(): solver={solver!r} is not supported by the bridge. "
            f"Omit the argument to use the default (Clarabel)."
        )

    code = f"""
    using DynamicPolynomials
    @polyvar {var_str}
    x = [{", ".join(variables)}]
    f = {dyn_str}
    g = {constr_str}
    h = {inv_str}
    sys = LyapunovToolbox.PolynomialSystem(x, f, "custom", zeros(length(x)),
                                           g, h, String[], "")
    opts = LyapunovToolbox.SOSOptions(ε = {epsilon},
                                       V_degree = {V_degree},
                                       multiplier_degree = {multiplier_degree},
                                       verbose = {'true' if verbose else 'false'})

    try
        cert = LyapunovToolbox.synthesize_lyapunov(sys; opts = opts)
        vc = LyapunovToolbox.verify_certificate(cert, sys)
        (cert, vc.algebraic.ok && vc.sampling.ok && vc.dynamical.ok)
    catch e
        # v_s_iteration can fail on numerically fragile instances (γ-step /
        # V-step infeasibility due to solver tolerance).  Structural failures
        # (e.g. V_degree=1, which can never be SOS positive-definite) are NOT
        # caught here — they must propagate so the caller sees a genuine error.
        msg = sprint(showerror, e)
        if !occursin("v_s_iteration", msg) && !occursin("γ-step", msg) && !occursin("V-step", msg)
            rethrow()
        end
        # Fall back to the cheaper linearisation-based estimate.
        V0 = LyapunovToolbox.lyapunov_from_linearization(sys)
        cert = LyapunovToolbox.estimate_roa(sys, V0; opts = opts)
        vc = LyapunovToolbox.verify_certificate(cert, sys)
        (cert, vc.algebraic.ok && vc.sampling.ok && vc.dynamical.ok)
    end
    """
    try:
        jl_cert, is_verified = session.eval(code)
    except Exception as exc:
        raise JuliaBackendError(
            f"ROA certification failed. Dynamics: {dynamics}, "
            f"variables: {variables}. Error: {exc}"
        ) from exc

    # Marshal to Python.
    exponents, coefficients = _marshal_polynomial(jl_cert.V, variables)
    return RoaCertificate(
        exponents=exponents,
        coefficients=coefficients,
        rho=float(jl_cert.ρ),
        variables=tuple(variables),
        epsilon=float(epsilon),
        system_name=str(jl_cert.system_name),
        solver_status=str(jl_cert.solver_status),
        solve_time=float(jl_cert.solve_time),
        sdp_size=int(jl_cert.sdp_size),
        bounded=bool(jl_cert.bounded),
        boundedness_evidence=str(jl_cert.boundedness_evidence),
        verified=bool(is_verified),
    )


def synthesize(
    system: str,
    *,
    V_degree: int = 4,
    multiplier_degree: int = 2,
    vs_iterations: int = 10,
    epsilon: float = 1e-6,
    shape: str | None = None,
    verbose: bool = False,
    solver: str | None = None,
) -> tuple[RoaCertificate, dict[str, Any]]:
    """Run the full V–s iteration on one of the built-in benchmark systems.

    Keyword-only, with an explicit list so a typo like ``epslion=1e-5`` fails
    loudly instead of being silently swallowed by ``**kwargs``. ``epsilon`` is
    threaded into the Julia ``SOSOptions(ε = ...)`` — before, it was recorded on
    the certificate but never passed to the solver, so the certificate claimed a
    margin the SDP did not use.
    """
    if shape is not None:
        raise NotImplementedError(
            "synthesize(): the 'shape' argument is not implemented yet. "
            "v_s_iteration defaults to the quadratic shape Σ(x_i - x_eq_i)², "
            "which is the only shape currently supported. Remove 'shape=' from "
            "the call to use the default."
        )
    if solver is not None:
        raise NotImplementedError(
            "synthesize(): the 'solver' argument is not implemented yet. "
            "Use the default solver (Clarabel)."
        )
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
                                       vs_iterations = {vs_iterations},
                                       ε = {epsilon},
                                       verbose = {'true' if verbose else 'false'})
    result = LyapunovToolbox.v_s_iteration(sys; opts = opts)
    vc = LyapunovToolbox.verify_certificate(result.certificate, sys)
    (result.certificate, result,
     vc.algebraic.ok && vc.sampling.ok && vc.dynamical.ok)
    """
    try:
        jl_cert, jl_result, is_verified = session.eval(code)
    except Exception as exc:
        raise JuliaBackendError(f"Synthesis for '{system}' failed: {exc}") from exc

    variables = [str(v) for v in session.eval("sys.x")]
    exponents, coefficients = _marshal_polynomial(jl_cert.V, variables)
    cert = RoaCertificate(
        exponents=exponents,
        coefficients=coefficients,
        rho=float(jl_cert.ρ),
        variables=tuple(variables),
        epsilon=float(epsilon),
        system_name=str(jl_cert.system_name),
        solver_status=str(jl_cert.solver_status),
        solve_time=float(jl_cert.solve_time),
        sdp_size=int(jl_cert.sdp_size),
        bounded=bool(jl_cert.bounded),
        boundedness_evidence=str(jl_cert.boundedness_evidence),
        verified=bool(is_verified),
    )
    history = {
        "rho_history": [float(v) for v in jl_result.ρ_history],
        "beta_history": [float(v) for v in jl_result.β_history],
        "volume_history": [float(v) for v in jl_result.volume_history],
        "iterations": int(jl_result.iterations),
        "converged": bool(jl_result.converged),
        "n_sdps": int(jl_result.n_sdps),
        "total_solve_time": float(jl_result.total_solve_time),
    }
    return cert, history


def load_certificate(path: str | Path) -> RoaCertificate:
    """Load a certificate saved by :meth:`RoaCertificate.save`. Needs no Julia."""
    import json

    path = Path(path)
    data = json.loads(path.read_text())
    return RoaCertificate.from_dict(data)


def _marshal_polynomial(
    julia_poly: Any,
    variables: Sequence[str],
    session: JuliaSession | None = None,
) -> tuple[NDArray[np.int64], Array]:
    """Convert a Julia polynomial into ``(exponents, coefficients)``.

    The polynomial arrives as a backend-wrapped Julia object (a juliacall
    ``JlValue`` or a pyjulia ``julia.Any``), which must be passed to Julia
    functions as an argument rather than string-interpolated into Julia
    source (a wrapped object is a Python-side reference, not a Julia
    binding). All extraction runs in one Julia helper that returns plain
    arrays. Works with whichever backend the active session uses.
    """
    if session is None:
        session = get_session()

    helper_source = """
    function _py_marshal_poly(p, varnames)
        vs = DynamicPolynomials.variables(p)
        idx = [findfirst(v -> string(v) == vn, vs) for vn in varnames]
        ms = DynamicPolynomials.monomials(p)
        cs = DynamicPolynomials.coefficients(p)
        # A variable absent from the polynomial is a constant-direction column:
        # exponent 0 for every term. (A polynomial in fewer variables than the
        # state vector is legal — e.g. V depending on x1 only — and erroring on
        # it made marshalling fragile for no benefit.)
        E = [i === nothing ? 0 : degree(m, vs[i]) for m in ms, i in idx]
        return (E, cs)
    end
    """

    n_vars = len(variables)

    if session._backend == "juliacall":
        import juliacall

        jl = juliacall.Main
        jl.seval("using DynamicPolynomials")
        if not hasattr(jl, "_py_marshal_poly"):
            jl.seval(helper_source)
        exponents, coefficients = jl.seval("_py_marshal_poly")(
            julia_poly, list(variables)
        )
    else:
        from julia import Main as jl_main

        jl_main.eval("using DynamicPolynomials")
        if not hasattr(jl_main, "_py_marshal_poly"):
            jl_main.eval(helper_source)
        exponents, coefficients = jl_main._py_marshal_poly(
            julia_poly, list(variables)
        )

    exponents = np.asarray(exponents, dtype=np.int64)
    coefficients = np.asarray(coefficients, dtype=np.float64)
    if exponents.size == 0:
        exponents = np.zeros((1, n_vars), dtype=np.int64)
        coefficients = np.zeros(1, dtype=np.float64)
    return exponents, coefficients


_VALID_VARIABLE_RE = re.compile(r'^[A-Za-z][A-Za-z0-9_]*$')

def _validate_variable_names(variables: Sequence[str]) -> None:
    """Reject variable names that are not valid Julia identifiers (R1-17).

    Variable names are interpolated directly into Julia f-strings; an
    unvalidated name could inject Julia code (e.g. a name containing
    ``; using evil_mod; #``).
    """
    for v in variables:
        if not _VALID_VARIABLE_RE.match(v):
            raise ValueError(
                f"Variable name {v!r} is not a valid identifier. "
                f"Must match pattern /^[A-Za-z][A-Za-z0-9_]*$/."
            )

def _check_expression_safety(expr: str, variables: Sequence[str]) -> None:
    """Reject a dynamics/constraint string containing anything but arithmetic.

    Uses a **token** allowlist (R1-18), not a character-class filter.
    The character-class approach allows through strings like ``)(`` or
    ``;`` which are not valid in a Julia polynomial expression. This
    tokeniser validates that every token is either:
      - a variable name (from ``variables``)
      - a numeric literal (integer or float)
      - an operator: ``+``, ``-``, ``*``, ``/``, ``^``, ``(``, ``)``
    """
    tokens = re.findall(
        r'[A-Za-z_][A-Za-z0-9_]*|\d+\.?\d*(?:[eE][+-]?\d+)?|[+\-*/^()]',
        expr,
    )
    reconstructed = "".join(tokens)
    if reconstructed != re.sub(r'\s+', '', expr):
        raise ValueError(
            f"Expression {expr!r} contains characters outside the allowed set. "
            f"Only arithmetic (+, -, *, /, ^) on {variables} is allowed."
        )
    var_set = set(variables)
    for token in tokens:
        # tokens: identifiers, numbers, operators
        if token in var_set:
            continue
        if token in {"+", "-", "*", "/", "^", "(", ")"}:
            continue
        if re.match(r'^\d+\.?\d*(?:[eE][+-]?\d+)?$', token):
            continue
        raise ValueError(
            f"Disallowed token {token!r} in expression {expr!r}. "
            f"Only arithmetic (+, -, *, /, ^) on {variables} is allowed."
        )


def _vector_literal(
    exprs: Sequence[str] | str, variables: Sequence[str]
) -> str:
    """Validate and normalise a vector of expressions into a Julia literal.

    Accepts either a sequence of component strings (``["-x2", "x1 - x2"]``)
    or a single string already containing the whole Julia vector literal
    (``"[-x2, x1 - x2]"``).
    """
    if isinstance(exprs, str):
        inner = exprs.strip()
        if inner.startswith("[") and inner.endswith("]"):
            inner = inner[1:-1]
        for component in inner.split(","):
            _check_expression_safety(component.strip(), variables)
        return exprs
    components = list(exprs)
    for component in components:
        _check_expression_safety(component, variables)
    return "[" + ", ".join(components) + "]"
