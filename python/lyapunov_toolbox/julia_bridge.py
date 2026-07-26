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
        """Validate shapes, non-negative exponents, ``rho > 0``, and ``bounded``.

        Raise ``ValueError`` on any violation. This runs on every certificate
        that enters Python, including ones loaded from disk, so a hand-edited
        JSON file cannot smuggle in an invalid certificate.
        """
        raise NotImplementedError

    # ------------------------------------------------------------------ shape
    @property
    def n_vars(self) -> int:
        """Number of state variables."""
        raise NotImplementedError

    @property
    def degree(self) -> int:
        """Total degree of ``V``."""
        raise NotImplementedError

    # ------------------------------------------------------------- evaluation
    def evaluate(self, x: Array) -> Array:
        """Evaluate ``V`` at ``x``.

        Accepts ``(n_vars,)`` or ``(k, n_vars)`` and returns a scalar or ``(k,)``.
        Implemented in NumPy — ``np.prod(x[:, None, :] ** exponents, axis=-1) @ coefficients``
        — so evaluation costs nothing and needs no live Julia session. This is
        what makes plotting a saved certificate possible on a machine with no
        Julia at all.
        """
        raise NotImplementedError

    def contains(self, x: Array, tol: float = 0.0) -> NDArray[np.bool_] | bool:
        """Membership in ``{V <= rho + tol}``."""
        raise NotImplementedError

    def boundary_2d(
        self,
        axes: tuple[int, int] = (0, 1),
        fixed: dict[int, float] | None = None,
        n: int = 400,
    ) -> Array:
        """Points on ``{V = rho}`` in a 2-D slice, shape ``(n, 2)``, angularly ordered.

        Mirrors ``roa_boundary_2d`` in Julia but runs in NumPy by bisecting along
        rays. The two implementations must agree to ``1e-6``; the bridge test
        asserts it, because a plotting-side reimplementation that drifts from the
        certified object is how a figure ends up showing a set that was never
        certified.

        ``fixed`` pins the non-plotted variables by index; every figure produced
        from a slice must state the pinned values in its title.
        """
        raise NotImplementedError

    # ---------------------------------------------------------- serialisation
    def to_dict(self) -> dict[str, Any]:
        """JSON-serialisable dict with every field, plus a ``"format_version"`` key."""
        raise NotImplementedError

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> RoaCertificate:
        """Inverse of :meth:`to_dict`; validates ``format_version``."""
        raise NotImplementedError

    def save(self, path: str | Path) -> None:
        """Write :meth:`to_dict` as indented JSON."""
        raise NotImplementedError


class JuliaSession:
    """A started Julia runtime with ``LyapunovToolbox`` loaded.

    Starting Julia costs seconds and loading SumOfSquares costs more, so the
    session is a process-level singleton obtained through :func:`get_session`.
    Creating a second one with different settings raises rather than silently
    returning the first — two sessions that disagree about which project is
    active produce failures that make no sense at the call site.

    Parameters
    ----------
    backend:
        ``"juliacall"`` or ``"pyjulia"``. ``None`` reads
        ``LYAPUNOV_TOOLBOX_BRIDGE``, defaulting to ``"juliacall"``.
    project:
        Julia project directory. ``None`` reads ``LYAPUNOV_TOOLBOX_JULIA_PROJECT``,
        defaulting to :data:`DEFAULT_JULIA_PROJECT`.
    instantiate:
        Run ``Pkg.instantiate()`` on startup. Default ``False``: instantiating
        implicitly can trigger a multi-minute download from inside what looked
        like a plotting call. The error raised when the project is not
        instantiated must print the exact command to fix it.
    """

    def __init__(
        self,
        backend: Backend | None = None,
        project: str | Path | None = None,
        instantiate: bool = False,
    ) -> None:
        raise NotImplementedError

    # --------------------------------------------------------------- startup
    def _start_juliacall(self) -> Any:
        """Start via ``juliacall``.

        Set ``PYTHON_JULIACALL_PROJECT`` before importing ``juliacall`` — it reads
        the environment at import time and ignores it afterwards, which is the
        usual reason the wrong project ends up active.
        """
        raise NotImplementedError

    def _start_pyjulia(self) -> Any:
        """Start via ``julia`` (PyJulia).

        Use ``julia.api.Julia(runtime=..., compiled_modules=False)`` when the
        interpreter is statically linked. Detect that case and explain it in the
        exception rather than letting the process die.
        """
        raise NotImplementedError

    # ------------------------------------------------------------------ calls
    def eval(self, code: str) -> Any:
        """Evaluate Julia source in the ``Main`` module. Escape hatch, used sparingly."""
        raise NotImplementedError

    def call(self, fname: str, *args: Any, **kwargs: Any) -> Any:
        """Call ``LyapunovToolbox.<fname>`` with marshalled arguments.

        Wrap every Julia exception in :class:`JuliaBackendError`, preserving the
        Julia stack trace in the message — a bare ``JuliaError: MethodError`` with
        no trace is not actionable.
        """
        raise NotImplementedError

    @property
    def version(self) -> str:
        """Julia version string, for provenance lines in logs and benchmark tables."""
        raise NotImplementedError

    def close(self) -> None:
        """Release the session. Julia cannot truly be shut down in-process; this
        drops the reference and clears the singleton so tests can rebuild it."""
        raise NotImplementedError


def get_session(
    backend: Backend | None = None,
    project: str | Path | None = None,
    instantiate: bool = False,
) -> JuliaSession:
    """Return the process-wide :class:`JuliaSession`, starting it on first use."""
    raise NotImplementedError


def julia_available(backend: Backend | None = None) -> bool:
    """Can the backend be started and ``LyapunovToolbox`` loaded?

    Never raises. Used to skip tests and to choose an error message. It does
    start Julia if it can, so it is not free — call it once and cache.
    """
    raise NotImplementedError


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
    """Certify an ROA for a polynomial vector field given as strings.

    ``dynamics`` is a sequence of Julia-parsable polynomial expressions in
    ``variables``, e.g. ``["-x2", "x1 + (x1^2 - 1)*x2"]``. Strings rather than a
    Python AST because the alternative is writing a Python-to-Julia expression
    compiler, which is a project of its own and not this one.

    Security note: these strings are evaluated as Julia code. Restrict the parse
    to arithmetic on the declared variables and reject anything containing
    ``using``, ``import``, ``include``, ``run``, ``ccall``, or a semicolon, with a
    ``ValueError`` naming the offending token. This function is called from
    notebooks and from other repositories in the portfolio; it must not be a way
    to execute arbitrary code by passing a "dynamics" string.

    Returns a fully verified :class:`RoaCertificate`: ``verify_certificate`` is
    run on the Julia side before marshalling, and a certificate that fails it is
    never returned. Raises :class:`JuliaBackendError` if it does.
    """
    raise NotImplementedError


def synthesize(
    system: str,
    *,
    V_degree: int = 4,
    multiplier_degree: int = 2,
    vs_iterations: int = 10,
    shape: str | None = None,
    **kwargs: Any,
) -> tuple[RoaCertificate, dict[str, Any]]:
    """Run the full V–s iteration on one of the built-in benchmark systems.

    ``system`` is one of ``"van_der_pol"``, ``"inverted_pendulum"``,
    ``"unicycle_with_obstacle"``, ``"consensus_network"`` — the names in
    ``julia/src/systems.jl``. Returns the certificate and a history dict with
    ``rho_history``, ``beta_history``, ``volume_history``, ``converged``,
    ``iterations``, ``n_sdps``, ``total_solve_time``, suitable for
    :func:`~lyapunov_toolbox.plotting.plot_vs_convergence`.

    ``converged`` must be surfaced to the caller and must not be dropped: the
    iteration is a heuristic and a result that ran out of budget has to be
    reported as such.
    """
    raise NotImplementedError


def load_certificate(path: str | Path) -> RoaCertificate:
    """Load a certificate saved by :meth:`RoaCertificate.save`. Needs no Julia."""
    raise NotImplementedError


def _marshal_polynomial(julia_poly: Any, variables: Sequence[str]) -> tuple[NDArray[np.int64], Array]:
    """Convert a Julia polynomial into ``(exponents, coefficients)``.

    Column order of ``exponents`` follows ``variables`` as given, **not**
    ``MultivariatePolynomials.variables(p)``, whose order is the library's own.
    Getting this wrong yields a ``V`` that evaluates correctly in Julia and
    incorrectly in Python — a plot that disagrees with a certificate for no
    visible reason. Assert the round trip on the Julia side before returning.
    """
    raise NotImplementedError


def _check_expression_safety(expr: str, variables: Sequence[str]) -> None:
    """Reject a dynamics/constraint string containing anything but arithmetic.

    Allowed: the declared variable names, decimal literals, ``+ - * / ^``,
    parentheses, whitespace. Everything else raises ``ValueError`` naming the
    offending token. Implement as an allowlist over a tokenised scan, never as a
    denylist of dangerous words — a denylist here is not a security control.
    """
    raise NotImplementedError
