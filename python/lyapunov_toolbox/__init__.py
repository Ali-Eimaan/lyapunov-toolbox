"""lyapunov-toolbox — Python front end for SOS Lyapunov synthesis.

The certified computation happens in Julia (``julia/src``). This package is a
thin bridge plus the plotting layer, so that the rest of the portfolio —
``mpc-cbf``, ``cbf-safety-filter``, ``transition-viable-swarm`` — can obtain a
region-of-attraction certificate from Python without reimplementing any of the
semidefinite programming.

Public API
----------
Bridge
    :class:`JuliaSession`, :func:`julia_available`, :func:`certify_roa`,
    :func:`synthesize`, :class:`RoaCertificate`
Plotting
    :func:`use_style`, :func:`plot_roa_2d`, :func:`plot_phase_portrait`,
    :func:`plot_vdot_field`, :func:`plot_roa_comparison`,
    :func:`plot_vs_convergence`, :func:`save_figure`

Design note
-----------
Nothing here re-implements a Lyapunov computation, and nothing here silently
degrades to a Python approximation when Julia is missing. If the backend is
unavailable the import still succeeds — so that ``plotting`` remains usable on a
saved certificate — but any call that would need Julia raises
:class:`JuliaBackendError` with the actual cause attached. A bridge that
quietly returns a plausible number when the certified backend is absent is the
single worst thing this package could do.
"""

from __future__ import annotations

__version__ = "0.1.0"

__all__ = [
    "__version__",
    # julia_bridge
    "JuliaSession",
    "JuliaBackendError",
    "RoaCertificate",
    "julia_available",
    "certify_roa",
    "synthesize",
    "load_certificate",
    # plotting
    "use_style",
    "plot_roa_2d",
    "plot_phase_portrait",
    "plot_vdot_field",
    "plot_roa_comparison",
    "plot_vs_convergence",
    "save_figure",
]

from lyapunov_toolbox.julia_bridge import (
    JuliaBackendError,
    JuliaSession,
    RoaCertificate,
    certify_roa,
    julia_available,
    load_certificate,
    synthesize,
)
from lyapunov_toolbox.plotting import (
    plot_phase_portrait,
    plot_roa_2d,
    plot_roa_comparison,
    plot_vdot_field,
    plot_vs_convergence,
    save_figure,
    use_style,
)
