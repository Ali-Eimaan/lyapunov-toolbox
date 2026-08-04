"""Figures for ROA certificates.

Four of these plots go in the README and are the first thing a reader sees, so
they carry more weight than their line count suggests. The rules below are not
stylistic preferences — each one exists because its absence makes a figure
overstate what was proved.

Rules
-----
1. **Every function takes an optional ``ax`` and returns the ``Axes``** (or the
   ``Figure`` for multi-panel layouts). None of them call ``plt.show()`` or
   ``plt.savefig()``; the caller decides.
2. **A certified set and an estimated set are never drawn the same way.** The SOS
   sublevel set is a filled region with a solid edge; the empirically estimated
   true ROA boundary is a dashed line. A reader must be able to tell the proved
   thing from the simulated thing without reading the legend.
3. **A 2-D slice of a higher-dimensional set always states its pinned values in
   the title.** A slice presented as the whole set is the most common way an ROA
   figure misleads.
4. **Axis labels carry units.** They come from ``PolynomialSystem.state_labels``.
5. No global state beyond ``rcParams``.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]


def use_style() -> None:
    """Set the repository's Matplotlib defaults.

    ``rcParams`` only — no ``plt.style.use`` of a third-party sheet, so that a
    reader reproducing a figure gets the same result without installing
    anything. Figure size ``(7, 6)`` (ROA plots are square-ish, unlike the
    time-series figures in ``mpc-foundations``), dpi 130, ``mathtext.fontset="cm"``,
    grid on at alpha 0.25, top and right spines off, and an explicit colour cycle
    so certified/estimated/trajectory colours are stable across every figure in
    the README.
    """
    import matplotlib as mpl

    mpl.rcParams.update({
        "figure.figsize": (7, 6),
        "figure.dpi": 130,
        "mathtext.fontset": "cm",
        "axes.grid": True,
        "grid.alpha": 0.25,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.prop_cycle": mpl.cycler(  # type: ignore[attr-defined]
            color=["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"]
        ),
    })


def plot_roa_2d(
    certificate: Any,
    *,
    ax: Any = None,
    axes: tuple[int, int] = (0, 1),
    fixed: dict[int, float] | None = None,
    true_boundary: Array | None = None,
    labels: Sequence[str] | None = None,
    fill_alpha: float = 0.25,
    show_level_curves: bool = False,
) -> Any:
    """The core figure: certified sublevel set, optionally against ground truth."""
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()

    # Rule 3 of the module docstring: a slice of an n > 2 certificate must state
    # its pinned values. `boundary_2d` defaults unplotted coordinates to 0, which
    # for e.g. the unicycle (x_eq = (0,0,0,1)) is a point off the invariant
    # manifold entirely — an unlabelled slice through c = 0, not a picture of the
    # ROA. Refuse to draw it rather than mislabel it.
    if certificate.n_vars > 2:
        non_axis = set(range(certificate.n_vars)) - set(axes)
        unpinned = sorted(non_axis - set((fixed or {}).keys()))
        if unpinned:
            pinned_hint = ", ".join(f"{k}: ..." for k in unpinned)
            raise ValueError(
                f"plot_roa_2d: certificate has {certificate.n_vars} variables but "
                f"coordinate(s) {unpinned} are not pinned. Pass "
                f"fixed={{{pinned_hint}}} to state the slice, or plot a 2-D certificate."
            )

    boundary = certificate.boundary_2d(axes=axes, fixed=fixed)
    bx, by = boundary[:, 0], boundary[:, 1]

    # Fill certified set.
    ax.fill(bx, by, alpha=fill_alpha, color="C0", label="Certified ROA")
    ax.plot(bx, by, color="C0", linewidth=1.5)

    if show_level_curves:
        # Interior contours at ρ * [0.25, 0.5, 0.75]. These are genuine level
        # sets {V = ρ·frac}, re-bisected along each ray — NOT the boundary scaled
        # radially by sqrt(frac), which reproduces {V = ρ·frac} only for a
        # homogeneous quadratic V and silently under-draws every interior contour
        # of a degree-4/6 candidate.
        for frac in [0.25, 0.5, 0.75]:
            inner = certificate.boundary_2d(
                axes=axes, fixed=fixed, level=certificate.rho * frac
            )
            ax.plot(inner[:, 0], inner[:, 1], color="C0",
                    linewidth=0.5, linestyle="--", alpha=0.5)

    if true_boundary is not None:
        tx, ty = true_boundary[:, 0], true_boundary[:, 1]
        ax.plot(tx, ty, color="C1", linewidth=2, linestyle="--",
                label="Estimated true ROA")

        # Detect leakage: certified points outside true boundary.
        from matplotlib.path import Path as MplPath
        true_path = MplPath(np.column_stack([tx, ty]))
        inside = true_path.contains_points(np.column_stack([bx, by]))
        if not np.all(inside):
            leaking = np.column_stack([bx[~inside], by[~inside]])
            ax.scatter(leaking[:, 0], leaking[:, 1], color="red", s=10, zorder=5,
                       label="Leakage (FALSE CERTIFICATE)")
            import warnings
            warnings.warn("Certified set extends outside true ROA boundary — "
                          "this is a false certificate!", stacklevel=2)

    ax.set_xlabel(labels[axes[0]] if labels else f"$x_{axes[0]+1}$")
    ax.set_ylabel(labels[axes[1]] if labels else f"$x_{axes[1]+1}$")

    if fixed:
        pinned = ", ".join(f"$x_{{{k+1}}}$={v}" for k, v in sorted(fixed.items()))
        ax.set_title(f"Slice at {pinned}")
    else:
        ax.set_title("Region of Attraction")

    ax.set_aspect("equal")
    ax.legend()
    return ax


def plot_phase_portrait(
    dynamics: Any,
    *,
    ax: Any = None,
    axes: tuple[int, int] = (0, 1),
    xlim: tuple[float, float] = (-3.0, 3.0),
    ylim: tuple[float, float] = (-3.0, 3.0),
    density: float = 1.2,
    trajectories: Sequence[Array] = (),
    fixed: dict[int, float] | None = None,
    n_vars: int | None = None,
) -> Any:
    """Streamplot of the vector field, with optional overlaid trajectories.

    ``n_vars`` is the state dimension; it is *not* derived from ``len(fixed) + 2``
    because a caller may pin only some of the non-axis coordinates (and the
    ``dynamics`` callable then receives a wrong-length vector). Pass it
    explicitly; it defaults to ``len(fixed) + 2`` only for convenience when every
    non-axis coordinate is pinned.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()

    if fixed is None:
        fixed = {}
    if n_vars is None:
        n_vars = len(fixed) + 2
    else:
        # Guard against the exact bug this argument exists to fix: a caller who
        # pins only some coordinates and forgets n_vars.
        unpinned = [k for k in range(n_vars) if k not in set(axes) and k not in fixed]
        if unpinned:
            raise ValueError(
                f"plot_phase_portrait: n_vars={n_vars} but coordinate(s) {unpinned} "
                f"are neither pinned in `fixed` nor among `axes`. The dynamics "
                f"callable would receive a vector with zeros in those entries. "
                f"Pin them or pass the true n_vars."
            )

    x_vals = np.linspace(xlim[0], xlim[1], int(30 * density))
    y_vals = np.linspace(ylim[0], ylim[1], int(30 * density))
    X, Y = np.meshgrid(x_vals, y_vals)

    # Evaluate dynamics on grid.
    U = np.zeros_like(X)
    Vf = np.zeros_like(Y)
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            pt = np.zeros(n_vars)
            pt[axes[0]] = X[i, j]
            pt[axes[1]] = Y[i, j]
            for k, v in fixed.items():
                pt[k] = v
            dx = dynamics(pt)
            U[i, j] = dx[axes[0]] if hasattr(dx, "__getitem__") else dx[0]
            Vf[i, j] = dx[axes[1]] if hasattr(dx, "__getitem__") else dx[1]

    speed = np.sqrt(U**2 + Vf**2)
    ax.streamplot(X, Y, U, Vf, color="grey", linewidth=0.8 * speed / (speed.max() + 1e-8),
                  density=density, arrowsize=0.5)

    for traj in trajectories:
        ax.plot(traj[:, 0], traj[:, 1], color="C3", linewidth=1.5, alpha=0.7)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel(f"$x_{axes[0]+1}$")
    ax.set_ylabel(f"$x_{axes[1]+1}$")
    ax.set_aspect("equal")
    return ax


def plot_vdot_field(
    certificate: Any,
    dynamics: Any,
    *,
    ax: Any = None,
    axes: tuple[int, int] = (0, 1),
    fixed: dict[int, float] | None = None,
    n: int = 300,
) -> Any:
    """Heatmap of V̇ with the V̇ = 0 contour and the certified boundary."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    if ax is None:
        _, ax = plt.subplots()

    # Estimate bounds from the certificate boundary.
    boundary = certificate.boundary_2d(axes=axes, fixed=fixed)
    margin = 0.2
    x_min, x_max = boundary[:, 0].min() - margin, boundary[:, 0].max() + margin
    y_min, y_max = boundary[:, 1].min() - margin, boundary[:, 1].max() + margin

    xs = np.linspace(x_min, x_max, n)
    ys = np.linspace(y_min, y_max, n)
    X, Y = np.meshgrid(xs, ys)
    Vdot = np.zeros_like(X)

    # Evaluate V̇ on the grid.
    for i in range(n):
        for j in range(n):
            pt = np.zeros(certificate.n_vars)
            pt[axes[0]] = X[i, j]
            pt[axes[1]] = Y[i, j]
            if fixed:
                for k, v in fixed.items():
                    pt[k] = v
            # V̇ = ∇V · f
            eps = 1e-6
            v_plus = certificate.evaluate(pt + eps * np.eye(certificate.n_vars)[axes[0]])
            v_minus = certificate.evaluate(pt - eps * np.eye(certificate.n_vars)[axes[0]])
            dv_dx = (v_plus - v_minus) / (2 * eps)
            v_plus = certificate.evaluate(pt + eps * np.eye(certificate.n_vars)[axes[1]])
            v_minus = certificate.evaluate(pt - eps * np.eye(certificate.n_vars)[axes[1]])
            dv_dy = (v_plus - v_minus) / (2 * eps)
            dx = dynamics(pt)
            Vdot[i, j] = dv_dx * dx[axes[0]] + dv_dy * dx[axes[1]]

    vmax = max(abs(Vdot.min()), abs(Vdot.max()))
    norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
    im = ax.pcolormesh(X, Y, Vdot, cmap="RdBu_r", norm=norm, shading="auto")
    plt.colorbar(im, ax=ax, label=r"$\dot{V}$")

    # V̇ = 0 contour.
    ax.contour(X, Y, Vdot, levels=[0], colors="black", linewidths=1.5, linestyles="--")

    # Certified boundary.
    ax.plot(boundary[:, 0], boundary[:, 1], color="C0", linewidth=2, label="Certified boundary")
    ax.set_xlabel(f"$x_{axes[0]+1}$")
    ax.set_ylabel(f"$x_{axes[1]+1}$")
    ax.set_aspect("equal")
    ax.legend()
    return ax


def plot_roa_comparison(
    certificates: Sequence[Any],
    *,
    labels: Sequence[str],
    ax: Any = None,
    true_boundary: Array | None = None,
    axes: tuple[int, int] = (0, 1),
    fixed: dict[int, float] | None = None,
) -> Any:
    """Several certified sets on one axis — the degree-2/4/6 comparison."""
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()

    if true_boundary is not None:
        ax.plot(true_boundary[:, 0], true_boundary[:, 1], color="black",
                linewidth=2, linestyle="--", label="Estimated true ROA")

    for i, cert in enumerate(certificates):
        boundary = cert.boundary_2d(axes=axes, fixed=fixed)
        bx, by = boundary[:, 0], boundary[:, 1]
        color = f"C{i}"
        if i == len(certificates) - 1:
            ax.fill(bx, by, alpha=0.15, color=color)
        ax.plot(bx, by, color=color, linewidth=1.5, label=labels[i])

    ax.set_xlabel(f"$x_{axes[0]+1}$")
    ax.set_ylabel(f"$x_{axes[1]+1}$")
    ax.set_aspect("equal")
    ax.legend()
    return ax


def plot_vs_convergence(history: dict[str, Any], *, fig: Any = None) -> Any:
    """Two-panel V–s iteration history: level/shape level, and estimated volume."""
    import matplotlib.pyplot as plt

    if fig is None:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    else:
        ax1, ax2 = fig.subplots(1, 2)

    iters = np.arange(len(history["rho_history"]))
    ax1.plot(iters, history["rho_history"], "o-", color="C0", label=r"$\rho$")
    ax1.plot(iters, history["beta_history"], "s-", color="C1", label=r"$\beta$")
    ax1.set_xlabel("Iteration")
    ax1.set_ylabel("Level")
    ax1.legend()

    # Mark convergence point.
    if history.get("converged", False):
        ax1.axvline(history["iterations"] - 1, color="green", linestyle=":", alpha=0.5)
    else:
        ax1.annotate("budget exhausted", xy=(iters[-1], history["rho_history"][-1]),
                     xytext=(20, 20), textcoords="offset points",
                     arrowprops=dict(arrowstyle="->"))

    ax2.plot(iters, history["volume_history"], "o-", color="C2")
    ax2.set_xlabel("Iteration")
    ax2.set_ylabel("Estimated volume")
    if history.get("converged", False):
        ax2.axvline(history["iterations"] - 1, color="green", linestyle=":", alpha=0.5)

    fig.tight_layout()
    return fig


def plot_consensus_tube(
    certificate: Any,
    trajectories: Sequence[Array],
    *,
    fig: Any = None,
    U: Array | None = None,
) -> Any:
    """Consensus figure: agent states over time plus the certified disagreement set."""
    import matplotlib.pyplot as plt

    if fig is None:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    else:
        ax1, ax2 = fig.subplots(1, 2)

    # Left panel: agent states over time.
    for i, traj in enumerate(trajectories):
        t_vals = np.arange(len(traj))
        ax1.plot(t_vals, traj, linewidth=1, alpha=0.7, label=f"Agent {i+1}")
    ax1.set_xlabel("Time step")
    ax1.set_ylabel("State")
    ax1.set_title("Agent states converging to consensus")
    ax1.legend()

    # Right panel: disagreement coordinates with certified set.
    if U is not None:
        boundary = certificate.boundary_2d()
        ax2.fill(boundary[:, 0], boundary[:, 1], alpha=0.2, color="C0",
                 label="Certified disagreement set")
        ax2.plot(boundary[:, 0], boundary[:, 1], color="C0", linewidth=1.5)

        for traj in trajectories:
            z_traj = U.T @ traj.T  # Project into disagreement coordinates
            ax2.plot(z_traj[0, :], z_traj[1, :], linewidth=0.5, alpha=0.6)

    ax2.set_xlabel(r"$z_1$")
    ax2.set_ylabel(r"$z_2$")
    ax2.set_title("Disagreement coordinates")
    ax2.set_aspect("equal")
    ax2.legend()

    fig.tight_layout()
    return fig


def save_figure(
    fig: Any,
    name: str,
    *,
    directory: str | Path | None = None,
    dpi: int = 200,
) -> Path:
    """Save to ``docs/figures/<name>.png`` (default) and return the path."""
    if directory is None:
        directory = Path(__file__).resolve().parents[2] / "docs" / "figures"
    else:
        directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{name}.png"
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    return path
