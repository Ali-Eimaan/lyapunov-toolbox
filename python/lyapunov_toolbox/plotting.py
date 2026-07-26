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
    raise NotImplementedError


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
    """The core figure: certified sublevel set, optionally against ground truth.

    Draws ``{V <= rho}`` filled with a solid edge. When ``true_boundary`` is
    given (an ``(m, 2)`` array from ``true_roa_boundary_2d``) it is drawn dashed
    and unfilled, and the legend reports the area ratio.

    ``show_level_curves`` adds a few interior contours of ``V`` at
    ``rho * [0.25, 0.5, 0.75]``, which is how a reader sees the shape of the
    Lyapunov function rather than just its one certified level.

    If any part of the certified set falls outside ``true_boundary``, that is a
    false certificate — draw the offending region in red and emit a warning
    rather than quietly plotting it. A figure is the last place that error should
    be allowed to look normal.
    """
    raise NotImplementedError


def plot_phase_portrait(
    dynamics: Any,
    *,
    ax: Any = None,
    xlim: tuple[float, float] = (-3.0, 3.0),
    ylim: tuple[float, float] = (-3.0, 3.0),
    density: float = 1.2,
    trajectories: Sequence[Array] = (),
    fixed: dict[int, float] | None = None,
) -> Any:
    """Streamplot of the vector field, with optional overlaid trajectories.

    Background layer for :func:`plot_roa_2d`. Streamlines are drawn thin and grey
    so the certified set stays the subject of the figure.

    ``dynamics`` is a callable ``x -> ẋ`` on ``(k, n)`` arrays. Accepting a
    callable rather than a certificate keeps this usable for the systems in the
    other portfolio repositories.
    """
    raise NotImplementedError


def plot_vdot_field(
    certificate: Any,
    dynamics: Any,
    *,
    ax: Any = None,
    axes: tuple[int, int] = (0, 1),
    fixed: dict[int, float] | None = None,
    n: int = 300,
) -> Any:
    """Heatmap of ``V̇`` with the ``V̇ = 0`` contour and the certified boundary.

    The diagnostic plot: the certified boundary must lie strictly inside the
    region where ``V̇ < 0``, and seeing the two curves nearly touch is what tells
    the reader the bisection actually converged on the binding constraint rather
    than stopping early.

    Use a diverging colormap centred at zero (``TwoSlopeNorm(vcenter=0)``) —
    with a sequential map the sign change, which is the entire content of the
    figure, is invisible.
    """
    raise NotImplementedError


def plot_roa_comparison(
    certificates: Sequence[Any],
    *,
    labels: Sequence[str],
    ax: Any = None,
    true_boundary: Array | None = None,
    axes: tuple[int, int] = (0, 1),
    fixed: dict[int, float] | None = None,
) -> Any:
    """Several certified sets on one axis — the degree-2/4/6 comparison.

    Nested boundaries with increasing degree, no fill except on the largest, so
    the nesting stays legible. The legend reports each certified area and its
    fraction of the true ROA. That table is the quantitative claim the figure
    exists to support; without the numbers it is three coloured curves.
    """
    raise NotImplementedError


def plot_vs_convergence(history: dict[str, Any], *, fig: Any = None) -> Any:
    """Two-panel V–s iteration history: level/shape level, and estimated volume.

    Mark the iteration at which the improvement fell below tolerance. When
    ``history["converged"]`` is ``False``, annotate the final point with
    "budget exhausted" — a convergence plot that ends at an iteration cap and
    does not say so reads as convergence.
    """
    raise NotImplementedError


def plot_consensus_tube(
    certificate: Any,
    trajectories: Sequence[Array],
    *,
    fig: Any = None,
    U: Array | None = None,
) -> Any:
    """Consensus figure: agent states over time plus the certified disagreement set.

    Two panels. Left: ``xᵢ(t)`` for each agent, converging to a common value.
    Right: the disagreement coordinates ``z = Uᵀx`` with the certified
    ``{V(z) <= rho}`` drawn behind the trajectories.

    The right panel is the one that carries the argument — the left panel alone
    shows agents agreeing, which any stable consensus protocol does, certified or
    not.
    """
    raise NotImplementedError


def save_figure(fig: Any, name: str, *, directory: str | Path | None = None, dpi: int = 200) -> Path:
    """Save to ``docs/figures/<name>.png`` (default) and return the path.

    Resolves the default directory relative to this file, not the working
    directory, so a notebook run from anywhere writes to the same place. Creates
    the directory if needed. Saves with ``bbox_inches="tight"``.

    The four README heroes are ``roa_van_der_pol``, ``roa_inverted_pendulum``,
    ``roa_unicycle_obstacle``, ``roa_consensus``; ``.gitignore`` ignores
    ``docs/figures/`` except for those four.
    """
    raise NotImplementedError
