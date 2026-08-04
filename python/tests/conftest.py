"""Shared fixtures for the Python test suite.

The ``julia`` marker is handled here rather than in each test: tests that need a
live backend skip with a reason naming the resolved project path, so that a
skipped run tells the reader what to fix.
"""

from __future__ import annotations

import pytest


def pytest_addoption(parser):
    """Declare ``--require-julia`` so ``getoption`` never falls back to a default.

    ``pytest_collection_modifyitems`` reads this option; without an
    ``addoption`` hook the flag cannot be passed on the command line and the
    "fail instead of skip" CI mode is unreachable.
    """
    parser.addoption(
        "--require-julia",
        action="store_true",
        default=False,
        help="Fail (instead of skip) when the Julia backend is unavailable",
    )


def pytest_collection_modifyitems(config, items):
    """Skip ``@pytest.mark.julia`` tests when the backend cannot start."""
    import os

    from lyapunov_toolbox.julia_bridge import julia_available

    require = (
        config.getoption("--require-julia")
        or os.environ.get("LYAPUNOV_TOOLBOX_REQUIRE_JULIA") == "1"
    )

    if not hasattr(config, "_lyapunov_julia_available"):
        try:
            config._lyapunov_julia_available = julia_available()
        except Exception:
            config._lyapunov_julia_available = False

    if require and not config._lyapunov_julia_available:
        raise RuntimeError(
            "Julia backend required but unavailable "
            "(--require-julia / LYAPUNOV_TOOLBOX_REQUIRE_JULIA=1). "
            "Instantiate the Julia project with "
            "`julia --project=julia -e 'import Pkg; Pkg.instantiate()'`."
        )

    skip_julia = pytest.mark.skip(
        reason="Julia backend not available — set LYAPUNOV_TOOLBOX_REQUIRE_JULIA=1 to fail instead"
    )
    for item in items:
        if "julia" in item.keywords:
            if not config._lyapunov_julia_available:
                item.add_marker(skip_julia)


@pytest.fixture(scope="session")
def julia_session():
    """The shared :class:`JuliaSession`, started once per test session."""
    from lyapunov_toolbox.julia_bridge import get_session
    return get_session()


@pytest.fixture(scope="session")
def vdp_certificate(julia_session):
    """A verified reversed-Van-der-Pol certificate, computed once and reused."""
    from lyapunov_toolbox.julia_bridge import certify_roa
    f_expr = "[-x2, x1 + (x1^2 - 1)*x2]"
    cert = certify_roa(f_expr, variables=("x1", "x2"), V_degree=4)
    if not cert.verified:
        # A skipped fixture makes the *dependent* tests pass vacuously; a
        # certificate that fails verification is a real bug in the SOS code and
        # the suite must say so.
        pytest.fail(
            "Van der Pol certificate failed verification: the SOS program solved "
            "but verify_certificate rejected the result (algebraic/sampling/"
            "dynamical). Check solver tolerances or the multiplier degree."
        )
    return cert


@pytest.fixture
def saved_certificate(tmp_path, vdp_certificate):
    """Path to ``vdp_certificate`` written to a temporary JSON file."""
    path = tmp_path / "certificate.json"
    vdp_certificate.save(str(path))
    return str(path)


@pytest.fixture
def agg_backend():
    """Force the Matplotlib ``Agg`` backend for the duration of a test."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    yield
    plt.close("all")
