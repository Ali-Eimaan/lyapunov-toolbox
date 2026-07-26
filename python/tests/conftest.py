"""Shared fixtures for the Python test suite.

The ``julia`` marker is handled here rather than in each test: tests that need a
live backend skip with a reason naming the resolved project path, so that a
skipped run tells the reader what to fix.
"""

from __future__ import annotations

import pytest


def pytest_collection_modifyitems(config, items):
    """Skip ``@pytest.mark.julia`` tests when the backend cannot start.

    Call ``julia_available()`` **once** here and cache the result on ``config`` —
    calling it per test starts and probes Julia dozens of times and turns a
    30-second suite into a five-minute one.

    In CI, ``LYAPUNOV_TOOLBOX_REQUIRE_JULIA=1`` is set: under that flag, an
    unavailable backend must *fail* collection rather than skip. A silently
    skipped bridge suite is the failure mode this whole file exists to prevent.
    """
    raise NotImplementedError


@pytest.fixture(scope="session")
def julia_session():
    """The shared :class:`JuliaSession`, started once per test session."""
    raise NotImplementedError


@pytest.fixture(scope="session")
def vdp_certificate(julia_session):
    """A verified reversed-Van-der-Pol certificate, computed once and reused.

    Session-scoped because it costs a real SDP solve. Tests must treat it as
    immutable — :class:`RoaCertificate` is a frozen dataclass, which is the other
    half of why that decision was made.
    """
    raise NotImplementedError


@pytest.fixture
def saved_certificate(tmp_path, vdp_certificate):
    """Path to ``vdp_certificate`` written to a temporary JSON file."""
    raise NotImplementedError


@pytest.fixture
def agg_backend():
    """Force the Matplotlib ``Agg`` backend for the duration of a test."""
    raise NotImplementedError
