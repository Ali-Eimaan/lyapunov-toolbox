"""Self-links in the README and package metadata use the current GitHub username.

The old account name only resolves through GitHub's rename redirect. These checks
look at URL path segments, so a display name that still uses the old spelling is
left for the owner to decide.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# The retired slug as its own GitHub path segment. Display names are not URLs
# and are outside this check.
OLD_USERNAME_URL = re.compile(
    r"https?://github\.com/ali-eiman(?![A-Za-z0-9])",
    re.IGNORECASE,
)


def test_readme_and_pyproject_have_no_old_username_self_links():
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    pyproject = (ROOT / "python" / "pyproject.toml").read_text(encoding="utf-8")
    urls = pyproject.split("[project.urls]", 1)[1].split("\n[", 1)[0]

    for label, text in (
        ("README.md", readme),
        ("python/pyproject.toml [project.urls]", urls),
    ):
        found = OLD_USERNAME_URL.findall(text)
        assert found == [], f"{label} still links the old GitHub username: {found}"
