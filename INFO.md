# GitHub Portfolio — Detailed Specification & PhD-Readiness Analysis

**Owner:** Pakistani robotics PhD candidate (MS Mechatronics, 5 yrs embedded/robotics, ROS 2 Jazzy + PX4 + Gazebo)
**Thesis target:** *Transition-Viable Distributed MPC-CBF for Aerial-Ground Swarms: A Hybrid-Systems Theory of Asynchronous Split, Merge, and Morph Events with Set-Valued Resets*
**Application target:** Fall 2027 PhD intake at hardware-equipped labs (Tier A: KTH, ETH, TU Delft, TU/e; Tier B: Waterloo, UTIAS, Polytechnique Montréal, NTNU; Tier C: HKUST, HKU, NTU, UNIST)
**Document purpose:** Define the full repository structure, file-level contents, and identify any gaps that need closing before December 2026 application submissions.

---

## How to read this document

For each repo I specify:
1. **Identity** — name, one-line tagline, primary language, lines-of-code estimate at v1.0 release.
2. **Purpose** — what it proves to a hardware-equipped advisor reading your GitHub.
3. **File-level structure** — every directory and the files in it, with one-sentence descriptions of what each non-boilerplate file contains.
4. **Demo deliverables** — videos, GIFs, plots that go in the README.
5. **CI / quality signals** — what GitHub Actions runs and what badges go on the README.
6. **Citation hook** — how this repo links back to the TVF/RRC-CBF/AHTD thesis.

At the end I do a hard portfolio gap analysis against what hardware-equipped PhD admissions committees actually look for.

# Repo 2 — `lyapunov-toolbox`

**Tagline:** Sum-of-Squares Lyapunov function synthesis and region-of-attraction estimation in Python (SumOfSquares.jl backend) and MATLAB (YALMIP+Mosek).
**Language:** Julia (40%) + Python (40%) + MATLAB (20%).
**Estimated size at v1.0:** ~1,500 LOC.
**Why this exists:** Proves you can do *certified* control, not just simulation. Sanfelice/Dimarogonas/Lygeros students all use SOS; if you can't do this you're locked out of the theoretical tier.

## Directory structure

```
lyapunov-toolbox/
├── .github/workflows/test.yml
├── julia/
│   ├── Project.toml
│   ├── Manifest.toml
│   └── src/
│       ├── LyapunovToolbox.jl
│       ├── sos_lyapunov.jl
│       ├── roa_estimation.jl
│       └── polynomial_utils.jl
├── python/
│   ├── lyapunov_toolbox/
│   │   ├── __init__.py
│   │   ├── julia_bridge.py
│   │   └── plotting.py
│   └── pyproject.toml
├── matlab/
│   ├── compute_roa_yalmip.m
│   ├── van_der_pol_roa.m
│   └── pendulum_roa.m
├── examples/
│   ├── van_der_pol/
│   │   ├── van_der_pol_demo.jl
│   │   └── van_der_pol_demo.ipynb
│   ├── inverted_pendulum/
│   ├── unicycle_with_obstacle/
│   └── multi_agent_consensus/
├── tests/
│   ├── test_sos_lyapunov.jl
│   └── test_roa_estimation.jl
├── docs/
│   ├── theory_notes.md
│   ├── derivations/
│   │   ├── sos_decomposition.tex
│   │   └── positivstellensatz.tex
│   └── benchmarks.md
├── LICENSE
└── README.md
```

## File-level descriptions

- **`julia/src/sos_lyapunov.jl`** — given a polynomial vector field `f(x)` and a candidate polynomial Lyapunov template, returns the largest sublevel set `{x : V(x) ≤ ρ}` certified by SOS to be inside the ROA. Returns the SOS certificate as a `Polynomial` object so downstream code can verify or re-use it.
- **`julia/src/roa_estimation.jl`** — bisection over `ρ` with the V-s iteration to find the largest invariant sublevel set.
- **`julia/src/polynomial_utils.jl`** — basis manipulation, monomial enumeration, polynomial arithmetic helpers.
- **`python/lyapunov_toolbox/julia_bridge.py`** — thin PyJulia wrapper; lets the rest of the portfolio call SOS routines from Python.
- **`matlab/van_der_pol_roa.m`** — reproduces the canonical VDP ROA result from Tedrake's MIT 6.832; a sanity check that the toolbox is correct.
- **`examples/multi_agent_consensus/`** — a 3-agent consensus example where you certify a quadratic Lyapunov function for the closed-loop network, **as a stepping stone to RRC-CBF compositions**.
- **`docs/derivations/sos_decomposition.tex`** — Parrilo's SOS decomposition theorem in your own notation, with one fully-worked example.
- **`docs/derivations/positivstellensatz.tex`** — Putinar's Positivstellensatz, the workhorse you'll cite in your TVF paper.

## Demo deliverables

- 4 ROA plots in README: VDP, pendulum, unicycle, 3-agent network.

## CI

- Julia test workflow on every commit using `runtests.jl`.
- Python integration test verifying PyJulia bridge.

## Citation hook

The SOS routines here are *exactly* what you call to produce inner approximations of TVFs in `transition-viable-swarm`. In your arXiv preprint you'll say "all SOS computations were performed using the open-source `lyapunov-toolbox` library available at github.com/<you>/lyapunov-toolbox at commit `abc123`."

---