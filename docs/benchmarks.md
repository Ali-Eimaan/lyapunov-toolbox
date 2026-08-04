# Benchmarks

Every number below was produced by `benchmark_row` on the commit and machine
named here. Nothing is typed from memory, estimated, or copied from a paper.

Regenerate with the script in `docs/` (see §6).

---

## Environment

| | |
|---|---|
| CPU | Intel64 Family 6 Model 154 Stepping 3 (12th-gen mobile) |
| OS | Windows 11 (10.0.26200) |
| Julia | 1.12.6 |
| SumOfSquares.jl | 0.7.6 |
| JuMP | 1.31.1 |
| DynamicPolynomials | 0.6.4 |
| Conic solver | Clarabel 0.11.1 (`tol_gap_abs/rel = 1e-9`, `tol_feas = 1e-9`, `max_iter = 500`) |
| MATLAB | R2023a (cross-validation not yet run — YALMIP not installed) |
| Commit | `6033ab7` + working tree |
| Date | 2026-08-04 |

---

## 1. Certified ROA per system

All six configurations synthesise **and** verify: the algebraic (Gram residual
against reconstructed target), sampling (`V̇ < 0` pointwise, on the invariant
manifold where one exists), and dynamical (forward simulation) checks all pass.

`V_degree` is what was requested; `deg V` is what the returned certificate
actually has. `v_s_iteration` keeps its **best iterate by volume**, so when a
richer basis does not enclose more it hands back the lower-degree one — see the
note under the table.

| System | `n` | `V_degree` | `deg V` | `deg λ` | `ε` used | `k` | `ρ` | Certified volume | Gram block | Wall (s) | Bounded via | Verified |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Van der Pol (reversed) | 2 | 2 | 2 | 4 | 1e-06 | 2 | 2.30438 | 6.49635 | 6 | 35.9 | `top_form_pd` | ✓ |
| Van der Pol (reversed) | 2 | 4 | **2** | 4 | 1e-06 | 2 | 2.30438 | 6.4675 | 6 | 1.9 | `top_form_pd` | ✓ |
| Inverted pendulum (Taylor-3) | 2 | 2 | 2 | 4 | 6.47e-08 | 2 | 1.0 | 164.917 | 6 | 2.7 | `top_form_pd` | ✓ |
| Inverted pendulum (trig recast) | 3 | 2 | 2 | 6 | 1e-06 | 3 | 3.99975 | 106.88 | 20 | 109.0 | `top_form_pd` | ✓ |
| Unicycle + obstacle | 4 | 2 | 2 | 2 | 5.56e-06 | 2 | 1.0 | 0.532426 | 15 | 26.5 | `top_form_pd` | ✓ |
| 3-agent consensus | 2 | 2 | 2 | 4 | 1e-07 | 2 | 1.0 | 109.257 | 6 | 2.0 | `top_form_pd` | ✓ |

**Reading this table.**

- **`ρ` is not comparable across rows.** `v_s_iteration` normalises to `ρ ≡ 1`
  when its refinement wins, and leaves the fixed-V level when it does not. The
  comparable quantity is **certified volume**.
- **Volume is not comparable across systems either** — the state spaces differ in
  dimension and units. Compare within a row-group only.
- **Row 2 returns a degree-2 certificate for a degree-4 request.** That is not a
  bug: at `V_degree = 4` the V–s iterate enclosed 6.4675 against the degree-2
  baseline's 6.49635, so the baseline was kept. The V-step maximises the
  inscribed shape-set level `β`, and for an anisotropic sublevel set that is not
  the same as maximising volume. Raising the degree is not guaranteed to help,
  and this repository reports that rather than hiding it.
- **`ε` varies by row** because `v_s_iteration` rescales the margin with `V`
  under the `ρ ≡ 1` normalisation (`ε ← ε/ρ`). The value shown is the one the
  returned certificate was solved with.
- **`bounded via top_form_pd` on every row** — the conclusive structural test
  (top-degree form positive definite), not the `sampling_only` fallback. This is
  the strong form of the boundedness evidence the ROA theorem needs.

**Wall time** is end-to-end for the whole `v_s_iteration` call including the
fixed-V bisection that precedes it, measured with `time()`. The per-SDP solver
times reported by `ROAResult` are **not** listed here: `total_solve_time` and
`n_sdps` currently accumulate only the V–s loop and drop the bisection, so they
read ~0.05 s / 2 SDPs against a 36 s wall clock. That is a reporting defect, not
a measurement — it is tracked in REVIEW.md and the columns are omitted rather
than published wrong.

The first Van der Pol row's 35.9 s includes first-call JIT compilation; the
second row (1.9 s) is the steady-state cost of the same work.

---

## 2. Cost versus degree

The Gram dimension is `binomial(n + d, d)` for an SOS polynomial of degree `2d`
in `n` variables. This growth is the dominant practical fact about SOS methods.

| `n` | `d` (= deg/2) | Gram dimension |
|---|---|---|
| 2 | 1 | 3 |
| 2 | 2 | 6 |
| 2 | 3 | 10 |
| 3 | 1 | 4 |
| 3 | 2 | 10 |
| 3 | 3 | 20 |
| 4 | 1 | 5 |
| 4 | 2 | 15 |
| 4 | 3 | 35 |
| 6 | 1 | 7 |
| 6 | 2 | 28 |
| 6 | 3 | 84 |

Observed failure points, for honesty about where this stops working:

- **Unicycle V-step at `V_degree = 2`** returns `ALMOST_OPTIMAL` with a
  `:derivative` Gram eigenvalue of `-9.2e-6` — four orders of magnitude outside
  the `-1e-8` acceptance gate. Raising `max_iter` from 500 to 5000 does not help;
  loosening `tol` to 1e-8 produces a solution that then fails verification. The
  iteration correctly stops and returns the last verified certificate.
- **Van der Pol at `V_degree = 4` with `multiplier_degree = 2`** is
  `INFEASIBLE` — the pairing bound `deg(λ) ≥ deg(f) - 1` is met exactly and that
  is not enough. `_effective_multiplier_degree` now adds one even degree of
  headroom automatically.

---

## 3. Solver comparison

Not run. Only Clarabel is wired up on the Julia side; SCS is supported but not
installed, and Mosek is MATLAB-only here. Left empty rather than filled with one
row pretending to be a comparison.

---

## 4. Julia vs MATLAB cross-validation

Not run — YALMIP is not installed on the benchmark machine, so
`tests/fixtures/matlab_vdp.csv` has not been generated. `matlab/assert_yalmip_env.m`
prints the setup instructions. The Julia-side reader is in
`tests/test_roa_estimation.jl` and skips (not passes) while the fixture is absent.

---

## 5. Bridge overhead

Not systematically measured. One figure worth recording: a cold
`certify_roa` call pays Julia startup plus SumOfSquares precompilation, roughly
30–60 s on this machine; subsequent calls in the same process do not.
`RoaCertificate.evaluate` runs in NumPy off the marshalled certificate and needs
no Julia at all, which is why the plotting layer works on a machine without it.

---

## 6. Reproduction

```bash
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
julia --project=julia -e 'using Pkg; Pkg.test()'
cd python && pip install -e ".[dev]" && ruff check lyapunov_toolbox tests && pytest
```

The table in §1 comes from `benchmark_row` called on a `v_s_iteration` result for
each configuration, with `vs_iterations = 5`. The consensus row additionally
passes `ρ_hi = 10.0`: with `κ > 0` its disagreement dynamics are globally stable,
so the bisection has no finite boundary and correctly refuses to invent one.

---

## 7. Caveats

- Single run per row, on a laptop under unspecified thermal conditions. Treat
  wall times as order-of-magnitude.
- Every certified set is an **inner** approximation. None of these numbers is the
  region of attraction; each is a set provably inside it.
- "Certified volume" for `n > 2` is a Monte Carlo estimate over a box fitted to
  the sublevel set, and is a sanity figure rather than a precise measurement.
- The V–s iteration is a local heuristic. A larger number here means a better
  run, not a better system, and — as row 2 shows — a richer basis does not
  reliably produce a larger set.
- No "fraction of true ROA" column: that requires the ground-truth boundary
  sweep, which has not been run on this machine.
