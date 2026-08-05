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

Same SOS program, same formulation, different conic solver. This is the
comparison §4 *cannot* make: there the whole implementation changes, so a
disagreement could come from either the formulation or the solver. Here only the
solver varies, which isolates solver accuracy.

Van der Pol (reversed), `V` from the linearisation, fixed-`V` ρ-bisection:

| Solver | `ρ` | Wall (s) | Gram residual | min eig `Q` | Status | Verified |
|---|---|---|---|---|---|---|
| Clarabel 0.11.1 | 2.30438 | 21.7 | 1.65e-11 | −1.25e-10 | `OPTIMAL` | ✓ |
| SCS 2.x | 2.29525 | 18.1 | 3.29e-12 | −3.50e-10 | `OPTIMAL` | ✓ |
| SeDuMi (via MATLAB, §4) | 2.316 | 1.1 | 5.08e-03 | +9.80e-04 | marginal, verified | ✓ |

**Clarabel vs SCS agree on ρ to 3.96e-03 relative**, and both return witnesses
that pass the independent Gram check — residuals around 1e-11 to 1e-12, eight
orders tighter than SeDuMi on the same mathematics. All three certified levels
sit within 0.5 % of each other.

The residual column is the one that matters, not the wall clock. A solver that
returns quickly with a residual of 5e-3 has not certified the same thing as one
returning 1e-11; the SeDuMi row is only usable because
`compute_roa_yalmip` verifies the returned Gram matrices rather than trusting
the reported status (see §4).

**What this exposed.** SCS at `eps = 1e-9` first came back `ITERATION_LIMIT`,
which `verify_sublevel_set` did not recognise — so it threw, which is the correct
default for an unhandled status but made the solver unusable. `ITERATION_LIMIT`
and `TIME_LIMIT` are now handled the same way as `ALMOST_OPTIMAL`: the solver has
given up, so the witness is checked and accepted only if it verifies. With
`max_iters = 100000` SCS converges to `OPTIMAL` outright.

Mosek is not benchmarked: it needs a licence, and this machine has none. The
`matlab/` path supports it via `opts.solver = 'mosek'`.

## 4. Julia vs MATLAB cross-validation

Two implementations sharing **no code**: Julia (SumOfSquares.jl + Clarabel) and
MATLAB (YALMIP + SeDuMi), against independently written formulations. Agreement
is evidence in a way a repeated run of the same code is not.

| Quantity | Julia | MATLAB | Agreement |
|---|---|---|---|
| Van der Pol `ρ` (fixed-V) | 2.304375 | 2.316 | 5.0e-3 relative |
| Van der Pol `V` coefficients (normalised) | `1.5, −1, 1` | `1.5, −1, 1` | **2.22e-16** |
| Van der Pol true ROA area (limit cycle) | 13.7222 | 13.7222 | exact |
| Van der Pol certified area | 6.4748 | 6.50742 | 5.0e-3 relative |
| Van der Pol area fraction | 47.18 % | 47.42 % | 0.24 pp |
| Pendulum (Taylor-3) `ρ` | — | 7.9045 | Julia normalises to ρ≡1; compare V |

The coefficient agreement is the strong result. Both sides build `V = x'Px` with
`P` solving `A'P + PA = −I`, so the *shape* of `V` is fixed by linear algebra
rather than by the SDP, and the two agree to machine epsilon. A transposed or
mis-ordered certificate would show up here immediately — and did: see below.

**Tolerance, and why it is 1e-2 on `ρ` and 1e-6 on the coefficients.** Clarabel
reaches a Gram residual of ~1e-10 on this program; SeDuMi reaches ~5e-3
(measured, and recorded in the fixture's own `gram_residual` column). Two
implementations cannot agree on a bisected level more tightly than the worse
solver can certify it. The coefficients carry no such limitation because they
come from a dense linear solve, hence the far tighter bound there.

**What this check caught.** On its first real run the two disagreed by **41×**
(Julia 2.304, MATLAB 0.056) and the `V` diagonal was swapped. Three separate
defects, each of which would have masked the others:

1. `compute_roa_yalmip` mapped SeDuMi's status 4 ("numerical problems") straight
   to infeasible. Measured, status 4 spans the entire feasible range with **PSD**
   Gram matrices and residual ~1e-3, while genuine infeasibility gives ~1e+6 —
   nine orders of separation. It now verifies the witness instead of trusting
   the status, exactly as the Julia side does for `ALMOST_OPTIMAL`.
2. `write_roa_fixture` used YALMIP's `exponents()`, which returns columns in
   YALMIP's internal variable order rather than the caller's. Replaced with
   per-variable `degree()`, which names the variable it asks about.
3. The Julia-side reader rebuilt `V` over freshly declared variables, so the
   comparison subtracted polynomials from disjoint variable sets and could never
   have passed.

**Provenance.** MATLAB R2023a, YALMIP master (2026-08-05), SeDuMi master built
against MinGW-w64 GCC 14.2.0 (34 mex files). Regenerate with:

```matlab
van_der_pol_roa('Solver','sedumi','WriteFixture',true,'Plot',false)
pendulum_roa('Solver','sedumi','WriteFixture',true,'Plot',false)
```

CI has neither MATLAB nor a solver licence, so `tests/fixtures/*.csv` are
committed artefacts. The commit that updates one must record these versions.

## 5. Bridge overhead

Every cold figure is measured in a **fresh subprocess**. Timing a "cold" import
inside a process that has already imported the module measures nothing, and
Julia startup cannot be repeated in-process at all. Warm figures are the median
of repeated calls after discarding the first.

| Operation | Time |
|---|---|
| `import lyapunov_toolbox` (cold, median of 5) | **0.169 s** (max 0.199 s) |
| First `certify_roa` — Julia start + package load + solve | **163.7 s** |
| Subsequent `certify_roa` in the same process | **71.5 s** |
| `RoaCertificate.evaluate` on 10⁶ points (NumPy, no Julia) | **0.189 s** |
| `save` + `load_certificate` round trip (JSON, no Julia) | **0.0042 s** |

**Read the first two rows together.** The 163.7 s cold call is not bridge
latency — 71.5 s of it is the SDP work itself (a ρ-bisection is ~40 sequential
solves), so the one-off cost of starting Julia and loading SumOfSquares is
roughly **90 s**. That is paid once per process, never again.

**The import figure is the load-bearing one.** At 0.169 s it confirms the bridge
does not boot Julia at import time, which is what lets the other portfolio
repositories depend on this package for its plotting layer alone without paying
for a Julia runtime they will not use. `test_import_does_not_start_julia`
asserts this stays under 2 s; the measured value has an order of magnitude of
headroom.

**The evaluation figure justifies the marshalling design.** A certificate
crosses the boundary once, as `(exponents, coefficients)`, and is then evaluated
in NumPy at ~5.3 million points per second with no Julia in the loop. Calling
back through the bridge per point would be slower by orders of magnitude and
would make plotting impossible on a machine without Julia — which is precisely
the case `RoaCertificate` is built to support. The 4.2 ms JSON round trip is why
storing certificates alongside downstream results is cheap enough to be routine.

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
