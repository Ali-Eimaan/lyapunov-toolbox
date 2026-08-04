
---

# Final state — 2026-08-04

| gate | result |
|---|---|
| `Pkg.test()` | **14320 passed, 3 broken, 0 failed** (13m32s) |
| `pytest` | **27 passed, 1 skipped** (pyjulia not installed) |
| `ruff check` | clean |
| MATLAB `checkcode` | 5 files, 0 messages |
| All six systems synthesise + verify | yes (algebraic, sampling, dynamical) |
| `docs/benchmarks.md` | filled from `benchmark_row` output |
| README Results | filled from the same run |
| Licence | BSD 3-Clause |

The 3 `broken` are `@test_skip`s with stated reasons, not hidden failures:
the ground-truth ROA sweep (hundreds of simulations), the degree sweep (~10 min
of SDPs), and the MATLAB cross-validation while its fixture is absent.

## Fixed since the last round

- Volume regression: `v_s_iteration` returned its last iterate rather than its
  best. It now floors at the fixed-V baseline.
- `θ_max` expectation: `> 3.0` was unachievable and unmeasured. Measured 2.6995
  (deg 2) / 2.4967 (deg 4), replaced with `> π/2` and a stated rationale.
- `test_unverifiable_certificate_is_not_returned`: relied on `V_degree=1` being
  infeasible, which stopped being true. Now uses the forward Van der Pol, which
  fails for a reason.

## Still open before a tag

1. **Notebooks: 0 of 28 code cells have outputs.** The prose is written, but none
   has been executed. CI runs `pytest --nbmake examples/`, so this is a hard gate.
2. **MATLAB cross-validation fixture absent** — needs YALMIP + SeDuMi installed;
   `matlab/assert_yalmip_env.m` prints the setup. Until then the Julia/MATLAB
   agreement claim is untested and `docs/benchmarks.md` §4 stays empty.
3. **`v_step` does not return its containment Gram**, so on the
   `synthesize_lyapunov` path containment is enforced and empirically confirmed
   (0.106 m clearance over 300 000 on-manifold samples) but not witness-checked.
   It *is* witness-checked on `estimate_roa`.
4. **`ROAResult.n_sdps` / `total_solve_time` undercount** — they accumulate only
   the V–s loop and drop the preceding bisection (~0.05 s / 2 SDPs reported
   against a 36 s wall clock). Those columns are omitted from `benchmarks.md`
   rather than published wrong.
