# REVIEW.md — implementation review of `lyapunov-toolbox`

**Review round 1.** Reviewed at commit `b9692e4`. Every fix in this document is
numbered `R1-n`; a later review pass becomes `R2-n`, and so on.

The implementation is substantial and several parts are genuinely good — the
Putinar path in `verify_sublevel_set`, the `ALMOST_OPTIMAL` Gram gate, the
tangent-space Lyapunov seed, the `PolynomialSystem` constructor validation, and
the MATLAB `problem == 4 → infeasible` mapping are all correct and were not
obvious.

The problems cluster in three places, and the ordering below reflects that:

1. **Things that crash or are simply absent.** Two functions reference struct
   fields that do not exist; the flagship derivation is a byte-identical copy of
   the other one; the notebooks were never written or executed; the benchmark
   tables are empty.
2. **Soundness.** The certificate verifier is a tautology, the certificate
   carries Gram matrices belonging to a *different* `V`, and the single most
   important containment test asserts nothing. Together these mean **no claim in
   this repository is currently checked by anything**. That matters more here
   than in an ordinary library, because verifiability is the product.
3. **Contract drift.** Several routines silently do something weaker than their
   own docstring says — returning a sentinel where the doc says "throw", accepting
   a sampling heuristic as proof, recording an `epsilon` that was never used.

Nothing below is a matter of taste. Each item is either a crash, a wrong number, a
claim the code does not support, or an explicit contract in the file's own
docstring that the body violates.

**`README.md` and `matlab/` are not covered here — both have already been fixed.**

- **README.** The unsupported result numbers, the broken quickstart, the wrong CI
  badge target, the missing hero images, the wrong install environment variable,
  and the links to uncommitted PDFs were corrected directly. Two follow-ups it
  left behind are R1-43 and R1-44.
- **MATLAB.** `compute_roa_yalmip.m` was rewritten: the duplicate function
  declaration and stray `end` that made the file unparseable are gone, the
  bisection now returns only a level a solver accepted (and re-solves at it), the
  boundedness and (P1) positivity checks the help block promised are implemented,
  the solver is resolved by name up front instead of silently substituted, and
  Gram matrices and residuals are extracted into `info`. `van_der_pol_roa.m` now
  extracts exactly one period of the limit cycle via a Poincaré section — the
  previous multi-period trace made both `polyarea` and `inpolygon` meaningless —
  and tests containment rather than distance. `pendulum_roa.m` gained
  `WriteFixture` and a bounded sampling loop. A shared `write_roa_fixture.m` emits
  the documented schema.

  **These MATLAB changes are only partly verified.** All four files now *parse* —
  MATLAB parses a whole file before executing any of it, and a later run reached
  the first `sdpvar` line in each script, which is proof the duplicate-declaration
  and stray-`end` class of error is gone. Nothing past that has run, because
  YALMIP is not installed on the author's machine. The numerical behaviour, and
  in particular the four-output `solvesos` call and the `exponents` fallback in
  `write_roa_fixture.m`, remain unexercised. Running all three to completion is
  R1-44's sibling and must happen before the cross-validation fixture is
  regenerated.

  That first run also exposed a real defect, now fixed: the environment check
  lived inside `compute_roa_yalmip`, but both scripts call `sdpvar` before they
  ever reach it, so a missing YALMIP surfaced as a bare
  `Unrecognized function or variable 'sdpvar'`. The check now lives in
  [matlab/assert_yalmip_env.m](matlab/assert_yalmip_env.m) and runs at the top of
  every entry point.

## How to use this document

The fixes are grouped into four severity tiers (A–D). **Work the tiers in order** —
tier B items are hard to verify while tier A items still crash, and tier C items
are hard to see while tier B items mask them. Numbering is continuous across
tiers, so `R1-20` is unambiguous on its own.

Each fix carries a model assignment:

- **`deepseek-v4-pro`** — needs mathematical reasoning, SDP formulation, numerical
  analysis, or a judgement call about what the certificate means. Do not delegate
  these down; a plausible-looking wrong answer here is the failure mode the whole
  repository exists to avoid.
- **`deepseek-v4-flash`** — mechanical: a wrong field name, a duplicated
  declaration, a docstring that disagrees with its body. The correct answer is
  fully determined by the surrounding code.

Line numbers are as of `b9692e4` and will drift as fixes land; the anchors are the
function names.

### What has actually been executed

Findings were reached by reading; a subset has since been **reproduced by running
code**, and those are marked **`[verified]`** in their heading. What ran:

| Suite | Result |
|---|---|
| `ruff check lyapunov_toolbox tests` | clean |
| `pytest` (Python) | 18 passed, 10 skipped — every skip is a `@pytest.mark.julia` test, skipped because Julia is absent |
| `checkcode` on all five `matlab/*.m` | 0 messages, 0 parse errors (MATLAB R2023a) |
| Julia `Pkg.test()` | **not run — Julia is not installed on this machine** |
| MATLAB execution | **not run — YALMIP is not installed** |

The Julia suite is where every tier-B item lives, so the most important part of
this review remains unexecuted. Treat tier B as read-and-reasoned, not measured.

Note what the passing Python suite does *not* mean: R1-17, R1-18, R1-21 and R1-33
were each confirmed as real defects **while the suite was green**, because no test
covers them. A green suite here is evidence about the tests, not about the code.

### Verification gate

```bash
julia --project=julia -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
cd python && pytest && cd ..
pytest --nbmake examples/ --nbmake-timeout=1200
cd docs/derivations && latexmk -pdf sos_decomposition.tex && latexmk -pdf positivstellensatz.tex
```

Do not move to the next tier until this is green for the current one.

---

# Tier A — Blockers (R1-1 … R1-5)

Everything here either raises on the first call or is a deliverable that is
missing. None of it requires understanding the mathematics.

---

### R1-1 · `Base.show` for `LyapunovCertificate` reads five fields that do not exist **[verified]**

**File:** [julia/src/sos_lyapunov.jl:1043–1059](julia/src/sos_lyapunov.jl#L1043)
**Model:** `deepseek-v4-flash`
**Severity:** blocker

The struct is declared at [sos_lyapunov.jl:114–127](julia/src/sos_lyapunov.jl#L114)
with fields `V, ρ, x, multipliers, grams, ε, system_name, solver_status,
solve_time, sdp_size, bounded, boundedness_evidence`.

`show` reads `cert.V_degree` (:1045, :1056), `cert.is_bounded` (:1048),
`cert.iters` (:1049), and `cert.verified` (:1058). None exist. Displaying a
certificate in the REPL, in a notebook, or in any `@info` that interpolates one
raises `type LyapunovCertificate has no field V_degree`.

`length(cert.V)` at :1052 is also wrong — `V` is a polynomial, and `length` of a
polynomial is not its term count.

**Reproduced** (Julia 1.12.6), constructing a certificate by hand and calling
`sprint(show, MIME"text/plain"(), cert)`:

```
FieldError: type LyapunovCertificate has no field `V_degree`, available fields:
`V`, `ρ`, `x`, `multipliers`, `grams`, `ε`, `system_name`, `solver_status`,
`solve_time`, `sdp_size`, `bounded`, `boundedness_evidence`
```

**Fix.** Rewrite the body against the real fields:

- `V degree` → `maximum(degree, monomials(cert.V))`
- `bounded` → `cert.bounded`, and also print `cert.boundedness_evidence` (the
  docstring at :1035 asks for it, and `:sampling_only` vs `:top_form_pd` is
  exactly the distinction a reader needs — see R1-14)
- drop `iters` and `verified` — they belong to `ROAResult` and to the Python
  `RoaCertificate`, not here — or add them to the struct, but then update every
  construction site
- `sdp_size`, `solve_time`, `solver_status`, `system_name` are all real; print them
- term count → `length(terms(prune_terms(cert.V)))`

Add a test that constructs a certificate and calls
`sprint(show, MIME"text/plain"(), cert)`. There is currently no test that displays
one, which is why this survived.

---

### R1-2 · `benchmark_row` reads four fields that do not exist **[verified]**

**File:** [julia/src/roa_estimation.jl:617–648](julia/src/roa_estimation.jl#L617)
**Model:** `deepseek-v4-flash`
**Severity:** blocker

- :628 `cert.gram_dim` → the field is `sdp_size`
- :631 `cert.solver` → no such field; `LyapunovCertificate` never records the
  solver
- :629 `result.total_sdps` → the field is `n_sdps`
- :630 `result.total_time` → the field is `total_solve_time`

`benchmark_row` therefore raises on every call, which is why `docs/benchmarks.md`
is empty (R1-5) and why the README could not quote a reproducible number.

**Reproduced** (Julia 1.12.6), on a hand-built `ROAResult`:

```
FieldError: type LyapunovCertificate has no field `gram_dim`
```

It fails on the first of the four; the other three are behind it.

**Fix.** Correct the four names. For the solver, read it from `result.opts.solver`
rather than adding a struct field — it needs no change to `LyapunovCertificate`.
`fraction_true_roa = NaN` with the comment "must be set externally" is acceptable,
but make the caller's obligation explicit by returning `missing` instead of `NaN`,
so a `NaN` landing in a published table is not mistaken for a computed value.

---

### R1-3 · `positivstellensatz.tex` is a byte-identical copy of `sos_decomposition.tex`

**Files:** [docs/derivations/positivstellensatz.tex](docs/derivations/positivstellensatz.tex), [docs/derivations/sos_decomposition.tex](docs/derivations/sos_decomposition.tex)
**Model:** `deepseek-v4-pro`
**Severity:** blocker

`diff` reports the two files are identical. `positivstellensatz.tex` opens with the
comment `% sos_decomposition.tex — Parrilo's SOS decomposition` and carries
`\title{Sum-of-Squares Decomposition and the Gram Matrix Parameterisation}`. Both
compile to the same PDF.

The Positivstellensatz document does not exist. This is the most damaging single
item in the review: that document is the one cited in the thesis, it is the one a
theory-side advisor will actually open, and it is the stated reason the repository
exists. The README now says so explicitly rather than linking a PDF that would
have shown the wrong paper.

**Fix.** Write it. The section skeleton should be recovered from git:

```bash
git show b1167db:docs/derivations/positivstellensatz.tex > docs/derivations/positivstellensatz.tex
```

The two sections that carry the weight:

- **§5, the ROA theorem.** State it with both hypotheses and prove it in four
  steps. **Name the boundedness hypothesis explicitly at the step of the proof
  that uses it** — it is the hypothesis most often dropped, the proof is false
  without it, and `is_sublevel_set_bounded` exists precisely to check it. R1-14
  shows that function does not currently deliver what the proof needs; write the
  section and fix the function together.
- **§8, transition-viable funnels.** Write the composition condition formally —
  the reachable set under the reset map must land inside the next mode's certified
  set — and say precisely which piece is open for a *set-valued* reset. This is
  the paragraph that makes the portfolio and the research proposal one project.

Also state, in §3, that Putinar gives **no degree bound** on the multipliers, and
that this is why an infeasible result at a fixed degree proves nothing about the
system. That sentence is load-bearing for every "infeasible" the toolbox reports.

---

### R1-4 · All four notebooks are unwritten and unexecuted

**Files:** `examples/*/*.ipynb` (all four)
**Model:** `deepseek-v4-pro`
**Severity:** blocker

Measured state:

| notebook | code cells | with outputs | unfilled `**Write:**` |
|---|---|---|---|
| `van_der_pol_demo.ipynb` | 7 | 0 | 8 |
| `inverted_pendulum_demo.ipynb` | 7 | 0 | 7 |
| `unicycle_demo.ipynb` | 7 | 0 | 7 |
| `consensus_demo.ipynb` | 7 | 0 | 7 |

Every prose block is still the literal skeleton instruction `**Write:** ...`, and
no notebook has ever been executed. The code cells are implemented; only the prose
and the outputs are missing.

A reader on GitHub sees rendered markdown and stored outputs. Right now they see
29 instructions addressed to the implementer and no figures.

**Fix.** Replace every `**Write:**` block with 3–8 sentences of actual prose, first
person, plain and specific, with the numbers taken from the cell above it. Then
execute each notebook top to bottom in a fresh kernel and **commit with outputs**.

The blocks that carry the most weight, and must not be softened:

- `van_der_pol_demo.ipynb` §7 "Where this stops working" — report the degree at
  which solve time became unreasonable, on named hardware, and every
  `ALMOST_OPTIMAL` you hit and what you did about it.
- `inverted_pendulum_demo.ipynb` §4 — the honest cross-check. If the converged
  fraction is 1, report the largest \|θ\| actually tested. "No counterexample below
  \|θ\| = 0.9 rad" is a usable claim; "the Taylor certificate is valid" is not.
- `consensus_demo.ipynb` §6 — the switching counterexample. If the adversarial
  search finds nothing, say so with the number of sequences searched. Do not
  delete the section for producing an inconvenient answer.

The CI job `examples` already runs `pytest --nbmake examples/`, so a notebook that
does not execute will fail the build once R1-1 and R1-2 stop raising.

---

### R1-5 · `docs/benchmarks.md` is entirely empty

**File:** [docs/benchmarks.md](docs/benchmarks.md)
**Model:** `deepseek-v4-pro`
**Severity:** blocker

Every cell of the environment block (:19–29) and every data cell of all six tables
is blank, because `benchmark_row` raises on every call (R1-2).

Two structural facts to get right when filling it, both of which the README's
previous (now removed) result table got wrong:

- `unicycle_with_obstacle` builds a **4**-state system `[x, y, s, c]`
  ([systems.jl:305](julia/src/systems.jl#L305)), not 3.
- `consensus_network` raises `ArgumentError` for any κ ≤ 0
  ([systems.jl:384](julia/src/systems.jl#L384)), so no negative-κ row can be
  produced until R1-25 lands.

**Fix.** Fix R1-2, then generate every table from `benchmark_row` output and fill
the environment block — a benchmark without a named machine and solver version is
not a benchmark. Then populate the README "Results" table from it; that table is
currently a row of em-dashes with a pointer here, deliberately.

Fill the "first configuration that failed" row in §2 too (out of memory, over ten
minutes, or a non-`OPTIMAL` status). It is the most informative line in the
document and the one a reviewer looks for.

---

# Tier B — Soundness (R1-6 … R1-14)

Every item here allows a wrong certificate to be reported as correct, or removes a
check that was supposed to catch one. Fix them before anything cosmetic.

---

### R1-6 · `verify_certificate` check 1 is a tautology and verifies nothing **[verified]**

**File:** [julia/src/sos_lyapunov.jl:907–916](julia/src/sos_lyapunov.jl#L907) and [:937–945](julia/src/sos_lyapunov.jl#L937)
**Model:** `deepseek-v4-pro`
**Severity:** critical

```julia
p_recon = polynomial(basis' * Q * basis)
r = gram_residual(p_recon, basis, Q)
```

The polynomial is *reconstructed from* `(basis, Q)` and then checked against
`(basis, Q)`. `r` is identically zero up to floating-point noise, for every input,
including a `Q` that has nothing to do with the certificate. The only surviving
content is the `eigmin` test.

The whole point of check 1 was to confirm the Gram matrix reproduces its **target**
polynomial — for `:derivative` that is
`-(V̇ + ε·margin) + λ(V - ρ) - Σσⱼgⱼ - Σμₖhₖ`; for `:positivity` it is
`V - ε·margin`; for `:λ` and `:σ_j` it is the multiplier itself. As written, a
certificate whose Gram matrices belong to a completely different program passes.

**Reproduced** (Julia 1.12.6). A certificate was built whose stored Gram pairs are
pure fiction — a *random* PSD matrix over the basis `[1, x₁, x₂]`, unrelated to the
system, to `V`, or to `ρ`:

```
== verify_certificate on FABRICATED Gram matrices ==
  algebraic.ok = true
    positivity: gram_residual = 0.0, eigmin = 1.258481
    derivative: gram_residual = 0.0, eigmin = 1.258481

== rho inflated x1000, same witnesses ==
  algebraic.ok = true   (unchanged: rho is never consulted)
```

`gram_residual` is **exactly** `0.0`, as it must be: the polynomial is
reconstructed from the matrix and then compared with the matrix it came from.
Multiplying `ρ` by 1000 changes nothing, because `ρ` never enters the check at all.

The verifier accepts any PSD matrix of any size over any basis. It is not a
verifier.

`_witness_grams_ok` at [:596–614](julia/src/sos_lyapunov.jl#L596) already does this
correctly, via `_derivative_target`. The logic exists; `verify_certificate` just
does not use it.

**Fix.** Reconstruct each target from the certificate's own data and compare
against that:

1. Reconstruct the targets inside the two-argument
   `verify_certificate(cert, sys)` from `sys.f`, `cert.V`, `cert.ρ`,
   `cert.multipliers`, and `cert.ε`. Prefer reconstruction over storing them — a
   stored target could itself be stale. You will need the margin exponent `k`
   (R1-31) on the certificate to rebuild `:derivative` exactly.
2. For each `(key, (basis, Q))`, compute `gram_residual(target[key], basis, Q)` and
   require `< tol`.
3. Keep the `eigmin > -tol` test.

Then add the regression test that does not currently exist: build a valid
certificate, replace one stored `Q` with a random PSD matrix of the same size, and
assert verification **fails**. The present implementation passes that test, which
is how you will know the fix landed.

**Expect this fix to turn the suite red, and expect that to be correct.** A
prototype of it was written and run during review (then reverted — the Julia side
is not the reviewer's to change). With a real algebraic check in place, the full
suite went from `17148 pass / 0 fail` to `17145 pass / 1 fail / 1 error`:

- `verify_certificate / passes on a genuine certificate` **fails** — because the
  certificate it is handed carries the γ-step's Gram matrices for the V-step's
  `V` (that is R1-7), and a real check notices. Fix R1-7 and R1-6 together, in
  that order, or this test will keep failing for the right reason.
- `verify_certificate / calls no solver` **errors** — it calls the one-argument
  `verify_certificate(cert)`, which has to start raising once check 1 needs
  `sys`. Rewrite it against the two-argument form.

Both are tests that were only passing because the thing they tested did nothing.

The one-argument `verify_certificate(cert)` at [:901](julia/src/sos_lyapunov.jl#L901)
also returns `sampling = (ok = false, ...)` and `dynamical = (ok = false, ...)`
unconditionally. A caller doing `all([r.algebraic.ok, r.sampling.ok,
r.dynamical.ok])` gets `false` for a perfectly good certificate. Either make the
one-argument form raise `ArgumentError("verify_certificate needs the system; call
verify_certificate(cert, sys)")`, or delete it. Returning a struct whose fields say
`ok = false` for "not applicable" is a trap.

---

### R1-7 · The certificate from `v_s_iteration` carries Gram matrices for a different `V`

**File:** [julia/src/roa_estimation.jl:287–294](julia/src/roa_estimation.jl#L287)
**Model:** `deepseek-v4-pro`
**Severity:** critical

```julia
cert = LyapunovCertificate(
    V_current, ρ_current, x,
    w_γ.multipliers, w_γ.grams,     # ← witness from the γ-step
    ...)
```

`w_γ` is the witness returned by `multiplier_step(sys, V_scaled, 1.0)` at
[:240](julia/src/roa_estimation.jl#L240) — computed for `V_scaled`, the **previous**
iterate. `V_current` at :271 is `V_new`, the output of the V-step. The certificate
therefore pairs a new `V` with a witness proving something about the old one.

The V-step's own witness (`w_v.grams`, containing `:positivity`, `:derivative`,
`:shape_containment` for the correct `V`) is computed at
[sos_lyapunov.jl:728–732](julia/src/sos_lyapunov.jl#L728) and then discarded.

Because R1-6 makes check 1 a no-op, nothing catches this. Once R1-6 is fixed this
surfaces immediately as a residual failure — fix both together.

**Fix.** Build the certificate from `w_v`, not `w_γ`. The V-step also needs to
return its multipliers: `λ` was fixed and is known, but the `σ` and `μ` it creates
at [sos_lyapunov.jl:689–698](julia/src/sos_lyapunov.jl#L689) are created
anonymously inside the loop and never read back. Capture them into a `multipliers`
dict the way `_gather_witness` does, and return them in `w_v`.

---

### R1-8 · The Van der Pol containment test asserts nothing

**File:** [tests/test_roa_estimation.jl:72–78](tests/test_roa_estimation.jl#L72)
**Model:** `deepseek-v4-pro`
**Severity:** critical

```julia
min_dist = minimum(norm([bx - cx, by - cy]) for (cx, cy) in ...)
@test min_dist > 0 || error("boundary point ($bx, $by) is on or outside the limit cycle")
```

`min_dist` is the distance from a certified-boundary point to the *nearest point of
the limit cycle*. It is positive for any point not lying exactly on the cycle —
including every point **outside** it. The test passes for a certified set ten times
larger than the true ROA.

This is the reference soundness check of the whole toolbox, named as such in the
file's own header comment, and it is vacuous.

**Fix.** Test containment, not distance. The limit cycle is a closed curve, so use
a point-in-polygon test — `_point_in_polygon` already exists at
[roa_estimation.jl:590](julia/src/roa_estimation.jl#L590):

```julia
@test all(_point_in_polygon((bx, by), cycle_pts) for (bx, by) in boundary)
```

and additionally report the **minimum clearance** — the smallest distance from a
boundary point to the cycle — as a number, since a certificate that only just fits
is worth knowing about.

Two supporting problems in the same block:

- The "limit cycle" at :68–70 is `X_fwd[half:end, :]` — the second half of a 25 s
  integration, which for μ = 1 is roughly *four* periods, not one. A polygon that
  wraps four times is not a valid argument to any containment test. Extract exactly
  one period: detect a Poincaré crossing (e.g. `x₂ = 0` with `x₁ > 0`) and take the
  arc between two consecutive crossings. The MATLAB side had the identical bug and
  has already been fixed — copy the method from `limit_cycle_one_period` in
  [matlab/van_der_pol_roa.m](matlab/van_der_pol_roa.m), including its
  settled-orbit and star-shapedness checks, so the two sides of the
  cross-validation measure the same object.
- `@test cond || error(...)` converts a failing assertion into a thrown exception,
  aborting the rest of the suite instead of recording one failure. See R1-40.

---

### R1-9 · The manifold Lyapunov seed falls back to a ball about the origin for a non-origin equilibrium

**File:** [julia/src/sos_lyapunov.jl:443–445](julia/src/sos_lyapunov.jl#L443)
**Model:** `deepseek-v4-pro`
**Severity:** critical

```julia
if any(abs.(real(evt)) .< 1e-10)
    return sum(x[i]^2 for i in 1:n)
end
```

The docstring of this very function, at
[:402–406](julia/src/sos_lyapunov.jl#L402), describes exactly this as the bug it
was written to fix:

> "Without the manifold seed the old `P = I` fallback seeded a ball about the
> *origin*, which is not a Lyapunov function for a non-origin equilibrium on the
> manifold, and every invariant system produced a degenerate (unverifiable)
> certificate."

The fallback path still does it. And this is the path the unicycle takes. Working
the Jacobian at the goal `(0, 0, 0, 1)` by hand:

```
A = [-k_ρ   0   0   0
       0    0   0   0
       0  -k_α  0   0
       0    0   0   0]     eigenvalues {-k_ρ, 0, 0, 0}
```

Row 2 (`ẏ = v·s`, a product of two first-order-small quantities) is identically
zero, so the reduced dynamics on the tangent space have a zero eigenvalue, the
branch triggers, and `V = x² + y² + s² + c²` is returned. That `V` equals **1** at
the equilibrium, not 0. It is a ball centred on the origin of the ambient space,
which is not a point of interest for this system.

Everything downstream inherits it: the bisection certifies sublevel sets of a
function that is not a Lyapunov function for the goal, `roa_volume` measures the
wrong set, and the notebook plots it as "the certified ROA".

**Fix.** This needs a real decision, not a patch:

1. **At minimum,** never return a seed that fails `V(x_eq) = 0`. Centre the
   fallback: `sum((x[i] - eq[i])^2 for i in 1:n)`. That alone makes the object a
   candidate Lyapunov function again.
2. **The zero eigenvalues are the real problem.** A non-hyperbolic linearisation
   means no quadratic seed exists from the linearisation, full stop. Options, in
   order of preference:
   - Change the unicycle feedback so the closed loop *is* exponentially stable at
     the goal in the recast coordinates. The current
     `v = -k_ρ(x·c + y·s)`, `ω = -k_α(y·c - x·s)` cannot stabilise `y` at first
     order. A feedback with a nonzero `∂ẏ/∂y` at the goal fixes it and is the
     honest fix — the example is supposed to demonstrate a constrained
     certificate, not a centre manifold.
   - Or accept the non-hyperbolic case explicitly: raise `ArgumentError` naming
     the zero eigenvalues and require the caller to supply `V0`. Refusing is
     defensible; silently substituting a wrong seed is not.
3. Whichever you pick, add an assertion at the top of `estimate_roa` and
   `v_s_iteration`: `abs(V(x_eq)) < 1e-12`, raising otherwise. That single check
   would have caught this.

---

### R1-10 · The unicycle Hurwitz check is weakened to the point of being vacuous

**File:** [julia/src/systems.jl:323–327](julia/src/systems.jl#L323)
**Model:** `deepseek-v4-pro`
**Severity:** high

```julia
n_stable = count(real(ev) .< -1e-9)
n_stable >= 1 || throw(ArgumentError("unicycle closed-loop is not locally stabilising: ..."))
```

The function's own docstring at [:290–292](julia/src/systems.jl#L290) requires:

> "the linearisation at the goal must be Hurwitz — verify that before running any
> SOS program, because an unstable closed loop makes every downstream
> infeasibility uninformative."

"At least 1 of 4 eigenvalues has negative real part" is not that check. As computed
in R1-9, the actual spectrum is `{-k_ρ, 0, 0, 0}` — the check passes with three
zero eigenvalues. The comment at :320–322 rationalises one zero eigenvalue as the
invariant `s² + c² = 1` direction, which is legitimate, but there are three.

**Fix.** Check what the mathematics needs: project onto the tangent space of the
invariant manifold and require **all** eigenvalues of the reduced Jacobian to be
strictly in the left half-plane, exactly as `_lyapunov_on_manifold` does at
[sos_lyapunov.jl:428](julia/src/sos_lyapunov.jl#L428). One zero eigenvalue per
invariant is expected and should be excluded by construction, not by lowering the
bar. With the current feedback this check will fail — that failure is correct and
is the same finding as R1-9; fix the feedback, not the check.

---

### R1-11 · `bisect_level_set` returns the doubling cap instead of throwing

**File:** [julia/src/roa_estimation.jl:93–97](julia/src/roa_estimation.jl#L93)
**Model:** `deepseek-v4-pro`
**Severity:** high

```julia
if !found_infeasible && opts.verbose
    @warn("... the system appears globally stable — returning doubling cap " *
          "ρ_hi = $ρ_hi as a generous over-estimate.")
end
```

Three problems, in increasing order of seriousness:

1. The docstring at [:59–63](julia/src/roa_estimation.jl#L59) says "throw if no
   infeasible `ρ` is found … Do not silently return the cap." The body does not throw.
2. The warning only fires when `opts.verbose` is set, which defaults to `false`.
   The default behaviour is silent.
3. The message calls the result "a generous over-estimate". An over-estimate of a
   certified ROA is a **false certificate**. The returned level is in fact still
   verified feasible, so the object is sound — but a message telling a reader the
   toolbox knowingly reports over-estimates is worse than the bug, and someone will
   quote it.

**Fix.** Throw, as documented:

```julia
found_infeasible || throw(ErrorException(
    "bisect_level_set: no infeasible ρ found after 30 doublings (ρ = $ρ_hi). " *
    "For a globally stable system there is no finite boundary; the caller must " *
    "supply an explicit ρ_hi or handle the global case."))
```

`test_sos_lyapunov.jl`'s "globally stable linear system / no finite level is ever
infeasible" testset expects exactly this. Also initialise `best` to a verified value
rather than to `ρ_lo = 0.0` ([:100](julia/src/roa_estimation.jl#L100)): if no
bisection midpoint is feasible, the function returns 0.0 and then asserts
feasibility at 0.0, which is degenerate rather than informative.

Replace the `@assert` at [:122](julia/src/roa_estimation.jl#L122) with an explicit
`throw` — `@assert` can be elided under `--check-bounds=no`, and this is the last
line of defence between the bisection and a reported number.

---

### R1-12 · `beta_step` returns a `1e10` sentinel instead of failing

**File:** [julia/src/sos_lyapunov.jl:803–809](julia/src/sos_lyapunov.jl#L803)
**Model:** `deepseek-v4-pro`
**Severity:** high

```julia
opts.verbose && @warn("beta_step: no infeasible β found — V dominates shape; returning large β")
return 1e10
```

`1e10` then flows into `v_s_iteration` as `β_new`, into the monotonicity check at
[roa_estimation.jl:253](julia/src/roa_estimation.jl#L253) (which it trivially
passes), into `β_history`, into the convergence test at :298 (`1e10 - β_old <= tol *
max(1, β_old)` is false, so the loop keeps running), and into `v_step` as the shape
level to certify. The plotted convergence history then shows a jump to 1e10.

Same pattern as R1-11: warn-only, verbose-gated, magic number.

**Fix.** Throw with a message naming the actual condition — the shape set `p` is too
small relative to `V`, which usually means `V` was rescaled by a large `ρ₀` and the
default `p = ‖x - x_eq‖²` is no longer commensurate. Then handle it at the call
site: either rescale the shape set alongside `V`, or raise to the user. A sentinel
that survives three downstream consumers is not error handling.

---

### R1-13 · Six tests assert nothing

**Files:** [tests/test_roa_estimation.jl:85](tests/test_roa_estimation.jl#L85), [:91](tests/test_roa_estimation.jl#L91), [:258](tests/test_roa_estimation.jl#L258), [:302](tests/test_roa_estimation.jl#L302), [:406](tests/test_roa_estimation.jl#L406); [tests/test_sos_lyapunov.jl:239](tests/test_sos_lyapunov.jl#L239)
**Model:** `deepseek-v4-pro`
**Severity:** high

| line | testset | body |
|---|---|---|
| roa:85 | "certified area is a nontrivial fraction of the true ROA" | `@test_skip false` |
| roa:91 | "degree-4 V beats degree-2, degree-6 beats degree-4" | `@test_skip false` |
| roa:258 | "shaped shape-set recovers more volume than x'x" | `@test true` |
| roa:302 | "attempting to certify in x coordinates is infeasible" | `@test true` |
| roa:406 | MATLAB cross-validation, fixture absent | `@test true` |
| sos:239 | "throws on an unrecognised solver status" | `@test true` |

`@test true` is a passing test that checks nothing; it is strictly worse than a
deleted test, because the suite reports it as green.

The two that matter most:

- **roa:91** — "degree-4 beats degree-2, degree-6 beats degree-4" is the diagnostic
  that tells you the V-step is doing anything at all. Skipped for runtime. Put it
  behind `if get(ENV, "LYAPUNOV_SLOW_TESTS", "0") == "1"` and **run it in CI**,
  rather than skipping it everywhere.
- **roa:302** — "certifying in x coordinates is infeasible" is the trap the entire
  consensus example exists to demonstrate. Build the raw-x system (it is ten lines:
  the `f_x` construction already exists at
  [systems.jl:421–435](julia/src/systems.jl#L421), just do not project) and assert
  `verify_sublevel_set` reports infeasible.

**Fix.** Implement all six. For sos:239, inject a stub optimizer whose
`termination_status` returns `MOI.OTHER_ERROR` and assert the `ErrorException` —
JuMP's `MOI.Utilities.MockOptimizer` does this in a few lines. For roa:406, `@info`
plus `@test true` on a missing fixture is a silent pass; use `@test_skip` with the
reason so the summary line distinguishes it from a real pass.

---

### R1-14 · `is_sublevel_set_bounded` returns `:sampling_only` and every caller treats it as proof

**File:** [julia/src/polynomial_utils.jl:315–349](julia/src/polynomial_utils.jl#L315); callers at [roa_estimation.jl:141](julia/src/roa_estimation.jl#L141), [:281](julia/src/roa_estimation.jl#L281)
**Model:** `deepseek-v4-pro`
**Severity:** high

The function correctly distinguishes `:top_form_pd` (conclusive) from
`:sampling_only` (not proof) — its own docstring says "Never report
`:sampling_only` as proof." Both callers then do:

```julia
bounded, evidence = is_sublevel_set_bounded(...)
bounded || throw(...)
```

`bounded` is `true` in the `:sampling_only` case, so the evidence is recorded on the
certificate and otherwise ignored.

The sampling itself is also weak: it draws points on a **single** sphere of radius
50 ([:339–347](julia/src/polynomial_utils.jl#L339)). A sublevel set that is bounded
but extends past radius 50 is falsely reported unbounded (conservative, tolerable);
an unbounded set whose escape cone is narrow is falsely reported bounded (a false
certificate). In 4-D, 20 000 uniform directions cover a narrow cone poorly.

Boundedness is the hypothesis that turns "V̇ < 0 on Ω_ρ" into "Ω_ρ ⊆ ROA" — see
R1-3 §5. A certificate resting on a sampling heuristic should not be presented the
same way as one resting on a proof.

**Fix.**

1. Make the structural test work at every even degree. It already handles `deg == 2`
   exactly (eigenvalue test) and `deg > 2` via `is_sos(V_top - ε‖x‖^deg)`. Verify the
   SOS path is actually reached and passing for the degree-4 and degree-6 Van der
   Pol candidates — if it is not, that is the real gap.
2. Make `:sampling_only` a **refusal by default**. Add
   `require_proof::Bool = true` to `estimate_roa` / `v_s_iteration`; under the
   default, throw unless `evidence === :top_form_pd`. A caller who knowingly wants
   the weaker result passes `require_proof = false` and gets it recorded on the
   certificate.
3. Sample a **range** of radii (a geometric sweep from the estimated boundary out to
   `radius`), not one sphere.

---

# Tier C — Correctness and contract drift (R1-15 … R1-31)

Wrong numbers, silently dropped arguments, and docstrings the body does not honour.

---

### R1-15 · `taylor_recast` computes Taylor coefficients by finite differences and is unusable above order 3

**File:** [julia/src/polynomial_utils.jl:391–399](julia/src/polynomial_utils.jl#L391)
**Model:** `deepseek-v4-pro`
**Severity:** high

`_taylor_coefficient` builds the k-th derivative from a central-difference stencil
divided by `h^k` with `h = 1e-3`, plus one Richardson step. The stencil algebra is
right; the conditioning is fatal.

Relative error is roughly `eps / h^k`:

| k | `h^k` | error floor |
|---|---|---|
| 3 | 1e-9 | ~1e-7 |
| 5 | 1e-15 | ~10 |
| 7 | 1e-21 | ~1e5 |

At order 7 — a configuration the pendulum notebook and `pendulum_roa.m` both sweep
over — the field would be built from numerical noise. And `jacobian_at`'s docstring
([:99–105](julia/src/polynomial_utils.jl#L99)) states the package's own position:
"Obtained by symbolic differentiation … not by finite differences — `f` is
polynomial, so there is no reason to accept truncation error here."

Note that `inverted_pendulum` does **not** call `taylor_recast` — it builds the sine
series in closed form at [systems.jl:225](julia/src/systems.jl#L225), which is
correct. So the damage is confined to direct callers. That makes the function a trap
rather than an active bug, which is arguably worse: it looks available and is
silently wrong.

**Fix.** Replace the finite-difference machinery with exact coefficients. Two
workable routes:

- **Preferred:** restrict `expr_fn` to a symbolic form and differentiate exactly.
  Since every use in this repository is a truncated `sin`/`cos`, accept a
  `Val{:sin}` / `Val{:cos}` dispatch or a small table of known series and raise
  `ArgumentError` for anything else. Honest and exact.
- **Or:** automatic differentiation via `TaylorSeries.jl` — but that is a new
  dependency, which the repository's dependency policy discourages.

Either way, delete `_central_difference_coefficients` and `_taylor_coefficient` once
unused, and keep the `truncation_error` metadata (a good idea) computed from the
exact next coefficient.

---

### R1-16 · Boundary and box routines are centred on the origin and ignore `xlim` / `ylim` / `fixed`

**Files:** [julia/src/roa_estimation.jl:418–457](julia/src/roa_estimation.jl#L418) (`roa_boundary_2d`), [:478–515](julia/src/roa_estimation.jl#L478) (`true_roa_boundary_2d`), [:376–400](julia/src/roa_estimation.jl#L376) (`_estimate_box`)
**Model:** `deepseek-v4-pro`
**Severity:** high

`roa_boundary_2d` accepts `xlim`, `ylim`, and `fixed` and uses none of them. It
evaluates `V(x => [r*d[1], r*d[2]])` — a hard-coded 2-vector — so:

- it raises or misbehaves for any system with `n > 2`, i.e. the trig pendulum
  (n = 3), the unicycle (n = 4), and any consensus network with n > 3 agents;
- it marches rays from the **origin**, but the trig pendulum's equilibrium is
  `(0, 1, 0)` and the unicycle's is `(0, 0, 0, 1)`. Rays from the origin do not
  describe those sublevel sets at all;
- the documented `fixed` slice support does not exist, so the 2-D slice figures the
  notebooks promise cannot be produced correctly.

`true_roa_boundary_2d` has the same shape problem (`x0 = [x_eq + …, y_eq + …]`,
2 entries only) and also ignores `xlim`/`ylim`/`fixed`.

`_estimate_box` marches from the origin too ([:387](julia/src/roa_estimation.jl#L387)),
so `roa_volume` measures the wrong region for every non-origin-equilibrium system —
which is exactly the set of systems whose volumes go in the benchmark table.

**Fix.**

1. Centre every ray march on `x_eq`, not `0`. The bisection point becomes
   `eq .+ r .* d_full` where `d_full` is the full-dimension direction.
2. Honour `fixed`: build the full-length point with the pinned coordinates and vary
   only the two axis coordinates. Add an `axes::Tuple{Int,Int} = (1, 2)` keyword to
   match the Python `boundary_2d` signature.
3. Use `xlim`/`ylim` for the bisection bracket instead of the hard-coded `100.0`, or
   remove them from the signature. A parameter that is accepted and ignored is worse
   than one that is absent.
4. Pass an explicit `equilibrium` to `_estimate_box`.

Add a test that plots a slice of the 4-state unicycle certificate and asserts the
returned points satisfy `V = ρ` to `1e-6` — the existing "every returned point
satisfies V = rho" test only exercises 2-D systems.

---

### R1-17 · `certify_roa` interpolates unvalidated `variables` directly into Julia source **[verified]**

**File:** [python/lyapunov_toolbox/julia_bridge.py:423–441](python/lyapunov_toolbox/julia_bridge.py#L423)
**Model:** `deepseek-v4-pro`
**Severity:** high

```python
var_str = " ".join(variables)
code = f"""
@polyvar {var_str}
x = [{", ".join(variables)}]
...
"""
```

`_check_expression_safety` is applied to `dynamics`, `constraints`, and `invariants`
via `_vector_literal` — but **never to `variables`**. A caller passing

```python
variables=["x1", "x2); run(`whoami`); @polyvar (y"]
```

gets that text spliced into the `@polyvar` line and executed. `certify_roa` is
called from notebooks and is the entry point other repositories in the portfolio are
meant to use, so this is a genuine path from user data to code execution.

**Fix.** Validate every variable name against a strict identifier pattern before it
touches the f-string:

```python
_IDENT = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
for v in variables:
    if not _IDENT.match(v):
        raise ValueError(f"variable name {v!r} is not a plain identifier")
if len(set(variables)) != len(variables):
    raise ValueError(f"duplicate variable names: {variables}")
```

Add the corresponding cases to `TestExpressionSafety` in
[python/tests/test_julia_bridge.py](python/tests/test_julia_bridge.py) — the existing
parametrised list covers `dynamics` only.

---

### R1-18 · `_check_expression_safety` is a character-class filter, not a token allowlist **[verified]**

**File:** [python/lyapunov_toolbox/julia_bridge.py:602–623](python/lyapunov_toolbox/julia_bridge.py#L602)
**Model:** `deepseek-v4-pro`
**Severity:** medium

```python
allowed_pattern = re.compile(
    r'^[\d\s\+\-\*\/\^\(\)\.' + ''.join(re.escape(v) for v in variables) + r']+$'
)
```

The variable names are concatenated **into a character class**, so what is permitted
is the set of *letters* appearing in any variable name, not the names themselves.
With `variables = ("x1", "x2")` the class is `[0-9\s+\-*/^().x12]`, and any string
built from those characters passes regardless of whether it is a well-formed
expression. It happens to be safe for the current variable sets because no dangerous
Julia identifier can be spelled from `{x,1,2}` — but that safety is accidental and
depends on the caller's naming.

The `forbidden` set at :613 is a denylist. It only runs after the allowlist has
already rejected the string, so it is used purely to improve the error message; that
is fine, but the comment should say so, because a denylist that looks like a control
invites someone to rely on it.

**Reproduced.** With `variables = ("x1", "x2")`, every one of these is *accepted*
by `_check_expression_safety`, and none is a well-formed expression:

```
'x1x2x1'   '1x1..2'   '(((('   'x1^^^x2'   '2..5*x1'
```

**Fix.** Tokenise and check tokens:

```python
_TOKEN = re.compile(r"[A-Za-z_][A-Za-z0-9_]*|\d+\.?\d*|[+\-*/^()]|\s+")
```

Scan the expression, require full coverage, and require every identifier token to be
in `variables`. Reject on the first token that is not. Keep the message naming the
offending token — `test_rejects_undeclared_variable` depends on it.

---

### R1-19 · The trig pendulum silently changes the controller, then claims to certify the real pendulum

**File:** [julia/src/systems.jl:247–249](julia/src/systems.jl#L247)
**Model:** `deepseek-v4-pro`
**Severity:** high

```julia
# In the recast: sin θ → s, and the θ in u is replaced by s (since s ≈ θ near 0).
θdot_expr = (g/l)*s - (b/(m*l^2))*θdot + (1.0/(m*l^2))*(-(K[1]*s + K[2]*θdot))
```

The plant is recast exactly, but the **feedback** is changed from `u = -K₁θ - K₂θ̇`
to `u = -K₁ sin θ - K₂θ̇`. That is a different closed-loop system, not a different
representation of the same one.

The `notes` string at [:254–260](julia/src/systems.jl#L254) nevertheless says "The
certificate is about the true pendulum … unlike the Taylor recast", and the whole
point of the pendulum example is a like-for-like comparison between the two
recastings. Comparing a certificate for `u = -K₁ sin θ` against one for `u = -K₁θ`
is not that comparison, and the notebook and tests both present it as if it were.

`u = -K₁ sin θ` is a perfectly reasonable controller (a real implementation
measuring `sin θ` would produce it), so the fix is not necessarily to change the
dynamics.

**Fix.** Decide and say so:

- **If you keep it:** rewrite `notes` to state that the *controller* uses `sin θ`,
  not `θ`, and that the trig and Taylor systems therefore differ in both plant and
  feedback. Update the notebook §2 table and the test in
  [tests/test_roa_estimation.jl](tests/test_roa_estimation.jl) that compares them —
  it currently claims a comparison "on the SAME physical system".
- **If you want the honest comparison:** carry `θ` as a fourth state with `θ̇` as its
  derivative, so `u = -K₁θ - K₂θ̇` stays polynomial. This costs one state, keeps the
  invariant, and is the version that supports the claim the example makes.

---

### R1-20 · `synthesize` drops `epsilon` and `shape`, then records the `epsilon` it did not use

**File:** [python/lyapunov_toolbox/julia_bridge.py:468–520](python/lyapunov_toolbox/julia_bridge.py#L468)
**Model:** `deepseek-v4-flash`
**Severity:** medium

The generated Julia at :492–494 sets only `V_degree`, `multiplier_degree`, and
`vs_iterations`. `epsilon` is never passed — Julia uses its default `1e-6`. But
:512 records

```python
epsilon=float(kwargs.get("epsilon", 1e-6)),
```

so a caller who passes `epsilon=1e-4` gets a certificate whose `epsilon` field says
`1e-4` while the SDP was solved with `1e-6`. That is a false provenance record, on
the one object whose purpose is provenance.

`shape` is accepted at :474 and never used. Everything else in `**kwargs` is
swallowed silently, which is how a caller's typo becomes a wrong result rather than
an error.

**Fix.** Pass `epsilon` into `SOSOptions(ε = ...)`. Either implement `shape`
(interpolate a validated Julia expression, via the same safety check as `dynamics`)
or remove the parameter. Replace `**kwargs` with explicit keywords and raise
`TypeError` on unknown ones.

---

### R1-21 · `get_session` silently ignores a conflicting backend or project **[verified]**

**File:** [python/lyapunov_toolbox/julia_bridge.py:379–387](python/lyapunov_toolbox/julia_bridge.py#L379)
**Model:** `deepseek-v4-flash`
**Severity:** medium

```python
if JuliaSession._instance is not None:
    return JuliaSession._instance
```

The arguments are discarded. `JuliaSession.__init__` does raise on a second
construction ([:267–271](python/lyapunov_toolbox/julia_bridge.py#L267)), but
`get_session` never reaches it.

`test_conflicting_session_config_raises` in
[python/tests/test_julia_bridge.py](python/tests/test_julia_bridge.py) asserts the
opposite of the current behaviour.

**Fix.** Compare the requested `backend` and `project` against the live session's and
raise `RuntimeError` naming both values when they differ. Return the existing session
only when the request is compatible (including the common case where both arguments
are `None`).

---

### R1-22 · The `vdp_certificate` fixture converts a verification failure into a skip

**File:** [python/tests/conftest.py:60–61](python/tests/conftest.py#L60)
**Model:** `deepseek-v4-flash`
**Severity:** medium

```python
if not cert.verified:
    pytest.skip("Van der Pol certificate failed verification")
```

A certificate that fails verification is the most important failure this suite can
detect. Skipping on it makes the suite green, and because the fixture is
session-scoped every downstream test skips too — a silent, cascading pass.

**Fix.** `pytest.fail(...)` with the verification detail in the message. Reserve
`skip` for "the backend is not available", which `pytest_collection_modifyitems`
already handles.

---

### R1-23 · `_marshal_polynomial` raises when `V` does not contain every variable

**File:** [python/lyapunov_toolbox/julia_bridge.py:559–569](python/lyapunov_toolbox/julia_bridge.py#L559)
**Model:** `deepseek-v4-flash`
**Severity:** medium

```julia
idx = [findfirst(v -> string(v) == vn, vs) for vn in varnames]
any(isnothing, idx) && error("marshalling: a variable was not found in the polynomial")
```

`vs = variables(p)` lists only the variables that actually appear. A legitimate
sparse candidate — `V = x₁²` in a two-state system, or any `V` whose V-step drove a
coefficient to zero — makes `findfirst` return `nothing` and the marshalling fails.

The variable-order handling here is otherwise correct (it maps the caller's order
onto the polynomial's, which is the right thing) and should be kept.

**Fix.** Treat a missing variable as an all-zero exponent column rather than an
error:

```julia
E = [i === nothing ? 0 : degree(m, vs[i]) for m in ms, i in idx]
```

Add a test with a deliberately sparse `V` and assert Python and Julia evaluations
agree.

---

### R1-24 · `containment_fraction` uses the global RNG

**File:** [julia/src/roa_estimation.jl:554–584](julia/src/roa_estimation.jl#L554)
**Model:** `deepseek-v4-flash`
**Severity:** medium

Both sampling loops call bare `rand()`, and the function takes no `rng`. Every other
sampler in the package threads an explicit `Random.AbstractRNG`. This one produces a
different answer on every call, and it is the function that computes the "fraction of
true ROA" number destined for `docs/benchmarks.md`.

**Fix.** Add `rng::Random.AbstractRNG = Random.default_rng()` and thread it through
both loops and the `_estimate_box` call at :572. Same treatment for any other bare
`rand()` in the file.

---

### R1-25 · `consensus_network` rejects κ ≤ 0, so the documented negative-κ case cannot run

**File:** [julia/src/systems.jl:384](julia/src/systems.jl#L384)
**Model:** `deepseek-v4-flash`
**Severity:** medium

`κ > 0 || throw(ArgumentError(...))` blocks the case the consensus notebook §4 and
`consensus_demo.jl` step 4 both ask for: *"include one negative κ to show the failure
is reported as a raised ArgumentError rather than a tiny ρ."*

The intended behaviour is that a destabilising κ makes the *linearisation*
non-Hurwitz and `lyapunov_from_linearization` raises — a meaningful failure at the
right place. The constructor pre-empts it with a different, less informative one.

**Fix.** Allow any real κ. Keep a guard against a κ that makes the polynomial
degenerate if there is one, but let the Hurwitz check in
`lyapunov_from_linearization` be what reports a destabilising coupling.

---

### R1-26 · `roa_volume` silently falls back to Monte Carlo when `:grid` is requested with n > 2

**File:** [julia/src/roa_estimation.jl:333](julia/src/roa_estimation.jl#L333)
**Model:** `deepseek-v4-flash`
**Severity:** medium

```julia
if method === :grid && n <= 2
    ...
else
    # Monte Carlo
```

A caller asking for `:grid` on a 3-D or 4-D system gets Monte Carlo without being
told. The docstring at [:323–324](julia/src/roa_estimation.jl#L323) says `:grid` is
"reproducible to the last digit, so it is the one the tests use" — a test that asks
for `:grid` and silently receives a stochastic estimate is not reproducible, and the
`se = 0.0` returned by the grid branch means the caller cannot tell which they got.

The `if n > 2 error(...)` at :347–350 is dead code inside a branch already guarded by
`n <= 2`.

**Fix.** Raise `ArgumentError` for `:grid` with `n > 2`, naming the dimension and
suggesting `:montecarlo`. Delete the dead branch. Consider returning the method
actually used as a third element so a benchmark row cannot misattribute the number.

---

### R1-27 · The Julia fixture reader parses only `ρ` and ignores the rest of the schema

**File:** [tests/test_roa_estimation.jl:410–411](tests/test_roa_estimation.jl#L410)
**Model:** `deepseek-v4-flash`
**Severity:** medium

```julia
data = readdlm(fixture_path, ',')
ρ_matlab = data[1, 1]
```

The reader takes one cell and stops. It has no header handling, so with the fixture
writer now emitting a header row it would read the literal string `"rho"` as a
number.

The MATLAB side of this has been fixed: [matlab/write_roa_fixture.m](matlab/write_roa_fixture.m)
now emits the schema documented in
[tests/fixtures/README.md](tests/fixtures/README.md) — long format, header row, one
row per monomial of `V`, scalar fields repeated per row:

```
rho,monomial,coeff,degree,epsilon,solver,solver_status,gram_residual,min_eig_Q,timestamp
2.9700...e+00,"[2,0]",1.5...e+00,2,1e-06,mosek,feasible,3.2e-10,-1.1e-11,2026-08-04T13:59:00
```

`monomial` is a quoted exponent tuple in the variable order of `x`.

**Fix.** Rewrite the reader against that layout: skip the header, parse the exponent
tuples, reconstruct `V` as a polynomial, and compare **both** `ρ` and the
coefficients against the Julia result. Comparing the coefficients is the point —
two implementations that agree on `ρ` while disagreeing on `V` have not
cross-validated anything.

Also fix the surrounding testset (R1-13, roa:406): `@info` plus `@test true` on a
missing fixture is a silent pass.

---

### R1-28 · The V–s convergence test at iteration 1 compares two incommensurable β values

**File:** [julia/src/roa_estimation.jl:226](julia/src/roa_estimation.jl#L226), [:298](julia/src/roa_estimation.jl#L298)
**Model:** `deepseek-v4-pro`
**Severity:** medium

`β0` is computed from the **unnormalised** seed `V` at :226. Inside the loop, `β_new`
is computed from `V_scaled = V_current / ρ_current` at :249. The monotonicity check
at :253 correctly exempts iteration 1 for this reason — but the convergence check at
:298 does not:

```julia
if β_new - β_old <= opts.vs_improvement_tol * max(1.0, β_old)
```

At `iter == 1`, `β_old` is still `β0`. Since `β` scales with `V`, the two numbers
differ by roughly `ρ0`, so the comparison is arbitrary: for `ρ0 > 1` it can declare
convergence on the first iteration and exit with `converged = true` after doing no
work at all.

**Fix.** Compute `β0` on the normalised seed — `beta_step(sys, V / ρ0, shape)` — so
every recorded `β` lives on the same scale. The monotonicity exemption at :253 can
then be removed too, which is the stronger check. Record both the raw and normalised
β in the history if the unnormalised value is wanted for reporting.

---

### R1-29 · The V-step's `ALMOST_OPTIMAL` gate omits the residual check

**File:** [julia/src/sos_lyapunov.jl:742–751](julia/src/sos_lyapunov.jl#L742)
**Model:** `deepseek-v4-pro`
**Severity:** medium

`_v_step_grams_ok` checks only `eigmin(Q) > -1e-8` on the three constraints.
`_gram_certificate_ok` ([:209](julia/src/sos_lyapunov.jl#L209)) and
`_witness_grams_ok` ([:596](julia/src/sos_lyapunov.jl#L596)) both additionally
require `gram_residual < 1e-7` against the target polynomial, which is the half that
detects "the solver returned a PSD matrix for a different problem".

A PSD Gram matrix that does not reproduce its target is not a certificate. This is
the accept path for a marginal solve on the step that produces the `V` everything
downstream is built on.

**Fix.** Reconstruct the three targets — `V - ε·margin`, the derivative expression
with the fixed `λ` substituted, and `(1 - V) - σ_β(β - p)` — and apply
`gram_residual` to each, as `_witness_grams_ok` does. Factor the shared logic into
one helper so the three copies cannot drift.

---

### R1-30 · `v_s_iteration` breaks with a verbose-only warning where the docstring says throw

**File:** [julia/src/roa_estimation.jl:243–246](julia/src/roa_estimation.jl#L243), [:265–268](julia/src/roa_estimation.jl#L265)
**Model:** `deepseek-v4-flash`
**Severity:** medium

The docstring at [:189–190](julia/src/roa_estimation.jl#L189) states:

> "**Feasibility is preserved.** The previous `V` is always feasible for the V-step,
> so an infeasible V-step is a bug, not a hard instance. Throw."

Both handlers `break` after an `opts.verbose &&` warning. With the default
`verbose = false` the iteration silently returns the previous certificate and reports
`converged = false`, indistinguishable from budget exhaustion.

**Fix.** Throw `ErrorException` on an infeasible V-step, naming the iteration and the
solver status. For the γ-step, throwing is also right for the same reason (the
previous `(V, ρ)` was feasible by construction). If you decide a γ-step failure is
recoverable, record it in `ROAResult` as a distinct terminal state — not as
`converged = false`, which already means something else.

---

### R1-31 · The positivity margin is `ε‖x‖^{2k}`, and the comment claims it is `ε‖x‖²`

**File:** [julia/src/sos_lyapunov.jl:505–512](julia/src/sos_lyapunov.jl#L505), same block in [`v_step` :679–682](julia/src/sos_lyapunov.jl#L679)
**Model:** `deepseek-v4-flash`
**Severity:** medium

```julia
k = max(1, cld(dVd, 2))
margin = sum((x .- sys.equilibrium) .^ 2)^k
# "For equilibrium-at-origin systems this is exactly the guide's ε‖x‖²."
```

For the reversed Van der Pol, `deg V̇ = 4`, so `k = 2` and the margin is `ε‖x‖⁴`, not
`ε‖x‖²`. The comment is wrong for every system in the repository except a linear one.

The certificate remains **sound** — `V̇ ≤ -ε‖x - x_eq‖^{2k} < 0` away from the
equilibrium is still strict negativity — so this is a documentation and
reproducibility defect, not a correctness one. But the margin is much weaker near the
equilibrium than the stated one, and `ε` is a recorded benchmark parameter, so a
reader comparing against the literature will misinterpret it.

Centring the margin on `sys.equilibrium` is correct and good — keep that.

**Fix.** Correct the comment to state the actual margin, explain why the degree is
matched to `deg V̇` (so the SOS degrees line up without inflating the basis), and note
that `ε‖x‖^{2k}` is a weaker but still strict margin. Record `k` on the certificate
so `verify_certificate` can reconstruct the target (needed by R1-6), and report it in
the benchmark row alongside `ε`.

---

# Tier D — Polish (R1-32 … R1-44)

Real defects, low blast radius. Batch them.

---

### R1-32 · `van_der_pol(reversed = false)` is labelled "(reversed)" **[verified]**

**File:** [julia/src/systems.jl:117](julia/src/systems.jl#L117) · **Model:** `deepseek-v4-flash`

The name is hard-coded regardless of the flag, so the forward oscillator — used as
the negative test and as the limit-cycle ground truth — appears in every log line,
plot title, and benchmark row as "Van der Pol (reversed)". The `notes` field is
correctly branched two lines above; branch the name too.

---

### R1-33 · `plot_roa_2d` interior level curves are only correct for a quadratic `V` **[verified]**

**File:** [python/lyapunov_toolbox/plotting.py:87–92](python/lyapunov_toolbox/plotting.py#L87) · **Model:** `deepseek-v4-flash`

`ax.plot(bx * np.sqrt(frac), by * np.sqrt(frac))` scales the boundary radially, which
reproduces `{V = ρ·frac}` only when `V` is a homogeneous quadratic. For the degree-4
and degree-6 candidates the drawn contours are not level sets of `V` at all.

**Reproduced.** For `V = x₁⁴ + x₂⁴`, `ρ = 1`, the contour at `ρ/4` crosses the
`x₁` axis at `(1/4)^{1/4} = 0.7071`. The `sqrt(frac)` shortcut draws it at
`0.5` — **29.3% too small**. The error grows with degree, and every curve in the
figure is wrong in the same direction, which makes it look like a plausible
nesting rather than a bug.

Re-bisect along each ray at the reduced level instead — call
`certificate.boundary_2d` with a temporary `rho * frac`, or add a `level` argument to
`boundary_2d`.

---

### R1-34 · A slice of an n > 2 certificate with `fixed=None` is titled as if it were the whole set

**File:** [python/lyapunov_toolbox/plotting.py:114–118](python/lyapunov_toolbox/plotting.py#L114) · **Model:** `deepseek-v4-flash`

The title says "Region of Attraction" whenever `fixed` is falsy. But `boundary_2d`
defaults every unplotted coordinate to **0** — not to the equilibrium — so for the
unicycle (`x_eq = (0,0,0,1)`) the figure is an unlabelled slice through `c = 0`, a
point not on the invariant manifold at all.

Rule 3 in this module's own docstring exists for this case. Raise `ValueError` when
`certificate.n_vars > 2` and `fixed` does not pin every non-axis coordinate.
`test_slice_title_states_pinned_values` covers the labelling but not this.

---

### R1-35 · `plot_phase_portrait` infers the state dimension from `len(fixed)`

**File:** [python/lyapunov_toolbox/plotting.py:144](python/lyapunov_toolbox/plotting.py#L144) · **Model:** `deepseek-v4-flash`

`n_vars = len(fixed) + 2` is right only when the caller pins every non-axis
coordinate. Pin one of three extra coordinates and `dynamics` receives a wrong-length
vector. Take `n_vars` as an explicit argument, or derive it from the certificate the
caller already has.

---

### R1-36 · Python `boundary_2d` has no multiple-crossing detection, so it silently disagrees with Julia

**File:** [python/lyapunov_toolbox/julia_bridge.py:178–204](python/lyapunov_toolbox/julia_bridge.py#L178) · **Model:** `deepseek-v4-flash`

`roa_boundary_2d` in Julia raises when a ray crosses `{V = ρ}` more than once
([roa_estimation.jl:438](julia/src/roa_estimation.jl#L438)); the Python
implementation returns the innermost crossing without comment.
`test_boundary_2d_matches_julia` asserts they agree to `1e-6`, which will not hold on
a non-star-shaped set — one raises, the other returns a wrong curve.

Port the coarse-scan crossing count and raise `ValueError` to match.

---

### R1-37 · `simulate` reports a time grid that does not match the step it integrated with

**File:** [julia/src/systems.jl:501–502](julia/src/systems.jl#L501) · **Model:** `deepseek-v4-flash`

`nsteps = ceil((t1 - t0)/dt)` but `t = range(t0, t1; length = nsteps + 1)`, so the
reported spacing is `(t1-t0)/nsteps` while RK4 advances by `dt`. Whenever `(t1 - t0)`
is not an exact multiple of `dt`, the returned `t` mislabels every sample and the
trajectory overshoots `t1`. Either integrate with the adjusted step
`(t1-t0)/nsteps`, or build `t` as `t0 .+ dt .* (0:nsteps)` and clamp the last step.

Separately, the returned status is `:ok`, while the skeleton documented
`:converged` / `:diverged` / `:timeout`. Callers ignore it today
([sos_lyapunov.jl:1009](julia/src/sos_lyapunov.jl#L1009) discards it) — the dynamical
check could report `:diverged` directly instead of re-deriving it.

---

### R1-38 · Dead parameters and dead code

**Model:** `deepseek-v4-flash`

- [sos_lyapunov.jl:573–575](julia/src/sos_lyapunov.jl#L573)
  `_poly_variable_polynomial(μ, basis)` ignores `basis` entirely. Drop the parameter
  and the call sites that pass it.
- [roa_estimation.jl:347–350](julia/src/roa_estimation.jl#L347) unreachable
  `error("grid method only supported for n ≤ 2")` inside a branch guarded by
  `n <= 2` (see R1-26).
- [julia_bridge.py:410](python/lyapunov_toolbox/julia_bridge.py#L410)
  `certify_roa(solver=...)` is accepted and never used — either thread it into
  `SOSOptions` or remove it.
- [polynomial_utils.jl:379–399](julia/src/polynomial_utils.jl#L379)
  `_central_difference_coefficients` / `_taylor_coefficient` become dead once R1-15
  lands.

---

### R1-39 · Julia ray sweeps duplicate the endpoint

**Model:** `deepseek-v4-flash`

`range(0, 2π; length = n)` at [roa_estimation.jl:420](julia/src/roa_estimation.jl#L420)
and :484 emits both `0` and `2π` — a duplicated vertex. Harmless for `polyarea`, but it makes the point count off by one
and adds a zero-length polygon edge. Python's `boundary_2d` already gets this right
with `endpoint=False`, and the MATLAB scripts have already been switched to
`(0:n-1)*(2*pi/n)`; match them.

---

### R1-40 · `@test cond || error(...)` turns a failed assertion into an aborted suite

**Files:** [tests/test_roa_estimation.jl:76](tests/test_roa_estimation.jl#L76), [:297](tests/test_roa_estimation.jl#L297), [:415](tests/test_roa_estimation.jl#L415) and elsewhere · **Model:** `deepseek-v4-flash`

The macro sees `@test (cond || error(msg))`. When `cond` is false, `error` throws
before `@test` can record a failure, so the run aborts at the first bad assertion
instead of reporting all of them. Use `@test cond` and put the diagnostic in the
message form Test.jl supports, or compute the message into a variable first.

`grep -n '@test .*|| error(' tests/` finds them all.

---

### R1-41 · `inverted_pendulum` docstring states LQR weights the code does not use

**File:** [julia/src/systems.jl:141–143](julia/src/systems.jl#L141) vs [:204–205](julia/src/systems.jl#L204) · **Model:** `deepseek-v4-flash`

Docstring: "compute the LQR gain for the linearisation with `Q = I`, `R = 1`". Code:
`Q = Diagonal([10.0, 1.0])`, with an inline comment acknowledging the mismatch. The
code is right; correct the docstring rather than leaving a comment that points at it.

---

### R1-42 · `--require-julia` is read but never registered

**File:** [python/tests/conftest.py:20](python/tests/conftest.py#L20) · **Model:** `deepseek-v4-flash`

`config.getoption("--require-julia", default=False)` works only because of the
`default` fallback; there is no `pytest_addoption` hook, so `pytest --require-julia`
fails with "unrecognized arguments". Either add the hook or drop the option and rely
on `LYAPUNOV_TOOLBOX_REQUIRE_JULIA`, which CI already sets.

---

### R1-43 · The repository owner slug is spelled three different ways

**Files:** [README.md](README.md) (`ali-eiman`), [julia/Project.toml:4](julia/Project.toml#L4) (`Ali-Eiman`), [python/pyproject.toml](python/pyproject.toml) `[project.urls]` (`Ali-Eimaan`) · **Model:** `deepseek-v4-flash`

The README was left on `ali-eiman` during this review because that is what the author
wrote and the correct GitHub handle could not be verified from the repository. Two of
the three must be wrong.

**Fix.** Confirm the actual GitHub handle, then make all three agree — including the
citation block in the README, which is the one that ends up in a paper. A citation
hook pointing at a 404 is worse than none.

---

### R1-44 · Verify the rewritten README quickstart actually runs

**File:** [README.md](README.md) "Quickstart" · **Model:** `deepseek-v4-flash`

The quickstart was rewritten during this review against the verified signatures of
`certify_roa` and `RoaCertificate` (checked by introspection: parameter names, return
type, attribute names). It was **not executed**, because no Julia installation was
available in the review environment.

**Fix.** Run it in a clean checkout after the tier A and tier B fixes land, and paste
the real output into the block as a comment. The standing rule for that file is that
every line was executed before it was committed.

---

# Summary

| Tier | Items | `pro` | `flash` |
|---|---|---|---|
| A — Blockers | R1-1 … R1-5 | 3 (R1-3, R1-4, R1-5) | 2 |
| B — Soundness | R1-6 … R1-14 | 9 (all) | 0 |
| C — Correctness & contract | R1-15 … R1-31 | 7 (R1-15…19, R1-28, R1-29) | 10 |
| D — Polish | R1-32 … R1-44 | 0 | 13 |
| **Total** | **44** | **19** | **25** |

`README.md` and `matlab/` are excluded — both were fixed in this pass rather than
written up. The MATLAB changes are unrun; see the note at the top.

## Suggested sequencing

1. **R1-1 + R1-2** (`flash`) first, in one pass — two wrong-name fixes that unblock
   every downstream gate.
2. **R1-6 + R1-7 together** (`pro`). They are one bug seen from two sides: the
   verifier does not check the target, so nobody noticed the certificate carries the
   wrong witness. Fixing either alone will look like it made things worse.
3. **R1-9 + R1-10 together** (`pro`) — both are the unicycle's non-hyperbolic
   linearisation. Decide the feedback question once.
4. **R1-8** (`pro`), porting the limit-cycle extraction already written for the
   MATLAB side so the two halves of the cross-validation measure the same object.
5. **R1-2 → R1-5 → README "Results"** last, once the numbers exist.
6. Everything else in numeric order.

## The one thing to hold onto

The recurring pattern across R1-6, R1-8, R1-11, R1-12, R1-13, R1-14, R1-22, and
R1-30 is a check that was specified as strict, implemented as permissive, and left
with a docstring still describing the strict version. In isolation each is minor.
Together they mean the repository currently reports success in every case where it
was designed to report failure — and this is a repository whose entire claim is that
the author verifies their own work.

Fix the verification path before the features. A toolbox that certifies a little less
than it could is a good toolbox; one that cannot tell you when it is wrong is not.
