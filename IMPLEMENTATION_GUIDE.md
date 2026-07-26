# IMPLEMENTATION_GUIDE.md

Implementation instructions for `lyapunov-toolbox`. The repository skeleton is
complete: every file, module, function, and test exists with its signature and
docstring. Your job is to replace every `error("not implemented — ...")`, every
`raise NotImplementedError`, and every `TODO(impl)` with working code, and to fill
the prose sections of the notebooks and the `.tex` derivations.

**Read this whole document before writing code.** Section 3 (Global conventions)
and Section 4 (The mathematics, once) contain decisions that, if violated,
produce failures far from where the mistake was made — and, worse, produce
*plausible wrong answers* rather than errors.

---

## 0. Executive summary

| | |
|---|---|
| **What this is** | A Sum-of-Squares toolbox that certifies regions of attraction for polynomial vector fields, in Julia (SumOfSquares.jl + Clarabel) with an independent MATLAB implementation (YALMIP + Mosek) and a Python bridge. |
| **Who reads it** | Robotics/controls PhD admissions committees at theory-heavy labs. Correctness and honesty outrank features. |
| **Languages** | Julia ~40%, Python ~40%, MATLAB ~20%. |
| **Target size** | ~1,500 LOC of implementation, excluding notebooks. |
| **Definition of done** | §12. |
| **Non-negotiable** | An ROA estimate that is **too large is a false certificate**. Every design decision in this repository trades away tightness for soundness. Never make the opposite trade. |

### The one thing to understand before anything else

This repository produces *proofs*. A number that is slightly too conservative is
a worse result. A number that is slightly too large is a **wrong** result, and a
reader who catches one will discount every other claim in the portfolio. Every
tolerance, every solver-status mapping, and every sampling check in this guide
exists to make the second failure mode impossible. If you are ever tempted to
relax one of them to make something pass, you have found the exact spot where
this repository would stop being worth having.

---

## 1. Repository map

```
lyapunov-toolbox/
├── .github/workflows/test.yml   DONE — do not weaken a job to make it pass
├── .gitignore                   DONE
├── IMPLEMENTATION_GUIDE.md      this file
├── INFO.md                      portfolio spec — read for context, do not edit
├── README.md                    STUB — rewrite LAST (§11)
├── LICENSE                      DONE
├── julia/
│   ├── Project.toml             DONE except the solver dep (§2 Phase 0)
│   ├── Manifest.toml            EMPTY — generate in Phase 0
│   ├── src/
│   │   ├── LyapunovToolbox.jl   module + exports (§5.0)
│   │   ├── polynomial_utils.jl  primitives (§5.1)   ← START HERE
│   │   ├── systems.jl           the four plants (§5.2)
│   │   ├── sos_lyapunov.jl      the SDPs (§5.3)
│   │   └── roa_estimation.jl    the search loops (§5.4)
│   └── test/runtests.jl         DONE — thin includer for tests/
├── python/
│   ├── pyproject.toml           DONE
│   ├── lyapunov_toolbox/
│   │   ├── __init__.py          re-exports (§6.0)
│   │   ├── julia_bridge.py      the bridge (§6.1)
│   │   └── plotting.py          figures (§6.2)
│   └── tests/                   conftest + bridge integration tests (§8.2)
├── matlab/
│   ├── compute_roa_yalmip.m     the independent implementation (§7)
│   ├── van_der_pol_roa.m        script (§7)
│   └── pendulum_roa.m           script (§7)
├── examples/                    4 dirs, each a .jl script + a notebook (§9)
├── tests/
│   ├── fixtures/                committed MATLAB cross-validation output
│   ├── test_polynomial_utils.jl (§8.1)
│   ├── test_sos_lyapunov.jl     (§8.1)
│   └── test_roa_estimation.jl   (§8.1)
└── docs/
    ├── theory_notes.md          OUTLINE — fill every checkbox (§10.1)
    ├── benchmarks.md            SCHEMA — fill from code output (§10.2)
    └── derivations/
        ├── preamble.tex         DONE (shared macros)
        ├── sos_decomposition.tex        OUTLINE (§10.3)
        └── positivstellensatz.tex       OUTLINE — flagship (§10.3)
```

### Deviations from `INFO.md`, and why

Three files are not in the original specification. Each was added for a reason;
keep all three.

1. **`julia/src/systems.jl`** — the four benchmark systems appear in the tests,
   in four examples, in two MATLAB scripts, and in the benchmark table. Without
   one definition those copies drift, and a Julia/MATLAB agreement test that
   compares two *different* Van der Pol oscillators is worse than no test.
2. **`julia/test/runtests.jl`** — `Pkg.test()` requires this exact path. It
   contains no assertions; it includes the files in `tests/`, which stay where
   the specification puts them.
3. **`docs/derivations/preamble.tex`** — shared macros, matching the convention
   in the sibling `mpc-foundations` repository.

One substitution: the specification says the Python bridge is **PyJulia**. The
default is `juliacall` (PythonCall.jl) because PyJulia requires the `python-jl`
launcher on statically linked CPython, which is the usual reason such a bridge
works on one machine and not the next. **PyJulia is still supported** as a
selectable backend, is installable via the `pyjulia` extra, and both backends are
tested for agreement. If you disagree with the substitution, keep both anyway —
the cost is one branch in `JuliaSession._start_*`.

---

## 2. Build order

Follow this order. Each phase is testable on its own. Do not start a phase until
the previous one's tests pass, and do not implement the whole package and then
run the suite once — you will not be able to localise the failure.

| Phase | Files | Gate before moving on |
|---|---|---|
| **0** | `julia/Project.toml`, `julia/Manifest.toml` | `Pkg.instantiate()` succeeds; `using LyapunovToolbox` loads; `is_sos` on `x^2` returns true |
| **1** | `polynomial_utils.jl` | `test_polynomial_utils.jl` fully green — no solver involved, should take < 5 s |
| **2** | `systems.jl` (Van der Pol + `simulate` only) | the reversed VDP constructs; its Jacobian is `[0 -1; 1 -1]`; the forward VDP's limit cycle is computed by `simulate` |
| **3** | `sos_lyapunov.jl`: `default_solver`, `is_sos`, `sos_decomposition` | the `is_sos` and `sos_decomposition` testsets pass, **including Motzkin** |
| **4** | `lyapunov_from_linearization`, `verify_sublevel_set` (unconstrained, P2) | the `verify_sublevel_set` testset passes, including monotonicity |
| **5** | `roa_estimation.jl`: `bisect_level_set`, `estimate_roa`, `roa_volume`, `roa_boundary_2d` | Van der Pol reference result: certified set inside the limit cycle |
| **6** | `verify_certificate` | its testset passes, including the two tampering cases |
| **7** | `v_step`, `beta_step`, `v_s_iteration` | β monotone; V–s beats fixed-V bisection |
| **8** | `systems.jl` remainder: pendulum, unicycle, consensus; `trig_recast`; (P3) | their testsets pass |
| **9** | `python/lyapunov_toolbox/julia_bridge.py` | bridge integration tests pass with a real Julia |
| **10** | `python/lyapunov_toolbox/plotting.py` | four hero figures render and save |
| **11** | `examples/*/*.jl` | all four exit 0 under CI |
| **12** | `examples/*/*.ipynb` | all four execute under `pytest --nbmake` |
| **13** | `matlab/*.m` | VDP agrees with Julia to 1e-3; fixture committed |
| **14** | `docs/*`, `README.md` | PDFs compile; every README number traceable to code |

### Phase 0 in detail

```bash
julia --project=julia -e 'using Pkg; Pkg.add("Clarabel"); Pkg.resolve(); Pkg.instantiate()'
julia --project=julia -e 'using LyapunovToolbox'
```

If `instantiate` complains that a package is not registered, a UUID in
`Project.toml` is wrong. **Do not hand-fix it.** Delete the offending `[deps]`
line and run `Pkg.add("ThatPackage")`, which writes the authoritative UUID. Then
add `Clarabel` to `[extras]` and to the `test` target, commit the generated
`Manifest.toml`, and move on.

---

## 3. Global conventions

These are load-bearing. Violating any one of them produces a bug that looks like
it lives somewhere else.

### 3.1 Sign conventions

- **Continuous time throughout.** `ẋ = f(x)`. There is no discrete-time path in
  this repository; do not add one, and do not write code that half-assumes one.
- **`V̇ = ∇V · f`**, and the certificate wants `V̇ < 0`. Every SOS constraint is
  written on `-(V̇ + ε‖x‖²)`, i.e. on the *negative* of the Lie derivative. There
  is no place in this repository where `+V̇` appears inside an SOS constraint. If
  you find yourself writing one, you have a sign error upstream.
- **LQR feedback is `u = -Kx`.** The pendulum and the unicycle both depend on it.
- **The Lyapunov equation is `AᵀP + PA = -Q`** (continuous time). MATLAB's `lyap`
  uses this convention; verify the residual anyway rather than trusting it.

### 3.2 Variable ordering

`x` is a `Vector` of variables in **physical order**, fixed by the system that
created it. Never obtain the ordering from `MultivariatePolynomials.variables(p)`
— that is the library's sort order, and it differs. This single rule prevents:

- a Lie derivative that pairs `∂V/∂x₁` with `f₂`;
- an exponent matrix that evaluates correctly in Julia and incorrectly in Python;
- a Gram matrix that cannot be re-read against its own basis.

Every function that takes a polynomial and a variable vector uses the vector.

### 3.3 Degrees

- `deg V` and `deg λ` are **even**. Validate and raise otherwise.
- **The pairing rule:** `deg(λ) + deg(V) ≥ deg(V̇)`. `deg(V̇) = deg(V) - 1 + deg(f)`,
  so `deg(λ) ≥ deg(f) - 1`, rounded up to even. If this is violated, the
  S-procedure term cannot cancel the top-degree part of `V̇` and the program is
  infeasible **for structural reasons that look exactly like numerical ones**.
  Check it in `SOSOptions` construction and raise with both degrees in the
  message.
- An SOS polynomial of degree `2d` is parameterised by a basis of degree `d`.
  `sos_basis(x, d)` returns the degree-`d` basis. Halving errors here are the
  single most common bug in SOS code.

### 3.4 Numerical tolerances

| Purpose | Value |
|---|---|
| `ε` (positive-definiteness margin) | `1e-6` |
| Gram residual, accept | `1e-7` |
| Gram minimum eigenvalue, accept | `> -1e-8` |
| Lyapunov / Riccati equation residual | `1e-9` relative |
| ρ-bisection relative stopping width | `1e-4` |
| Julia ↔ MATLAB agreement on ρ | `1e-3` relative |
| Julia ↔ Python evaluation of `V` | `1e-10` |
| `boundary_2d` cross-implementation | `1e-6` |
| Equilibrium residual `‖f(x_eq)‖∞` | `1e-12` |
| Trig recast exactness | `1e-12` |
| Simulation convergence, `‖x‖` | `1e-3` |

Do not tighten a tolerance to make a test pass, and do not loosen one either. If
a test fails at `1e-7`, the code is wrong.

### 3.5 Failure behaviour

Raise; do not warn; never return a plausible-looking wrong answer.

- `ArgumentError` for bad user input: shape mismatch, odd degree, non-Hurwitz
  linearisation, `f(x_eq) ≠ 0`, a degree pairing violation.
- `ErrorException` (Julia) / `RuntimeError` (Python) for numerical failure: a
  Lyapunov residual over tolerance, a bisection that found no infeasible upper
  bound, a Gram matrix that does not reconstruct its polynomial.
- **Infeasibility of an SOS program is not an exception.** It is
  `feasible = false` on the returned tuple. Only the top-level `synthesize_lyapunov`
  raises, and only when *nothing* could be certified.
- **An unrecognised solver status is an exception.** Mapping it to "infeasible"
  silently shrinks every reported ROA; mapping it to "feasible" produces a false
  certificate. Neither is acceptable. Throw.
- `ALMOST_OPTIMAL` / `SLOW_PROGRESS` is accepted **only** after the returned Gram
  matrices pass `gram_residual < 1e-7` and `eigmin > -1e-8`. This is the single
  most important line of code in the repository. Write it once, in
  `verify_sublevel_set`, and route everything through it.

### 3.6 Determinism

Every random draw goes through an explicitly passed `Random.AbstractRNG` (Julia)
or `numpy.random.Generator` (Python). No `Random.seed!`, no `np.random.seed`, no
module-level global RNG. When a randomized test fails, the assertion message must
contain the seed that reproduces it.

Conic solvers are deterministic given identical input. If two runs with the same
seed disagree, something is iterating over an unordered collection — usually a
`Dict` of multipliers. Use ordered containers or sort the keys.

### 3.7 Style

- Julia: 4-space indent, ~92-column lines, `snake_case` functions, `CamelCase`
  types. Docstrings on every exported symbol. No `@warn` in library code.
- Python: type hints on every public signature — they are already written, keep
  them. `ruff check lyapunov_toolbox tests` clean, line length 100. No `print` in
  library code.
- MATLAB: full help block on every function, ALL-CAPS function name in the
  summary line, `inputParser` for options.
- **Keep the existing docstrings.** Extend them if you learn something while
  implementing; do not delete the rationale paragraphs — they are the point of
  the file. Comment density: match what is there. Explain *why*, never *what*.

### 3.8 Dependencies

Julia: `SumOfSquares`, `JuMP`, `DynamicPolynomials`, `MultivariatePolynomials`,
`SemialgebraicSets`, `Clarabel` (or `SCS`), plus stdlib. **Do not** add
`DifferentialEquations`, `Plots`, `ControlSystems`, or `DynamicPolynomials`
extensions — `simulate` is a hand-written RK4 for a reason (§5.2), plotting lives
in Python, and the Riccati/Lyapunov solves are small enough to write directly.

Python: `numpy`, `matplotlib`, `juliacall`. Nothing else at runtime.

MATLAB: YALMIP + one SDP solver. Nothing else.

---

## 4. The mathematics, once

Everything below is implemented somewhere in the package. Read this section
before §5; the code comments assume it.

### 4.1 The certificate

Given `ẋ = f(x)` polynomial with `f(0) = 0`, and `V` polynomial with `V(0) = 0`:

```
(P1)  V - ε‖x‖²                        ∈ Σ[x]          V is positive definite
(P2)  λ                                ∈ Σ[x]
      -(V̇ + ε‖x‖²) + λ·(V - ρ)         ∈ Σ[x]
```

**Why (P2) works.** On `Ω_ρ = {V ≤ ρ}` we have `V - ρ ≤ 0` and `λ ≥ 0`, so
`λ(V - ρ) ≤ 0`. The SOS condition gives
`-(V̇ + ε‖x‖²) + λ(V - ρ) ≥ 0`, hence `V̇ ≤ -ε‖x‖² - λ(V-ρ)·(-1) ≤ -ε‖x‖² < 0`
for `x ≠ 0` on `Ω_ρ`. Three lines; write them out in the `.tex`.

**The hypothesis everyone drops:** `Ω_ρ` must be **bounded**. Without it, "V̇ < 0
on `Ω_ρ`" does not give invariance or convergence. `V = x₁²` in a two-state
system is the counterexample to keep in mind. `is_sublevel_set_bounded` exists
for this and every certificate-producing entry point calls it.

**Why `ε > 0`.** SOS gives `≥ 0`, not `> 0`. Without the margin, `V̇ = 0` on a set
of measure zero inside `Ω_ρ` is permitted and the convergence argument fails.

### 4.2 Constrained version (Putinar)

With inequality constraints `gⱼ(x) ≥ 0` and equality invariants `hₖ(x) = 0`:

```
(P3)  λ, σⱼ ∈ Σ[x],     μₖ free (NOT SOS)
      -(V̇ + ε‖x‖²) + λ·(V - ρ) - Σⱼ σⱼ·gⱼ - Σₖ μₖ·hₖ  ∈ Σ[x]
```

- `σⱼ` must be SOS: `gⱼ ≥ 0` on the set, so `-σⱼgⱼ ≤ 0` there only if `σⱼ ≥ 0`.
- `μₖ` must be **free-sign**: `hₖ = 0` on the set, so `μₖhₖ = 0` regardless.
  Constraining `μₖ` to be SOS is valid but strictly weaker, and the difference is
  measurable — a test measures it.
- Putinar's theorem needs the **Archimedean** condition. Add a redundant ball
  constraint `N - ‖x‖² ≥ 0` to any compact constraint set to obtain it. Do this
  for the unicycle, record the `N` used, and say so in the `.tex`.

### 4.3 The search

`ρ` enters (P2) bilinearly with `λ`, so the joint problem is not an SDP. Two
facts break it:

1. **Feasibility is monotone in `ρ`.** If `Ω_ρ` admits a certificate then so does
   every `Ω_ρ'` with `ρ' < ρ`, because it is a subset. Therefore **bisect**.
   This one line is the entire justification for the method; state it in a
   comment at the top of `bisect_level_set`.
2. For fixed `ρ`, (P2) is **linear in `λ`** — an ordinary SDP.

Searching over `V` too gives the V–s iteration, with `ρ` normalised to 1:

```
(P4) V-step, λ and β fixed:
     V - ε‖x‖²                        ∈ Σ[x],   V(0) = 0
     -(V̇ + ε‖x‖²) + λ·(V - 1)         ∈ Σ[x]
     (1 - V) - σ_β·(β - p)            ∈ Σ[x],   σ_β ∈ Σ[x]
     maximise β
```

- **`ρ ≡ 1` is a normalisation, not a restriction:** `{V ≤ ρ} = {V/ρ ≤ 1}`.
  Without it the V-step is degenerate and reports rescalings as improvements.
- **The shape set.** `max ρ` is meaningless once `V` is free (rescale `V`). The
  objective is `max β` subject to `{p ≤ β} ⊆ {V ≤ 1}` with `p` a *fixed* positive
  definite polynomial, usually `p = xᵀx`. That containment is
  `(1 - V) - σ_β(β - p) ∈ Σ[x]` — check the S-procedure direction: on `{p ≤ β}`,
  `β - p ≥ 0`, so `1 - V ≥ σ_β(β-p) ≥ 0`, i.e. `V ≤ 1`. ✓
- `(β, σ_β)` is bilinear too. Bisect on `β`.

**The iteration is monotone in `β` by construction** (the previous iterate is
feasible for the next step) **and has no optimality guarantee.** Both halves of
that sentence go in the docstring, the notebooks, and the README.

### 4.4 Recasting

- **Taylor:** `sin θ ≈ θ - θ³/6`. This is a *different vector field*. A
  certificate for it is not a certificate for the original system without a
  truncation-error bound, which this repository does not compute. The pendulum
  example measures the gap empirically instead, and says that is what it is doing.
- **Trigonometric:** `s = sin θ`, `c = cos θ`, `ṡ = c θ̇`, `ċ = -s θ̇`, invariant
  `s² + c² - 1 = 0`. Exact. One extra state per angle, plus a free-sign multiplier.
- **Sampling on a manifold.** For a trig-recast system, a point drawn from the
  ambient box is not a state of the system. Testing `V̇` there proves nothing and
  can make an invalid certificate look valid. Sample `θ` and set `(s, c)`.

---

## 5. Julia implementation

### 5.0 `LyapunovToolbox.jl`

Nothing to implement beyond keeping the include order (`polynomial_utils` →
`systems` → `sos_lyapunov` → `roa_estimation`) and the export list in sync.
Verify with:

```julia
julia --project=julia -e 'using LyapunovToolbox; for n in names(LyapunovToolbox) @assert isdefined(LyapunovToolbox, n) n end'
```

### 5.1 `polynomial_utils.jl`

**`monomial_basis(x, degree; mindegree=0)`** — use
`DynamicPolynomials.monomials(x, mindegree:degree)`. Assert the cardinality
equals `binomial(n+degree, degree) - binomial(n+mindegree-1, mindegree-1)`; that
one assertion catches every enumeration bug.

**`sos_basis(x, d)`** = `monomial_basis(x, d)`. It exists as a separate name so
that call sites read correctly and the halving is explicit.

**`gram_dimension(n, d)`** = `binomial(n + d, d)`.

**`lie_derivative(V, f, x)`** — `sum(differentiate(V, x[i]) * f[i] for i in eachindex(x))`.
Assert `length(f) == length(x)` first and raise `ArgumentError` with both lengths.
Do **not** write `differentiate(V, variables(V))` (§3.2).

**`jacobian_at(f, x, x0)`** — `differentiate(f, x)` gives a matrix of
polynomials; substitute `x => x0` and convert to `Float64`. Symbolic, not finite
differences: `f` is polynomial, so there is no reason to accept truncation error.

**`taylor_recast(expr_fn, x, order; about)`** — build the multivariate Taylor
polynomial by repeated symbolic differentiation of `expr_fn` applied to
`ForwardDiff`-free scalar arguments. For the systems here, only *univariate*
expansion of `sin`/`cos` in one angle is needed; implement the univariate case
properly, and raise `ArgumentError` for the multivariate one rather than shipping
a half-correct general version. Include the truncation-error estimate in the
returned metadata so the pendulum example can report it.

**`trig_recast(f_trig, x, θ_index)`** — returns `(f_poly, x_new, invariants)`:
replace `θ` by two variables `s, c`; substitute `sin θ → s`, `cos θ → c`; append
`ṡ = c·θ̇` and `ċ = -s·θ̇` to the field in the position where `θ̇` sat; return
`invariants = [s^2 + c^2 - 1]`. Order the new state vector so the substitution is
unambiguous and document the order in the returned system's `notes`.

**`prune_terms(p; tol)`** — drop terms with
`abs(coef) < tol * maximum(abs, coefficients(p))`. **Display and storage only.**
Never verify a certificate against a pruned polynomial: pruning can break the SOS
inequality. Verify first, prune second.

**`evaluate_on_grid(p, x, xs, ys; fixed)`** — returns `(length(ys), length(xs))`,
`[row = y, col = x]`, so it goes straight into a contour plot with no transpose.
Test with an asymmetric polynomial (`p = x₁`), because a symmetric one passes
whether or not you transposed.

**`gram_residual(p, basis, Q)`** — build `basis' * Q * basis`, subtract `p`,
return the Frobenius norm of the coefficient vector of the difference. Must not
call a solver and must not ask the modelling layer whether it is satisfied.

**`is_sublevel_set_bounded(V, ρ, x; ...)`** — two checks:

1. Extract the top-degree homogeneous part `V_top`. If `V_top - ε·(xᵀx)^(deg/2)`
   is SOS, every sublevel set is bounded → `(true, :top_form_pd)`. This is
   conclusive.
2. Otherwise sample the sphere of radius `radius`; if any sample has `V ≤ ρ`,
   return `(false, :unbounded_witness)`; if none does, return
   `(true, :sampling_only)`.

`:sampling_only` is **not proof** and callers must record which they got. Note
that check 1 calls a solver, so pass the options through.

### 5.2 `systems.jl`

**`PolynomialSystem` inner constructor** — validate `length(f) == length(x)`,
that every `f[i]` is polynomial in `x` only, and `‖f(equilibrium)‖∞ < 1e-12`.
Put the actual residual in the `ArgumentError` message.

**`van_der_pol(; reversed=true, μ=1.0)`**

```
reversed:  ẋ₁ = -x₂,   ẋ₂ = x₁ + μ(x₁² - 1)x₂
forward:   ẋ₁ =  x₂,   ẋ₂ = -x₁ - μ(x₁² - 1)x₂
```

Jacobian of the reversed system at the origin is `[0 -1; 1 -μ]` — Hurwitz for
`μ > 0`. Verify that in the constructor. State labels `("x₁", "x₂")`, no units
(it is a normalised oscillator; do not invent units).

**`inverted_pendulum(...)`**

```
θ̈ = (g/l)·sin θ - (b/(m l²))·θ̇ + u/(m l²),    u = -K·[θ, θ̇]
```

Linearisation about upright: `A = [0 1; g/l  -b/(m l²)]`, `B = [0; 1/(m l²)]`.
Note `A[2,1] = +g/l` — upright is **unstable** open-loop. If you get a Hurwitz `A`
before feedback, the sign is wrong.

When `K === nothing`, compute the LQR gain with `Q = diag([10, 1])`, `R = [1]`.
Write the continuous Riccati solve directly (Hamiltonian matrix, take the stable
invariant subspace via `schur` + reordering) rather than adding a dependency; the
system is 2×2 and the code is fifteen lines. Verify `eig(A - BK)` is Hurwitz and
that the Riccati residual is under `1e-9` before returning.

`recast = :taylor` → replace `sin θ` with the odd Taylor polynomial of order
`taylor_order` (assert odd order). `recast = :trig` → `trig_recast`, states
`(s, c, θ̇)`, invariant `s² + c² - 1`.

State labels `("θ [rad]", "θ̇ [rad/s]")`.

**`unicycle_with_obstacle(...)`** — `(ẋ, ẏ, θ̇) = (v cos θ, v sin θ, ω)` with a
smooth polynomial stabilising feedback in the recast states, goal at the origin,
constraint `g(x,y) = (x - x_o)² + (y - y_o)² - r² ≥ 0`, plus the Archimedean ball
constraint (§4.2). Recast the heading exactly — states `(x, y, s, c)`, invariant
`s² + c² - 1`.

Choosing the feedback is the only genuinely open design decision in this file.
Requirements: polynomial in the recast states, `f(0) = 0`, and a Hurwitz
linearisation at the goal. A polynomial approximation of the standard
`(k_ρ, k_α)` polar controller works; so does a simpler
`v = -k₁(x·c + y·s)`, `ω = -k₂(y·c - x·s) - k₃·(something)` form. **Verify the
Hurwitz condition numerically before running any SOS program** — an infeasibility
caused by a non-stabilising gain is indistinguishable from a real one, and you
will lose an afternoon.

**`consensus_network(; n=3, laplacian, coupling, κ)`** — dynamics
`ẋᵢ = -Σⱼ aᵢⱼ φ(xᵢ - xⱼ)`, `φ(z) = z + κz³`. Default `L` is the 3-node path graph.

**Return the system in disagreement coordinates.** Build `U ∈ ℝ^{n×(n-1)}` with
orthonormal columns spanning `𝟙^⊥` (QR of `[𝟙 | I]`, drop the first column, or
the normalised difference basis — either is fine, but store `U` and verify
`Uᵀ𝟙 = 0` and `UᵀU = I`). Then `z = Uᵀx`, `x = Uz` on `𝟙^⊥`, and
`ż = Uᵀ f(Uz)`. Verify `f(Uz)` has no `𝟙` component so the reduction is exact.
Store `U` in `notes` and expose it — the example needs it to map back.

The certificate must be invariant to the average state; verify that numerically
in the example rather than asserting it.

**`simulate(sys, x0, tspan; dt, max_norm)`** — fixed-step RK4, hand-written.
Returns `(t, X, status)` with `status ∈ (:converged, :diverged, :timeout)`.
Bail out as soon as `‖x‖∞ > max_norm`. `X` is `(length(t), n)` — time first,
matching the Python convention in the sibling repositories.

### 5.3 `sos_lyapunov.jl`

**`default_solver()`** — `Clarabel.Optimizer` with tolerances set here and only
here:

```julia
optimizer_with_attributes(Clarabel.Optimizer,
    "tol_gap_abs" => 1e-9, "tol_gap_rel" => 1e-9,
    "tol_feas" => 1e-9, "max_iter" => 500, "verbose" => false)
```

If Clarabel is unavailable, fall back to `SCS` with `eps_abs = eps_rel = 1e-9` and
**log which one is in use** — a benchmark row that does not name its solver is not
reproducible.

**`is_sos(p, x; solver, ε)`** — short-circuit `false` for odd `maxdegree(p)`
without building a model. Otherwise:

```julia
model = SOSModel(solver)
@constraint(model, p - ε * sum(x.^2)^k >= 0)   # k chosen so degrees match, or ε = 0
optimize!(model)
```

Map `OPTIMAL` → `true`, `INFEASIBLE`/`DUAL_INFEASIBLE`/`INFEASIBLE_OR_UNBOUNDED`
→ `false`, `ALMOST_OPTIMAL` → verify the Gram matrix and decide, anything else →
throw. This is the first thing to get working (Phase 3): it is the smallest
program that exercises the whole stack.

**`sos_decomposition(p, x; solver)`** — build the model, `optimize!`, extract the
Gram matrix with `gram_matrix(constraint_ref)` (from SumOfSquares), symmetrise,
and check `eigmin`. Clip eigenvalues in `[-1e-8, 0)` to zero and record the
largest clipped magnitude. If `eigmin < -1e-7`, **do not clip** — return a failure
status. Also verify `gram_residual < 1e-7` before returning success.

**`lyapunov_from_linearization(sys; Q)`** — `A = jacobian_at(sys.f, sys.x)`;
raise `ArgumentError` naming the offending eigenvalues if `maximum(real(eigvals(A))) ≥ -1e-9`.
Solve `AᵀP + PA = -Q` (default `Q = I`) via the Kronecker form
`(I⊗Aᵀ + Aᵀ⊗I) vec(P) = -vec(Q)` — `n ≤ 6` here, so the `n²×n²` solve is free —
then symmetrise, verify the residual under `1e-9·max(1, ‖P‖_F)`, verify `P ≻ 0`,
and return `V = x'Px` as a polynomial.

**`verify_sublevel_set(sys, V, ρ; opts)`** — the workhorse. Build (P2), or (P3)
when `sys.constraints` or `sys.invariants` is non-empty:

```julia
model = SOSModel(opts.solver)
Vdot  = lie_derivative(V, sys.f, sys.x)
@variable(model, λ, SOSPoly(sos_basis(sys.x, opts.multiplier_degree ÷ 2)))
expr = -(Vdot + opts.ε * sum(sys.x.^2)) + λ * (V - ρ)
for g in sys.constraints
    @variable(model, σ, SOSPoly(sos_basis(sys.x, opts.constraint_multiplier_degree ÷ 2)))
    expr -= σ * g
end
for h in sys.invariants
    @variable(model, μ, Poly(monomial_basis(sys.x, opts.constraint_multiplier_degree)))
    expr -= μ * h          # FREE-SIGN — Poly, not SOSPoly. §4.2.
end
@constraint(model, cref, expr >= 0)
optimize!(model)
```

Return `(feasible, witness)` with the multipliers, every Gram pair, the raw
status string, and `solve_time`. Do not mutate `opts`.

**`multiplier_step`** is an alias for this — implement one, have the other call it.

**`v_step(sys, λ, shape, β; opts)`** — (P4). Parameterise `V` over
`monomial_basis(x, V_degree; mindegree=2)` **and** add `V(0) = 0` as an explicit
linear equality on the coefficients. `mindegree=2` alone is not enough: without
the equality, the solver can place mass on terms that make `∇V(0) ≠ 0`. Objective
`Max β` with `β` a JuMP variable only if you also fix `σ_β` — otherwise bisect on
`β` and solve a feasibility problem. Bisecting is simpler and is what the skeleton
assumes.

**`beta_step(sys, V, shape; opts)`** — bisection on `β` over
`(1 - V) - σ_β(β - p) ∈ Σ[x]`, `σ_β ∈ Σ[x]`. Same bisection machinery as
`bisect_level_set`; factor the bisection into a private helper taking a
`ρ -> Bool` oracle and reuse it for both.

**`verify_certificate(cert; tol, samples, rng)`** — **must not call a solver.**
Three checks, reported separately:

1. For each stored `(basis, Q)`: `gram_residual < tol` and `eigmin(Symmetric(Q)) > -tol`.
2. Rejection-sample `samples` points from `{V ≤ ρ}`; assert `V̇(x) < 0` for every
   one with `‖x‖ > 1e-6`. Report the worst `V̇` and its location. For a system
   with invariants, **sample on the manifold** (§4.4).
3. Integrate from ~100 interior points; assert the trajectory stays in `{V ≤ ρ}`
   and reaches `‖x‖ < 1e-3`.

Check 1 is the proof. Checks 2 and 3 catch the case where the program solved is
not the program intended — a sign error, the wrong `f`, a stale `V` — all of
which check 1 would pass happily.

Return a `NamedTuple` with each check's result and the diagnostics, not a single
`Bool`.

### 5.4 `roa_estimation.jl`

**`bisect_level_set`** — bracket, then bisect.

```
if ρ_hi === nothing
    ρ_hi = 1e-3
    for _ in 1:30
        verify_sublevel_set(sys, V, ρ_hi) infeasible && break
        ρ_hi *= 2
    end
    still feasible after 30 doublings  →  throw
end
best = ρ_lo
while (ρ_hi - ρ_lo) > opts.bisection_tol * max(1, ρ_lo)
    mid = (ρ_lo + ρ_hi) / 2
    if feasible(mid);  best = mid;  ρ_lo = mid  else  ρ_hi = mid  end
end
return best   # a VERIFIED value, never the midpoint of the final bracket
```

The `throw` on 30 doublings is not a limitation to work around. For a globally
stable polynomial system there is no finite boundary, and returning the cap as if
it were one is a silent lie. Say so in the exception message.

**`estimate_roa`** — `bisect_level_set`, then `is_sublevel_set_bounded`, then
package a `LyapunovCertificate`. Throw if not bounded, with the evidence symbol
in the message.

**`v_s_iteration`** — the loop in the docstring. Three assertions, all required:

- `β_new ≥ β_old - 1e-8`, else throw with both values (§4.3 says why it cannot
  decrease, so a decrease is a bug).
- the V-step is never infeasible (the previous `V` is feasible for it), else throw.
- renormalise `V ← V/ρ` after each V-step so `ρ ≡ 1` holds on entry to the next
  γ-step.

Record `ρ`, `β`, and volume per iteration. Set `converged = false` when the
budget runs out, and never present that as convergence.

**`roa_volume`** — `:grid` for tests (deterministic), `:montecarlo` for reporting.
Monte Carlo returns `(volume, standard_error)`; report the error, or a 3%
"improvement" that is noise ends up in the README. For `n > 3`, say in the docs
that the number is a sanity check and that `β` is the comparison that matters.

**`roa_boundary_2d`** — march along `n` rays from the origin, bisect on radius for
`V(r·d) = ρ`. Detect multiple crossings along a ray (sample coarsely first, count
sign changes of `V - ρ`) and **raise** rather than returning the innermost one.
Angular order comes free from the ray construction.

**`true_roa_boundary_2d`** — bisect on radius per direction, using `simulate` and
a convergence test. It is an empirical estimate at a stated `tspan` and `dt`;
label it as such everywhere it appears. For Van der Pol, prefer the limit-cycle
route (§9.1) and say why it is stronger.

**`containment_fraction`** — return `(fraction_of_true_captured, fraction_of_certified_outside_true)`
in that order. The second **must be zero**. If it is not, the certificate is
wrong and no volume number matters.

**`benchmark_row`** — a `NamedTuple` with every column in `docs/benchmarks.md` §1.
Generate the table from these; do not type it.

---

## 6. Python implementation

### 6.0 `__init__.py`

Uncomment the re-export block; keep `julia_bridge`'s Julia startup lazy. Verify:

```bash
python -c "import lyapunov_toolbox as m; print(sorted(set(m.__all__) - set(dir(m))))"
```

Must print `[]`. And `import lyapunov_toolbox` must take under 2 seconds — there
is a test for it, because downstream repositories import this package for its
plotting alone.

### 6.1 `julia_bridge.py`

**Startup.** For `juliacall`, set `PYTHON_JULIACALL_PROJECT` **before** importing
`juliacall` — it reads the environment at import time and ignores it afterwards,
which is the usual reason the wrong project ends up active. For PyJulia, use
`julia.api.Julia(runtime=..., compiled_modules=False)` when
`sysconfig.get_config_var("Py_ENABLE_SHARED")` is falsy, and explain that in the
exception rather than letting the process die.

**Never instantiate implicitly.** `instantiate=False` by default. The error when
the project is not instantiated must contain the resolved path and the literal
command to fix it. There is a test asserting on that message content; that is
deliberate, because it is the error every new user hits first.

**Marshalling.** Polynomials cross as `(exponents, coefficients)`:

```julia
# Julia side
exps = [Int(degree(m, v)) for m in monomials(p), v in x]   # note the ORDER: x, not variables(p)
coefs = coefficients(p)
```

Assert the round trip on the Julia side before returning. A permuted exponent
matrix evaluates correctly on each side in isolation and disagrees only when
compared — which is why the test compares.

**`RoaCertificate.evaluate`** in NumPy:

```python
np.prod(x[:, None, :] ** self.exponents[None, :, :], axis=-1) @ self.coefficients
```

This is why plotting needs no Julia. Keep it exact to `1e-10` against Julia; a
test enforces it.

**`_check_expression_safety`** — `certify_roa` evaluates its argument as Julia
source, and it is called from notebooks and from other repositories. Implement an
**allowlist** over a tokenised scan: declared variable names, decimal literals,
`+ - * / ^`, parentheses, whitespace. Everything else raises `ValueError` naming
the offending token. A denylist of dangerous words is not a security control and
must not be what ships here.

**`certify_roa`** runs `verify_certificate` on the Julia side before marshalling
and raises `JuliaBackendError` if it fails. A certificate with `verified=False`
is never returned.

### 6.2 `plotting.py`

`use_style` sets `rcParams` only: figsize `(7, 6)`, dpi 130, grid alpha 0.25,
`mathtext.fontset="cm"`, top/right spines off, and an explicit colour cycle so
"certified" and "estimated" are the same colours in every figure.

Every function takes an optional `ax` and returns the `Axes` (or `Figure`). None
call `plt.show()` or `plt.savefig()`.

Three rules that are correctness properties, not style:

1. **Certified sets are filled with a solid edge; estimated true boundaries are
   dashed and unfilled.** A reader must distinguish the proved thing from the
   simulated thing without reading the legend.
2. **Every 2-D slice states its pinned values in the title.** A slice presented
   as the whole set is how an ROA figure misleads. There is a test asserting on
   the title text.
3. **If any part of the certified set falls outside the true boundary, draw it in
   red and warn.** That is a false certificate; the figure is the last place it
   should be allowed to look normal.

`plot_vdot_field` uses `TwoSlopeNorm(vcenter=0)` with a diverging colormap — with
a sequential map the sign change, which is the entire content of the figure, is
invisible.

`save_figure` resolves `docs/figures/` relative to the module file, not the CWD.

---

## 7. MATLAB implementation

The MATLAB side is an **independent** implementation of the same mathematics
against a different modelling layer and a different solver. That independence is
the point: agreement between two implementations sharing no code is evidence in a
way that a repeated run of the same code is not.

**Do not refactor the two toward a shared formulation.** If they disagree, find
out which is wrong.

`compute_roa_yalmip` — the outline is in the file's `TODO(impl)` block. The
critical points:

- `Vdot = jacobian(V0, x) * f;` — YALMIP's `jacobian` returns a row vector, so
  this is row × column. A transpose error gives a scalar of the wrong sign and a
  certificate that looks plausible.
- Treat `sol.problem == 4` (numerical problems) as **infeasible**, and record it.
  Accepting it inflates `ρ`, which is the one direction of error this file must
  not make.
- Any other nonzero `sol.problem` that is not `1` is an error: rethrow with
  `yalmiperror(sol.problem)`.
- Error, with an actionable message, when neither Mosek nor SeDuMi is on the
  path. Do not fall through to a default.
- Extract the Gram matrices via `sosd`/`getbase` and store them in `info`, so the
  MATLAB result can be re-verified without MATLAB.

`van_der_pol_roa` and `pendulum_roa` — outlines are in the files. Both must
**error**, not warn, on their soundness assertions (certified set inside the limit
cycle; certified set excluding `θ = ±π`).

CI has neither MATLAB nor Mosek. The fixture in `tests/fixtures/` is a committed
artefact and its commit message must record the MATLAB, YALMIP, and solver
versions. That message is the only provenance there is.

---

## 8. Tests

Every test body is a stub with a docstring or comment stating precisely what it
asserts. Implement the assertion described — do not weaken it, do not re-scope
it, and do not delete a test because it is inconvenient. If a test looks wrong,
implement it as written first; if it then fails and you are confident the test is
at fault, say so in your final report rather than editing it silently.

### 8.1 Julia

The suite is organised so that the cheap, solver-free tests
(`test_polynomial_utils.jl`) run first and localise most failures. Budget: the
default suite under three minutes. Guard anything slower behind
`get(ENV, "LYAPUNOV_SLOW_TESTS", "0") == "1"`, and run the full set in CI.

The tests that matter most, in order:

1. **`rejects the Motzkin polynomial`** — nonnegative but not SOS. This is the
   test that shows the author knows what SOS actually is.
2. **`certified set really has V̇ < 0 throughout`** — the anti-tautology test.
   Rejection-sample the certified set and check pointwise, with no solver. If the
   SDP is right and this fails, the Lie derivative or a sign is wrong.
3. **`feasibility is monotone in rho`** — the property that justifies bisection.
   Report the whole boolean sequence on failure.
4. **`certified set excludes theta = pi`** — cheapest end-to-end soundness check
   in the repository.
5. **`fails when rho is inflated after the fact`** — takes a valid certificate,
   multiplies `ρ` by 5 without re-solving, and requires verification to fail. A
   verifier that only re-reads its own stored Gram matrices would pass this and
   is therefore not a verifier.
6. **`attempting to certify in x coordinates is infeasible`** (consensus) —
   documents the trap that motivates the whole example.
7. **`no finite level is ever infeasible`** (globally stable linear system) —
   requires `bisect_level_set` to throw rather than report the doubling cap.

**Do not hard-code a `ρ` from a paper or from this guide.** The value depends on
`ε`, the multiplier degree, and the solver. Assert *properties* — containment in
the true ROA, monotonicity in degree, monotonicity in `κ` — and record the
numbers you obtain in `docs/benchmarks.md`. An assertion against a borrowed
constant tests nothing about this code.

Target: line coverage above 85% on `julia/src`, excluding the `show` method.

### 8.2 Python

Bridge tests are **integration** tests. A mocked bridge passes on a machine with
no Julia, which is precisely the machine where the bridge is broken. Tests
needing Julia carry `@pytest.mark.julia`; `conftest.py` skips them locally with a
reason naming the resolved project path, and **fails** them in CI, where
`LYAPUNOV_TOOLBOX_REQUIRE_JULIA=1` is set. A skipped bridge suite in CI is a
failed bridge suite.

Call `julia_available()` **once** in `pytest_collection_modifyitems` and cache it.
Per-test probing turns a 30-second suite into a five-minute one.

---

## 9. Examples

Each directory has a `.jl` script (run by CI on every commit, must exit non-zero
on any failed assertion) and a notebook (executed by `pytest --nbmake`).

### 9.1 Division of labour

The `.jl` scripts do the computation and write CSV to `docs/figures/data/`. The
notebooks read the CSV and do the plotting and the prose. This keeps the Julia
scripts free of a plotting dependency, keeps CI cheap, and means a notebook can
be re-run to restyle a figure without re-solving an SDP. Be consistent about it
across all four.

Van der Pol ground truth: integrate the **forward** oscillator from a point
outside the limit cycle, discard the transient, keep one period. That closed
curve is the exact ROA boundary of the reversed system. Prefer it to
`true_roa_boundary_2d`'s convergence sweep and say why in the output: the limit
cycle is exact up to integration error; the sweep is empirical at a finite horizon.

### 9.2 Notebook rules

1. They must execute top to bottom in a fresh kernel. CI runs
   `pytest --nbmake examples/`, 1200 s timeout per notebook. Budget accordingly.
2. Commit them **with outputs**. A reader on GitHub sees the rendered result; an
   output-stripped notebook shows nothing.
3. Every plot gets axis labels with units and a legend. Every claim in prose gets
   a number from a cell above it.
4. Replace each `**Write:**` block with actual prose. 3–8 sentences per block,
   first person, plain and specific.
5. Every notebook ends with a Takeaways cell — three or four bullets with
   numbers. No filler.
6. Where a notebook asks you to report a failure, a limitation, or a sensitivity
   — the pendulum's Taylor-model cross-check, the unicycle's slice caveats, the
   consensus switching counterexample — **report it honestly with the numbers.**
   Those sections establish credibility. A notebook where everything works and
   nothing is qualified reads as untested.
7. The consensus notebook's section 6 is the most important prose in the
   repository: it is the negative result that motivates the thesis. If the
   adversarial search finds no counterexample, say so with the number of
   sequences searched. Do not delete the section for producing an inconvenient
   answer.

---

## 10. Documentation

### 10.1 `docs/theory_notes.md`

Fill every checkbox. Sections 1 (notation) and 6 (computing the ROA) must be
airtight — the notation table has to match the code identifiers exactly, because
that is what makes the repository navigable.

### 10.2 `docs/benchmarks.md`

Every number comes from `benchmark_row` or from a printed script summary. **Fill
the environment block first**; a benchmark without one is not a benchmark.
Include the first configuration that *failed* — ran out of memory, exceeded ten
minutes, or returned a non-`OPTIMAL` status. The failure row is the most
informative line in the table.

### 10.3 `.tex` derivations

```bash
cd docs/derivations && latexmk -pdf positivstellensatz.tex
```

Order of importance: `positivstellensatz.tex` first — it is the flagship and the
document cited in the thesis — then `sos_decomposition.tex`.

Write the proofs in your own words. Where a result is quoted rather than proved
(Putinar's theorem, Hilbert's classification, Khalil Thm 4.1), **say so
explicitly and cite it**. An honest citation is stronger than a hand-waved proof.

Two sections carry disproportionate weight:

- `sos_decomposition.tex` §4, the fully worked example. Do it by hand, exhibit the
  free parameter that makes the Gram matrix non-unique, then show the solver's
  `Q` and explain why they differ while both being valid. That explanation is what
  distinguishes understanding SOS from calling a library.
- `positivstellensatz.tex` §5, the ROA theorem. Call out the boundedness
  hypothesis **explicitly at the step of the proof that uses it**. It is the
  hypothesis most often dropped, and the proof is false without it.

---

## 11. README.md

Rewrite last, once the numbers exist. Order:

1. Tagline (keep the existing one), CI + coverage badges.
2. The four hero figures as a 2×2 table, each captioned with its certified `ρ`
   and its fraction of the true ROA.
3. What is implemented, as a short Julia/MATLAB/Python table.
4. Install: three commands, no prose.
5. Quickstart: ten lines that certify the Van der Pol ROA and print `ρ`. Test it
   in a clean checkout before committing.
6. Results: the abridged benchmark table. Quote what you measured. If it is 78%,
   write 78%.
7. Example index, one line each, linking script and notebook.
8. Theory: links to `theory_notes.md` and the two compiled PDFs.
9. **"What this does not do"** — inner approximations only; floating-point
   certificates, not exact ones; no robustness; the V–s iteration is a local
   heuristic; the measured degree/dimension wall; no discrete time. Naming the
   limits is what makes everything above credible.
10. The citation hook, in the exact commit-pinned form to paste into a paper.
11. References, matching `theory_notes.md` §10.

---

## 12. Definition of done

```bash
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
julia --project=julia -e 'using Pkg; Pkg.test()'
LYAPUNOV_SLOW_TESTS=1 julia --project=julia -e 'using Pkg; Pkg.test()'
cd python && pip install -e ".[dev]" && ruff check lyapunov_toolbox tests && pytest && cd ..
for f in examples/*/*.jl; do julia --project=julia "$f"; done
pytest --nbmake examples/ --nbmake-timeout=1200
cd docs/derivations && latexmk -pdf sos_decomposition.tex && latexmk -pdf positivstellensatz.tex
```

MATLAB, manually:

```matlab
addpath('matlab'); van_der_pol_roa('WriteFixture', true); pendulum_roa
```

Plus, by inspection:

- [ ] No `not implemented` error remains in `julia/` or `tests/`.
- [ ] No `NotImplementedError` remains in `python/`.
- [ ] No `TODO(impl)` remains anywhere.
- [ ] No `**Write:**` block remains unfilled in the notebooks.
- [ ] Notebooks are committed with outputs.
- [ ] `julia/Manifest.toml` is generated and committed.
- [ ] `docs/benchmarks.md` has a filled environment block and no empty cells.
- [ ] `tests/fixtures/matlab_vdp.csv` is committed with version provenance in the
      commit message.
- [ ] Every number in `README.md` was produced by code in this repository.
- [ ] For every system: `containment_fraction`'s **second** return value is zero.

That last item is the one that matters. Everything else in this list is
housekeeping; that one is the difference between a toolbox and a false claim.

---

## 13. If you get stuck

Report the blocker; do not route around it.

- **Do not** widen a certified set to make a comparison look better.
- **Do not** accept an `ALMOST_OPTIMAL` solve without the Gram check.
- **Do not** loosen a tolerance, delete a test, or weaken a CI job.
- **Do not** hard-code a `ρ` from a paper so a test passes.
- **Do not** add a dependency that implements the thing you were asked to write.
- **Do** implement everything that is not blocked, and state precisely what you
  left undone and why.

A partial implementation with an honest gap list is useful. A complete-looking
one with a silently weakened test is not — and in a repository whose entire
purpose is to demonstrate that the author verifies their own claims, it is the
worst possible outcome.
