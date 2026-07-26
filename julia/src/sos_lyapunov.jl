# sos_lyapunov.jl — the semidefinite programs.
#
# Everything that talks to JuMP/SumOfSquares lives here. The outer loops
# (bisection, V–s iteration) live in roa_estimation.jl and call into this file;
# keeping the split clean is what makes it possible to unit-test the programs
# against known-feasible and known-infeasible instances without also testing the
# search heuristics.
#
# ------------------------------------------------------------------------------
# The programs, written out once so the code below can refer to them by name.
#
# (P1) Positivity of the candidate
#          V - ε‖x‖²  is SOS,        V(0) = 0
#
# (P2) S-procedure certificate of the sublevel set, ρ fixed
#          λ  is SOS
#          -(V̇ + ε‖x‖²) + λ·(V - ρ)  is SOS
#      On {V ≤ ρ} the term λ·(V - ρ) ≤ 0, hence V̇ ≤ -ε‖x‖² < 0 for x ≠ 0.
#      Note this is LINEAR in λ for fixed ρ and V, which is what makes the
#      ρ-bisection a sequence of ordinary SDPs rather than a bilinear mess.
#
# (P3) Constrained version, constraints g ≥ 0 and equalities h = 0 (Putinar)
#          λ, σⱼ  SOS,   μₖ  free
#          -(V̇ + ε‖x‖²) + λ·(V - ρ) - Σⱼ σⱼ·gⱼ - Σₖ μₖ·hₖ  is SOS
#      σⱼ SOS and gⱼ ≥ 0 on the set means -σⱼgⱼ ≤ 0 there; the μₖ are free-sign
#      because hₖ = 0 makes the sign irrelevant. Using an SOS multiplier on an
#      equality is not wrong but it is strictly weaker — do not do it by
#      accident. See docs/derivations/positivstellensatz.tex.
#
# (P4) V-step of the V–s iteration, λ fixed, ρ normalised to 1
#          V - ε‖x‖²  is SOS,   V(0) = 0
#          -(V̇ + ε‖x‖²) + λ·(V - 1)  is SOS
#          (1 - V) - σ_β·(β - p)  is SOS,   σ_β SOS      [shape containment]
#      maximise β. The normalisation ρ ≡ 1 removes the scaling degeneracy
#      V ↦ V/ρ; without it the V-step is unbounded and the solver will happily
#      report a "better" ρ that is pure rescaling.
# ------------------------------------------------------------------------------

"""
    SOSOptions

Every knob that changes the numerical answer, in one place, so that a benchmark
row can record the full configuration that produced it.

# Fields
- `solver` — a JuMP optimizer factory. `nothing` means [`default_solver`](@ref).
- `ε::Float64 = 1e-6` — the positive-definiteness margin in (P1)/(P2). It buys
  strictness: SOS gives `≥ 0`, and `V̇ ≤ -ε‖x‖²` is what actually rules out
  `V̇ = 0` on a set of measure zero inside the sublevel set. Do not set it to 0.
  Do not set it large enough to dominate the dynamics either — report it.
- `V_degree::Int = 2` — total degree of the Lyapunov candidate. Must be even.
- `multiplier_degree::Int = 2` — total degree of `λ`. Must be even. The pairing
  rule is in the guide: `deg(λ) + deg(V) ≥ deg(V̇)` or the S-procedure term
  cannot cancel the top-degree part of `V̇` and the program is infeasible for
  structural reasons that look like numerical ones.
- `constraint_multiplier_degree::Int = 2` — degree of the `σⱼ` in (P3).
- `bisection_tol::Float64 = 1e-4` — relative width at which ρ-bisection stops.
- `max_bisection_iters::Int = 40`
- `vs_iterations::Int = 10` — outer V–s rounds.
- `vs_improvement_tol::Float64 = 1e-4` — relative gain below which the V–s loop
  is declared converged.
- `verbose::Bool = false` — per-iteration logging. Library code prints nothing
  unless this is set.
- `solver_silent::Bool = true`
- `feasibility_tol::Float64 = 1e-7` — threshold on solver primal/dual residuals
  below which a `MOI.OPTIMAL` result is accepted.

Construct with keywords; all fields have the defaults above.
"""
Base.@kwdef struct SOSOptions
    solver = nothing
    ε::Float64 = 1e-6
    V_degree::Int = 2
    multiplier_degree::Int = 2
    constraint_multiplier_degree::Int = 2
    bisection_tol::Float64 = 1e-4
    max_bisection_iters::Int = 40
    vs_iterations::Int = 10
    vs_improvement_tol::Float64 = 1e-4
    verbose::Bool = false
    solver_silent::Bool = true
    feasibility_tol::Float64 = 1e-7
end

"""
    LyapunovCertificate

A certificate that `{x : V(x) ≤ ρ}` is inside the region of attraction, carrying
enough information to be re-verified without re-solving anything.

Returning the certificate rather than a bare `ρ` is the whole point of the
package: downstream code in `transition-viable-swarm` composes these objects,
and a number with no witness cannot be composed or audited.

# Fields
- `V` — the Lyapunov polynomial, **unpruned** (see [`prune_terms`](@ref)).
- `ρ::Float64` — certified level.
- `x::Vector` — the variables `V` is expressed in, in order.
- `multipliers::Dict{Symbol,Any}` — `:λ`, `:σ` (vector, one per inequality
  constraint), `:μ` (vector, per equality). The witnesses for (P2)/(P3).
- `grams::Dict{Symbol,Tuple{Vector,Matrix{Float64}}}` — `(basis, Q)` pairs for
  every SOS constraint in the program, keyed by the same symbols plus
  `:positivity` and `:derivative`. This is what [`verify_certificate`](@ref)
  checks; without the Gram matrices the certificate is an assertion, not a proof.
- `ε::Float64`
- `system_name::String`
- `solver_status::String` — the raw `termination_status` string.
- `solve_time::Float64` — seconds, summed over every SDP that produced this.
- `sdp_size::Int` — largest Gram dimension, for the benchmark table.
- `bounded::Bool`, `boundedness_evidence::Symbol` — from
  [`is_sublevel_set_bounded`](@ref). A certificate with `bounded = false` is not
  a certificate; the constructors refuse to return one.
"""
struct LyapunovCertificate
    V::Any
    ρ::Float64
    x::Vector
    multipliers::Dict{Symbol,Any}
    grams::Dict{Symbol,Tuple{Vector,Matrix{Float64}}}
    ε::Float64
    system_name::String
    solver_status::String
    solve_time::Float64
    sdp_size::Int
    bounded::Bool
    boundedness_evidence::Symbol
end

"""
    default_solver()

The conic solver used when `SOSOptions.solver` is `nothing`.

Clarabel is preferred: pure Julia, MIT-licensed, no binary artifact, and stable
on the ill-conditioned SDPs that degree-6 Lyapunov candidates produce. SCS is the
documented fallback. Mosek is *not* used from Julia — it appears only in
`matlab/`, which is where the licence question belongs.

Return an optimizer factory suitable for `SOSModel(...)`. Set solver tolerances
here, once, so every program in the package is solved to the same accuracy;
tolerances scattered across call sites are how two "identical" runs end up
disagreeing in the third digit.
"""
function default_solver()
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.3")
end

"""
    is_sos(p, x; solver = nothing, ε = 0.0)

Is the polynomial `p` a sum of squares? Returns `Bool`.

The elementary entry point, and the one to implement and test first: it is the
smallest program that exercises the whole SumOfSquares/JuMP/solver stack, so
when it works, later infeasibilities are about the mathematics rather than the
installation.

Must return `false` — not throw, not hang — for a polynomial that is nonnegative
but not SOS. The Motzkin polynomial `x²y⁴ + x⁴y² - 3x²y²z² + z⁶` is the standard
witness and is in the test suite for exactly this reason.
"""
function is_sos(p, x::AbstractVector; solver = nothing, ε::Real = 0.0)
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.3")
end

"""
    sos_decomposition(p, x; solver = nothing)

Return `(basis, Q, status)` such that `p ≈ basis' * Q * basis` with `Q ⪰ 0`.

`Q` is symmetrised and projected to its PSD part before return (clip eigenvalues
in `[-tol, 0)` to 0 and record the largest clipped magnitude in the log when
`verbose`). Do **not** silently clip a genuinely negative eigenvalue — if the
most negative eigenvalue exceeds `1e-7` in magnitude, the decomposition failed
and the function must say so through `status`, not return a matrix that looks fine.

The caller can turn `(basis, Q)` into an explicit sum of squares via a Cholesky
of `Q`; a worked instance of that appears in `docs/derivations/sos_decomposition.tex`.
"""
function sos_decomposition(p, x::AbstractVector; solver = nothing)
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.3")
end

"""
    lyapunov_from_linearization(sys; Q = I)

Quadratic seed `V(x) = x' P x` where `P` solves the continuous Lyapunov equation
`A'P + PA = -Q` for `A = ∂f/∂x(0)`.

This is the starting candidate for every ROA computation in the package. It is
guaranteed to be a valid Lyapunov function on *some* neighbourhood of the origin
whenever `A` is Hurwitz, which is what makes the subsequent ρ-bisection
well-posed: there is always a feasible `ρ > 0` to find.

Raise `ArgumentError`, with the offending eigenvalues in the message, if `A` is
not Hurwitz. That situation is not a numerical difficulty to be worked around —
it means the origin is not locally asymptotically stable and there is nothing to
certify. Solve the Lyapunov equation directly (Bartels–Stewart via
`LinearAlgebra`, or a Kronecker-product solve for the small `n` here); symmetrise
`P` afterwards, and verify the residual `‖A'P + PA + Q‖_F < 1e-9·max(1, ‖P‖_F)`.
"""
function lyapunov_from_linearization(sys::PolynomialSystem; Q = nothing)
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.3")
end

"""
    verify_sublevel_set(sys, V, ρ; opts = SOSOptions())

Solve (P2) — or (P3) when `sys` carries constraints or invariants — for the
given fixed `V` and `ρ`.

Returns `(feasible::Bool, witness::NamedTuple)` where `witness` holds the
multiplier polynomials, the Gram pairs, the raw solver status, and the solve
time. This is the primitive the bisection calls; it must be cheap to call
repeatedly and must never mutate `opts`.

Feasibility mapping, and it must be exactly this:

| `termination_status` | result |
|---|---|
| `OPTIMAL` | feasible |
| `ALMOST_OPTIMAL` / `SLOW_PROGRESS` | feasible **only if** the returned Gram matrices pass [`gram_residual`](@ref) and a minimum-eigenvalue check |
| `INFEASIBLE`, `DUAL_INFEASIBLE`, `INFEASIBLE_OR_UNBOUNDED` | infeasible |
| anything else | throw |

The middle row is where correctness is won or lost. A conic solver reporting
`ALMOST_OPTIMAL` on a badly scaled SOS program is common and frequently
meaningless; accepting it blindly is how a toolbox certifies a `ρ` that is too
large. The `throw` on the last row is deliberate — an unrecognised status is a
bug or a solver failure, and treating it as "infeasible" would silently shrink
every reported ROA.
"""
function verify_sublevel_set(sys::PolynomialSystem, V, ρ::Real; opts::SOSOptions = SOSOptions())
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.3")
end

"""
    v_step(sys, λ, shape, β; opts = SOSOptions())

The V-step of the V–s iteration: (P4) with the multiplier `λ` fixed, `ρ ≡ 1`, and
the shape-set level `β` fixed, searching over the coefficients of `V` and the
containment multiplier `σ_β`.

Returns `(feasible, V_new, witness)`.

Two constraints exist only to keep the problem non-degenerate, and both must be
present:

1. `V(0) = 0` as a linear equality on the coefficients — not merely
   `mindegree = 1` on the basis, because the solver can and will place mass on
   linear terms otherwise, giving a `V` with a nonzero gradient at the origin.
2. A normalisation. With `ρ ≡ 1` the scaling freedom is already fixed, but the
   *shape* freedom is not: without the containment constraint the V-step can
   inflate `V` in directions where `V̇` is comfortably negative and gain nothing.
   That is why the objective is `max β` on the shape set rather than anything
   defined in terms of `V` alone.
"""
function v_step(sys::PolynomialSystem, λ, shape, β::Real; opts::SOSOptions = SOSOptions())
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.3")
end

"""
    multiplier_step(sys, V, ρ; opts = SOSOptions())

The λ-step: `V` and `ρ` fixed, solve (P2)/(P3) for the multipliers alone.

Identical program to [`verify_sublevel_set`](@ref) — this is a named alias used
by [`v_s_iteration`](@ref) for readability. Implement one and have the other call
it; two copies of an S-procedure that must stay in sync is a guaranteed bug.
"""
function multiplier_step(sys::PolynomialSystem, V, ρ::Real; opts::SOSOptions = SOSOptions())
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.3")
end

"""
    beta_step(sys, V, shape; opts = SOSOptions())

Largest `β` with `{p(x) ≤ β} ⊆ {V(x) ≤ 1}`, via

```
(1 - V) - σ_β·(β - p)  is SOS,   σ_β  SOS.
```

Bilinear in `(β, σ_β)`, so bisect on `β` exactly as [`bisect_level_set`](@ref)
bisects on `ρ`. `shape` is normally `p(x) = x'x`, but any positive definite
polynomial works; when the ROA is strongly anisotropic (the unicycle) a shaped
`p` matched to the geometry recovers substantially more volume, and the example
is expected to show that comparison rather than assert it.
"""
function beta_step(sys::PolynomialSystem, V, shape; opts::SOSOptions = SOSOptions())
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.3")
end

"""
    synthesize_lyapunov(sys; opts = SOSOptions(), shape = nothing, V0 = nothing)

Top-level synthesis: seed from the linearisation (or `V0`), then run
[`v_s_iteration`](@ref) to convergence. Returns a [`LyapunovCertificate`](@ref).

Refuses to return a certificate whose sublevel set is not bounded — throws with
the boundedness evidence in the message.
"""
function synthesize_lyapunov(sys::PolynomialSystem; opts::SOSOptions = SOSOptions(),
                             shape = nothing, V0 = nothing)
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.3")
end

"""
    verify_certificate(cert::LyapunovCertificate; tol = 1e-7, samples = 50_000, rng = Random.default_rng())

Independent re-verification of a certificate. **Must not call a solver.**

Three checks, reported individually in the returned `NamedTuple` so a partial
failure is visible rather than collapsed into one `false`:

1. **Algebraic.** For each stored `(basis, Q)`: `gram_residual < tol` and
   `minimum(eigvals(Symmetric(Q))) > -tol`. This re-derives the SOS identity
   from the witness, independently of whatever the solver believed.
2. **Numerical sampling.** Draw points uniformly in `{V ≤ ρ}` (rejection-sample
   from a bounding box) and assert `V̇(x) < 0` for every one with `‖x‖ > 1e-6`.
   Report the worst `V̇` found and where.
3. **Dynamical.** Integrate `sys` from a sample of interior points and assert the
   trajectory neither leaves `{V ≤ ρ}` nor fails to converge to the origin.

Check 1 is the proof; checks 2 and 3 catch the case where the program that was
solved is not the program that was intended — a sign error in
[`lie_derivative`](@ref), the wrong `f`, a stale `V`. Check 1 alone would pass
happily in all three of those cases, which is precisely why all three are here.
"""
function verify_certificate(cert::LyapunovCertificate; tol::Real = 1e-7, samples::Integer = 50_000,
                            rng::Random.AbstractRNG = Random.default_rng())
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.3")
end

"""
    Base.show(io::IO, ::MIME"text/plain", cert::LyapunovCertificate)

Human-readable summary: system name, `deg V`, `ρ`, `ε`, largest SDP block,
solve time, solver status, and the boundedness evidence. Print the pruned `V`
when it has at most twelve terms and the degree and term count otherwise.

Worth writing properly — this is what appears in every notebook cell and every
console session, and a certificate that prints its own provenance is much harder
to misquote later.
"""
function Base.show(io::IO, ::MIME"text/plain", cert::LyapunovCertificate)
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.3")
end
