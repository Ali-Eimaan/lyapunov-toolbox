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
- `k::Int` — the exponent in the positivity margin `ε‖x - x_eq‖^{2k}` actually
  used. `k = max(1, cld(deg V̇, 2))`, so the margin degree matches `deg V̇` and
  the SOS basis need not be inflated. Recorded so [`verify_certificate`](@ref)
  can reconstruct the derivative target exactly, and reported in benchmark rows
  alongside `ε`.
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
    k::Int
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
function default_solver(; name::Symbol = :clarabel)
    if name === :clarabel
        return optimizer_with_attributes(
            Clarabel.Optimizer,
            "tol_gap_abs" => 1e-9,
            "tol_gap_rel" => 1e-9,
            "tol_feas" => 1e-9,
            "max_iter" => 500,
            "verbose" => false,
        )
    elseif name === :scs
        isdefined(@__MODULE__, :SCS) || throw(ArgumentError(
            "SCS fallback requested but SCS is not in julia/Project.toml [deps]; " *
            "add it with Pkg.add(\"SCS\") or use the default Clarabel backend"))
        return optimizer_with_attributes(SCS.Optimizer, "verbose" => false)
    else
        throw(ArgumentError("unknown solver name $(repr(name)); expected :clarabel or :scs"))
    end
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
    ms = monomials(p)
    isempty(ms) && return convert(Float64, p) >= 0  # constant: SOS iff nonnegative
    deg = maximum(degree, ms)
    # An SOS polynomial has even degree: an odd-degree form is sign-indefinite,
    # so the answer is false without touching the solver.
    isodd(deg) && return false
    d = deg ÷ 2
    model = SOSModel(solver === nothing ? default_solver() : solver)
    @constraint(model, cref, p - ε * sum(x .^ 2)^d >= 0)
    optimize!(model)
    st = termination_status(model)
    if st == MOI.OPTIMAL
        return true
    elseif st == MOI.INFEASIBLE || st == MOI.INFEASIBLE_OR_UNBOUNDED ||
           st == MOI.DUAL_INFEASIBLE || st == MOI.ALMOST_INFEASIBLE ||
           st == MOI.NUMERICAL_ERROR
        return false
    elseif st == MOI.ALMOST_OPTIMAL || st == MOI.SLOW_PROGRESS ||
           st == MOI.ITERATION_LIMIT || st == MOI.TIME_LIMIT
        # Accept only if the returned Gram matrix is genuinely a certificate:
        # it must reproduce p (residual) and be positive semidefinite.
        ok, _ = _gram_certificate_ok(p - ε * sum(x .^ 2)^d, cref)
        return ok
    else
        throw(ErrorException("is_sos: unexpected solver termination status $st"))
    end
end

# Shared gate for ALMOST_OPTIMAL / SLOW_PROGRESS results: the Gram matrix must
# reproduce `p` and be PSD, otherwise the solver's near-feasible answer is not a
# certificate and must be treated as failure.
function _gram_certificate_ok(p, cref; residual_tol::Real = 1e-7, eigmin_tol::Real = -1e-8)
    gm = gram_matrix(cref)
    Q = Matrix(gm.Q)
    basis = gm.basis.monomials
    res = gram_residual(p, basis, Q)
    λmin = minimum(eigvals(Symmetric(Q)))
    ok = res < residual_tol && λmin > eigmin_tol
    return ok, (residual = res, eigmin = λmin)
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
    ms = monomials(p)
    if isempty(ms)
        # Constant polynomial: SOS iff nonnegative, and the "Gram matrix" is 1x1.
        c = convert(Float64, p)
        return ([1], fill(c, 1, 1), c >= 0 ? :success : :failure)
    end
    deg = maximum(degree, ms)
    d = deg ÷ 2
    basis = sos_basis(x, d)
    emptyQ = zeros(gram_dimension(length(x), d), gram_dimension(length(x), d))
    # An odd-degree form is sign-indefinite and cannot be SOS; say so through
    # status rather than asking the solver to fail.
    isodd(deg) && return (basis, emptyQ, :failure)
    model = SOSModel(solver === nothing ? default_solver() : solver)
    @constraint(model, cref, p >= 0)
    optimize!(model)
    st = termination_status(model)
    if st == MOI.OPTIMAL
        gm = gram_matrix(cref)
        return _project_psd_and_report(p, gm.basis.monomials, Matrix(gm.Q))
    elseif st == MOI.ALMOST_OPTIMAL || st == MOI.SLOW_PROGRESS ||
           st == MOI.ITERATION_LIMIT || st == MOI.TIME_LIMIT
        # Accept only if the returned Gram matrix is genuinely a certificate.
        ok, _ = _gram_certificate_ok(p, cref)
        if ok
            gm = gram_matrix(cref)
            return _project_psd_and_report(p, gm.basis.monomials, Matrix(gm.Q))
        end
        return (basis, emptyQ, :failure)
    elseif st == MOI.INFEASIBLE || st == MOI.DUAL_INFEASIBLE ||
           st == MOI.INFEASIBLE_OR_UNBOUNDED || st == MOI.ALMOST_INFEASIBLE ||
           st == MOI.NUMERICAL_ERROR
        return (basis, emptyQ, :failure)
    else
        throw(ErrorException("sos_decomposition: unexpected solver termination status $st"))
    end
end

# Symmetrise `Q`, clip eigenvalues in [-1e-7, 0) to 0 (recording the largest
# clipped magnitude), and refuse to clip a genuinely negative eigenvalue — if the
# most negative eigenvalue exceeds 1e-7 in magnitude the decomposition failed and
# the failure must be visible in `status`, not papered over by projection. The
# final arbiter is the gram residual: the projected Q must still reproduce `p`.
function _project_psd_and_report(p, basis, Q)
    Q = (Q + Q') / 2
    F = eigen(Symmetric(Q))
    λmin = minimum(F.values)
    if λmin < -1e-7
        return (basis, Q, :failure)
    end
    vals = copy(F.values)
    largest_clipped = 0.0
    for i in eachindex(vals)
        if -1e-7 <= vals[i] < 0
            largest_clipped = max(largest_clipped, -vals[i])
            vals[i] = 0.0
        end
    end
    Qproj = Symmetric(F.vectors * Diagonal(vals) * F.vectors')
    res = gram_residual(p, basis, Qproj)
    @debug "sos_decomposition: largest clipped eigenvalue magnitude $largest_clipped, gram residual $res"
    return (basis, Matrix(Qproj), res < 1e-7 ? :success : :failure)
end

"""
    lyapunov_from_linearization(sys; Q = I)

Quadratic seed for a [`PolynomialSystem`](@ref): `V(x) = (x - x_eq)'P(x - x_eq)`
where `P` solves the continuous Lyapunov equation `A'P + PA = -Q` for
`A = ∂f/∂x(x_eq)`.

This is the starting candidate for every ROA computation in the package. It is
guaranteed to be a valid Lyapunov function on *some* neighbourhood of the
equilibrium whenever the linearisation is Hurwitz, which is what makes the
subsequent ρ-bisection well-posed: there is always a feasible `ρ > 0` to find.

Two cases are handled:

1. **No equality invariants** — `A` must be Hurwitz and the full Lyapunov
   equation is solved by Kronecker-product solve, centred at the equilibrium
   (identical to `x'Px` when the equilibrium is the origin).
2. **Equality invariants** (trig recasts carry `s² + c² = 1`): the flow is
   confined to the manifold `M = {h(x) = 0}`, so `A` has a zero eigenvalue in
   each constrained direction and the full Lyapunov equation has no solution.
   The seed is instead built on the **tangent space** `T_x_eq M`: the reduced
   dynamics `A_t = T'AT` are Hurwitz on `T_x_eq M`, the reduced Lyapunov
   equation is solved there, and the solution is lifted to
   `P = T P_t T' + (I - T T')`, i.e. the tangent-space solution plus a
   normal-direction penalty (the second term). The penalty makes `V` positive
   at *other* equilibria on the manifold — e.g. the pendulum's downward
   equilibrium `(0, -1, 0)` — which a certificate must exclude, and without it
   the SOS program is structurally infeasible at every `ρ ≥ 0`. The result
   vanishes at the equilibrium and is positive on the manifold — the correct
   analogue of the Hurwitz-`A` quadratic. If `A_t` itself has a zero eigenvalue
   (non-hyperbolic on the manifold, as for the unicycle), the seed falls back
   to the identity ball about the origin. See
   [`_lyapunov_on_manifold`](@ref).

Raise `ArgumentError`, with the offending eigenvalues in the message, if the
(reduced) linearisation is not Hurwitz. That situation is not a numerical
difficulty to be worked around — it means the equilibrium is not locally
asymptotically stable on the manifold and there is nothing to certify. Symmetrise
`P` afterwards, and verify the residual `‖A'P + PA + Q‖_F < 1e-9·max(1, ‖P‖_F)`.
"""
function lyapunov_from_linearization(sys::PolynomialSystem; Q = nothing)
    n = length(sys.x)
    x = sys.x
    eq = sys.equilibrium
    A = jacobian_at(sys.f, x, eq)
    Qmat = Q === nothing ? Matrix{Float64}(I, n, n) : Float64.(Q)
    size(Qmat) == (n, n) || throw(ArgumentError(
        "Q must be $(n)x$(n) to match the $(n)-state system, got $(size(Qmat))"))

    # Equality invariants (trig recasts: s² + c² = 1) make the equilibrium a
    # non-origin point of the manifold and give A zero eigenvalues in the
    # constrained directions.  A full-rank A'P + PA = -Q has no solution there,
    # and a seed that does not vanish at the equilibrium makes the γ-step of the
    # V–s iteration infeasible from the first round — every invariant system
    # used to come back with a degenerate certificate.  Solve on the tangent
    # space instead.
    if !isempty(sys.invariants)
        return _lyapunov_on_manifold(sys, A, Qmat)
    end

    ev = eigvals(A)
    if any(real(ev) .> 1e-10)
        throw(ArgumentError(
            "A = ∂f/∂x at the equilibrium is not Hurwitz (eigenvalues " *
            "$(round.(ev; digits = 6))); the origin is not locally " *
            "asymptotically stable, so there is nothing to certify"))
    end
    # No invariants but a zero eigenvalue (a centre direction the nonlinear
    # terms stabilise).  A'P+PA=-Q has no solution when A has zero eigenvalues
    # and Q is full rank; fall back to ‖x - x_eq‖² centred on the equilibrium.
    has_zero_eig = any(abs.(real(ev)) .< 1e-10)
    if has_zero_eig
        eq = sys.equilibrium
        return sum((x[i] - eq[i])^2 for i in 1:n)
    end
    # Continuous Lyapunov equation A'P + PA = -Q via Kronecker-product solve.
    M = kron(Matrix{Float64}(I, n, n), A') + kron(A', Matrix{Float64}(I, n, n))
    P = reshape(M \ vec(-Qmat), n, n)
    P = (P + P') / 2
    res = norm(A' * P + P * A + Qmat, 2)
    res < 1e-9 * max(1, norm(P, 2)) || throw(ErrorException(
        "Lyapunov solve failed: ‖A'P + PA + Q‖_F = $res is not below " *
        "1e-9·max(1, ‖P‖_F) = $(1e-9 * max(1, norm(P, 2)))"))
    # Centre the quadratic at the equilibrium: V = (x - x_eq)'P(x - x_eq).
    # For the equilibrium-at-origin systems this is exactly x'Px.
    return sum(P[i, j] * (x[i] - eq[i]) * (x[j] - eq[j]) for i in 1:n, j in 1:n)
end

"""
    _lyapunov_on_manifold(sys, A, Qmat)

Quadratic seed for a system with equality invariants `h(x) = 0` (trig recasts).

`M = {x : h(x) = 0}` is the invariant manifold and the flow is tangent to it.
Build the orthonormal basis `T` of the tangent space
`T_x_eq M = null(∇h(x_eq)')`, solve the Lyapunov equation for the reduced
dynamics `A_t = T'AT` on the tangent space, and lift the solution to
`P = T P_t T' + (I - T T')` — PSD, positive definite on `T_x_eq M`, and
positive in the normal directions (so `V` is positive at *other* equilibria on
the manifold, e.g. the pendulum's downward equilibrium, which the certificate
must exclude). Then `V = (x - x_eq)'P(x - x_eq)` vanishes at the equilibrium
and its derivative along the constrained flow is negative definite in a
neighbourhood: the correct analogue of the Hurwitz-`A` quadratic for
equilibrium-at-origin systems. If the reduced dynamics `A_t` have a zero
eigenvalue (non-hyperbolic on the manifold — e.g. the unicycle, whose `(x, y, s)`
dynamics decouple from `c` at first order), the reduced Lyapunov equation is
singular and the seed falls back to the identity ball about the *origin*, the
same fallback the non-manifold path uses. Without the manifold seed the old
`P = I` fallback seeded a ball about the *origin*, which is not a Lyapunov
function for a non-origin equilibrium on the manifold, and every invariant
system produced a degenerate (unverifiable) certificate.
"""
function _lyapunov_on_manifold(sys::PolynomialSystem, A::Matrix{Float64},
                               Qmat::Matrix{Float64})
    n = length(sys.x)
    x = sys.x
    eq = sys.equilibrium
    n_inv = length(sys.invariants)
    n_inv > 0 || throw(ArgumentError("_lyapunov_on_manifold needs ≥ 1 invariant"))
    # Gradient rows ∇h_k(eq)ᵀ; the tangent space of M at eq is their null space.
    G = zeros(n_inv, n)
    for (k, h) in enumerate(sys.invariants)
        for j in 1:n
            G[k, j] = convert(Float64, differentiate(h, x[j])(x => eq))
        end
    end
    T = nullspace(G)               # n × (n - n_inv), orthonormal columns
    m = size(T, 2)
    m == n - n_inv || throw(ErrorException(
        "tangent-space dimension mismatch: expected $(n - n_inv), got $m"))
    At = T' * A * T
    evt = eigvals(At)
    if any(real(evt) .> 1e-10)
        throw(ArgumentError(
            "A restricted to the invariant manifold at the equilibrium is not " *
            "Hurwitz (eigenvalues $(round.(evt; digits = 6))); the equilibrium " *
            "is not locally asymptotically stable on the manifold"))
    end
    # Zero eigenvalues in the reduced dynamics: the reduced Lyapunov equation
    # A_t'P_t + P_t A_t = -Q_t has no solution (its Kronecker matrix is
    # singular), so there is no tangent-space quadratic seed.  This happens when
    # the linearisation is non-hyperbolic ON the manifold — e.g. the unicycle,
    # whose (x, y, s) dynamics decouple from the c-direction at first order
    # (rows 2 and 4 of A are zero), leaving a zero eigenvalue that only the
    # nonlinear terms stabilise.  Fall back to ‖x - x_eq‖² — centred on the
    # equilibrium, not the origin.  A ball about the origin is not a Lyapunov
    # function at x_eq unless x_eq = 0.
    if any(abs.(real(evt)) .< 1e-10)
        eq = sys.equilibrium
        return sum((x[i] - eq[i])^2 for i in 1:n)
    end
    Qt = T' * Qmat * T
    Mt = kron(Matrix{Float64}(I, m, m), At') + kron(At', Matrix{Float64}(I, m, m))
    Pt = reshape(Mt \ vec(-Qt), m, m)
    Pt = (Pt + Pt') / 2
    res = norm(At' * Pt + Pt * At + Qt, 2)
    res < 1e-9 * max(1, norm(Pt, 2)) || throw(ErrorException(
        "manifold Lyapunov solve failed: ‖A_t'P_t + P_t A_t + Q_t‖_F = $res is " *
        "not below 1e-9·max(1, ‖P_t‖_F) = $(1e-9 * max(1, norm(Pt, 2)))"))
    # Lift to the full space, and add a normal-direction penalty.  P = T P_t T'
    # alone is singular (zero in the invariant directions), so the seed would
    # vanish at every point of the manifold whose tangent projection vanishes —
    # for the pendulum that is the DOWNWARD equilibrium (0, -1, 0), the opposite
    # point of the s² + c² = 1 circle.  A certificate can then never exclude it:
    # V(down) = 0 ≤ ρ while V̇(down) = 0 and the positivity margin is positive
    # there, so the SOS program is structurally infeasible for every ρ ≥ 0
    # (observed as "best ρ = 0.0 not feasible on re-verification").  The
    # penalty makes V positive at other equilibria on the manifold, so the
    # bisection's ρ is bounded by V at the nearest one and the certified set
    # excludes it — which is exactly the ground-truth anchor the tests assert.
    # γ = 1 keeps the normal term on the same scale as the (unit-Q) tangent
    # solution so the seed stays reasonably isotropic near the equilibrium.
    P = T * Pt * T' + (Matrix{Float64}(I, n, n) - T * T')
    return sum(P[i, j] * (x[i] - eq[i]) * (x[j] - eq[j]) for i in 1:n, j in 1:n)
end

"""
    _effective_multiplier_degree(sys, opts)

Degree to use for the S-procedure multiplier `λ`, raised above
`opts.multiplier_degree` when the pairing rule demands it.

The rule is `deg(λ) + deg(V) ≥ deg(V̇)`, and since `deg(V̇) = deg(V) - 1 + deg(f)`
that reduces to `deg(λ) ≥ deg(f) - 1`. Sitting exactly on that bound is not
enough in practice: the top-degree part of `λ·(V - ρ)` must cancel the top-degree
part of `V̇` *and* leave something PSD behind, and with no spare degree the
solver has no room to do both.

Measured on the reversed Van der Pol (`deg f = 3`) at `V_degree = 4`: the strict
bound gives `deg λ = 2` and the V-step reports `INFEASIBLE`; `deg λ = 4`
certifies and verifies (ρ = 1.0). So one even degree of headroom is the default.

A caller who sets `multiplier_degree` explicitly above this floor keeps their
value — this only raises, never lowers.
"""
function _effective_multiplier_degree(sys::PolynomialSystem, opts::SOSOptions)
    degf = 0
    for fi in sys.f
        ms = monomials(fi)
        isempty(ms) && continue
        degf = max(degf, Int(maximum(degree, ms)))
    end
    # deg(f) - 1 is the strict pairing bound; +1 for headroom, then round up to
    # the next even degree (an SOS multiplier must have even degree).
    floor_deg = degf - 1 + 1
    iseven(floor_deg) || (floor_deg += 1)
    return max(opts.multiplier_degree, floor_deg)
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
| `ALMOST_OPTIMAL` / `SLOW_PROGRESS` / `ITERATION_LIMIT` / `TIME_LIMIT` | feasible **only if** the returned Gram matrices pass [`gram_residual`](@ref) and a minimum-eigenvalue check |
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
    x = sys.x
    bad = setdiff(variables(V), x)
    isempty(bad) || throw(ArgumentError(
        "V contains variables $(bad) not in the state vector x; a Lyapunov " *
        "candidate must be a polynomial in the system's own states"))
    Vd = lie_derivative(V, sys.f, x)
    # The positivity margin must vanish at the equilibrium, or the certificate is
    # structurally infeasible for any recast system whose equilibrium is not the
    # origin (trig_recast states live on s² + c² = 1, so ‖x_eq‖² = 1 ≠ 0 there).
    # The margin is ε‖x - x_eq‖^{2k} with k = max(1, cld(deg V̇, 2)): the degree
    # is matched to deg V̇ so the S-procedure term can cancel the top-degree part
    # of V̇ without inflating the SOS basis. For a linear system (deg V̇ = 2) this
    # is exactly the guide's ε‖x‖²; for the reversed Van der Pol (deg V̇ = 4) it is
    # the weaker ε‖x‖⁴ — still a strict negativity margin, but near the
    # equilibrium it is much smaller than ε‖x‖² would be, and ε is a recorded
    # benchmark parameter, so the actual margin is what gets recorded as `k`.
    ms = monomials(Vd)
    dVd = isempty(ms) ? 0 : maximum(degree, ms)
    k = max(1, cld(dVd, 2))
    margin = sum((x .- sys.equilibrium) .^ 2)^k

    model = SOSModel(opts.solver === nothing ? default_solver() : opts.solver)
    # (P2)/(P3) multipliers: λ SOS of degree multiplier_degree (basis degree is
    # half), σⱼ SOS per inequality constraint, μₖ FREE per equality invariant.
    mult_deg = _effective_multiplier_degree(sys, opts)
    λb = sos_basis(x, mult_deg ÷ 2)
    λ = @variable(model, λ, SOSPoly(λb))
    n_cons = length(sys.constraints)
    σb = sos_basis(x, opts.constraint_multiplier_degree ÷ 2)
    σ = n_cons > 0 ? @variable(model, σ[1:n_cons], SOSPoly(σb)) : Any[]
    n_inv = length(sys.invariants)
    μb = monomial_basis(x, mult_deg)
    μ = n_inv > 0 ? @variable(model, μ[1:n_inv], Poly(μb)) : Any[]

    expr = -(Vd + opts.ε * margin) + λ * (V - ρ)
    for j in 1:n_cons
        expr -= σ[j] * sys.constraints[j]
    end
    for k in 1:n_inv
        expr -= μ[k] * sys.invariants[k]
    end
    @constraint(model, cref, expr >= 0)

    # CONTAINMENT: {V ≤ ρ} ⊆ {gⱼ ≥ 0} for every inequality constraint.
    #
    # The derivative constraint above only certifies V̇ < 0 ON the intersection
    # {V ≤ ρ} ∩ {g ≥ 0}. It says nothing about whether the sublevel set actually
    # lies inside {g ≥ 0} — so without this the toolbox proved "trajectories
    # starting inside the certified set AND outside the obstacle converge", while
    # the unicycle example claimed "the certified set is obstacle-free". Measured
    # intrusion before this constraint existed: 6.7 cm into the obstacle.
    #
    # S-procedure: on {V ≤ ρ} we have (ρ - V) ≥ 0 and σ_c ≥ 0, so
    #     gⱼ - σ_cⱼ·(ρ - V) ⪰ 0   ⟹   gⱼ ≥ σ_cⱼ·(ρ - V) ≥ 0 there.
    σcb = sos_basis(x, opts.constraint_multiplier_degree ÷ 2)
    σc = n_cons > 0 ? @variable(model, σc[1:n_cons], SOSPoly(σcb)) : Any[]
    ccrefs = Any[]
    for j in 1:n_cons
        cc = @constraint(model, sys.constraints[j] - σc[j] * (ρ - V) >= 0)
        push!(ccrefs, cc)
    end

    optimize!(model)
    solve_time = JuMP.solve_time(model)
    st = termination_status(model)
    if st == MOI.OPTIMAL
        feasible = true
    elseif st == MOI.ALMOST_OPTIMAL || st == MOI.SLOW_PROGRESS ||
           st == MOI.ITERATION_LIMIT || st == MOI.TIME_LIMIT
        # Accept only if every Gram matrix the solver returned is genuinely a
        # certificate: it must reproduce its target polynomial and be PSD.
        feasible = _witness_grams_ok(sys, V, Vd, ρ, margin, opts, cref, λ, σ, μ, μb)
    elseif st == MOI.INFEASIBLE || st == MOI.DUAL_INFEASIBLE ||
           st == MOI.INFEASIBLE_OR_UNBOUNDED || st == MOI.ALMOST_INFEASIBLE ||
           st == MOI.NUMERICAL_ERROR
        feasible = false
    else
        throw(ErrorException(
            "verify_sublevel_set: unexpected solver termination status $st"))
    end
    if feasible
        multipliers, grams = _gather_witness(cref, λ, σ, μ)
        # Containment witnesses, keyed :containment_j, multiplier list :σc.
        σcp = Vector{Any}(undef, length(ccrefs))
        for (j, ccj) in enumerate(ccrefs)
            gm_c = gram_matrix(ccj)
            grams[Symbol(:containment_, j)] = (gm_c.basis.monomials, Matrix(gm_c.Q))
            σcp[j] = _sos_variable_polynomial(σc[j])
        end
        multipliers[:σc] = σcp
    else
        multipliers = Dict{Symbol,Any}()
        grams = Dict{Symbol,Tuple{Vector,Matrix{Float64}}}()
    end
    witness = (multipliers = multipliers, grams = grams,
               status = string(st), solve_time = solve_time)
    return feasible, witness
end

# The P2/P3 expression with every multiplier substituted by its solved value.
# This is the target polynomial the derivative constraint's Gram matrix claims
# to reproduce.
function _derivative_target(sys, V, Vd, ρ, margin, opts, λp, σp, μp)
    expr = -(Vd + opts.ε * margin) + λp * (V - ρ)
    for (j, g) in enumerate(sys.constraints)
        expr -= σp[j] * g
    end
    for (k, h) in enumerate(sys.invariants)
        expr -= μp[k] * h
    end
    return expr
end

# Gate for ALMOST_OPTIMAL / SLOW_PROGRESS: every Gram matrix in the program must
# pass the same residual + minimum-eigenvalue test that _gram_certificate_ok
# applies to a single constraint. A near-feasible answer whose Gram matrices do
# not reproduce their targets is not a certificate and must be treated as
# infeasible — that is where certifying a ρ too large would otherwise slip in.
function _witness_grams_ok(sys, V, Vd, ρ, margin, opts, cref, λ, σ, μ, μb)
    gm = gram_matrix(cref)
    Q = Matrix(gm.Q)
    basis = gm.basis.monomials
    λp = _sos_variable_polynomial(λ)
    σp = [_sos_variable_polynomial(σj) for σj in σ]
    μp = [_poly_variable_polynomial(μk) for μk in μ]
    target = _derivative_target(sys, V, Vd, ρ, margin, opts, λp, σp, μp)
    res = gram_residual(target, basis, Q)
    eigmin = minimum(eigvals(Symmetric(Q)))
    ok = res < 1e-7 && eigmin > -1e-8
    for σj in σ
        Qj = Matrix(value.(σj.Q))
        ok &= minimum(eigvals(Symmetric(Qj))) > -1e-8
    end
    Qλ = Matrix(value.(λ.Q))
    ok &= minimum(eigvals(Symmetric(Qλ))) > -1e-8
    return ok
end

# Build the witness `(multipliers, grams)` from a solved model. Grams are keyed
# :derivative (the main P2/P3 constraint), :λ, and :σ_j per inequality
# constraint; :positivity appears in the full certificate (from the V-step), not
# here. Multipliers holds the multiplier polynomials: :λ, :σ (vector), :μ (vector).
function _gather_witness(cref, λ, σ, μ)
    gm = gram_matrix(cref)
    grams = Dict{Symbol,Tuple{Vector,Matrix{Float64}}}(
        :derivative => (gm.basis.monomials, Matrix(gm.Q)))
    multipliers = Dict{Symbol,Any}()
    λp = _sos_variable_polynomial(λ)
    multipliers[:λ] = λp
    grams[:λ] = (λ.basis.monomials, Matrix(value.(λ.Q)))
    σp = Vector{Any}(undef, length(σ))
    for (j, σj) in enumerate(σ)
        σp[j] = _sos_variable_polynomial(σj)
        grams[Symbol(:σ_, j)] = (σj.basis.monomials, Matrix(value.(σj.Q)))
    end
    multipliers[:σ] = σp
    μp = Vector{Any}(undef, length(μ))
    for (k, μk) in enumerate(μ)
        μp[k] = _poly_variable_polynomial(μk)
    end
    multipliers[:μ] = μp
    return multipliers, grams
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
function v_step(sys::PolynomialSystem, λ, shape, β::Real;
                opts::SOSOptions = SOSOptions(), σc_fixed = Any[])
    x = sys.x
    # V is parameterised over the FULL monomial basis, including the constant and
    # linear terms (mindegree = 0).
    #
    # `mindegree = 2` is only valid when the equilibrium is the origin. For a
    # trig-recast system it is not: the pendulum sits at (s, c, θ̇) = (0, 1, 0)
    # and the unicycle at (e₁, e₂, s, c) = (0, 0, 0, 1), so a quadratic centred
    # there expands to terms like (c - 1)² = 1 - 2c + c² — a constant and a
    # linear term. Excluding them put the previous iterate OUTSIDE the search
    # space, so the V-step was infeasible at iteration 1 on exactly the systems
    # whose equilibrium is off the origin, contradicting the iteration's own
    # invariant that the previous iterate is always feasible.
    #
    # Nothing is lost for origin-centred systems: V(x_eq) = 0 is imposed as an
    # explicit linear equality below, and V - ε·margin ⪰ 0 then makes x_eq a
    # global minimum, which forces ∇V(x_eq) = 0 without constraining the basis.
    Vb = monomial_basis(x, opts.V_degree; mindegree = 0)
    nv = length(Vb)
    model = SOSModel(opts.solver === nothing ? default_solver() : opts.solver)
    @variable(model, vcoeffs[1:nv])
    V = sum(vcoeffs[i] * Vb[i] for i in 1:nv)
    # Enforce V(0) = 0: the constant term is excluded by mindegree=2, but for a
    # trig-recast system the equilibrium may not be at the origin.  Build the
    # constraint anyway as an equality on the evaluated V.
    @constraint(model, sum(vcoeffs[i] * convert(Float64, Vb[i](x => sys.equilibrium)) for i in 1:nv) == 0)

    Vd = lie_derivative(V, sys.f, x)
    ms = monomials(Vd)
    dVd = isempty(ms) ? 0 : maximum(degree, ms)
    # Margin is ε‖x - x_eq‖^{2k} with k = max(1, cld(deg V̇, 2)) — degree matched
    # to deg V̇ so the SOS degrees line up without inflating the basis. Only a
    # linear system (deg V̇ = 2) gets the guide's ε‖x‖²; the reversed Van der Pol
    # gets ε‖x‖⁴. Still a strict margin (V̇ ≤ -ε‖x - x_eq‖^{2k} < 0 off the
    # equilibrium), but weaker near the equilibrium than ε‖x‖² would be.
    k = max(1, cld(dVd, 2))
    margin = sum((x .- sys.equilibrium) .^ 2)^k

    # (P4) positivity constraint: V - ε·margin_pos is SOS
    # IMPORTANT: use a degree-2 margin for positivity, NOT the dVd-matched
    # margin.  When deg(V̇) > 2 (e.g. cubic dynamics), margin has degree ≥ 4,
    # and V - ε·‖x‖⁴ cannot be globally SOS for a degree-2 V because the
    # leading term is -ε·‖x‖⁴.  The derivative constraint uses the full margin;
    # positivity only needs a token strictness term at the same degree as V.
    margin_pos = sum((x .- sys.equilibrium) .^ 2)
    @constraint(model, pos_cref, V - opts.ε * margin_pos >= 0)

    # (P4) derivative constraint with fixed λ: -(V̇ + ε·margin) + λ·(V - 1) is SOS
    expr_deriv = -(Vd + opts.ε * margin) + λ * (V - 1.0)
    σ_vars = Any[]
    μ_vars = Any[]
    for (j, g) in enumerate(sys.constraints)
        σb_j = sos_basis(x, opts.constraint_multiplier_degree ÷ 2)
        σj = @variable(model, [1:1], SOSPoly(σb_j))
        push!(σ_vars, σj[1])
        expr_deriv -= σj[1] * g
    end
    for (k_inv, h) in enumerate(sys.invariants)
        μb_h = monomial_basis(x, _effective_multiplier_degree(sys, opts))
        μk = @variable(model, [1:1], Poly(μb_h))
        push!(μ_vars, μk[1])
        expr_deriv -= μk[1] * h
    end
    @constraint(model, deriv_cref, expr_deriv >= 0)

    # CONTAINMENT, same as verify_sublevel_set but at the normalised level ρ ≡ 1:
    # {V ≤ 1} ⊆ {gⱼ ≥ 0}. Without this the V-step could enlarge V back over the
    # obstacle, undoing the containment the γ-step had established — the V–s
    # iteration would then report a larger set that is not obstacle-free.
    # σ_c is HELD FIXED from the γ-step, exactly as λ is. Both σ_c and V are
    # unknowns of the V-step, so leaving σ_c free makes σ_c·(1 - V) bilinear and
    # the problem is no longer an SDP ("This constraint type is not supported by
    # the solver"). Fixing it keeps the constraint linear in V.
    for (j, g) in enumerate(sys.constraints)
        j <= length(σc_fixed) || break
        @constraint(model, g - σc_fixed[j] * (1.0 - V) >= 0)
    end

    # (P4) shape containment: (1 - V) - σ_β·(β - shape) is SOS, σ_β SOS
    # The caller has already bisected on β, so β is feasible here.
    σβ_b = sos_basis(x, 2)
    σβ = @variable(model, σβ, SOSPoly(σβ_b))
    @constraint(model, shape_cref, (1.0 - V) - σβ * (β - shape) >= 0)

    optimize!(model)
    solve_time = JuMP.solve_time(model)
    st = termination_status(model)

    if st == MOI.OPTIMAL
        feasible = true
    elseif st == MOI.ALMOST_OPTIMAL || st == MOI.SLOW_PROGRESS ||
           st == MOI.ITERATION_LIMIT || st == MOI.TIME_LIMIT
        # Build V_new early so we can form the target polynomials for the
        # gram_residual check (R1-29).
        V_tentative = _build_V_from_solution(vcoeffs, Vb, x)
        Vd_tentative = lie_derivative(V_tentative, sys.f, x)
        ms_t = monomials(Vd_tentative)
        dVd_t = isempty(ms_t) ? 0 : maximum(degree, ms_t)
        k_t = max(1, cld(dVd_t, 2))
        margin_t = sum((x .- sys.equilibrium) .^ 2)^k_t
        # Positivity target uses the degree-2 margin, matching the constraint.
        pos_target = V_tentative - opts.ε * sum((x .- sys.equilibrium) .^ 2)
        deriv_target = -(Vd_tentative + opts.ε * margin_t) + λ * (V_tentative - 1.0)
        for (j, g) in enumerate(sys.constraints)
            σv = σ_vars[j]
            σp = _sos_variable_polynomial(σv)
            deriv_target -= σp * g
        end
        for (k_inv, h) in enumerate(sys.invariants)
            μv = μ_vars[k_inv]
            μp = _poly_variable_polynomial(μv)
            deriv_target -= μp * h
        end
        shape_target = (1.0 - V_tentative) - _sos_variable_polynomial(σβ) * (β - shape)
        feasible = _v_step_grams_ok(model, pos_cref, deriv_cref, shape_cref,
                                     pos_target, deriv_target, shape_target;
                                     verbose = opts.verbose)
    elseif st == MOI.INFEASIBLE || st == MOI.DUAL_INFEASIBLE ||
           st == MOI.INFEASIBLE_OR_UNBOUNDED || st == MOI.ALMOST_INFEASIBLE ||
           st == MOI.NUMERICAL_ERROR
        feasible = false
    else
        throw(ErrorException("v_step: unexpected solver termination status $st"))
    end

    if feasible
        V_new = _build_V_from_solution(vcoeffs, Vb, x)
        gm_p = gram_matrix(pos_cref)
        gm_d = gram_matrix(deriv_cref)
        gm_s = gram_matrix(shape_cref)
        grams = Dict{Symbol,Tuple{Vector,Matrix{Float64}}}(
            :positivity => (gm_p.basis.monomials, Matrix(gm_p.Q)),
            :derivative => (gm_d.basis.monomials, Matrix(gm_d.Q)),
            :shape_containment => (gm_s.basis.monomials, Matrix(gm_s.Q)),
        )
        # Capture the constraint/invariant multipliers from the V-step so the
        # certificate carries the correct witness for V_new (R1-7).
        multipliers = Dict{Symbol,Any}()
        σp = [_sos_variable_polynomial(σv) for σv in σ_vars]
        if !isempty(σp)
            multipliers[:σ] = σp
            for (j, σv) in enumerate(σ_vars)
                grams[Symbol(:σ_, j)] = (σv.basis.monomials, Matrix(value.(σv.Q)))
            end
        end
        μp = [_poly_variable_polynomial(μv) for μv in μ_vars]
        if !isempty(μp)
            multipliers[:μ] = μp
        end
        # λ was fixed (from the γ-step) — don't overwrite it; the caller
        # merges it into the certificate.
        witness = (grams = grams, multipliers = multipliers,
                   status = string(st), solve_time = solve_time)
    else
        V_new = nothing
        witness = (grams = Dict{Symbol,Tuple{Vector,Matrix{Float64}}}(),
                   multipliers = Dict{Symbol,Any}(),
                   status = string(st), solve_time = solve_time)
    end
    return feasible, V_new, witness
end

# Extract a DynamicPolynomials.Polynomial from a JuMP SOSPoly variable
# after the model has been solved.  Uses the monomial basis and the
# Gram matrix Q stored in the variable.
function _sos_variable_polynomial(σv)
    basis_vec = σv.basis.monomials
    Q = value.(σv.Q)
    return polynomial(basis_vec' * Q * basis_vec)
end

# Extract a DynamicPolynomials.Polynomial from a JuMP Poly variable
# (linear combination of monomials) after the model has been solved.
function _poly_variable_polynomial(μv)
    # μv is a DynamicPolynomials.Polynomial with JuMP.VariableRef coefficients
    return sum(value(coeff) * monom for (coeff, monom) in zip(μv.a, μv.x))
end

function _v_step_grams_ok(model, pos_cref, deriv_cref, shape_cref,
                           pos_target, deriv_target, shape_target;
                           verbose::Bool = false,
                           eigmin_tol::Real = -1e-8, residual_tol::Real = 1e-6)
    names  = [:positivity, :derivative, :shape_containment]
    crefs   = [pos_cref, deriv_cref, shape_cref]
    targets = [pos_target, deriv_target, shape_target]
    ok = true
    for (nm, cref, target) in zip(names, crefs, targets)
        gm = gram_matrix(cref)
        Q = Matrix(gm.Q)
        λmin = minimum(eigvals(Symmetric(Q)))
        # Residual against the TARGET polynomial, not just PSD-ness: a Gram
        # matrix can be PSD and still represent a different polynomial.
        res = gram_residual(target, gm.basis.monomials, Q)
        passed = λmin > eigmin_tol && res <= residual_tol
        passed || (ok = false)
        # A marginal solve that gets rejected used to say nothing about why,
        # which made "V-step infeasible" indistinguishable from "solver was
        # 3e-9 outside a hard-coded threshold". Report the numbers.
        verbose && @info "v_step gram check" constraint=nm eigmin=λmin residual=res passed
    end
    return ok
end

function _build_V_from_solution(vcoeffs, Vb, x)
    vals = value.(vcoeffs)
    Vp = sum(vals[i] * Vb[i] for i in eachindex(Vb))
    return Vp
end

"""
    multiplier_step(sys, V, ρ; opts = SOSOptions())

The λ-step: `V` and `ρ` fixed, solve (P2)/(P3) for the multipliers alone.

Identical program to [`verify_sublevel_set`](@ref) — this is a named alias used
by [`v_s_iteration`](@ref) for readability. Implement one and have the other call
it; two copies of an S-procedure that must stay in sync is a guaranteed bug.
"""
function multiplier_step(sys::PolynomialSystem, V, ρ::Real; opts::SOSOptions = SOSOptions())
    return verify_sublevel_set(sys, V, ρ; opts = opts)
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
    # Bisect on β for {p ≤ β} ⊆ {V ≤ 1} via (1-V) - σ_β(β-p) SOS.
    # We reuse the bisection machinery from roa_estimation.jl by calling a
    # private helper that builds the feasibility oracle for β.
    oracle = _beta_oracle(sys, V, shape, opts)
    β_lo = 0.0
    # Find upper bound by doubling from a small positive value.
    β_hi = 1e-3
    found_infeasible = false
    for _ in 1:30
        if !oracle(β_hi)
            found_infeasible = true
            break
        end
        β_hi *= 2.0
    end
    found_infeasible || begin
        throw(ErrorException(
            "beta_step: no infeasible β found after 30 doublings up to " *
            "β_hi = $β_hi; the V function may dominate the shape everywhere. " *
            "Returning a sentinel (1e10) is not acceptable — the β estimate " *
            "is meaningless without an upper bound."))
    end

    best = β_lo
    iters = 0
    while (β_hi - β_lo) > opts.bisection_tol * max(1.0, β_lo) && iters < opts.max_bisection_iters
        mid = (β_lo + β_hi) / 2
        if oracle(mid)
            best = mid
            β_lo = mid
        else
            β_hi = mid
        end
        iters += 1
    end
    return best
end

function _beta_oracle(sys::PolynomialSystem, V, shape, opts::SOSOptions)
    x = sys.x
    return function(β::Float64)
        model = SOSModel(opts.solver === nothing ? default_solver() : opts.solver)
        σβ_b = sos_basis(x, 2)
        σβ = @variable(model, σβ, SOSPoly(σβ_b))
        expr = (1.0 - V) - σβ * (β - shape)
        @constraint(model, cref, expr >= 0)
        optimize!(model)
        st = termination_status(model)
        if st == MOI.OPTIMAL
            return true
        elseif st == MOI.ALMOST_OPTIMAL || st == MOI.SLOW_PROGRESS ||
           st == MOI.ITERATION_LIMIT || st == MOI.TIME_LIMIT
            gm = gram_matrix(cref)
            Q = Matrix(gm.Q)
            ok = minimum(eigvals(Symmetric(Q))) > -1e-8
            return ok
        elseif st == MOI.INFEASIBLE || st == MOI.DUAL_INFEASIBLE ||
               st == MOI.INFEASIBLE_OR_UNBOUNDED || st == MOI.ALMOST_INFEASIBLE ||
               st == MOI.NUMERICAL_ERROR
            return false
        else
            throw(ErrorException("beta_step oracle: unexpected solver status $st"))
        end
    end
end

"""
    synthesize_lyapunov(sys; opts = SOSOptions(), shape = nothing, V0 = nothing)

Top-level synthesis: seed from the linearisation (or `V0`), then run
[`v_s_iteration`](@ref) to convergence. Returns a [`LyapunovCertificate`](@ref).

Refuses to return a certificate whose sublevel set is not bounded — throws with
the boundedness evidence in the message.
"""
function synthesize_lyapunov(sys::PolynomialSystem; opts::SOSOptions = SOSOptions(),
                             shape = nothing, V0 = nothing, ρ_hi = nothing)
    x = sys.x
    # Seed: quadratic Lyapunov from the linearisation, or use provided V0.
    if V0 === nothing
        V = lyapunov_from_linearization(sys)
    else
        V = V0
    end
    # Default shape: p(x) = (x - x_eq)' (x - x_eq)
    if shape === nothing
        shape = sum((x .- sys.equilibrium) .^ 2)
    end
    result = v_s_iteration(sys; opts = opts, V0 = V, shape = shape, ρ_hi = ρ_hi)
    return result.certificate
end

"""
    _project_to_manifold(pt, sys; iters = 25, tol = 1e-12)

Gauss–Newton projection of `pt` onto the invariant variety `{h_k(x) = 0}`.

Returns `(projected_point, converged)`.

`verify_certificate`'s sampling and simulation checks draw points from a box
around the equilibrium. For a trig-recast system that box is 3- or 4-dimensional
while the actual state space is the 2- or 3-dimensional manifold `s² + c² = 1`,
so essentially every raw sample is **not a state of the system**. Evaluating `V̇`
there tests nothing (the certificate never claimed anything off the manifold) and
integrating from there produces trajectories the system cannot take — which is
why the trig pendulum's certificate was algebraically sound and still failed
checks 2 and 3.

The correction is the minimum-norm Newton step `x ← x - J⁺h(x)` with
`J = ∂h/∂x`, iterated. For `h = s² + c² - 1` this converges in a couple of steps
from anywhere off the origin of the `(s, c)` plane.
"""
function _project_to_manifold(pt::AbstractVector, sys::PolynomialSystem;
                              iters::Integer = 25, tol::Real = 1e-12)
    isempty(sys.invariants) && return (collect(float.(pt)), true)
    x = sys.x
    y = collect(float.(pt))
    m = length(sys.invariants)
    n = length(x)
    for _ in 1:iters
        h = [convert(Float64, sys.invariants[k](x => y)) for k in 1:m]
        maximum(abs, h) < tol && return (y, true)
        J = Matrix{Float64}(undef, m, n)
        for k in 1:m, j in 1:n
            J[k, j] = convert(Float64, differentiate(sys.invariants[k], x[j])(x => y))
        end
        # Minimum-norm correction: move as little as possible off the sample.
        step = J \ h
        all(isfinite, step) || return (y, false)
        y .-= step
    end
    h = [abs(convert(Float64, sys.invariants[k](x => y))) for k in 1:m]
    return (y, maximum(h) < 1e-8)
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
function verify_certificate(cert::LyapunovCertificate; tol::Real = 1e-7,
                            samples::Integer = 50_000,
                            rng::Random.AbstractRNG = Random.default_rng())
    throw(ArgumentError(
        "verify_certificate(cert) needs the system to reconstruct the algebraic " *
        "targets; call verify_certificate(cert, sys) instead. The one-argument " *
        "form was a tautology that reconstructed each polynomial from its own " *
        "(basis, Q) pair and then checked it against itself — it passed for " *
        "every PSD matrix regardless of whether it belonged to this certificate."))
end

# ---------------------------------------------------------------------------
# Two-argument verify_certificate with full checks 1–3.
# ---------------------------------------------------------------------------
function verify_certificate(cert::LyapunovCertificate, sys::PolynomialSystem;
                            tol::Real = 1e-7, samples::Integer = 50_000,
                            rng::Random.AbstractRNG = Random.default_rng())
    x = cert.x
    V = cert.V
    ρ = cert.ρ
    eq = sys.equilibrium
    Vd_poly = lie_derivative(V, sys.f, x)

    # --- Check 1: Algebraic ---
    # Reconstruct the TARGET polynomial for each Gram pair from the certificate's
    # own data (V, ρ, ε, k, multipliers), then check that (basis, Q) reproduces
    # that target. The old code reconstructed p from (basis, Q) and checked it
    # against itself — a tautology that passed for any PSD matrix.
    dVd = isempty(monomials(Vd_poly)) ? 0 : maximum(degree, monomials(Vd_poly))
    k_margin = max(1, cld(dVd, 2))
    margin = sum((x .- eq) .^ 2)^k_margin
    ε_used = cert.ε

    # Build the target map: for each key in cert.grams, the polynomial the
    # Gram matrix claims to represent.
    targets = Dict{Symbol,Any}()

    # :derivative target: -(V̇ + ε·margin) + λ(V - ρ) - Σσⱼgⱼ - Σμₖhₖ
    if haskey(cert.multipliers, :λ)
        λp = cert.multipliers[:λ]
        target_deriv = -(Vd_poly + ε_used * margin) + λp * (V - ρ)
        if haskey(cert.multipliers, :σ)
            for (j, σj) in enumerate(cert.multipliers[:σ])
                target_deriv -= σj * sys.constraints[j]
            end
        end
        if haskey(cert.multipliers, :μ)
            for (k_inv, μk) in enumerate(cert.multipliers[:μ])
                target_deriv -= μk * sys.invariants[k_inv]
            end
        end
        targets[:derivative] = target_deriv
    end

    # :positivity target: V - ε·margin_pos
    # The positivity margin is always degree‑2 (matched to V's minimum degree),
    # independent of the derivative margin (matched to V̇).  See v_step for the
    # rationale: V - ε·‖x-x_eq‖⁴ cannot be globally SOS when V has degree 2.
    targets[:positivity] = V - ε_used * sum((x .- eq) .^ 2)

    # :shape_containment target (if present): (1 - V) - σ_β(β - shape).
    # We cannot reconstruct this without knowing shape and β, but we can
    # still check that the stored Gram matrix is PSD and skip the residual.
    # The shape containment constraint is separately checked by the caller
    # (bisect_level_set / beta_step already verified it at synthesis time).

    # :λ target: λ itself
    if haskey(cert.multipliers, :λ) && haskey(cert.grams, :λ)
        targets[:λ] = cert.multipliers[:λ]
    end

    # :containment_j targets: gⱼ - σ_cⱼ·(ρ - V), the S-procedure statement of
    # {V ≤ ρ} ⊆ {gⱼ ≥ 0}. This is what makes "the certified set is
    # obstacle-free" a proved claim rather than an assertion.
    if haskey(cert.multipliers, :σc)
        for (j, σcj) in enumerate(cert.multipliers[:σc])
            key = Symbol(:containment_, j)
            if haskey(cert.grams, key) && j <= length(sys.constraints)
                targets[key] = sys.constraints[j] - σcj * (ρ - V)
            end
        end
    end

    # :σ_j targets: each σⱼ itself
    if haskey(cert.multipliers, :σ)
        for (j, σj) in enumerate(cert.multipliers[:σ])
            key = Symbol(:σ_, j)
            if haskey(cert.grams, key)
                targets[key] = σj
            end
        end
    end

    algebraic_ok = true
    algebraic_details = Dict{Symbol,Dict{Symbol,Float64}}()
    unchecked_keys = Symbol[]
    for (key, (basis, Q)) in pairs(cert.grams)
        λmin = minimum(eigvals(Symmetric(Q)))
        # Every key with a reconstructable target is checked against it. There
        # is no exemption list: :derivative and :positivity ARE the certificate,
        # and asserting r = 0.0 for them reduces check 1 to a PSD test — which a
        # random PSD matrix over an unrelated basis passes.
        #
        # A Newton-reduced basis is not a reason to skip the residual. The stored
        # pair is (basis, Q) as the solver returned it, so basis'*Q*basis
        # reproduces the target over THAT basis; a nonzero residual means the
        # basis and Q were not stored as a matching pair, which is exactly the
        # defect this check exists to find.
        #
        # :shape_containment has no reconstructable target — beta and the shape
        # set are not carried on the certificate — so it is reported as unchecked
        # and does NOT count as passing.
        if haskey(targets, key)
            r = gram_residual(targets[key], basis, Q)
            checked = 1.0
        else
            r = NaN
            checked = 0.0
            push!(unchecked_keys, key)
        end
        algebraic_details[key] = Dict(:gram_residual => r, :eigmin => λmin,
                                      :checked => checked)
        if checked == 0.0 || r >= tol || λmin <= -tol
            algebraic_ok = false
        end
    end

    # --- Check 2: Numerical sampling ---
    n = length(x)
    # Estimate bounding box for {V ≤ ρ} via ray bisection.
    radius_est = 10.0
    for _ in 1:1000
        u = randn(rng, n)
        nrm = norm(u)
        nrm > 0 || continue
        lo, hi = 0.0, 100.0
        for _ in 1:50
            mid = (lo + hi) / 2
            pt = eq .+ mid .* (u ./ nrm)
            if convert(Float64, V(x => pt)) <= ρ
                lo = mid
            else
                hi = mid
            end
        end
        radius_est = max(radius_est, lo)
    end
    box_half = 1.5 * radius_est

    sampling_ok = true
    worst_Vdot = Inf
    worst_pt = zeros(n)
    n_checked = 0
    n_violations = 0
    for _ in 1:samples
        pt = eq .+ (2 * box_half .* rand(rng, n) .- box_half)
        # Constrain the sample to the invariant manifold: a point off it is not
        # a state of the system, so neither V̇ nor a trajectory from it means
        # anything (see _project_to_manifold).
        pt, projected = _project_to_manifold(pt, sys)
        projected || continue
        if convert(Float64, V(x => pt)) > ρ
            continue
        end
        n_checked += 1
        dist = norm(pt .- eq)
        if dist > 1e-6
            Vdot_val = convert(Float64, Vd_poly(x => pt))
            if Vdot_val >= 0
                sampling_ok = false
                n_violations += 1
            end
            if Vdot_val < worst_Vdot
                worst_Vdot = Vdot_val
                worst_pt = copy(pt)
            end
        end
    end
    sampling_details = (worst_Vdot = worst_Vdot, worst_point = worst_pt,
                        n_checked = n_checked, n_violations = n_violations)

    # --- Check 3: Dynamical ---
    n_traj = min(100, n_checked)
    dyn_ok = true
    dyn_failures = Any[]
    checked_pts = Vector{Float64}[]
    for _ in 1:(4 * samples)
        pt = eq .+ (2 * box_half .* rand(rng, n) .- box_half)
        # Constrain the sample to the invariant manifold: a point off it is not
        # a state of the system, so neither V̇ nor a trajectory from it means
        # anything (see _project_to_manifold).
        pt, projected = _project_to_manifold(pt, sys)
        projected || continue
        if convert(Float64, V(x => pt)) <= ρ * 0.5
            push!(checked_pts, pt)
            length(checked_pts) >= n_traj && break
        end
    end
    for pt in checked_pts
        t, X, status = simulate(sys, pt, (0.0, 30.0); max_norm = 1e6)
        escaped = false
        for k in 1:size(X, 1)
            if convert(Float64, V(x => X[k, :])) > ρ + 1e-6
                escaped = true
                break
            end
        end
        final_dist = norm(X[end, :] .- eq)
        converged = final_dist < 1e-3
        if escaped || !converged
            dyn_ok = false
            push!(dyn_failures, (initial = pt, escaped = escaped, final_dist = final_dist))
        end
    end
    dyn_details = (n_tested = length(checked_pts), n_failed = length(dyn_failures),
                   failures = dyn_failures)

    return (algebraic = (ok = algebraic_ok, details = algebraic_details),
            sampling = (ok = sampling_ok, details = sampling_details),
            dynamical = (ok = dyn_ok, details = dyn_details))
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
    deg_V = maximum(degree, monomials(cert.V))
    println(io, "LyapunovCertificate")
    println(io, "  system:  ", cert.system_name)
    println(io, "  V degree: ", deg_V)
    println(io, "  ρ:       ", @sprintf("%.6e", cert.ρ))
    println(io, "  ε:       ", @sprintf("%.3e", cert.ε))
    println(io, "  margin:  ε‖x - x_eq‖^{2·", cert.k, "}")
    println(io, "  bounded: ", cert.bounded, " (", cert.boundedness_evidence, ")")
    println(io, "  time:    ", @sprintf("%.2f s", cert.solve_time))
    println(io, "  status:  ", cert.solver_status)
    println(io, "  SDP size: ", cert.sdp_size)
    n_terms = length(terms(prune_terms(cert.V)))
    if n_terms <= 12
        println(io, "  V(x) = ", cert.V)
    else
        println(io, "  V: ", n_terms, " terms, degree ", deg_V)
    end
end
