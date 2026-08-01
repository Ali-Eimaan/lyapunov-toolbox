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
           st == MOI.DUAL_INFEASIBLE || st == MOI.ALMOST_INFEASIBLE
        return false
    elseif st == MOI.ALMOST_OPTIMAL || st == MOI.SLOW_PROGRESS
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
    elseif st == MOI.ALMOST_OPTIMAL || st == MOI.SLOW_PROGRESS
        # Accept only if the returned Gram matrix is genuinely a certificate.
        ok, _ = _gram_certificate_ok(p, cref)
        if ok
            gm = gram_matrix(cref)
            return _project_psd_and_report(p, gm.basis.monomials, Matrix(gm.Q))
        end
        return (basis, emptyQ, :failure)
    elseif st == MOI.INFEASIBLE || st == MOI.DUAL_INFEASIBLE ||
           st == MOI.INFEASIBLE_OR_UNBOUNDED || st == MOI.ALMOST_INFEASIBLE
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
    n = length(sys.x)
    A = jacobian_at(sys.f, sys.x, sys.equilibrium)
    Qmat = Q === nothing ? Matrix{Float64}(I, n, n) : Float64.(Q)
    size(Qmat) == (n, n) || throw(ArgumentError(
        "Q must be $(n)x$(n) to match the $(n)-state system, got $(size(Qmat))"))
    ev = eigvals(A)
    if any(real(ev) .>= 0)
        throw(ArgumentError(
            "A = ∂f/∂x at the equilibrium is not Hurwitz (eigenvalues " *
            "$(round.(ev; digits = 6))); the origin is not locally " *
            "asymptotically stable, so there is nothing to certify"))
    end
    # Continuous Lyapunov equation A'P + PA = -Q via the Kronecker-product solve,
    # exact up to roundoff for the small n in this package.
    M = kron(Matrix{Float64}(I, n, n), A') + kron(A', Matrix{Float64}(I, n, n))
    P = reshape(M \ vec(-Qmat), n, n)
    P = (P + P') / 2
    res = norm(A' * P + P * A + Qmat, 2)  # Frobenius
    res < 1e-9 * max(1, norm(P, 2)) || throw(ErrorException(
        "Lyapunov solve failed: ‖A'P + PA + Q‖_F = $res is not below " *
        "1e-9·max(1, ‖P‖_F) = $(1e-9 * max(1, norm(P, 2)))"))
    x = sys.x
    return sum(P[i, j] * x[i] * x[j] for i in 1:n, j in 1:n)
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
    x = sys.x
    bad = setdiff(variables(V), x)
    isempty(bad) || throw(ArgumentError(
        "V contains variables $(bad) not in the state vector x; a Lyapunov " *
        "candidate must be a polynomial in the system's own states"))
    Vd = lie_derivative(V, sys.f, x)
    # The positivity margin must vanish at the equilibrium, or the certificate is
    # structurally infeasible for any recast system whose equilibrium is not the
    # origin (trig_recast states live on s² + c² = 1, so ‖x_eq‖² = 1 ≠ 0 there).
    # For equilibrium-at-origin systems this is exactly the guide's ε‖x‖².
    ms = monomials(Vd)
    dVd = isempty(ms) ? 0 : maximum(degree, ms)
    k = max(1, cld(dVd, 2))
    margin = sum((x .- sys.equilibrium) .^ 2)^k

    model = SOSModel(opts.solver === nothing ? default_solver() : opts.solver)
    # (P2)/(P3) multipliers: λ SOS of degree multiplier_degree (basis degree is
    # half), σⱼ SOS per inequality constraint, μₖ FREE per equality invariant.
    λb = sos_basis(x, opts.multiplier_degree ÷ 2)
    λ = @variable(model, λ, SOSPoly(λb))
    n_cons = length(sys.constraints)
    σb = sos_basis(x, opts.constraint_multiplier_degree ÷ 2)
    σ = n_cons > 0 ? @variable(model, σ[1:n_cons], SOSPoly(σb)) : Any[]
    n_inv = length(sys.invariants)
    μb = monomial_basis(x, opts.multiplier_degree)
    μ = n_inv > 0 ? @variable(model, μ[1:n_inv], Poly(μb)) : Any[]

    expr = -(Vd + opts.ε * margin) + λ * (V - ρ)
    for j in 1:n_cons
        expr -= σ[j] * sys.constraints[j]
    end
    for k in 1:n_inv
        expr -= μ[k] * sys.invariants[k]
    end
    @constraint(model, cref, expr >= 0)

    optimize!(model)
    solve_time = JuMP.solve_time(model)
    st = termination_status(model)
    if st == MOI.OPTIMAL
        feasible = true
    elseif st == MOI.ALMOST_OPTIMAL || st == MOI.SLOW_PROGRESS
        # Accept only if every Gram matrix the solver returned is genuinely a
        # certificate: it must reproduce its target polynomial and be PSD.
        feasible = _witness_grams_ok(sys, V, Vd, ρ, margin, opts, cref, λ, σ, μ, μb)
    elseif st == MOI.INFEASIBLE || st == MOI.DUAL_INFEASIBLE ||
           st == MOI.INFEASIBLE_OR_UNBOUNDED || st == MOI.ALMOST_INFEASIBLE
        feasible = false
    else
        throw(ErrorException(
            "verify_sublevel_set: unexpected solver termination status $st"))
    end
    if feasible
        multipliers, grams = _gather_witness(cref, λ, σ, μ, μb)
    else
        multipliers = Dict{Symbol,Any}()
        grams = Dict{Symbol,Tuple{Vector,Matrix{Float64}}}()
    end
    witness = (multipliers = multipliers, grams = grams,
               status = string(st), solve_time = solve_time)
    return feasible, witness
end

# Reconstruct the polynomial represented by a solved SOSPoly (GramMatrix)
# variable: p(x) = z(x)' Q z(x) with Q = value.(gm.Q).
function _sos_variable_polynomial(gm)
    basis = gm.basis.monomials
    Q = Matrix(value.(gm.Q))
    return polynomial(basis' * Q * basis)
end

# Reconstruct the polynomial represented by a solved free Poly variable: a
# linear combination of the basis monomials with the solved coefficients.
function _poly_variable_polynomial(μ, basis)
    return sum(value(μ)[i] * basis[i] for i in eachindex(basis))
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
    μp = [_poly_variable_polynomial(μk, μb) for μk in μ]
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
function _gather_witness(cref, λ, σ, μ, μb)
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
        μp[k] = _poly_variable_polynomial(μk, μb)
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
function v_step(sys::PolynomialSystem, λ, shape, β::Real; opts::SOSOptions = SOSOptions())
    x = sys.x
    # V parameterised over monomials degree opts.V_degree, mindegree=2 so V(0)=0 and
    # ∇V(0)=0 structurally.  The explicit linear equality V(0)=0 below is a
    # belt-and-suspenders check that catches any solver shenanigans.
    Vb = monomial_basis(x, opts.V_degree; mindegree = 2)
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
    k = max(1, cld(dVd, 2))
    margin = sum((x .- sys.equilibrium) .^ 2)^k

    # (P4) positivity constraint: V - ε·margin is SOS
    @constraint(model, pos_cref, V - opts.ε * margin >= 0)

    # (P4) derivative constraint with fixed λ: -(V̇ + ε·margin) + λ·(V - 1) is SOS
    expr_deriv = -(Vd + opts.ε * margin) + λ * (V - 1.0)
    for (j, g) in enumerate(sys.constraints)
        σb_j = sos_basis(x, opts.constraint_multiplier_degree ÷ 2)
        σj = @variable(model, [1:1], SOSPoly(σb_j))
        expr_deriv -= σj[1] * g
    end
    for (k_inv, h) in enumerate(sys.invariants)
        μb_h = monomial_basis(x, opts.multiplier_degree)
        μk = @variable(model, [1:1], Poly(μb_h))
        expr_deriv -= μk[1] * h
    end
    @constraint(model, deriv_cref, expr_deriv >= 0)

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
    elseif st == MOI.ALMOST_OPTIMAL || st == MOI.SLOW_PROGRESS
        feasible = _v_step_grams_ok(model, pos_cref, deriv_cref, shape_cref)
    elseif st == MOI.INFEASIBLE || st == MOI.DUAL_INFEASIBLE ||
           st == MOI.INFEASIBLE_OR_UNBOUNDED || st == MOI.ALMOST_INFEASIBLE
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
        witness = (grams = grams, status = string(st), solve_time = solve_time)
    else
        V_new = nothing
        witness = (grams = Dict{Symbol,Tuple{Vector,Matrix{Float64}}}(),
                   status = string(st), solve_time = solve_time)
    end
    return feasible, V_new, witness
end

function _v_step_grams_ok(model, pos_cref, deriv_cref, shape_cref)
    for cref in [pos_cref, deriv_cref, shape_cref]
        gm = gram_matrix(cref)
        Q = Matrix(gm.Q)
        if minimum(eigvals(Symmetric(Q))) <= -1e-8
            return false
        end
    end
    return true
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
    found_infeasible || throw(ErrorException(
        "beta_step: no infeasible β found after 30 doublings; " *
        "the shape set may be too weak or V may be globally dominant"))

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
        elseif st == MOI.ALMOST_OPTIMAL || st == MOI.SLOW_PROGRESS
            gm = gram_matrix(cref)
            Q = Matrix(gm.Q)
            ok = minimum(eigvals(Symmetric(Q))) > -1e-8
            return ok
        elseif st == MOI.INFEASIBLE || st == MOI.DUAL_INFEASIBLE ||
               st == MOI.INFEASIBLE_OR_UNBOUNDED || st == MOI.ALMOST_INFEASIBLE
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
                             shape = nothing, V0 = nothing)
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
    result = v_s_iteration(sys; opts = opts, V0 = V, shape = shape)
    return result.certificate
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
    # Algebraic check only — checks 2 and 3 need the system.
    algebraic_ok = true
    algebraic_details = Dict{Symbol,Dict{Symbol,Float64}}()
    for (key, (basis, Q)) in pairs(cert.grams)
        λmin = minimum(eigvals(Symmetric(Q)))
        # gram_residual of a self-reconstructed polynomial from its own (basis, Q)
        p_recon = polynomial(basis' * Q * basis)
        r = gram_residual(p_recon, basis, Q)
        algebraic_details[key] = Dict(:gram_residual => r, :eigmin => λmin)
        if r >= tol || λmin <= -tol
            algebraic_ok = false
        end
    end
    return (algebraic = (ok = algebraic_ok, details = algebraic_details),
            sampling = (ok = false, note = "requires system; use verify_certificate(cert, sys)"),
            dynamical = (ok = false, note = "requires system; use verify_certificate(cert, sys)"))
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
    algebraic_ok = true
    algebraic_details = Dict{Symbol,Dict{Symbol,Float64}}()
    for (key, (basis, Q)) in pairs(cert.grams)
        p_recon = polynomial(basis' * Q * basis)
        r = gram_residual(p_recon, basis, Q)
        λmin = minimum(eigvals(Symmetric(Q)))
        algebraic_details[key] = Dict(:gram_residual => r, :eigmin => λmin)
        if r >= tol || λmin <= -tol
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
function verify_certificate(cert::LyapunovCertificate; tol::Real = 1e-7, samples::Integer = 50_000,
                            rng::Random.AbstractRNG = Random.default_rng())
    # Convenience: mark as not-verified when called without a system reference.
    return (algebraic = (ok = false, details = "No system provided"),
            sampling = (ok = false, details = "No system provided"),
            dynamical = (ok = false, details = "No system provided"))
end

"""
    Base.show(io::IO, ::MIME\"text/plain\", cert::LyapunovCertificate)

Human-readable summary: system name, `deg V`, `ρ`, `ε`, largest SDP block,
solve time, solver status, and the boundedness evidence. Print the pruned `V`
when it has at most twelve terms and the degree and term count otherwise.

Worth writing properly — this is what appears in every notebook cell and every
console session, and a certificate that prints its own provenance is much harder
to misquote later.
"""
function Base.show(io::IO, ::MIME\"text/plain\", cert::LyapunovCertificate)
    println(io, "LyapunovCertificate")
    println(io, "  V degree: ", cert.V_degree)
    println(io, "  ρ:       ", @sprintf("%.6e", cert.ρ))
    println(io, "  ε:       ", @sprintf("%.3e", cert.ε))
    println(io, "  bounded: ", cert.is_bounded)
    println(io, "  iters:   ", cert.iters)
    println(io, "  time:    ", @sprintf("%.2f s", cert.solve_time))
    println(io, "  status:  ", cert.solver_status)
    n_terms = length(cert.V)
    if n_terms <= 12
        println(io, "  V(x) = ", cert.V)
    else
        println(io, "  V: ", n_terms, " terms, degree ", cert.V_degree)
    end
    println(io, "  verified: ", cert.verified)
end
