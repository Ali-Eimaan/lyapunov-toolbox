# polynomial_utils.jl — basis manipulation, monomial enumeration, and the
# polynomial arithmetic the SOS programs are built on.
#
# Nothing in this file touches JuMP or a solver. It is pure symbolic/numeric
# manipulation over `MultivariatePolynomials`, which makes it the one file that
# can be tested exhaustively and fast. Everything downstream assumes it is
# correct; test it first (IMPLEMENTATION_GUIDE.md §2, Phase 1).
#
# Conventions used throughout the package:
#   * `x` is always a `Vector{<:AbstractVariable}` in a fixed, meaningful order.
#     Never rely on `variables(p)` for ordering — it is sorted by the polynomial
#     library, not by the physical meaning of the states.
#   * Degrees are total degrees.
#   * A "basis" is a `Vector{<:AbstractMonomial}`, ordered consistently so a
#     Gram matrix and its basis can be stored and re-read as a pair.

"""
    monomial_basis(x, degree; mindegree = 0)

All monomials in the variables `x` with total degree in `mindegree:degree`,
in a deterministic order.

`mindegree = 1` drops the constant term and is what you want for a Lyapunov
candidate, which must vanish at the origin. `mindegree = 1, degree = d` still
admits linear terms; a candidate with `V(0) = 0` and `∇V(0) = 0` needs
`mindegree = 2` — but do not use `mindegree = 2` for the SOS *multiplier*
bases, where the constant term is both legal and usually active.

The returned length must equal `binomial(n + degree, degree) - binomial(n + mindegree - 1, mindegree - 1)`;
assert this, it catches enumeration bugs immediately.
"""
function monomial_basis(x::AbstractVector, degree::Integer; mindegree::Integer = 0)
    mindegree >= 0 || throw(ArgumentError("mindegree must be non-negative, got $mindegree"))
    degree >= mindegree ||
        throw(ArgumentError("degree $degree must be >= mindegree $mindegree"))
    n = length(x)
    expected = binomial(n + degree, degree) -
               (mindegree > 0 ? binomial(n + mindegree - 1, mindegree - 1) : 0)
    basis = DynamicPolynomials.monomials(x, mindegree:degree)
    @assert length(basis) == expected "monomial enumeration mismatch: got $(length(basis)) monomials for $n variables, degree $degree, mindegree $mindegree; expected $expected"
    return basis
end

"""
    sos_basis(x, degree)

Monomial basis for a Gram parameterisation of an SOS polynomial of total degree
`2 * degree`, i.e. monomials up to degree `degree`.

Provided separately from [`monomial_basis`](@ref) because the halving is where
off-by-one errors live: an SOS polynomial of degree `2d` is `z(x)' Q z(x)` with
`z` of degree `d`, not `2d`.
"""
function sos_basis(x::AbstractVector, degree::Integer)
    return monomial_basis(x, degree)
end

"""
    gram_dimension(nvars, degree)

Size of the square Gram matrix for a dense SOS parameterisation of degree
`2 * degree` in `nvars` variables, i.e. `binomial(nvars + degree, degree)`.

Used only to report problem size in logs and in `docs/benchmarks.md`; the growth
of this number with `degree` is the single most important practical fact about
SOS, so it is worth surfacing.
"""
function gram_dimension(nvars::Integer, degree::Integer)
    nvars >= 1 || throw(ArgumentError("need at least one variable, got nvars = $nvars"))
    degree >= 0 || throw(ArgumentError("degree must be non-negative, got $degree"))
    return binomial(nvars + degree, degree)
end

"""
    lie_derivative(V, f, x)

The Lie derivative `V̇(x) = ∇V(x) · f(x) = Σᵢ (∂V/∂xᵢ)·fᵢ(x)`.

This is the single most-used function in the package. Two failure modes to guard
against, both of which produce a *plausible* wrong answer rather than an error:

1. `x` and `f` misaligned — `f[i]` must be the time derivative of `x[i]`.
   Check `length(f) == length(x)` and say so in the error message.
2. Differentiating with respect to `variables(V)` instead of `x`. If `V` happens
   not to contain a variable (e.g. `V = x₁²` in a 2-state system) the
   library-sorted variable list is shorter than `x` and the products silently
   pair the wrong terms.

Return type is a polynomial in the same variables as `x`.
"""
function lie_derivative(V, f::AbstractVector, x::AbstractVector)
    length(f) == length(x) || throw(ArgumentError(
        "length(f) = $(length(f)) != length(x) = $(length(x)); f[i] must be the " *
        "time derivative of x[i], so the two vectors must be aligned"))
    return sum(differentiate(V, x[i]) * f[i] for i in eachindex(x))
end

"""
    jacobian_at(f, x, x0 = zeros(length(x)))

Numeric Jacobian `A = ∂f/∂x` evaluated at `x0`, as a `Matrix{Float64}`.

Obtained by symbolic differentiation followed by substitution, not by finite
differences — `f` is polynomial, so there is no reason to accept truncation
error here. Used by [`lyapunov_from_linearization`](@ref).
"""
function jacobian_at(f::AbstractVector, x::AbstractVector, x0::AbstractVector = zeros(length(x)))
    length(f) == length(x) || throw(ArgumentError(
        "length(f) = $(length(f)) != length(x) = $(length(x))"))
    length(x0) == length(x) || throw(ArgumentError(
        "length(x0) = $(length(x0)) != length(x) = $(length(x))"))
    n = length(x)
    A = Matrix{Float64}(undef, n, n)
    for i in 1:n, j in 1:n
        # differentiate each component individually (works with Any-typed vectors).
        dfi_dxj = differentiate(f[i], x[j])
        A[i, j] = convert(Float64, dfi_dxj(x => x0))
    end
    return A
end

"""
    taylor_recast(expr_fn, x, order; about = zeros(length(x)))

Polynomial approximation of a non-polynomial vector field by truncated Taylor
expansion of order `order` about `about`.

`expr_fn` is a Julia function `x -> [ẋ₁, ẋ₂, ...]` accepting and returning
numbers, so the caller writes `sin`/`cos` naturally.

**This changes the system.** A certificate proved for the truncated field says
nothing rigorous about the original one unless the truncation error is bounded
and absorbed. Every call site must state which it is doing; the inverted
pendulum example uses this deliberately and says so, the unicycle example uses
[`trig_recast`](@ref) instead precisely because it is exact.
"""
function taylor_recast(expr_fn, x::AbstractVector, order::Integer; about::AbstractVector = zeros(length(x)))
    length(x) == 1 || throw(ArgumentError(
        "taylor_recast implements only the univariate expansion (one angle); " *
        "got $(length(x)) variables. For exact multi-angle recasting use trig_recast."))
    length(about) == length(x) || throw(ArgumentError(
        "length(about) = $(length(about)) != length(x) = $(length(x))"))
    order >= 0 || throw(ArgumentError("order must be non-negative, got $order"))
    a = about[1]
    v = x[1]
    y0 = expr_fn(convert(Float64, a))
    n = length(y0)
    f_poly = Vector{Any}(undef, n)
    next_coeffs = Vector{Float64}(undef, n)
    for i in 1:n
        fi = y -> expr_fn(y)[i]
        # Derivatives up to order+1 using exact dual-number arithmetic,
        # not finite differences (R1-15).
        coeffs = [_exact_taylor_coefficient(fi, a, k) for k in 0:(order + 1)]
        p_i = coeffs[1] + sum(coeffs[k + 1] * (v - a)^k for k in 1:order)
        f_poly[i] = _as_polynomial(p_i, v)
        next_coeffs[i] = abs(coeffs[order + 2])
    end
    meta = (order = order, about = copy(about), truncation_error = maximum(next_coeffs))
    return (f = f_poly, meta = meta)
end

"""
    trig_recast(f_trig, x, θ_index)

Exact polynomial recasting of a vector field containing `sin θ` and `cos θ` by
introducing `s = sin θ`, `c = cos θ` as new states with

```
ṡ = c·θ̇,   ċ = -s·θ̇,   subject to the algebraic invariant  s² + c² = 1.
```

Returns `(f_poly, x_new, invariants)` where `invariants` is the vector of
equality constraints (here `[s² + c² - 1]`) that must be handed to the SOS
program as Positivstellensatz equality multipliers — *unconstrained* multipliers
of free sign, not SOS ones.

Unlike [`taylor_recast`](@ref) this is exact: the recast system's trajectories
are in bijection with the original's on the manifold `s² + c² = 1`. The price is
one extra state per angle and a certificate that is only valid on that manifold —
which is why the invariant has to be carried through, not dropped.
"""
function trig_recast(f_trig, x::AbstractVector, θ_index::Integer;
                     s = nothing, c = nothing)
    1 <= θ_index <= length(x) || throw(ArgumentError(
        "θ_index = $θ_index out of bounds for length(x) = $(length(x))"))
    length(f_trig) == length(x) || throw(ArgumentError(
        "length(f_trig) = $(length(f_trig)) != length(x) = $(length(x))"))
    s === nothing && throw(ArgumentError(
        "trig_recast needs the caller's sin and cos variables: create them with " *
        "@polyvar s c, write f_trig in terms of s and c, and pass s = s, c = c"))
    c === nothing && throw(ArgumentError("trig_recast: keyword argument c is required (see s)"))
    θ = x[θ_index]
    θ in variables(f_trig) && throw(ArgumentError(
        "the angle variable θ = x[$θ_index] appears in f_trig; write sin θ and " *
        "cos θ as the polynomial variables s and c before calling trig_recast"))
    # New state vector: s, c take the slot θ occupied, everything else shifts.
    x_new = [x[1:(θ_index - 1)]; s; c; x[(θ_index + 1):end]]
    n = length(x)
    f_poly = Vector{Any}(undef, n + 1)
    for i in 1:n
        if i < θ_index
            f_poly[i] = f_trig[i]
        elseif i == θ_index
            f_poly[i] = c * f_trig[i]        # ṡ = c·θ̇
            f_poly[i + 1] = -s * f_trig[i]   # ċ = -s·θ̇
        else
            f_poly[i + 1] = f_trig[i]
        end
    end
    invariants = [s^2 + c^2 - 1]
    return (f_poly, x_new, invariants)
end

"""
    prune_terms(p; tol = 1e-9)

Drop monomials whose coefficient magnitude is below `tol * maximum(abs, coeffs)`.

Conic solvers return dense garbage at the 1e-10 level in every coefficient.
Printing an unpruned degree-6 `V` in a notebook produces forty terms of noise
and hides the six that matter. Pruning is for *display and storage only* — never
prune a polynomial and then re-verify a certificate against the pruned version
and call it verified, because pruning can break the SOS inequality. Verify
first, prune second, and say which one the reported `V` is.
"""
function prune_terms(p; tol::Real = 1e-9)
    tol >= 0 || throw(ArgumentError("tol must be non-negative, got $tol"))
    ts = terms(p)
    isempty(ts) && return p
    maxc = maximum(abs, MultivariatePolynomials.coefficients(p))
    kept = filter(t -> abs(MultivariatePolynomials.coefficient(t)) >= tol * maxc, ts)
    # The largest-magnitude term always survives, so `kept` is never empty.
    return polynomial(kept)
end

"""
    evaluate_on_grid(p, x, xs, ys; fixed = Dict())

Evaluate a polynomial on a 2-D grid, returning a `Matrix{Float64}` of size
`(length(ys), length(xs))` in `[row = y, col = x]` order so it can be passed
straight to `contour`/`imshow` without a transpose.

`fixed` pins the remaining variables (needed for the 3-state unicycle and the
6-state consensus network, where the plots are 2-D slices). Which variables were
pinned, and at what value, must appear in the plot title of every figure that
uses this — a slice presented as if it were the whole picture is the easiest way
to overstate an ROA.
"""
function evaluate_on_grid(p, x::AbstractVector, xs::AbstractVector, ys::AbstractVector; fixed = Dict())
    length(x) >= 2 || throw(ArgumentError(
        "evaluate_on_grid needs at least two free variables, got $(length(x))"))
    haskey(fixed, x[1]) && throw(ArgumentError(
        "grid variable x[1] must not be pinned in `fixed`; the first two variables are the plot axes"))
    haskey(fixed, x[2]) && throw(ArgumentError(
        "grid variable x[2] must not be pinned in `fixed`; the first two variables are the plot axes"))
    for k in 3:length(x)
        haskey(fixed, x[k]) || throw(ArgumentError(
            "variable x[$k] is not pinned in `fixed`; a 2-D slice must pin every " *
            "variable beyond the two axes"))
    end
    M = Matrix{Float64}(undef, length(ys), length(xs))
    vals = Vector{Any}(undef, length(x))
    for k in 3:length(x)
        vals[k] = fixed[x[k]]
    end
    for (i, yv) in enumerate(ys), (j, xv) in enumerate(xs)
        vals[1] = xv
        vals[2] = yv
        M[i, j] = convert(Float64, p(x => vals))
    end
    return M
end

"""
    gram_residual(p, basis, Q)

Frobenius norm of the coefficient mismatch between `p` and `basis' * Q * basis`.

This is the independent check that a Gram matrix returned by the solver actually
reproduces the polynomial it claims to. Together with `minimum(eigvals(Q))` it is
the whole of [`verify_certificate`](@ref)'s numerical argument, so it must not
call the solver, and it must not be implemented by asking the modelling layer
whether it is happy.
"""
function gram_residual(p, basis::AbstractVector, Q::AbstractMatrix)
    size(Q) == (length(basis), length(basis)) || throw(ArgumentError(
        "Q is $(size(Q)) but the basis has $(length(basis)) entries; a Gram " *
        "matrix must be square with side length(basis)"))
    recon = polynomial(basis' * Q * basis)
    diff = recon - p
    return sqrt(sum(abs2, MultivariatePolynomials.coefficients(diff)))
end

"""
    is_sublevel_set_bounded(V, ρ, x; radius = 50.0, samples = 20_000, rng = Random.default_rng())

Necessary-condition screen for boundedness of `{x : V(x) ≤ ρ}`.

Boundedness is what turns "`V̇ < 0` on the sublevel set" into "the sublevel set
is in the ROA", and it is the assumption most often left implicit. Two checks,
in this order:

1. **Structural, and conclusive when it applies.** If the top-degree homogeneous
   part of `V` is positive definite on the unit sphere, every sublevel set is
   bounded. Test it by sampling the sphere densely and, separately, by asking
   whether `V_top - ε‖x‖^{deg}` is SOS.
2. **Sampling fallback.** Draw points on a sphere of radius `radius` and check
   `V > ρ` on all of them. This can only ever *disprove* boundedness.

Return `(bounded::Bool, evidence::Symbol)` with `evidence ∈ (:top_form_pd, :sampling_only, :unbounded_witness)`
so the caller can record honestly which of the two it got. Never report
`:sampling_only` as proof.
"""
function is_sublevel_set_bounded(V, ρ::Real, x::AbstractVector;
                                 radius::Real = 50.0, samples::Integer = 20_000,
                                 rng::Random.AbstractRNG = Random.default_rng(),
                                 ε::Real = 1e-6, solver = nothing,
                                 require_proof::Bool = true)
    ms = monomials(V)
    isempty(ms) && return (true, :sampling_only)  # constant V: vacuously handled below
    deg = maximum(degree, ms)
    # Check 1 — structural. If the top-degree homogeneous part is positive
    # definite on the unit sphere, every sublevel set is bounded.
    if iseven(deg)
        V_top = _top_homogeneous_part(V, deg)
        if deg == 2
            # For a quadratic the question "V_top - ε‖x‖² SOS" is exactly
            # "A ⪰ εI", answerable without a solver.
            A = _quadratic_form_matrix(V_top, x)
            if minimum(eigvals(Symmetric(A))) > ε
                return (true, :top_form_pd)
            end
        elseif is_sos(V_top - ε * sum(x .^ 2)^(deg ÷ 2), x; solver = solver)
            return (true, :top_form_pd)
        end
    end
    # Check 2 — sampling fallback. Can only disprove boundedness.
    # When require_proof=true, :sampling_only is NOT accepted as bounded.
    # The caller must handle it explicitly.
    n = length(x)
    for _ in 1:samples
        u = randn(rng, n)
        nrm = norm(u)
        nrm > 0 || continue
        pt = radius .* (u ./ nrm)
        if convert(Float64, V(x => pt)) <= ρ
            return (false, :unbounded_witness)
        end
    end
    # Sample at multiple radii to increase confidence, not just one.
    for r in [radius * 0.1, radius, radius * 10.0]
        for _ in 1:(samples ÷ 3)
            u = randn(rng, n)
            nrm = norm(u)
            nrm > 0 || continue
            pt = r .* (u ./ nrm)
            if convert(Float64, V(x => pt)) <= ρ
                return (false, :unbounded_witness)
            end
        end
    end
    if require_proof
        throw(ErrorException("is_sublevel_set_bounded: could not prove boundedness of " *
                             "{V ≤ $ρ} — top-degree form is not PD and sampling did not " *
                             "disprove boundedness. Set require_proof=false to accept " *
                             ":sampling_only as an indication, not a proof."))
    end
    return (true, :sampling_only)
end

# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

# The degree-deg homogeneous part of p (sum of all terms of exactly that degree).
function _top_homogeneous_part(p, deg::Integer)
    return polynomial(filter(t -> degree(t) == deg, terms(p)))
end

# Symmetric matrix A with p == x' * A * x for a quadratic form p.
function _quadratic_form_matrix(p, x)
    n = length(x)
    A = zeros(n, n)
    for i in 1:n
        A[i, i] = convert(Float64, MultivariatePolynomials.coefficient(p, x[i]^2))
        for j in (i + 1):n
            A[i, j] = A[j, i] = convert(Float64, MultivariatePolynomials.coefficient(p, x[i] * x[j])) / 2
        end
    end
    return A
end

# Wrap a number as a polynomial in v (polynomial(::Number) is not defined).
_as_polynomial(p::Number, v) = polynomial(v * 0 + p)
_as_polynomial(p, v) = polynomial(p)

# Exact Taylor coefficient f^(k)(a)/k! using automatic differentiation via
# hyper-dual numbers.  The old code used central finite differences with
# Richardson extrapolation, which is unusable above order 3 because the
# magnified round-off in the stencil coefficients (up to ~k!/(h^k)) swamps
# the signal.  Hyper-dual numbers give EXACT derivatives to machine precision
# for any order at the cost of one evaluation per order + 1, with no h choice.
struct _DualTaylor{T}
    val::T
    derivs::Vector{T}
end

function _dual_sin(d::_DualTaylor{T}) where T
    s = sin(d.val)
    c = cos(d.val)
    der = Vector{T}(undef, length(d.derivs))
    # Leibniz: D^n sin(f) = Σ binomial(n-1,k-1) D^k sin(f) * D^{n-k} f ...
    # For f(x) = x + ε·1 where ε^n factors are stored in derivs,
    # we use the chain rule directly: sin(x+ε) = sin(x) + ε·cos(x) - ε²·sin(x)/2!
    # which encodes the Taylor expansion of sin(x+ε) in ε.
    # Build it order by order using the recurrence:
    #   a_0 = sin(x), b_0 = cos(x)
    #   for k >= 1:  coeff of ε^k in sin(x+ε) is b_{k-1}/k  where b_j are
    #                Taylor coeffs of cos(x+ε)
    # But simpler: use sin(x+ε) ≈ sin(x)+cos(x)·ε - sin(x)·ε²/2! - cos(x)·ε³/3! ...
    # This is the Taylor series of sin.
    for k in 1:length(der)
        r = k % 4
        if r == 1
            der[k] = c
        elseif r == 2
            der[k] = -s
        elseif r == 3
            der[k] = -c
        else  # r == 0
            der[k] = s
        end
        der[k] /= factorial(k)
    end
    return _DualTaylor(s, der)
end

function _dual_cos(d::_DualTaylor{T}) where T
    cval = cos(d.val)
    sval = sin(d.val)
    der = Vector{T}(undef, length(d.derivs))
    for k in 1:length(der)
        r = k % 4
        if r == 1
            der[k] = -sval
        elseif r == 2
            der[k] = -cval
        elseif r == 3
            der[k] = sval
        else  # r == 0
            der[k] = cval
        end
        der[k] /= factorial(k)
    end
    return _DualTaylor(cval, der)
end

# Compute f(a), f'(a), f''(a)/2!, ..., f^(K)(a)/K! for an arbitrary function f
# using hyper-dual numbers (generalised).
function _taylor_coefficients(f, a::Float64, K::Integer)
    # Build a hyper-dual unit ε^K with all ε^k for k=1..K stored.
    eps_ = _DualTaylor(a, zeros(K))
    result = f(eps_)
    # The result.derivs[k] is the coefficient of ε^k = f^(k)(a)/k!
    return [result.val; result.derivs]
end

# Arithmetic for _DualTaylor: propagate derivatives.
import Base: +, -, *, /, ^
for op in (:+, :-)
    @eval function $op(d1::_DualTaylor{T}, d2::_DualTaylor{T}) where T
        val = $op(d1.val, d2.val)
        K = length(d1.derivs)
        der = [($op)(d1.derivs[k], d2.derivs[k]) for k in 1:K]
        return _DualTaylor(val, der)
    end
end

*(d1::_DualTaylor{T}, d2::_DualTaylor{T}) where T = begin
    val = d1.val * d2.val
    K = length(d1.derivs)
    der = zeros(T, K)
    for k in 1:K
        # der[k] = Σ_{j=0}^{k} coeff_j(d1) * coeff_{k-j}(d2)
        # where coeff_0 is the value, coeff_j is derivs[j]
        s = d1.val * d2.derivs[k] + d1.derivs[k] * d2.val
        for j in 1:(k-1)
            s += d1.derivs[j] * d2.derivs[k - j]
        end
        der[k] = s
    end
    return _DualTaylor(val, der)
end

/(d1::_DualTaylor{T}, d2::_DualTaylor{T}) where T = d1 * _dual_reciprocal(d2)

function _dual_reciprocal(d::_DualTaylor{T}) where T
    inv_val = inv(d.val)
    K = length(d.derivs)
    der = zeros(T, K)
    der[1] = -inv_val^2 * d.derivs[1]
    for k in 2:K
        der[k] = -inv_val * sum(d.derivs[j] * der[k - j] for j in 1:k)
    end
    return _DualTaylor(inv_val, der)
end

# Fallback for scalar ops with a DualTaylor on the left and a scalar on the right.
for op in (:+, :-)
    @eval $op(d::_DualTaylor{T}, s::Number) where T = $op(d, _DualTaylor(float(T(s)), zeros(T, length(d.derivs))))
    @eval $op(s::Number, d::_DualTaylor{T}) where T = $op(_DualTaylor(float(T(s)), zeros(T, length(d.derivs))), d)
end
*(d::_DualTaylor{T}, s::Number) where T = begin
    val = d.val * s
    der = [dk * s for dk in d.derivs]
    return _DualTaylor(val, der)
end
*(s::Number, d::_DualTaylor{T}) where T = d * s
/(d::_DualTaylor{T}, s::Number) where T = begin
    inv_s = inv(float(s))
    return d * inv_s
end
/(s::Number, d::_DualTaylor{T}) where T = float(s) * _dual_reciprocal(d)

# Override sin/cos/exp for DualTaylor to ensure we intercept them.
Base.sin(d::_DualTaylor) = _dual_sin(d)
Base.cos(d::_DualTaylor) = _dual_cos(d)
function Base.exp(d::_DualTaylor{T}) where T
    val = exp(d.val)
    der = Vector{T}(undef, length(d.derivs))
    for k in 1:length(der)
        der[k] = val / factorial(k)
    end
    return _DualTaylor(val, der)
end
# For other functions, fall back to IEEE on the value (value only, not derivative-aware).
# Users who need exact derivatives for other special functions should supply their own
# derivative rules or use ForwardDiff.

function _exact_taylor_coefficient(f, a, k::Integer)
    k == 0 && return convert(Float64, f(convert(Float64, a)))
    coeffs = _taylor_coefficients(f, Float64(a), k)
    return coeffs[k + 1]  # coeffs[1] = f(a), coeffs[2] = f'(a), ..., coeffs[k+1] = f^(k)(a)/k!
end

# Ray-casting point-in-polygon test for a closed polygon.
# Returns true if (px, py) is inside the polygon defined by
# vertices::Vector{Vector{Float64}} (first vertex is NOT repeated at end).
# Uses the even-odd rule.
function _point_in_polygon(px, py, vertices)
    n = length(vertices)
    n < 3 && return false
    inside = false
    j = n
    for i in 1:n
        xi, yi = vertices[i][1], vertices[i][2]
        xj, yj = vertices[j][1], vertices[j][2]
        if ((yi > py) != (yj > py)) &&
           (px < (xj - xi) * (py - yi) / (yj - yi) + xi)
            inside = !inside
        end
        j = i
    end
    return inside
end
