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
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.1")
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
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.1")
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
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.1")
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
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.1")
end

"""
    jacobian_at(f, x, x0 = zeros(length(x)))

Numeric Jacobian `A = ∂f/∂x` evaluated at `x0`, as a `Matrix{Float64}`.

Obtained by symbolic differentiation followed by substitution, not by finite
differences — `f` is polynomial, so there is no reason to accept truncation
error here. Used by [`lyapunov_from_linearization`](@ref).
"""
function jacobian_at(f::AbstractVector, x::AbstractVector, x0::AbstractVector = zeros(length(x)))
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.1")
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
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.1")
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
function trig_recast(f_trig, x::AbstractVector, θ_index::Integer)
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.1")
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
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.1")
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
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.1")
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
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.1")
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
                                 rng::Random.AbstractRNG = Random.default_rng())
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.1")
end
