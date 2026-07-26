# Tests for julia/src/polynomial_utils.jl.
#
# This file is pure symbolic/numeric manipulation with no solver in the loop, so
# it should run in well under a second and should be the first thing that goes
# green. Every later failure is much easier to localise when these pass.

using Test
using LinearAlgebra
using Random
using DynamicPolynomials
using LyapunovToolbox

@testset "polynomial_utils" begin

    @testset "monomial_basis" begin
        @testset "cardinality matches the binomial count" begin
            # For n variables and total degree ≤ d the count is binomial(n+d, d);
            # with mindegree m subtract binomial(n+m-1, m-1). Check n ∈ 1:4, d ∈ 0:5.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "mindegree drops exactly the low-degree monomials" begin
            # mindegree=1 excludes the constant and nothing else; mindegree=2 also
            # excludes every linear monomial.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "ordering is deterministic across calls" begin
            # Two calls with identical arguments return identical vectors. A Gram
            # matrix is meaningless without a stable basis order.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "sos_basis halves the degree" begin
            # sos_basis(x, d) == monomial_basis(x, d); the polynomial it
            # parameterises has degree 2d. Assert the degree of
            # basis'*Q*basis for a random PSD Q is exactly 2d.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end
    end

    @testset "lie_derivative" begin
        @testset "matches hand-computed derivative on a known pair" begin
            # V = x1^2 + 2*x2^2, f = [-x1, -x2 + x1^2]
            # V̇ = 2*x1*(-x1) + 4*x2*(-x2 + x1^2) = -2x1^2 - 4x2^2 + 4x1^2*x2
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "matches finite differences at random points" begin
            # d/dt V(x(t)) along the flow, approximated by (V(x+h f(x)) - V(x))/h,
            # agrees with lie_derivative to O(h). Use h = 1e-6 and a 1e-4 tolerance.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "uses the caller's variable order, not the library's" begin
            # V that omits one state entirely (e.g. V = x2^2 in a 2-state system)
            # must still differentiate against the full x. This is the regression
            # test for the misalignment failure described in the docstring.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "rejects mismatched lengths" begin
            # length(f) != length(x) must throw ArgumentError, not zip-truncate.
            # Silent truncation here produces a V̇ that is missing a term and a
            # certificate that is simply wrong, with nothing to indicate it.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end
    end

    @testset "jacobian_at" begin
        @testset "reproduces the analytic Jacobian of the reversed Van der Pol" begin
            # A = [0 -1; 1 -1] at the origin for μ = 1.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "agrees with central differences away from the origin" begin
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end
    end

    @testset "recasting" begin
        @testset "taylor_recast reproduces the truncated series" begin
            # sin(x) at order 3 gives x - x^3/6 exactly, coefficient by coefficient.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "taylor_recast error grows as the stated order" begin
            # ‖f_true - f_poly‖ on |x| ≤ r scales like r^(order+1); check the ratio
            # between r and r/2 is within 20% of 2^(order+1).
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "trig_recast is exact on the invariant manifold" begin
            # Sample θ, set s = sin θ, c = cos θ, and assert the recast field equals
            # the original field to 1e-12 at every sample. Not approximately — exactly.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "trig_recast returns the s^2 + c^2 - 1 invariant" begin
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end
    end

    @testset "gram_residual" begin
        @testset "is zero for an exactly reconstructed polynomial" begin
            # Build p = basis' * Q * basis for random PSD Q, then assert the
            # residual is below 1e-12.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "is nonzero when a coefficient is perturbed" begin
            # Perturb one entry of Q by 1e-3 and assert the residual exceeds 1e-6.
            # A residual function that always returns ~0 would pass the test above.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end
    end

    @testset "prune_terms" begin
        @testset "removes only sub-tolerance coefficients" begin
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "tolerance is relative to the largest coefficient" begin
            # Scaling the whole polynomial by 1e6 must not change which terms survive.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end
    end

    @testset "evaluate_on_grid" begin
        @testset "row/column order is [y, x]" begin
            # Use an asymmetric polynomial (p = x1) so a transpose is detectable,
            # and assert grid[i, j] == p(xs[j], ys[i]).
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "fixed variables are substituted, not ignored" begin
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end
    end

    @testset "is_sublevel_set_bounded" begin
        @testset "positive definite quadratic is bounded with structural evidence" begin
            # V = x1^2 + x2^2 → (true, :top_form_pd).
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "detects an unbounded sublevel set" begin
            # V = x1^2 has {V ≤ 1} unbounded in x2 → (false, :unbounded_witness).
            # This is the case that silently invalidates an ROA claim.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "indefinite quartic top form is not reported as bounded" begin
            # V = x1^4 - x1^2*x2^2 + x2^4 is bounded; V = x1^4 - x2^4 is not.
            # Assert the second never returns evidence stronger than :sampling_only.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end
    end
end
