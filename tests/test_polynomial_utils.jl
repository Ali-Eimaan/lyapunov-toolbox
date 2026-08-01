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
            for n in 1:4, d in 0:5, m in 0:2
                m <= d || continue
                @polyvar vars[1:n]
                b = monomial_basis(vars, d; mindegree = m)
                expected = binomial(n + d, d) - (m > 0 ? binomial(n + m - 1, m - 1) : 0)
                @test length(b) == expected
            end
        end

        @testset "mindegree drops exactly the low-degree monomials" begin
            @polyvar vars[1:3]
            b0 = monomial_basis(vars, 3)
            b1 = monomial_basis(vars, 3; mindegree = 1)
            b2 = monomial_basis(vars, 3; mindegree = 2)
            @test any(m -> m == 1, b0)
            @test !any(m -> m == 1, b1)
            @test all(m -> sum(exponents(m)) >= 1, b1)
            @test all(m -> sum(exponents(m)) >= 2, b2)
            # degree-2 and degree-3 monomials are present in all three
            @test length(b0) == length(b1) + 1
            @test length(b1) == length(b2) + 3
        end

        @testset "ordering is deterministic across calls" begin
            @polyvar vars[1:4]
            @test monomial_basis(vars, 3) == monomial_basis(vars, 3)
            @test monomial_basis(vars, 2; mindegree = 1) == monomial_basis(vars, 2; mindegree = 1)
        end

        @testset "sos_basis halves the degree" begin
            @polyvar vars[1:2]
            for d in 0:3
                @test sos_basis(vars, d) == monomial_basis(vars, d)
                basis = sos_basis(vars, d)
                n = length(basis)
                Q = randn(MersenneTwister(7 + d), n, n)
                Q = Q' * Q + I  # PSD
                p = polynomial(basis' * Q * basis)
                @test maximum(degree, monomials(p)) == 2d
            end
        end
    end

    @testset "lie_derivative" begin
        @testset "matches hand-computed derivative on a known pair" begin
            @polyvar x[1:2]
            V = x[1]^2 + 2 * x[2]^2
            f = [-x[1], -x[2] + x[1]^2]
            Vdot = lie_derivative(V, f, x)
            expected = -2 * x[1]^2 - 4 * x[2]^2 + 4 * x[1]^2 * x[2]
            @test Vdot == expected
        end

        @testset "matches finite differences at random points" begin
            @polyvar x[1:2]
            V = x[1]^2 + 2 * x[2]^2
            f = [-x[1], -x[2] + x[1]^2]
            Vdot = lie_derivative(V, f, x)
            rng = MersenneTwister(3)
            h = 1e-6
            for _ in 1:20
                x0 = randn(rng, 2)
                fv = [convert(Float64, f[1](x => x0)), convert(Float64, f[2](x => x0))]
                x1 = x0 .+ h .* fv
                fd = (convert(Float64, V(x => x1)) - convert(Float64, V(x => x0))) / h
                @test abs(fd - convert(Float64, Vdot(x => x0))) < 1e-4
            end
        end

        @testset "uses the caller's variable order, not the library's" begin
            @polyvar x[1:2]
            V = x[2]^2  # omits x[1] entirely
            f = [-x[1], -x[2]]
            Vdot = lie_derivative(V, f, x)
            @test Vdot == -2 * x[2]^2
        end

        @testset "rejects mismatched lengths" begin
            @polyvar x[1:2]
            V = x[1]^2 + x[2]^2
            @test_throws ArgumentError lie_derivative(V, [-x[1]], x)
            @test_throws ArgumentError lie_derivative(V, [-x[1], -x[2], x[1]], x)
        end
    end

    @testset "jacobian_at" begin
        @testset "reproduces the analytic Jacobian of the reversed Van der Pol" begin
            @polyvar x[1:2]
            μ = 1.0
            f = [-x[2], x[1] + μ * (x[1]^2 - 1) * x[2]]
            A = jacobian_at(f, x, [0.0, 0.0])
            @test size(A) == (2, 2)
            @test A ≈ [0.0 -1.0; 1.0 -1.0] atol = 1e-12
            # Hurwitz: eigenvalues have negative real part
            @test all(real(eigvals(A)) .< 0)
        end

        @testset "agrees with central differences away from the origin" begin
            @polyvar x[1:2]
            μ = 1.0
            f = [-x[2], x[1] + μ * (x[1]^2 - 1) * x[2]]
            rng = MersenneTwister(4)
            h = 1e-6
            for _ in 1:10
                x0 = 0.5 .* randn(rng, 2)
                A = jacobian_at(f, x, x0)
                Afd = zeros(2, 2)
                for j in 1:2
                    ep = zeros(2); ep[j] = h
                    fp = [convert(Float64, f[i](x => x0 .+ ep)) for i in 1:2]
                    fm = [convert(Float64, f[i](x => x0 .- ep)) for i in 1:2]
                    Afd[:, j] = (fp .- fm) ./ (2h)
                end
                @test A ≈ Afd atol = 1e-5
            end
        end
    end

    @testset "recasting" begin
        @testset "taylor_recast reproduces the truncated series" begin
            @polyvar v
            f_poly, meta = taylor_recast(y -> [sin(y)], [v], 3)
            p = f_poly[1]
            @test coefficient(p, v) ≈ 1.0 atol = 1e-9
            @test coefficient(p, v^3) ≈ -1 / 6 atol = 1e-9
            @test coefficient(p, v^2) ≈ 0.0 atol = 1e-9
            @test maximum(degree, monomials(p)) == 3
            @test meta.order == 3
        end

        @testset "taylor_recast error grows as the stated order" begin
            @polyvar v
            order = 3
            # exp, not sin: every coefficient of exp is nonzero, so the leading
            # truncation term really is the (order+1)-th one (sin's 4th coefficient
            # vanishes and the ratio would come out 2^(order+2)).
            f_poly, _ = taylor_recast(y -> [exp(y)], [v], order)
            p = f_poly[1]
            err(r) = maximum(abs(exp(y) - convert(Float64, p(v => y))) for y in range(-r, r; length = 2001))
            r = 0.4
            ratio = err(r) / err(r / 2)
            @test ratio ≈ 2^(order + 1) rtol = 0.2
        end

        @testset "trig_recast is exact on the invariant manifold" begin
            @polyvar x1 x2 θ s c
            x = [x1, x2, θ]
            # original field: ẋ1 = cosθ·x1 + x2, ẋ2 = sinθ·x2, θ̇ = x1 + x2
            f_trig = [c * x1 + x2, s * x2, x1 + x2]
            f_poly, x_new, _ = trig_recast(f_trig, x, 3; s = s, c = c)
            @test x_new == [x1, x2, s, c]
            rng = MersenneTwister(5)
            for _ in 1:50
                θ0 = 2π * rand(rng)
                x10, x20 = randn(rng), randn(rng)
                s0, c0 = sin(θ0), cos(θ0)
                vals = [x10, x20, s0, c0]
                # recast components 1:2 are the original field with cosθ → c, sinθ → s
                @test convert(Float64, f_poly[1](x_new => vals)) ≈ c0 * x10 + x20 atol = 1e-12
                @test convert(Float64, f_poly[2](x_new => vals)) ≈ s0 * x20 atol = 1e-12
                # ṡ = c·θ̇ = cosθ·θ̇ and ċ = -s·θ̇ = -sinθ·θ̇: exactly d/dt(sinθ), d/dt(cosθ)
                θdot = x10 + x20
                @test convert(Float64, f_poly[3](x_new => vals)) ≈ c0 * θdot atol = 1e-12
                @test convert(Float64, f_poly[4](x_new => vals)) ≈ -s0 * θdot atol = 1e-12
            end
        end

        @testset "trig_recast returns the s^2 + c^2 - 1 invariant" begin
            @polyvar x1 x2 θ s c
            x = [x1, x2, θ]
            f_trig = [c * x1, s * x2, x1]
            _, _, invs = trig_recast(f_trig, x, 3; s = s, c = c)
            @test length(invs) == 1
            @test invs[1] == s^2 + c^2 - 1
        end
    end

    @testset "gram_residual" begin
        @testset "is zero for an exactly reconstructed polynomial" begin
            @polyvar x[1:2]
            basis = sos_basis(x, 2)
            n = length(basis)
            Q = randn(MersenneTwister(11), n, n)
            Q = Q' * Q + I
            p = polynomial(basis' * Q * basis)
            @test gram_residual(p, basis, Q) < 1e-12
        end

        @testset "is nonzero when a coefficient is perturbed" begin
            @polyvar x[1:2]
            basis = sos_basis(x, 2)
            n = length(basis)
            Q = randn(MersenneTwister(12), n, n)
            Q = Q' * Q + I
            p = polynomial(basis' * Q * basis)
            Qp = copy(Q)
            Qp[1, 2] += 1e-3
            @test gram_residual(p, basis, Qp) > 1e-6
        end
    end

    @testset "prune_terms" begin
        @testset "removes only sub-tolerance coefficients" begin
            @polyvar x[1:2]
            p = x[1]^2 + 0.5 * x[1] * x[2] + 1e-12 * x[2]^2
            q = prune_terms(p; tol = 1e-9)
            @test Set(monomials(q)) == Set([x[1]^2, x[1] * x[2]])
            @test coefficient(q, x[1]^2) ≈ 1.0
            @test coefficient(q, x[1] * x[2]) ≈ 0.5
        end

        @testset "tolerance is relative to the largest coefficient" begin
            @polyvar x[1:2]
            p = x[1]^2 + 1e-12 * x[2]^2
            q1 = prune_terms(p; tol = 1e-9)
            q2 = prune_terms(1e6 * p; tol = 1e-9)
            @test Set(monomials(q1)) == Set(monomials(q2))
            @test Set(monomials(q1)) == Set([x[1]^2])
            @test coefficient(q2, x[1]^2) ≈ 1e6
        end
    end

    @testset "evaluate_on_grid" begin
        @testset "row/column order is [y, x]" begin
            @polyvar x[1:2]
            p = x[1]  # asymmetric: a transpose would be detected immediately
            xs = range(-1.0, 1.0; length = 5)
            ys = range(-2.0, 2.0; length = 7)
            M = evaluate_on_grid(p, x, xs, ys)
            @test size(M) == (7, 5)
            for i in eachindex(ys), j in eachindex(xs)
                @test M[i, j] ≈ xs[j]
            end
        end

        @testset "fixed variables are substituted, not ignored" begin
            @polyvar x[1:3]
            p = x[1] + x[2] * x[3]
            xs = range(-1.0, 1.0; length = 4)
            ys = range(-1.0, 1.0; length = 6)
            M = evaluate_on_grid(p, x, xs, ys; fixed = Dict(x[3] => 0.5))
            @test size(M) == (6, 4)
            for i in eachindex(ys), j in eachindex(xs)
                @test M[i, j] ≈ xs[j] + 0.5 * ys[i]
            end
        end
    end

    @testset "is_sublevel_set_bounded" begin
        @testset "positive definite quadratic is bounded with structural evidence" begin
            @polyvar x[1:2]
            V = x[1]^2 + x[2]^2
            bounded, evidence = is_sublevel_set_bounded(V, 1.0, x; rng = MersenneTwister(21))
            @test bounded
            @test evidence == :top_form_pd
        end

        @testset "detects an unbounded sublevel set" begin
            @polyvar x[1:2]
            V = x[1]^2  # {V ≤ 1} is a vertical strip: unbounded in x2
            bounded, evidence = is_sublevel_set_bounded(V, 1.0, x; rng = MersenneTwister(22))
            @test !bounded
            @test evidence == :unbounded_witness
        end

        @testset "indefinite quartic top form is not reported as bounded" begin
            @polyvar x[1:2]
            Vb = x[1]^4 - x[1]^2 * x[2]^2 + x[2]^4   # genuinely bounded
            Vu = x[1]^4 - x[2]^4                      # genuinely unbounded
            bounded_b, _ = is_sublevel_set_bounded(Vb, 1.0, x; rng = MersenneTwister(23))
            @test bounded_b
            bounded_u, evidence_u = is_sublevel_set_bounded(Vu, 1.0, x; rng = MersenneTwister(23))
            @test !(bounded_u && evidence_u == :top_form_pd)
            @test evidence_u in (:sampling_only, :unbounded_witness)
        end
    end
end
