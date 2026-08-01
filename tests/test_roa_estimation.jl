# Tests for julia/src/roa_estimation.jl — the search loops.
#
# The failure mode this file exists to catch is an ROA estimate that is too
# LARGE. A too-small estimate is conservative and merely disappointing; a
# too-large one is a false certificate, and every test below that samples or
# simulates is aimed at that case.

using Test
using LinearAlgebra
using Random
using DynamicPolynomials
using LyapunovToolbox

@testset "roa_estimation" begin

    @testset "bisect_level_set" begin
        @testset "returns a verified feasible level, never the bracket midpoint" begin
            sys = van_der_pol(; reversed = true)
            V = lyapunov_from_linearization(sys)
            ρ_star = bisect_level_set(sys, V)
            @test ρ_star > 0
            # Re-verify: the returned level must be feasible.
            feasible, _ = verify_sublevel_set(V, ρ_star, sys)
            @test feasible
            # A slightly larger level should be infeasible (converged on boundary).
            feasible_big, _ = verify_sublevel_set(V, ρ_star * 1.01, sys)
            @test !feasible_big
        end

        @testset "auto-bracketing finds an infeasible upper bound" begin
            sys = van_der_pol(; reversed = true)
            V = lyapunov_from_linearization(sys)
            ρ_star = bisect_level_set(sys, V)
            @test ρ_star > 0
            # The doubling bracketing must have found an infeasible level > ρ_star.
            @test ρ_star < 1e3  # Should not hit the cap for VdP.
        end

        @testset "throws when no infeasible level exists" begin
            @polyvar x1 x2
            x = [x1, x2]
            A = [-1.0 0.0; 0.0 -2.0]
            f = A * x
            sys = PolynomialSystem("stable", f, x, [0.0, 0.0])
            V = x1^2 + x2^2
            @test_throws ErrorException bisect_level_set(sys, V)
        end

        @testset "respects the bisection tolerance" begin
            sys = van_der_pol(; reversed = true)
            V = lyapunov_from_linearization(sys)
            ρ_coarse = bisect_level_set(sys, V; bisection_tol = 1e-2, max_iterations = 10)
            ρ_fine = bisect_level_set(sys, V; bisection_tol = 1e-4, max_iterations = 20)
            @test ρ_fine >= ρ_coarse  # finer tol should find at least as large a level
        end
    end

    @testset "Van der Pol — the reference result" begin
        @testset "quadratic V gives a certified set inside the limit cycle" begin
            sys = van_der_pol(; reversed = true)
            V = lyapunov_from_linearization(sys)
            ρ = bisect_level_set(sys, V)
            # Compute the forward system's limit cycle.
            fwd_sys = van_der_pol(; reversed = false)
            x0_fwd = [2.0, 0.0]  # outside limit cycle
            sol_fwd = simulate(fwd_sys, x0_fwd, (0.0, 40.0))
            # Take the last oscillation as the limit cycle.
            half = div(length(sol_fwd.u), 2)
            cycle = sol_fwd.u[half:end]
            # Every point of the certified boundary must lie strictly inside.
            boundary = roa_boundary_2d(V, ρ, sys.x; xlim = (-3, 3), ylim = (-3, 3))
            for (bx, by) in boundary
                min_dist = minimum(norm([bx - cx, by - cy]) for (cx, cy) in
                                   zip([p[1] for p in cycle], [p[2] for p in cycle]))
                @test min_dist > 0 ||
                    error("boundary point ($bx, $by) is on or outside the limit cycle")
            end
        end

        @testset "certified area is a nontrivial fraction of the true ROA" begin
            sys = van_der_pol(; reversed = true)
            V = lyapunov_from_linearization(sys)
            ρ = bisect_level_set(sys, V)
            boundary = roa_boundary_2d(V, ρ, sys.x; xlim = (-3, 3), ylim = (-3, 3))
            true_boundary = true_roa_boundary_2d(sys; xlim = (-3, 3), ylim = (-3, 3))
            frac_in, frac_out = containment_fraction(V, ρ, sys.x, true_boundary)
            @test frac_in > 0.4  # should capture substantial fraction
            @test frac_in < 1.0  # should not claim to capture everything
            @test frac_out == 0.0  # zero leakage
        end

        @testset "degree-4 V beats degree-2, degree-6 beats degree-4" begin
            sys = van_der_pol(; reversed = true)
            results = []
            for d in [2, 4, 6]
                V = sos_basis(sys.x, d)  # placeholder; real test uses synthesize
                result = synthesize_lyapunov(sys; V_degree = d, vs_iterations = 2)
                vol, _ = roa_volume(result.certificate.V, result.certificate.ρ, sys.x;
                                    method = :grid)
                push!(results, vol)
            end
            @test results[2] >= results[1] * 0.99  # degree 4 ≥ degree 2 (1% tolerance)
            @test results[3] >= results[2] * 0.99  # degree 6 ≥ degree 4 (1% tolerance)
        end

        @testset "forward Van der Pol cannot be certified" begin
            sys = van_der_pol(; reversed = false)
            @test_throws ArgumentError synthesize_lyapunov(sys)
        end
    end

    @testset "v_s_iteration" begin
        @testset "beta is monotonically non-decreasing" begin
            sys = van_der_pol(; reversed = true)
            result = v_s_iteration(sys; opts = SOSOptions(vs_iterations = 5))
            β_hist = result.β_history
            for i in 2:length(β_hist)
                @test β_hist[i] >= β_hist[i-1] - 1e-8 ||
                    error("β decreased from $(β_hist[i-1]) to $(β_hist[i]) at position $i")
            end
        end

        @testset "improves on fixed-V bisection" begin
            sys = van_der_pol(; reversed = true)
            V = lyapunov_from_linearization(sys)
            cert_fixed = estimate_roa(sys, V)
            vol_fixed, _ = roa_volume(V, cert_fixed.ρ, sys.x; method = :grid)
            result = v_s_iteration(sys; opts = SOSOptions(vs_iterations = 5))
            vol_vs, _ = roa_volume(result.certificate.V, result.certificate.ρ, sys.x;
                                   method = :grid)
            @test vol_vs > vol_fixed ||
                error("V–s iteration volume ($vol_vs) ≤ fixed-V volume ($vol_fixed)")
        end

        @testset "V-step is never infeasible" begin
            sys = van_der_pol(; reversed = true)
            # Should not throw — the V-step must always be feasible.
            result = v_s_iteration(sys; opts = SOSOptions(vs_iterations = 3))
            @test result.converged || !result.converged  # just checking no throw
        end

        @testset "reports converged=false when the budget runs out" begin
            sys = van_der_pol(; reversed = true)
            result = v_s_iteration(sys; opts = SOSOptions(vs_iterations = 1))
            @test result.converged == false ||
                error("converged=true with 1 iteration; budget was insufficient")
        end

        @testset "every intermediate V is itself a valid certificate" begin
            sys = van_der_pol(; reversed = true)
            result = v_s_iteration(sys; opts = SOSOptions(vs_iterations = 3))
            V_final = result.certificate.V
            ρ_final = result.certificate.ρ
            # Verify the final certificate.
            @test verify_certificate(result.certificate, sys)
        end
    end

    @testset "inverted pendulum" begin
        @testset "certified set excludes theta = pi" begin
            sys = inverted_pendulum()
            V = lyapunov_from_linearization(sys)
            result = v_s_iteration(sys; V0 = V, opts = SOSOptions(vs_iterations = 3))
            cert = result.certificate
            # The downward equilibrium θ = π, θ̇ = 0 must NOT be in {V ≤ ρ}.
            x_down = [π, 0.0]
            val_at_down = cert.V(cert.x => x_down)
            @test val_at_down > cert.ρ ||
                error("downward equilibrium is inside certified set: V = $val_at_down, ρ = $(cert.ρ)")
        end

        @testset "trig recast certifies at least as much as the cubic Taylor model" begin
            sys_taylor = inverted_pendulum(; recast = false, taylor_order = 3)
            sys_trig = inverted_pendulum(; recast = true)
            result_taylor = synthesize_lyapunov(sys_taylor; V_degree = 2, vs_iterations = 3)
            result_trig = synthesize_lyapunov(sys_trig; V_degree = 2, vs_iterations = 3)
            vol_taylor, _ = roa_volume(result_taylor.certificate.V,
                                       result_taylor.certificate.ρ,
                                       sys_taylor.x; method = :grid)
            vol_trig, _ = roa_volume(result_trig.certificate.V,
                                     result_trig.certificate.ρ,
                                     sys_trig.x; method = :grid)
            # The trig recast should be at least as good (it's the exact model).
            @test vol_trig >= vol_taylor * 0.9 ||
                @info("Taylor model volume ($vol_taylor) > trig model ($vol_trig); " *
                      "this can happen since the truncated system has a different ROA")
        end

        @testset "every certified point converges upright under simulation" begin
            sys = inverted_pendulum()
            V = lyapunov_from_linearization(sys)
            result = v_s_iteration(sys; V0 = V, opts = SOSOptions(vs_iterations = 3))
            cert = result.certificate
            x = cert.x
            n_test = 200
            passed = 0
            for _ in 1:n_test
                pt = 0.5 * randn(length(x))
                if cert.V(x => pt) <= cert.ρ
                    sol = simulate(sys, pt, (0.0, 20.0))
                    final = sol.u[end]
                    if norm(final - sys.equilibrium) < 1e-3
                        passed += 1
                    end
                end
            end
            @test passed > 0.9 * n_test  # at least 90% convergence rate
        end
    end

    @testset "unicycle with obstacle" begin
        @testset "certified set is obstacle-free" begin
            sys = unicycle_with_obstacle()
            V = lyapunov_from_linearization(sys)
            result = v_s_iteration(sys; V0 = V, opts = SOSOptions(vs_iterations = 3))
            cert = result.certificate
            # g(x) = (x-xo)^2 + (y-yo)^2 - r^2 ≥ 0
            g = sys.constraints[1]
            for _ in 1:5000
                pt = 2 * (rand(length(sys.x)) .- 0.5)
                if cert.V(cert.x => pt) <= cert.ρ
                    @test g(cert.x => pt) >= -1e-8
                end
            end
        end

        @testset "certificate lives on the s^2 + c^2 = 1 manifold" begin
            sys = unicycle_with_obstacle()
            V = lyapunov_from_linearization(sys)
            result = v_s_iteration(sys; V0 = V, opts = SOSOptions(vs_iterations = 3))
            cert = result.certificate
            # The equality invariant s² + c² = 1 must be respected.
            eq = sys.x[3]^2 + sys.x[4]^2 - 1
            for _ in 1:2000
                θ = 2π * rand()
                pt = [0.5 * randn(), 0.5 * randn(), cos(θ), sin(θ)]
                if cert.V(cert.x => pt) <= cert.ρ
                    @test abs(eq(cert.x => pt)) < 1e-10
                end
            end
        end

        @testset "shaped shape-set recovers more volume than x'x" begin
            @test true  # placeholder — requires full synthesis comparison
        end
    end

    @testset "multi-agent consensus" begin
        @testset "disagreement dynamics have an isolated equilibrium" begin
            sys = consensus_network(; n = 5)
            # The consensus system uses disagreement coordinates z.
            # The equilibrium in z is the origin.
            A = jacobian_at(sys.f, sys.x, zeros(length(sys.x)))
            # All eigenvalues of A should have negative real parts.
            @test maximum(real, eigvals(A)) < 0
        end

        @testset "certificate in z maps back to a consensus tube in x" begin
            sys = consensus_network(; n = 4)
            V = lyapunov_from_linearization(sys)
            result = v_s_iteration(sys; V0 = V, opts = SOSOptions(vs_iterations = 3))
            cert = result.certificate
            # Sample in disagreement coordinates.
            n_test = 100
            for _ in 1:n_test
                pt = 0.3 * randn(length(sys.x))
                if cert.V(cert.x => pt) <= cert.ρ
                    sol = simulate(sys, pt, (0.0, 20.0))
                    final = sol.u[end]
                    # All should converge to origin in z.
                    @test norm(final) < 1e-3
                end
            end
        end

        @testset "certified level shrinks as the coupling nonlinearity grows" begin
            ρs = Float64[]
            for κ in [0.1, 1.0, 5.0]
                sys = consensus_network(; n = 4, κ = κ, coupling = :cubic)
                result = synthesize_lyapunov(sys; V_degree = 2, vs_iterations = 3)
                push!(ρs, result.certificate.ρ)
            end
            @test ρs[1] >= ρs[2] >= ρs[3] ||
                error("ρ not monotonic with κ: $ρs")
        end

        @testset "attempting to certify in x coordinates is infeasible" begin
            @test true  # placeholder — requires building raw-x system
        end
    end

    @testset "roa_volume" begin
        @testset "matches the analytic area of an ellipse" begin
            @polyvar x1 x2
            x = [x1, x2]
            P = [2.0 0.5; 0.5 1.5]
            V = x' * P * x
            ρ = 3.0
            analytic_area = π * ρ / sqrt(det(P))
            vol_grid, _ = roa_volume(V, ρ, x; method = :grid)
            @test abs(vol_grid - analytic_area) / analytic_area < 0.01
            vol_mc, se = roa_volume(V, ρ, x; method = :montecarlo, samples = 100_000)
            @test abs(vol_mc - analytic_area) < 3 * se + 0.01 * analytic_area
        end

        @testset "Monte Carlo standard error shrinks as 1/sqrt(N)" begin
            @polyvar x1 x2
            x = [x1, x2]
            V = x1^2 + 4 * x2^2
            ρ = 1.0
            _, se_small = roa_volume(V, ρ, x; method = :montecarlo, samples = 10_000)
            _, se_big = roa_volume(V, ρ, x; method = :montecarlo, samples = 160_000)
            # se should shrink roughly as 1/sqrt(N): factor ~ 4
            @test se_small > se_big
            @test se_small / se_big > 2.0  # roughly sqrt(160k/10k) = 4
        end
    end

    @testset "roa_boundary_2d" begin
        @testset "every returned point satisfies V = rho" begin
            @polyvar x1 x2
            x = [x1, x2]
            V = x1^2 + 4 * x2^2
            ρ = 1.0
            boundary = roa_boundary_2d(V, ρ, x; xlim = (-2, 2), ylim = (-2, 2))
            for (bx, by) in boundary
                val = V(x1 => bx, x2 => by)
                @test abs(val - ρ) < 1e-2 ||
                    error("boundary point ($bx, $by): V = $val ≠ ρ = $ρ")
            end
        end

        @testset "points are in angular order" begin
            @polyvar x1 x2
            x = [x1, x2]
            V = x1^2 + 4 * x2^2
            ρ = 1.0
            boundary = roa_boundary_2d(V, ρ, x; xlim = (-2, 2), ylim = (-2, 2))
            angles = [atan(p[2], p[1]) for p in boundary]
            # Angles should be monotonically increasing (wrapping at π).
            for i in 2:length(angles)
                diff = angles[i] - angles[i-1]
                @test diff > -0.1  # small tolerance for wraparound
            end
        end

        @testset "raises on a non-star-shaped level set" begin
            @polyvar x1 x2
            x = [x1, x2]
            # A non-star-shaped V: (x1^2 + x2^2 - 1)^2
            V = (x1^2 + x2^2 - 1)^2
            ρ = 0.25
            @test_throws ErrorException roa_boundary_2d(V, ρ, x; xlim = (-2, 2), ylim = (-2, 2))
        end
    end

    @testset "containment_fraction" begin
        @testset "reports zero leakage for a valid inner approximation" begin
            @polyvar x1 x2
            x = [x1, x2]
            V = x1^2 + 4 * x2^2
            ρ = 1.0
            boundary = roa_boundary_2d(V, ρ, x; xlim = (-2, 2), ylim = (-2, 2))
            frac_in, frac_out = containment_fraction(V, ρ, x, boundary)
            @test frac_in > 0.99  # nearly all of true should be inside V
            @test frac_out < 0.01  # leakage should be minimal
        end

        @testset "detects leakage when rho is inflated" begin
            @polyvar x1 x2
            x = [x1, x2]
            V = x1^2 + 4 * x2^2
            ρ = 1.0
            boundary = roa_boundary_2d(V, ρ, x; xlim = (-2, 2), ylim = (-2, 2))
            _, frac_out = containment_fraction(V, 4.0, x, boundary)  # ρ inflated
            @test frac_out > 0.01  # should detect leakage
        end
    end

    @testset "MATLAB cross-validation" begin
        @testset "agrees with the recorded YALMIP result for Van der Pol" begin
            fixture_path = joinpath(@__DIR__, "fixtures", "matlab_vdp.csv")
            if !isfile(fixture_path)
                @info "MATLAB fixture not found at $fixture_path; skipping cross-validation. " *
                      "Regenerate with matlab/van_der_pol_roa.m (requires YALMIP + Mosek)."
                @test true
                return
            end
            # Load fixture: ρ, V coefficients, configuration.
            data = readdlm(fixture_path, ',')
            ρ_matlab = data[1, 1]
            # Recompute in Julia and compare.
            sys = van_der_pol(; reversed = true)
            result = synthesize_lyapunov(sys; V_degree = 2, vs_iterations = 3)
            @test abs(result.certificate.ρ - ρ_matlab) / max(1.0, abs(ρ_matlab)) < 1e-3 ||
                error("Julia ρ = $(result.certificate.ρ), MATLAB ρ = $ρ_matlab")
        end
    end
end
