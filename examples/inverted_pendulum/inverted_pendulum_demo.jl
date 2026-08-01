# inverted_pendulum_demo.jl — certification under an approximated model.
#
# Run with:  julia --project=julia examples/inverted_pendulum/inverted_pendulum_demo.jl
#
# The point of this example is the gap between the model that is certified and
# the system that exists. Two recastings of sin(θ) are compared:
#
#   :taylor  — sin θ ≈ θ - θ³/6. A DIFFERENT vector field. A certificate proved
#              here is a statement about the truncated system, and transferring
#              it to the pendulum requires a bound on the truncation error that
#              this script does not prove and does not pretend to.
#   :trig    — exact, via s = sin θ, c = cos θ with s² + c² = 1. One more state,
#              a larger SDP, and a certificate about the actual pendulum.
#
# Running both, and reporting where the Taylor certificate stops being true of
# the real system, is the entire content of the example.

using LinearAlgebra
using Printf
using Random
using DelimitedFiles
using LyapunovToolbox

function main()
    # ------------------------------------------------------- 1. two formulations
    sys_taylor = inverted_pendulum(recast = :taylor, taylor_order = 3)
    sys_trig   = inverted_pendulum(recast = :trig)
    @printf("Taylor system: %d states, Gram dim @ deg 4 = %d\n",
            sys_taylor.n, gram_dimension(sys_taylor.n, 4))
    @printf("Trig system:   %d states, Gram dim @ deg 4 = %d\n",
            sys_trig.n, gram_dimension(sys_trig.n, 4))

    # ------------------------------------------------------- 2. certify each one
    opts = SOSOptions(; V_degree = 4)
    cert_taylor = synthesize_lyapunov(sys_taylor, opts)
    cert_trig   = synthesize_lyapunov(sys_trig, opts)
    @printf("Taylor: ρ = %.4f, #SDPs = %d, time = %.2fs\n",
            cert_taylor.ρ, cert_taylor.iters, cert_taylor.solve_time)
    @printf("Trig:   ρ = %.4f, #SDPs = %d, time = %.2fs\n",
            cert_trig.ρ, cert_trig.iters, cert_trig.solve_time)

    # ------------------------------------------ 3. map both back to (θ,θ̇) and compare
    # Taylor: direct (θ, θ̇)
    bnd_taylor = roa_boundary_2d(cert_taylor, 400)
    area_taylor = _polygon_area(bnd_taylor)

    # Trig: project from (s, c, θ̇) to (θ, θ̇) along s²+c²=1
    bnd_trig_4d = roa_boundary_2d(cert_trig, 400; slice_dims = [1, 3], slice_vals = [0.0, 1.0])
    # The trig certificate has states (s, c, θ̇) — project s=sinθ, c=cosθ → θ
    bnd_trig = zeros(size(bnd_trig_4d, 1), 2)
    for i in 1:size(bnd_trig_4d, 1)
        s, thd = bnd_trig_4d[i, 1], bnd_trig_4d[i, 2]
        th = atan(s, sqrt(max(0, 1 - s^2)))
        bnd_trig[i, :] = [th, thd]
    end
    area_trig = _polygon_area(bnd_trig)
    @printf("Projected areas — Taylor: %.4f, Trig: %.4f, ratio: %.2f%%\n",
            area_taylor, area_trig, area_taylor / area_trig * 100)

    # -------------------------------------------------- 4. the honest cross-check
    # True dynamics: θ̈ = (g/l)sinθ - (b/(ml²))θ̇ + u/(ml²), u = -K[θ;θ̇]
    # With defaults: m=1, l=1, b=0.1, g=9.81, K from LQR in system construction
    function true_rhs(y)
        th, thd = y[1], y[2]
        g, l, b, m = 9.81, 1.0, 0.1, 1.0
        K = sys_taylor.notes[:K]
        u = -K[1]*th - K[2]*thd
        thdd = (g/l)*sin(th) - (b/(m*l^2))*thd + u/(m*l^2)
        return [thd, thdd]
    end
    function simulate_true(y0, T)
        y = copy(y0)
        dt = 1e-3
        N = Int(T / dt)
        for _ in 1:N
            k1 = true_rhs(y)
            k2 = true_rhs(y .+ (dt/2).*k1)
            k3 = true_rhs(y .+ (dt/2).*k2)
            k4 = true_rhs(y .+ dt.*k3)
            y .+= (dt/6).*(k1 .+ 2k2 .+ 2k3 .+ k4)
            if any(abs.(y) .> 1e6); break; end
        end
        return y
    end

    # Sample Taylor interior
    rng = MersenneTwister(1)
    n_test = 200
    taylor_conv = 0; taylor_max_th = 0.0; taylor_min_fail_th = Inf
    for _ in 1:n_test
        th = 2*pi*(rand(rng) - 0.5)
        thd = 10*(2*rand(rng) - 1)
        if cert_taylor.V([th, thd]) <= cert_taylor.ρ
            yf = simulate_true([th, thd], 20.0)
            if norm(yf) < 1e-3
                taylor_conv += 1
                taylor_max_th = max(taylor_max_th, abs(th))
            else
                taylor_min_fail_th = min(taylor_min_fail_th, abs(th))
            end
        else
            n_test += 1  # resample
        end
    end
    @printf("Taylor cross-check: %d/%d converged", taylor_conv, n_test)
    if taylor_conv == n_test
        @printf(", max |θ| tested = %.3f rad\n", taylor_max_th)
    else
        @printf(", smallest failing |θ| = %.3f rad\n", taylor_min_fail_th)
    end

    # Trig interior (must be 100%)
    n_test_trig = 200
    trig_conv = 0
    for _ in 1:n_test_trig
        th = 2*pi*(rand(rng) - 0.5)
        thd = 10*(2*rand(rng) - 1)
        s, c = sin(th), cos(th)
        if cert_trig.V([s, c, thd]) <= cert_trig.ρ
            yf = simulate_true([th, thd], 20.0)
            if norm(yf) < 1e-3
                trig_conv += 1
            end
        else
            n_test_trig += 1
        end
    end
    @assert trig_conv == n_test_trig "Trig certificate failed cross-check: $trig_conv/$n_test_trig"
    @printf("Trig cross-check:   %d/%d converged ✓\n", trig_conv, n_test_trig)

    # ---------------------------------------------------- 5. saddle sanity check
    V_taylor_pi = cert_taylor.V([pi, 0.0])
    V_trig_pi   = cert_trig.V([sin(pi), cos(pi), 0.0])
    @assert V_taylor_pi > cert_taylor.ρ "Taylor certified set contains (±π, 0)!"
    @assert V_trig_pi   > cert_trig.ρ   "Trig certified set contains (±π, 0)!"
    println("  ✓ Neither certificate contains the saddle at (±π, 0)")

    # ------------------------------------------------------- 6. Taylor order sweep
    println("\n--- Taylor Order Sweep ---")
    println("| order | ρ      | area   | converged | min-fail-|θ| | time [s] |")
    println("|-------|--------|--------|-----------|---------------|----------|")
    for order in [3, 5, 7]
        sys_o = inverted_pendulum(recast = :taylor, taylor_order = order)
        cert_o = synthesize_lyapunov(sys_o, SOSOptions(; V_degree = 4))
        bnd_o = roa_boundary_2d(cert_o, 400)
        area_o = _polygon_area(bnd_o)
        # Cross-check
        conv = 0; max_th = 0.0; min_fail = Inf
        rng_o = MersenneTwister(1)
        n_o = 200
        for _ in 1:n_o
            th = 2*pi*(rand(rng_o) - 0.5)
            thd = 10*(2*rand(rng_o) - 1)
            if cert_o.V([th, thd]) <= cert_o.ρ
                yf = simulate_true([th, thd], 20.0)
                if norm(yf) < 1e-3
                    conv += 1; max_th = max(max_th, abs(th))
                else
                    min_fail = min(min_fail, abs(th))
                end
            else
                n_o += 1
            end
        end
        fail_str = conv == n_o ? "none" : @sprintf("%.3f", min_fail)
        @printf("| %5d | %6.4f | %6.4f | %3d/%-3d   | %13s | %8.2f |\n",
                order, cert_o.ρ, area_o, conv, n_o, fail_str, cert_o.solve_time)
    end

    # ------------------------------------------------------- 7. verify + export
    ok_t, w_t = verify_certificate(sys_taylor, cert_taylor)
    ok_g, w_g = verify_certificate(sys_trig, cert_trig)
    @assert ok_t && ok_g "Verification failed"
    @printf("Verify — Taylor: worst V̇ = %.3e, Trig: worst V̇ = %.3e\n", w_t, w_g)

    datadir = joinpath(dirname(dirname(@__DIR__)), "docs", "figures", "data")
    mkpath(datadir)
    open(joinpath(datadir, "inverted_pendulum.csv"), "w") do io
        writedlm(io, bnd_taylor, ',')
    end
    println("  Data written to docs/figures/data/inverted_pendulum.csv")
    return 0
end

_polygon_area(p::AbstractMatrix) = 0.5 * abs(sum(p[i, 1] * p[i+1, 2] - p[i+1, 1] * p[i, 2]
    for i in 1:(size(p, 1)-1)) + p[end, 1] * p[1, 2] - p[1, 1] * p[end, 2])

exit(main())
