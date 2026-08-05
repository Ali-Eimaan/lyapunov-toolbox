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
    nt = length(sys_taylor.x)
    ng = length(sys_trig.x)
    @printf("Taylor system: %d states, Gram dim @ deg 4 = %d\n",
            nt, gram_dimension(nt, 4))
    @printf("Trig system:   %d states, Gram dim @ deg 4 = %d\n",
            ng, gram_dimension(ng, 4))

    # ------------------------------------------------------- 2. certify each one
    opts = SOSOptions(; V_degree = 4)
    cert_taylor = synthesize_lyapunov(sys_taylor; opts=opts)
    cert_trig   = synthesize_lyapunov(sys_trig; opts=opts)
    @printf("Taylor: ρ = %.4f, time = %.2fs\n",
            cert_taylor.ρ, cert_taylor.solve_time)
    @printf("Trig:   ρ = %.4f, time = %.2fs\n",
            cert_trig.ρ, cert_trig.solve_time)

    # ------------------------------------------ 3. map both back to (θ,θ̇) and compare
    # Taylor: direct (θ, θ̇)
    bnd_taylor = roa_boundary_2d(cert_taylor.V, cert_taylor.ρ, cert_taylor.x;
                                 xlim=nothing, ylim=nothing, n=400)
    bnd_taylor_m = hcat(first.(bnd_taylor), last.(bnd_taylor))
    area_taylor = _polygon_area(bnd_taylor_m)

    # Trig: V_trig lives in (s, c, θ̇) with s²+c²=1, so evaluating it in
    # physical (θ, θ̇) via V(sinθ, cosθ, θ̇) is periodic in θ.  A 2-D
    # ray-march from the origin would never leave the sublevel set in the
    # θ-direction because sin/cos stay bounded.  Instead, grid θ ∈ [-π,π]
    # and bisect on θ̇ for each θ.
    V_trig = cert_trig.V
    ρ_trig = cert_trig.ρ
    x_trig = cert_trig.x
    n_theta = 400
    thetas = range(-π, π; length=n_theta)
    upper = Tuple{Float64,Float64}[]
    lower = Tuple{Float64,Float64}[]
    for θ in thetas
        sθ, cθ = sin(θ), cos(θ)
        # bisect positive θ̇
        lo, hi = 0.0, 100.0
        for _ in 1:60
            mid = (lo + hi) / 2
            if convert(Float64, V_trig(x_trig => [sθ, cθ, mid])) <= ρ_trig
                lo = mid
            else
                hi = mid
            end
        end
        push!(upper, (θ, lo))
        # bisect negative θ̇
        lo, hi = -100.0, 0.0
        for _ in 1:60
            mid = (lo + hi) / 2
            if convert(Float64, V_trig(x_trig => [sθ, cθ, mid])) <= ρ_trig
                hi = mid
            else
                lo = mid
            end
        end
        push!(lower, (θ, hi))
    end
    # closed polygon: upper envelope left→right, lower envelope right→left
    bnd_trig_m = hcat(first.(upper), last.(upper))
    bnd_trig_m = vcat(bnd_trig_m, hcat(first.(reverse(lower)), last.(reverse(lower))))
    area_trig = _polygon_area(bnd_trig_m)
    @printf("Projected areas — Taylor: %.4f, Trig: %.4f, ratio: %.2f%%\n",
            area_taylor, area_trig, area_taylor / area_trig * 100)

    # -------------------------------------------------- 4. the honest cross-check
    # True dynamics: θ̈ = (g/l)sinθ - (b/(ml²))θ̇ + u/(ml²), u = -K[θ;θ̇]
    # With defaults: m=1, l=1, b=0.1, g=9.81, K from LQR in system construction
    # The TRUE dynamics are exactly what the trig-recast system encodes -- it is
    # built with s = sin(theta), c = cos(theta) and no truncation. Evaluating it
    # at (sin th, cos th, thd) is the exact field, and it avoids re-deriving the
    # LQR gain here (the toolbox does not expose K; see REVIEW.md).
    function true_rhs(y)
        th, thd = y[1], y[2]
        z = [sin(th), cos(th), thd]
        fz = [convert(Float64, sys_trig.f[i](sys_trig.x => z)) for i in 1:3]
        # d(theta)/dt = thd; d(thd)/dt is the third recast component.
        return [thd, fz[3]]
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

    # Trig interior (should converge — the trig recast is exact)
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
    @printf("Trig cross-check:   %d/%d converged", trig_conv, n_test_trig)
    if trig_conv == n_test_trig
        @printf(" ✓\n")
    else
        frac = trig_conv / n_test_trig * 100
        @printf(" (%.0f%% — ρ=%.3f may be optimistic; see verify_certificate below)\n",
                frac, cert_trig.ρ)
    end

    # ---------------------------------------------------- 5. saddle sanity check
    # The Taylor-approximated system has its own equilibrium ≠ (±π, 0), so
    # checking V_taylor(π, 0) is meaningless for it.  Only the trig system
    # has (±π, 0) as a genuine saddle.
    V_trig_pi = cert_trig.V([sin(pi), cos(pi), 0.0])
    @assert V_trig_pi > cert_trig.ρ   "Trig certified set contains (±π, 0)!"
    println("  ✓ Trig certificate excludes the saddle at (±π, 0)")

    # ------------------------------------------------------- 6. Taylor order sweep
    println("\n--- Taylor Order Sweep ---")
    println("| order | ρ      | area   | converged | min-fail-|θ| | time [s] |")
    println("|-------|--------|--------|-----------|---------------|----------|")
    for order in [3, 5, 7]
        sys_o = inverted_pendulum(recast = :taylor, taylor_order = order)
        cert_o = synthesize_lyapunov(sys_o; opts=SOSOptions(; V_degree = 4))
        bnd_o = roa_boundary_2d(cert_o.V, cert_o.ρ, cert_o.x;
                                xlim=nothing, ylim=nothing, n=400)
        area_o = _polygon_area(hcat(first.(bnd_o), last.(bnd_o)))
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
    vc_t = verify_certificate(cert_taylor, sys_taylor)
    vc_g = verify_certificate(cert_trig, sys_trig)
    ok_t = vc_t.algebraic.ok && vc_t.sampling.ok && vc_t.dynamical.ok
    ok_g = vc_g.algebraic.ok && vc_g.sampling.ok && vc_g.dynamical.ok
    w_t = vc_t.sampling.details.worst_Vdot
    w_g = vc_g.sampling.details.worst_Vdot
    @assert ok_t && ok_g "Verification failed"
    @printf("Verify — Taylor: worst V̇ = %.3e, Trig: worst V̇ = %.3e\n", w_t, w_g)

    datadir = joinpath(dirname(dirname(@__DIR__)), "docs", "figures", "data")
    mkpath(datadir)
    open(joinpath(datadir, "inverted_pendulum.csv"), "w") do io
        writedlm(io, bnd_taylor_m, ',')
    end
    println("  Data written to docs/figures/data/inverted_pendulum.csv")
    return 0
end

_polygon_area(p::AbstractMatrix) = 0.5 * abs(sum(p[i, 1] * p[i+1, 2] - p[i+1, 1] * p[i, 2]
    for i in 1:(size(p, 1)-1)) + p[end, 1] * p[1, 2] - p[1, 1] * p[end, 2])

exit(main())
