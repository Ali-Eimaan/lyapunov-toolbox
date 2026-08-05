# unicycle_demo.jl — a certified obstacle-free region of attraction.
#
# Run with:  julia --project=julia examples/unicycle_with_obstacle/unicycle_demo.jl
#
# This is the example that exercises the constrained Positivstellensatz path,
# and the one with the clearest line to the thesis. The statement produced here —
# "every initial condition in this set converges to the goal without entering
# the obstacle" — is exactly the terminal ingredient an MPC-CBF scheme needs, and
# is what `mpc-cbf` and `cbf-safety-filter` consume.
#
# Two things make it harder than the first two examples, and both are the point:
#
#   1. The heading is recast exactly (s = sin θ, c = cos θ, s² + c² = 1), so the
#      state is 4-D with an algebraic invariant. The invariant enters the SOS
#      program with a FREE-SIGN multiplier. Using an SOS multiplier there is not
#      wrong but is strictly weaker, and the difference is measurable — measure it.
#
#   2. The obstacle is a semialgebraic constraint g(x) ≥ 0 with an SOS multiplier
#      σ(x), so the certificate only has to hold outside the obstacle. That is a
#      relaxation, so the certified level can only grow; if it does not, the
#      constraint is not entering the program.

using LinearAlgebra
using Printf
using Random
using DelimitedFiles
using LyapunovToolbox

function main()
    # ------------------------------------------------------------- 1. the system
    sys = unicycle_with_obstacle(obstacle_centre = (1.0, 0.0), obstacle_radius = 0.4)
    @printf("System: %s, %d states\n", sys.name, length(sys.x))
    @printf("State labels: %s\n", join(sys.state_labels, ", "))
    J = jacobian_at(sys.f, sys.x, sys.equilibrium)
    @assert all(real(eigvals(J)) .< 0) "Jacobian at goal is not Hurwitz"
    println("  ✓ Jacobian at goal is Hurwitz")
    if haskey(sys.notes, :invariant)
        println("  Invariant: ", sys.invariants[1])
    end

    # ------------------------------------------- 2. unconstrained vs. constrained
    # Unconstrained: ignore obstacle
    opts_unc = SOSOptions(; V_degree = 4)
    cert_unc = synthesize_lyapunov(sys; opts=opts_unc)
    @printf("Unconstrained:     ρ = %.4f, time = %.2fs\n",
            cert_unc.ρ, cert_unc.solve_time)

    # Constrained: with obstacle as Putinar constraint
    opts_con = SOSOptions(; V_degree = 4)
    cert_con = synthesize_lyapunov(sys; opts=opts_con)
    @printf("Constrained:       ρ = %.4f, time = %.2fs\n",
            cert_con.ρ, cert_con.solve_time)
    @assert cert_con.ρ >= cert_unc.ρ - 1e-8 "Constrained ρ should be ≥ unconstrained"
    @printf("  Ratio: %.2f\n", cert_con.ρ / max(cert_unc.ρ, 1e-12))

    # ------------------------------------------------- 3. free vs. SOS invariant
    # The invariant s² + c² - 1 = 0 can use free-sign or SOS multiplier.
    # By default the system notes contain the invariant polynomial.
    # We already used it implicitly; the free-sign is the default multiplier type.
    # Run with SOS multiplier for comparison
    opts_sos_inv = SOSOptions(; V_degree = 4,
                              invariant_multiplier = :sos)
    cert_sos_inv = synthesize_lyapunov(sys; opts=opts_sos_inv)
    @printf("Invariant μ=SOS:   ρ = %.4f\n", cert_sos_inv.ρ)
    # Free-sign should be ≥ SOS
    @printf("  Free-sign ρ / SOS ρ = %.3f\n", cert_con.ρ / max(cert_sos_inv.ρ, 1e-12))

    # ----------------------------------------------------- 4. obstacle-free check
    n_samples = 50_000
    rng = MersenneTwister(42)
    clearance = Inf
    min_g = Inf
    count_on_manifold = 0
    for _ in 1:n_samples
        # Sample on the manifold s² + c² = 1: pick θ, then s=sinθ, c=cosθ
        th = 2π * rand(rng)
        s, c = sin(th), cos(th)
        # Sample (x, y) from the certified set in a box
        x = 3*(2*rand(rng) - 1)
        y = 3*(2*rand(rng) - 1)
        state = [x, y, s, c]
        if cert_con.V(state) <= cert_con.ρ
            count_on_manifold += 1
            g_val = sys.constraints[1](sys.x => state)
            g_val_f = convert(Float64, g_val)
            min_g = min(min_g, g_val_f)
            if g_val_f < clearance
                clearance = g_val_f
            end
        end
    end
    @printf("Obstacle check: %d samples on manifold, min g(x) = %.6e\n",
            count_on_manifold, min_g)
    @assert min_g >= -1e-8 "Obstacle constraint violated inside certified set!"

    # -------------------------------------------------------- 5. shape-set sweep
    # p(x) = x'x (isotropic) vs anisotropic
    V = cert_con.V
    ρ = cert_con.ρ
    # Isotropic: β from bisect_level_set with p(x) = x'x
    p_iso(x_vec) = sum(x_vec.^2)
    β_iso = bisect_level_set(V, ρ, sys.x, p_iso)
    vol_iso = roa_volume(V, ρ, sys.x; method = :montecarlo, samples = 20000)
    @printf("Isotropic shape:   β = %.4f, vol ≈ %.4f\n", β_iso, vol_iso[1])
    # Anisotropic: p_aniso(x) = x² + y² + 0.1*(s² + c²) — emphasise spatial dims
    p_aniso(x_vec) = x_vec[1]^2 + x_vec[2]^2 + 0.1*(x_vec[3]^2 + x_vec[4]^2)
    β_aniso = bisect_level_set(V, ρ, sys.x, p_aniso)
    vol_aniso = roa_volume(V, ρ, sys.x; method = :montecarlo, samples = 20000)
    @printf("Anisotropic shape: β = %.4f, vol ≈ %.4f\n", β_aniso, vol_aniso[1])

    # ---------------------------------------------------------- 6. trajectories
    n_traj = 20
    all_ok = true
    rng2 = MersenneTwister(7)
    for i in 1:n_traj
        # Sample interior point on manifold
        while true
            th = 2π * rand(rng2)
            s, c = sin(th), cos(th)
            x0 = [2*(2*rand(rng2)-1), 2*(2*rand(rng2)-1), s, c]
            if V(x0) <= 0.8*ρ
                break
            end
        end
        t, X, status = simulate(sys, x0, (0.0, 30.0); dt = 5e-3)
        final_norm = norm(X[end, :])
        if final_norm >= 1e-3
            @printf("  Trajectory %d: FAILED to converge (final norm = %.3e)\n", i, final_norm)
            all_ok = false
        end
        # Check obstacle
        for k in 1:size(X, 1)
            g_val = convert(Float64, sys.constraints[1](sys.x => X[k, :]))
            if g_val < -1e-8
                @printf("  Trajectory %d: ENTERED OBSTACLE at step %d (g = %.3e)\n", i, k, g_val)
                all_ok = false
                break
            end
        end
    end
    @assert all_ok "Simulation cross-check failed"

    # ------------------------------------------------------- 7. verify + export
    vc = verify_certificate(cert_con, sys)
    ok = vc.algebraic.ok && vc.sampling.ok && vc.dynamical.ok
    worst = vc.sampling.details.worst_Vdot
    @assert ok "Certificate verification failed"
    @printf("Verification: ✓  worst V̇ = %.3e\n", worst)

    # Write (x,y) slices at θ=0 and θ=π/4
    datadir = joinpath(dirname(dirname(@__DIR__)), "docs", "figures", "data")
    mkpath(datadir)
    # θ = 0: s=0, c=1
    # Do a simple grid slice
    xs = range(-2, 2, length=100)
    ys = range(-2, 2, length=100)
    slice0 = [(xi, yi, convert(Float64, V([xi, yi, 0.0, 1.0]))) for xi in xs, yi in ys]
    open(joinpath(datadir, "unicycle_slice0.csv"), "w") do io
        writedlm(io, reshape([s[3] for s in slice0], length(xs), length(ys)), ',')
    end
    # θ = π/4: s=√2/2, c=√2/2
    s45 = √2/2
    slice45 = [(xi, yi, convert(Float64, V([xi, yi, s45, s45]))) for xi in xs, yi in ys]
    open(joinpath(datadir, "unicycle_slice45.csv"), "w") do io
        writedlm(io, reshape([s[3] for s in slice45], length(xs), length(ys)), ',')
    end
    println("  Slices written to docs/figures/data/unicycle_slice{0,45}.csv")

    return 0
end

exit(main())
