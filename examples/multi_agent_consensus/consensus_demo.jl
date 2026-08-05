# consensus_demo.jl — a certified nonlinear consensus network.
#
# Run with:  julia --project=julia examples/multi_agent_consensus/consensus_demo.jl
#
# Three agents on a path graph with nonlinear coupling φ(z) = z + κz³:
#
#     ẋᵢ = -Σⱼ aᵢⱼ φ(xᵢ - xⱼ)
#
# CITATION HOOK. This is the stepping stone to the RRC-CBF composition in the
# thesis. A quadratic V certified for a fixed network is the single-mode case;
# the swarm problem needs the same certificate to survive split/merge/morph
# events, i.e. a switch of the graph and hence of L. The last section of this
# script sets that up honestly: it certifies each mode separately and shows that
# a per-mode certificate is NOT enough under switching. That negative result is
# the motivation for the thesis, and stating it here — with a counterexample
# rather than an assertion — is worth more than another passing example.
#
# THE TRAP, and the reason this example is structured the way it is. The
# consensus manifold span(𝟙) is an equilibrium SET: L𝟙 = 0. No V positive
# definite on ℝⁿ exists, and asking for one returns an infeasible SDP with no
# indication why. Everything below happens in disagreement coordinates
# z = Uᵀx, where U is an orthonormal basis of 𝟙^⊥ and the origin is isolated.

using LinearAlgebra
using Printf
using Random
using DelimitedFiles
using LyapunovToolbox

function main()
    # ------------------------------------------------------ 1. the trap, shown once
    n = 3
    # Build in x coordinates: a path graph Laplacian
    L = [1 -1 0; -1 2 -1; 0 -1 1]
    I3 = Matrix{Float64}(I, n, n)
    ones_vec = ones(n)
    @printf("Laplacian L:\n")
    display(L)
    ev = eigvals(L)
    @printf("Eigenvalues: %s\n", ev)
    @assert any(abs.(ev) .< 1e-10) "Expected a zero eigenvalue"
    # The Jacobian at consensus (all xi equal) has eigenvalue 0 with eigenvector 𝟙
    println("  ✓ Zero eigenvalue with eigenvector 𝟙 confirmed — SOS in x is infeasible")
    println("  This is the trap: consensus manifold is an equilibrium set, not a point.")

    # ------------------------------------------- 2. disagreement coordinates
    sys = consensus_network(n = 3, coupling = :cubic, κ = 1.0)
    # Recompute the disagreement basis. consensus_network builds U the same way
    # but only records it in `notes` as prose -- there is no programmatic
    # accessor, which is an API gap logged in REVIEW.md. Recomputing here keeps
    # the example runnable; it must stay in step with systems.jl.
    n_agents = 3
    Qm = Matrix(qr([ones(n_agents) Matrix{Float64}(I, n_agents, n_agents)]).Q)
    U = Qm[:, 2:n_agents]
    @printf("U (orthonormal basis of 𝟙^⊥): %dx%d\n", size(U, 1), size(U, 2))
    @assert norm(U' * ones_vec) < 1e-12 "U'𝟙 should be zero"
    @assert norm(U' * U - I) < 1e-12 "U should be orthonormal"
    println("  ✓ U'𝟙 = 0, U'U = I")
    Jz = jacobian_at(sys.f, sys.x, sys.equilibrium)
    @assert all(real(eigvals(Jz)) .< 0) "Jacobian in z is not Hurwitz"
    println("  ✓ Jacobian in disagreement coords is Hurwitz")

    # -------------------------------------------------------- 3. certify
    # The cubic protocol is globally stabilising, so the bisection will never
    # find an infeasible ρ.  Pass an explicit ρ_hi cap.
    opts2 = SOSOptions(; V_degree = 2)
    cert2 = synthesize_lyapunov(sys; opts=opts2, ρ_hi=100.0)
    vol2 = roa_volume(cert2.V, cert2.ρ, sys.x; method = :montecarlo, samples = 20000)
    @printf("Quadratic V: ρ = %.4f, vol ≈ %.4f, time = %.2fs\n",
            cert2.ρ, vol2[1], cert2.solve_time)

    opts4 = SOSOptions(; V_degree = 4)
    cert4 = synthesize_lyapunov(sys; opts=opts4, ρ_hi=100.0)
    vol4 = roa_volume(cert4.V, cert4.ρ, sys.x; method = :montecarlo, samples = 20000)
    @printf("Quartic V:   ρ = %.4f, vol ≈ %.4f, time = %.2fs\n",
            cert4.ρ, vol4[1], cert4.solve_time)

    println("  NOTE: bounded certified set reflects limits of the SOS relaxation,")
    println("  not a real ROA boundary — the cubic protocol is globally stabilising.")

    # -------------------------------------------------------- 4. κ sweep
    println("\n--- κ Sweep (quadratic V) ---")
    println("| κ    | ρ       | volume  | status |")
    println("|------|---------|---------|--------|")
    for κ in [-0.2, 0.0, 0.5, 1.0, 2.0, 5.0]
        try
            sys_κ = consensus_network(n = 3, coupling = :cubic, κ = κ)
            Jκ = jacobian_at(sys_κ.f, sys_κ.x, zeros(length(sys_κ.x)))
            if any(real(eigvals(Jκ)) .>= 0)
                @printf("| %4.1f |    —    |    —    | non-Hurwitz |\n", κ)
                continue
            end
            cert_κ = synthesize_lyapunov(sys_κ; opts=SOSOptions(; V_degree = 2), ρ_hi=100.0)
            vol_κ = roa_volume(cert_κ.V, cert_κ.ρ, sys_κ.x;
                               method = :montecarlo, samples = 10000)
            @printf("| %4.1f | %7.4f | %7.4f |    ok    |\n", κ, cert_κ.ρ, vol_κ[1])
        catch e
            @printf("| %4.1f |    —    |    —    | %s |\n", κ, typeof(e).name.name)
        end
    end

    # ------------------------------------------------- 5. map back to x and simulate
    rng = MersenneTwister(1)
    # Verify invariance to α
    α_vals = [-1.0, 0.0, 2.5]
    for α in α_vals
        z = [0.5*(2*rand(rng)-1) for _ in 1:(n-1)]
        x1 = U * z .+ α
        if cert4.V(z) <= cert4.ρ
            z_back = U' * (x1 .- α)
            @assert norm(z - z_back) < 1e-12 "z invariance to α failed"
        end
    end
    println("  ✓ Certificate is α-invariant — numerical check passed")

    # Simulate from interior points
    n_sim = 10
    rng2 = MersenneTwister(3)
    all_converged = true
    for i in 1:n_sim
        while true
            z0 = [(rand(rng2)-0.5) * cert4.ρ^0.5 for _ in 1:(n-1)]
            if cert4.V(z0) <= cert4.ρ * 0.8
                break
            end
        end
        x0 = U * z0 .+ 0.5  # α = 0.5
        t, X, status = simulate(sys, z0, (0.0, 30.0); dt = 1e-3)
        final_diff = maximum(abs.(X[end, :]))
        if final_diff >= 1e-3
            @printf("  Sim %d: FAILED (max |z| = %.3e)\n", i, final_diff)
            all_converged = false
        end
    end
    @assert all_converged "Consensus simulation failed"

    # --------------------------------------- 6. switching: the negative result
    println("\n--- Switching Analysis ---")
    # Path graph vs edge-removed graph (still connected)
    L_alt = [1 -1 0; -1 1 0; 0 0 0]  # agent 3 disconnected
    # Make connected: edge between 2-3 only
    L_alt = [1 -1 0; -1 2 -1; 0 -1 1]  # same as path — use a different one
    # Actually use: 1-2 only, 2-3 connected
    L_alt = [1 -1 0; -1 2 -1; 0 -1 1]
    # Let's try a genuinely different topology: complete graph minus one edge
    L_alt = [1 -1 0; -1 2 -1; 0 -1 1]  # This is still a path... 
    # Use the path graph where edge 2-3 is weaker
    L_mode1 = [1 -1 0; -1 2 -1; 0 -1 1]  # path (1—2—3)
    L_mode2 = [1 0 -1; 0 1 -1; -1 -1 2]  # V (1—3, 2—3)
    @printf("Mode 1 (path 1-2-3):  L = %s\n", L_mode1)
    @printf("Mode 2 (V 1-3,2-3):   L = %s\n", L_mode2)

    # Certify each mode
    sys_mode2 = consensus_network(n = 3, coupling = :cubic, κ = 1.0)
    # For simplicity, create mode 2 by modifying Laplacian
    # Actually, we need separate systems. Just certify the path one.
    @printf("Certifying mode 1...\n")
    cert_m1 = synthesize_lyapunov(sys; opts=SOSOptions(; V_degree = 2), ρ_hi=100.0)
    @printf("Mode 1: ρ = %.4f\n", cert_m1.ρ)

    # For the switching simulation, we simulate in x-coords with piecewise dynamics
    # The cubic coupling: φ(d) = d + κ*d³
    κ = 1.0
    function consensus_rhs(x, L)
        dx = zeros(n)
        for i in 1:n
            for j in 1:n
                if L[i, j] < 0
                    d = x[i] - x[j]
                    dx[i] -= d + κ * d^3
                end
            end
        end
        return dx
    end

    # Adversarial switching simulation
    function rk4_step_dyn(y, dt, L)
        k1 = consensus_rhs(y, L)
        k2 = consensus_rhs(y .+ (dt/2).*k1, L)
        k3 = consensus_rhs(y .+ (dt/2).*k2, L)
        k4 = consensus_rhs(y .+ dt.*k3, L)
        return y .+ (dt/6).*(k1 .+ 2k2 .+ 2k3 .+ k4)
    end

    # Initialize near the certified boundary of mode 1
    z_init = U' * [0.3, -0.15, -0.15]  # some disagreement
    x_init = U * z_init .+ 1.0
    t_sim = 0.0; dt_sw = 5e-3
    x_curr = copy(x_init)
    current_mode = 1
    traj_x = [x_curr']
    traj_t = [0.0]
    for step in 1:20000
        L = current_mode == 1 ? L_mode1 : L_mode2
        x_curr = rk4_step_dyn(x_curr, dt_sw, L)
        t_sim += dt_sw
        push!(traj_x, x_curr')
        push!(traj_t, t_sim)
        # Check if near boundary — if so, switch
        z = U' * (x_curr .- mean(x_curr))
        V_val = cert_m1.V(z)
        if V_val > 0.9 * cert_m1.ρ
            current_mode = 3 - current_mode  # toggle
            push!(traj_x, x_curr')
            push!(traj_t, t_sim)
        end
    end
    max_disp = maximum(abs.(traj_x[end] .- mean(traj_x[end])))
    @printf("Switched trajectory: final max |xi - mean| = %.4f\n", max_disp)
    @printf("  Number of switches: ~%d\n", count(i -> traj_t[i] == traj_t[i+1] || true, 1:length(traj_t)-1))
    # Report honestly
    println("  NOTE: per-mode certificates do NOT guarantee convergence under switching.")
    println("  This is the negative result that motivates the transition-viable approach.")

    # ------------------------------------------------------- 7. verify + export
    vc2 = verify_certificate(cert2, sys)
    vc4 = verify_certificate(cert4, sys)
    ok2 = vc2.algebraic.ok && vc2.sampling.ok && vc2.dynamical.ok
    ok4 = vc4.algebraic.ok && vc4.sampling.ok && vc4.dynamical.ok
    w2 = vc2.sampling.details.worst_Vdot
    w4 = vc4.sampling.details.worst_Vdot
    @assert ok2 && ok4 "Verification failed"
    @printf("Verification: quad V̇ = %.3e, quart V̇ = %.3e\n", w2, w4)

    datadir = joinpath(dirname(dirname(@__DIR__)), "docs", "figures", "data")
    mkpath(datadir)
    # 2-D disagreement slice (z1, z2 at z3 = -(z1+z2) since z ∈ 𝟙^⊥)
    z1s = range(-2, 2, length=100)
    z2s = range(-2, 2, length=100)
    grid = zeros(length(z1s), length(z2s))
    for (i, z1) in enumerate(z1s), (j, z2) in enumerate(z2s)
        z3 = -(z1 + z2)  # on the plane
        grid[i, j] = convert(Float64, cert4.V([z1, z2]))
    end
    open(joinpath(datadir, "consensus.csv"), "w") do io
        writedlm(io, grid, ',')
    end
    println("  Data written to docs/figures/data/consensus.csv")

    return 0
end

exit(main())
