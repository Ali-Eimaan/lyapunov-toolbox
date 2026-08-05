# van_der_pol_demo.jl — the reference example.
#
# Run with:  julia --project=julia examples/van_der_pol/van_der_pol_demo.jl
#
# This script is executed by CI on every commit, so it must exit 0 on success
# and non-zero on any failed assertion. It is not a scratchpad: everything it
# prints ends up quoted in docs/benchmarks.md or in the README.
#
# The reversed Van der Pol is the reference case because its exact ROA boundary
# is computable — it is the forward oscillator's stable limit cycle — so every
# claim here can be checked against ground truth rather than against another SOS
# computation.

using LinearAlgebra
using Printf
using Random
using DelimitedFiles
using LyapunovToolbox

function main()
    # ---------------------------------------------------------------- 1. system
    sys = van_der_pol(μ = 1.0)
    @printf("System: %s\n", sys.name)
    @printf("State labels: %s\n", join(sys.state_labels, ", "))
    J = jacobian_at(sys.f, sys.x, sys.equilibrium)
    @printf("Jacobian at origin:\n")
    display(J)
    @assert all(real(eigvals(J)) .< 0) "Jacobian at origin is not Hurwitz"
    println("  ✓ Jacobian is Hurwitz")

    # ------------------------------------------------- 2. quadratic V, fixed-V ρ
    V0 = lyapunov_from_linearization(sys)
    cert2 = estimate_roa(sys, V0)
    @printf("Quadratic V: ρ = %.6f, time = %.2fs\n", cert2.ρ, cert2.solve_time)

    # ------------------------------------------------------ 3. ground truth ROA
    # Forward oscillator: ẋ = [x₂; -x₁ - μ(x₁²-1)x₂]
    function fwd_rhs(y)
        x1, x2 = y[1], y[2]
        return [x2, -x1 - (x1^2 - 1)*x2]
    end
    function rk4_step(y, dt)
        k1 = fwd_rhs(y)
        k2 = fwd_rhs(y .+ (dt/2) .* k1)
        k3 = fwd_rhs(y .+ (dt/2) .* k2)
        k4 = fwd_rhs(y .+ dt .* k3)
        return y .+ (dt/6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
    end
    y = [2.0, 0.0]
    dt = 1e-3
    n_steps = Int(40 / dt)
    traj = zeros(n_steps + 1, 2)
    traj[1, :] = y
    for k in 1:n_steps
        y = rk4_step(y, dt)
        traj[k + 1, :] = y
    end
    half = div(n_steps + 1, 2)
    cycle = traj[half:end, :]
    true_area = _polygon_area(cycle)
    @printf("True ROA area: %.6f (from limit cycle integration)\n", true_area)

    # ----------------------------------------------------- 4. containment check
    boundary2 = roa_boundary_2d(cert2, 400)
    cert_area2 = _polygon_area(boundary2)
    inside = all(_point_in_polygon(boundary2[i, 1], boundary2[i, 2], cycle) for i in 1:size(boundary2, 1))
    @assert inside "Quadratic certified boundary NOT contained in true ROA"
    @printf("  Area: %.6f, ratio: %.1f%%\n", cert_area2, cert_area2 / true_area * 100)

    # -------------------------------------------------- 5. degree sweep + V–s it.
    println("\n--- Degree Sweep ---")
    println("| d_V | ρ      | area   | area/true | #SDPs | time [s] |")
    println("|-----|--------|--------|-----------|-------|----------|")
    prev_area = -Inf
    certs = LyapunovCertificate[]
    for d_V in [2, 4, 6]
        opts = SOSOptions(; V_degree = d_V)
        cert = synthesize_lyapunov(sys, opts)
        push!(certs, cert)
        bnd = roa_boundary_2d(cert, 400)
        area = _polygon_area(bnd)
        ratio = area / true_area
        @assert area > prev_area "Area not increasing with degree (d=$d_V)"
        prev_area = area
        inside = all(_point_in_polygon(bnd[i, 1], bnd[i, 2], cycle) for i in 1:size(bnd, 1))
        @assert inside "Degree-$d_V certified set NOT contained in true ROA"
        @printf("| %3d | %6.4f | %6.4f | %8.1f%% | %5d | %8.2f |\n",
                d_V, cert.ρ, area, ratio * 100, cert.iters, cert.solve_time)
    end

    # ----------------------------------------------------------- 6. verification
    println("\n--- Verification ---")
    for cert in certs
        ok, worst = verify_certificate(sys, cert)
        @assert ok "Certificate verification failed for d_V=$(cert.V_degree)"
        @printf("  d_V=%-2d: ✓  worst V̇ = %.6e\n", cert.V_degree, worst)
    end

    # ------------------------------------------------------------- 7. the figure
    datadir = joinpath(dirname(dirname(@__DIR__)), "docs", "figures", "data")
    mkpath(datadir)
    open(joinpath(datadir, "van_der_pol.csv"), "w") do io
        writedlm(io, boundary2, ',')
    end
    println("\n  Figure data written to docs/figures/data/van_der_pol.csv")

    return 0
end

# --- Inline geometry helpers ---
_polygon_area(p::AbstractMatrix) = 0.5 * abs(sum(p[i, 1] * p[i+1, 2] - p[i+1, 1] * p[i, 2]
    for i in 1:(size(p, 1)-1)) + p[end, 1] * p[1, 2] - p[1, 1] * p[end, 2])
function _point_in_polygon(x, y, poly::AbstractMatrix)
    inside = false
    n = size(poly, 1)
    j = n
    for i in 1:n
        xi, yi = poly[i, 1], poly[i, 2]
        xj, yj = poly[j, 1], poly[j, 2]
        if ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
            inside = !inside
        end
        j = i
    end
    return inside
end

exit(main())
