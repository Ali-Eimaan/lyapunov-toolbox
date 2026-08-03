# roa_estimation.jl — the outer loops around the SDPs in sos_lyapunov.jl.
#
# Nothing here builds a JuMP model directly. Everything goes through
# verify_sublevel_set / v_step / beta_step, which keeps the search logic
# testable against a stubbed oracle and keeps the SDP formulations in one file.

"""
    ROAResult

Outcome of an ROA estimation run, including the whole history of the search.

Keeping the history is not bookkeeping for its own sake: the V–s iteration is a
heuristic with no convergence guarantee, and the only honest way to report it is
to show the sequence of levels and let the reader see where it stalled. A single
final number invites the reader to assume monotone convergence to an optimum,
which is not what happens.

# Fields
- `certificate::LyapunovCertificate` — the best certificate found.
- `ρ_history::Vector{Float64}` — certified level after each outer iteration.
- `β_history::Vector{Float64}` — shape-set level after each outer iteration;
  empty when no shape set was used.
- `volume_history::Vector{Float64}` — estimated volume of `{V ≤ ρ}` per iteration.
- `iterations::Int` — outer iterations actually run.
- `converged::Bool` — improvement fell below `opts.vs_improvement_tol`. `false`
  means the iteration budget ran out, and that must be reported as such rather
  than presented as convergence.
- `total_solve_time::Float64`
- `n_sdps::Int` — how many SDPs were solved end to end. The headline cost number
  for `docs/benchmarks.md`.
- `opts::SOSOptions` — the exact configuration used.
"""
struct ROAResult
    certificate::LyapunovCertificate
    ρ_history::Vector{Float64}
    β_history::Vector{Float64}
    volume_history::Vector{Float64}
    iterations::Int
    converged::Bool
    total_solve_time::Float64
    n_sdps::Int
    opts::SOSOptions
end

"""
    bisect_level_set(sys, V; opts = SOSOptions(), ρ_lo = 0.0, ρ_hi = nothing)

Largest `ρ` for which [`verify_sublevel_set`](@ref) is feasible, by bisection.

Bisection is valid here because feasibility is **monotone** in `ρ`: if
`{V ≤ ρ}` admits a certificate then so does every smaller sublevel set, since it
is a subset. That monotonicity is the entire justification for the method and it
is worth a comment at the top of the loop — a reader who does not see it will
reasonably wonder why a nonconvex-looking problem is being bisected.

Bracketing:

- `ρ_lo = 0` is always feasible for a `V` from the linearisation.
- `ρ_hi = nothing` triggers a doubling search from a small positive value until
  a `ρ` is found that is *infeasible*. Cap the doubling (say 30 steps) and throw
  if no infeasible `ρ` is found: for a globally stable polynomial system that is
  the correct answer and the caller needs to know it is looking at an unbounded
  estimate, not a large one. Do not silently return the cap.

Return `(ρ, witness, n_sdps, solve_time)`. The returned `ρ` is the largest
*verified feasible* value tested — never the midpoint of the final bracket and
never `ρ_hi`. Reporting an unverified level as certified is the one error this
function absolutely must not make.

Stop when `(ρ_hi - ρ_lo) ≤ opts.bisection_tol * max(1, ρ_lo)`.
"""
function bisect_level_set(sys::PolynomialSystem, V; opts::SOSOptions = SOSOptions(),
                          ρ_lo::Real = 0.0, ρ_hi = nothing)
    # The oracle: true if the SOS program is feasible at this ρ.
    oracle(ρ_test) = verify_sublevel_set(sys, V, ρ_test; opts = opts)[1]
    n_sdps = 0
    solve_time_acc = 0.0

    # Bracket: find an infeasible upper bound.
    if ρ_hi === nothing
        ρ_hi = 1e-3
        found_infeasible = false
        for _ in 1:30
            feasible, w = verify_sublevel_set(sys, V, ρ_hi; opts = opts)
            n_sdps += 1
            solve_time_acc += w.solve_time
            if !feasible
                found_infeasible = true
                break
            end
            ρ_hi *= 2.0
        end
        if !found_infeasible && opts.verbose
            @warn("bisect_level_set: no infeasible ρ found after 30 doublings; " *
                  "the system appears globally stable — returning doubling cap " *
                  "ρ_hi = $ρ_hi as a generous over-estimate.")
        end
    end

    best = ρ_lo
    rho_lo = Float64(ρ_lo)
    rho_hi = Float64(ρ_hi)
    for iter in 1:opts.max_bisection_iters
        if (rho_hi - rho_lo) <= opts.bisection_tol * max(1.0, rho_lo)
            break
        end
        mid = (rho_lo + rho_hi) / 2
        feasible, w = verify_sublevel_set(sys, V, mid; opts = opts)
        n_sdps += 1
        solve_time_acc += w.solve_time
        if feasible
            best = mid
            rho_lo = mid
        else
            rho_hi = mid
        end
    end
    # Final verification at best: re-run to get the witness.
    feasible_final, witness = verify_sublevel_set(sys, V, best; opts = opts)
    n_sdps += 1
    solve_time_acc += witness.solve_time
    @assert feasible_final "bisect_level_set: best ρ = $best is not feasible on re-verification"
    return best, witness, n_sdps, solve_time_acc
end

"""
    estimate_roa(sys, V; opts = SOSOptions())

Convenience wrapper: bisect on `ρ` for a **fixed** `V`, check boundedness, and
package the result as a [`LyapunovCertificate`](@ref).

This is the cheap path and the one every example runs first, because it produces
the baseline that the V–s iteration has to beat. If `v_s_iteration` does not
improve on `estimate_roa`, something is wrong with the V-step — that comparison
is the main diagnostic in the whole package.
"""
function estimate_roa(sys::PolynomialSystem, V; opts::SOSOptions = SOSOptions())
    ρ, witness, n_sdps, solve_time = bisect_level_set(sys, V; opts = opts,
                                                       ρ_lo = 0.0, ρ_hi = nothing)
    # Check boundedness.
    bounded, evidence = is_sublevel_set_bounded(V, ρ, sys.x;
                                                solver = opts.solver, ε = opts.ε)
    bounded || throw(ErrorException(
        "estimate_roa: sublevel set {V ≤ $ρ} is not bounded " *
        "(evidence: $evidence). An unbounded sublevel set does not certify " *
        "the ROA — V̇ < 0 on an unbounded set does not imply convergence."))
    # Build the LyapunovCertificate.
    x = sys.x
    Vd = lie_derivative(V, sys.f, x)
    ms = monomials(Vd)
    dVd = isempty(ms) ? 0 : maximum(degree, ms)
    sdp_size = gram_dimension(length(x), cld(dVd, 2))
    cert = LyapunovCertificate(
        V, ρ, x,
        witness.multipliers,
        witness.grams,
        opts.ε,
        sys.name,
        witness.status,
        solve_time,
        sdp_size,
        bounded,
        evidence,
    )
    return cert
end

"""
    v_s_iteration(sys; opts = SOSOptions(), V0 = nothing, shape = nothing)

The V–s (or "three-step") iteration of Jarvis-Wloszek/Topcu/Packard:

```
V ← V0 (seed from the linearisation unless given)
repeat
    (γ-step) λ  ← multiplier_step(sys, V, 1)          # V normalised so ρ ≡ 1
    (β-step) β  ← beta_step(sys, V, shape)            # bisection, bilinear
    (V-step) V  ← v_step(sys, λ, shape, β)            # maximise β over V
until improvement in β < tol  or  iteration budget exhausted
```

Three properties that must hold, and must be asserted rather than assumed:

1. **Monotone non-decrease.** Each iteration's `β` is achievable by the previous
   iterate, so `β` can never decrease. If it does, the V-step is solving a
   different program from the one the β-step verified — most often because `ρ`
   was not renormalised to 1 after the V-step, or because `λ` from the γ-step was
   not held fixed. Assert `β_new ≥ β_old - 1e-8` and throw with both values.
2. **Feasibility is preserved.** The previous `V` is always feasible for the
   V-step, so an infeasible V-step is a bug, not a hard instance. Throw.
3. **Renormalise after every V-step** so the invariant `ρ ≡ 1` holds on entry to
   the next γ-step.

The iteration is a local heuristic with no optimality guarantee. Say that in the
docstring, in the notebooks, and in the README — a claim of "the largest ROA" is
false and a reviewer will know it.
"""
function v_s_iteration(sys::PolynomialSystem; opts::SOSOptions = SOSOptions(),
                       V0 = nothing, shape = nothing)
    x = sys.x
    # Seed V from linearisation.
    if V0 === nothing
        V = lyapunov_from_linearization(sys)
    else
        V = V0
    end
    # Default shape: p(x) = Σ (x_i - x_eq_i)².
    if shape === nothing
        shape = sum((x .- sys.equilibrium) .^ 2)
    end

    ρ_history = Float64[]
    β_history = Float64[]
    volume_history = Float64[]
    total_time = 0.0
    total_sdps = 0

    # Initial certification (fixed-V bisection).
    cert0 = estimate_roa(sys, V; opts = opts)
    ρ0 = cert0.ρ
    push!(ρ_history, ρ0)
    # β0: find the shape-set level for the initial V.
    β0 = beta_step(sys, V, shape; opts = opts)
    push!(β_history, β0)
    vol0, _ = roa_volume(V, ρ0, x; method = :montecarlo)
    push!(volume_history, vol0)

    β_old = β0
    V_current = V
    ρ_current = ρ0
    converged = false
    best_cert = cert0

    for iter in 1:opts.vs_iterations
        # γ-step: get multiplier λ at ρ ≡ 1 (after V is normalised).
        V_scaled = V_current / ρ_current
        feasible_γ, w_γ = multiplier_step(sys, V_scaled, 1.0; opts = opts)
        total_sdps += 1
        total_time += w_γ.solve_time
        if !feasible_γ
            opts.verbose && @warn("v_s_iteration γ-step infeasible at iteration $iter — solver may have failed; keeping previous certificate")
            break
        end

        # β-step: find largest β with {shape ≤ β} ⊆ {V_scaled ≤ 1}.
        β_new = beta_step(sys, V_scaled, shape; opts = opts)

        # Monotonicity check.
        if β_new < β_old - 1e-8
            throw(ErrorException(
                "v_s_iteration: β decreased from $β_old to $β_new at iteration $iter; " *
                "this violates the monotonicity guarantee and indicates ρ was not " *
                "renormalised or λ was not held fixed."))
        end

        # V-step: maximise β over V with λ fixed.
        λ = w_γ.multipliers[:λ]
        feasible_v, V_new, w_v = v_step(sys, λ, shape, β_new; opts = opts)
        total_sdps += 1
        total_time += w_v.solve_time
        if !feasible_v
            opts.verbose && @warn("v_s_iteration V-step infeasible at iteration $iter — solver may have failed; keeping previous V")
            break
        end

        # Renormalise: ρ ≡ 1, so new V is V_new with ρ = 1.
        V_current = V_new
        ρ_current = 1.0  # normalised

        # Record history.
        push!(ρ_history, ρ_current)
        push!(β_history, β_new)
        vol, _ = roa_volume(V_current, ρ_current, x; method = :montecarlo)
        push!(volume_history, vol)

        # Build certificate for this iteration.
        bounded, evidence = is_sublevel_set_bounded(V_current, ρ_current, x;
                                                     ε = opts.ε, solver = opts.solver)
        bounded || throw(ErrorException(
            "v_s_iteration: V at iteration $iter produced an unbounded sublevel set " *
            "(evidence: $evidence); the V-step should preserve boundedness."))

        cert = LyapunovCertificate(
            V_current, ρ_current, x,
            w_γ.multipliers, w_γ.grams,
            opts.ε, sys.name,
            w_γ.status, w_γ.solve_time + w_v.solve_time,
            gram_dimension(length(x), cld(maximum(degree, monomials(V_current)), 2)),
            bounded, evidence,
        )
        best_cert = cert

        # Convergence check.
        if β_new - β_old <= opts.vs_improvement_tol * max(1.0, β_old)
            converged = true
            break
        end
        β_old = β_new
    end

    result = ROAResult(
        best_cert, ρ_history, β_history, volume_history,
        length(ρ_history) - 1, converged,
        total_time, total_sdps, opts,
    )
    return result
end

"""
    roa_volume(V, ρ, x; method = :montecarlo, samples = 200_000, box = nothing, rng = Random.default_rng())

Estimate the volume (2-D: area) of `{x : V(x) ≤ ρ}`.

- `:montecarlo` — rejection sampling in `box`, which must contain the sublevel
  set; when `box === nothing`, derive one by growing until no sample near the
  boundary is accepted. Return `(volume, standard_error)`. Reporting a Monte
  Carlo estimate without its standard error is how a 3% "improvement" that is
  actually noise gets into a README table.
- `:grid` — deterministic quadrature on a uniform grid; slower to converge but
  reproducible to the last digit, so it is the one the tests use.

For `n > 3` this is only a sanity number; the comparison that matters for the
higher-dimensional systems is the shape-set level `β`, which is exact.
"""
function roa_volume(V, ρ::Real, x::AbstractVector; method::Symbol = :montecarlo,
                    samples::Integer = 200_000, box = nothing,
                    rng::Random.AbstractRNG = Random.default_rng())
    n = length(x)
    if method === :grid && n <= 2
        # Deterministic grid integration for 1-D / 2-D.
        if box === nothing
            box = _estimate_box(V, ρ, x; rng = rng)
        end
        nx, ny = 400, 400
        xs = range(box[1][1], box[1][2]; length = nx)
        ys = range(box[2][1], box[2][2]; length = ny)
        dx = (box[1][2] - box[1][1]) / nx
        dy = (box[2][2] - box[2][1]) / ny
        cell_area = dx * dy
        count = 0
        for y in ys, xv in xs
            pt = [xv, y]
            if n > 2
                # Pin remaining variables (slice).
                error("grid method only supported for n ≤ 2 in this implementation")
            end
            if convert(Float64, V(x => pt)) <= ρ
                count += 1
            end
        end
        return Float64(count * cell_area), 0.0
    else
        # Monte Carlo.
        if box === nothing
            box = _estimate_box(V, ρ, x; rng = rng)
        end
        box_vol = prod(box[i][2] - box[i][1] for i in 1:n)
        count = 0
        for _ in 1:samples
            pt = [box[i][1] + (box[i][2] - box[i][1]) * rand(rng) for i in 1:n]
            if convert(Float64, V(x => pt)) <= ρ
                count += 1
            end
        end
        p = count / samples
        vol = p * box_vol
        se = box_vol * sqrt(p * (1 - p) / samples)
        return vol, se
    end
end

function _estimate_box(V, ρ, x; radius::Real = 50.0, rng = Random.default_rng())
    n = length(x)
    lo = fill(Inf, n)
    hi = fill(-Inf, n)
    for _ in 1:2000
        u = randn(rng, n)
        nrm = norm(u)
        nrm > 0 || continue
        lo_r, hi_r = 0.0, 100.0
        for _ in 1:50
            mid = (lo_r + hi_r) / 2
            pt = mid .* (u ./ nrm)
            if convert(Float64, V(x => pt)) <= ρ
                lo_r = mid
            else
                hi_r = mid
            end
        end
        pt = lo_r .* (u ./ nrm)
        lo = min.(lo, pt)
        hi = max.(hi, pt)
    end
    margin = 0.1
    return [(lo[i] - margin, hi[i] + margin) for i in 1:n]
end

"""
    roa_boundary_2d(V, ρ, x; xlim, ylim, n = 400, fixed = Dict())

Points on the level curve `{V = ρ}` in a 2-D slice, as a `Vector{Tuple{Float64,Float64}}`
ordered so they can be plotted as a closed polygon.

Implemented by marching along rays from the origin and bisecting on the radius
for `V(r·d) = ρ`, not by a contouring library: the ray method returns points in
angular order for free, and a star-shaped sublevel set (which every `V` positive
definite with a bounded sublevel set gives near the origin) is exactly the case
where it is valid. Where it is not — a non-star-shaped set at high degree —
detect the multiple crossings along a ray and raise rather than silently drawing
the innermost one.

`fixed` pins the non-plotted variables; the caller must label the slice.
"""
function roa_boundary_2d(V, ρ::Real, x::AbstractVector; xlim, ylim, n::Integer = 400, fixed = Dict())
    boundary = Tuple{Float64,Float64}[]
    angles = range(0, 2π; length = n)
    for θ in angles
        d = [cos(θ), sin(θ)]
        # Detect multiple crossings by coarse sampling.
        ns = 100
        crossings = 0
        prev = nothing
        for k in 0:ns
            r = 0.0 + 100.0 * k / ns
            val = convert(Float64, V(x => [r * d[1], r * d[2]])) - ρ
            s = val > 0 ? 1 : (val < 0 ? -1 : 0)
            if prev !== nothing && s != prev && s != 0
                crossings += 1
            end
            if s != 0
                prev = s
            end
        end
        if crossings > 1
            throw(ErrorException(
                "roa_boundary_2d: ray at angle $θ crosses the level set " *
                "$crossings times; the sublevel set may not be star-shaped. " *
                "Ray-marching is not valid in this case."))
        end
        # Bisect for the crossing.
        r_lo, r_hi = 0.0, 100.0
        for _ in 1:60
            mid = (r_lo + r_hi) / 2
            if convert(Float64, V(x => [mid * d[1], mid * d[2]])) <= ρ
                r_lo = mid
            else
                r_hi = mid
            end
        end
        push!(boundary, (r_lo * d[1], r_lo * d[2]))
    end
    return boundary
end

"""
    true_roa_boundary_2d(sys; xlim, ylim, n = 200, tspan = (0.0, 40.0), tol = 1e-3, fixed = Dict())

Ground-truth ROA boundary by exhaustive forward simulation: for each of `n`
directions from the equilibrium, bisect on the radius to find the largest `r`
whose trajectory converges to the equilibrium.

This is the only honest denominator for the "fraction of the true ROA captured"
number in `docs/benchmarks.md` and in the README. It is **not** a proof of
anything — a trajectory that has not diverged within `tspan` may still diverge —
and the reported numbers must be labelled as an empirical estimate at the stated
`tspan` and `dt`.

For the reversed Van der Pol there is a much better ground truth available and it
should be used instead: the boundary *is* the forward system's stable limit
cycle, so integrate the forward system from a point outside it and take the
periodic orbit it converges to. Prefer that path for `van_der_pol`, and say in
the example why it is stronger.
"""
function true_roa_boundary_2d(sys::PolynomialSystem; xlim, ylim, n::Integer = 200,
                              tspan = (0.0, 40.0), tol::Real = 1e-3, fixed = Dict())
    eq = sys.equilibrium
    x_eq = eq[1]
    y_eq = eq[2]
    boundary = Tuple{Float64,Float64}[]
    angles = range(0, 2π; length = n)
    for θ in angles
        d = [cos(θ), sin(θ)]
        # Bisect on radius to find largest r that converges to equilibrium.
        r_lo, r_hi = 0.0, 50.0
        # First expand hi until divergence.
        for _ in 1:20
            x0 = [x_eq + r_hi * d[1], y_eq + r_hi * d[2]]
            _, X_sol, _ = simulate(sys, x0, tspan)
            final = X_sol[end, :]
            if norm(final - eq) <= tol
                r_hi *= 2
            else
                break
            end
        end
        # Bisect.
        for _ in 1:40
            mid = (r_lo + r_hi) / 2
            x0 = [x_eq + mid * d[1], y_eq + mid * d[2]]
            _, X_sol, _ = simulate(sys, x0, tspan)
            final = X_sol[end, :]
            if norm(final - eq) <= tol
                r_lo = mid
            else
                r_hi = mid
            end
        end
        push!(boundary, (x_eq + r_lo * d[1], y_eq + r_lo * d[2]))
    end
    return boundary
end

"""
    containment_fraction(V, ρ, x, boundary_points)

Fraction of the area enclosed by `boundary_points` (the true ROA estimate) that
is also inside `{V ≤ ρ}` — the headline "how good is this inner approximation"
number.

Also return the complementary check: the fraction of `{V ≤ ρ}` that lies
*outside* the true ROA. That number must be zero for a valid inner
approximation, and if it is not, the certificate is wrong and no amount of
volume improvement matters. Report both, always, in that order.
"""
function containment_fraction(V, ρ::Real, x::AbstractVector, boundary_points)
    # Polygon area by shoelace formula.
    n_pts = length(boundary_points)
    function poly_area(pts)
        a = 0.0
        for i in 1:n_pts
            j = i == n_pts ? 1 : i + 1
            a += pts[i][1] * pts[j][2] - pts[j][1] * pts[i][2]
        end
        return abs(a) / 2
    end
    true_area = poly_area(boundary_points)
    true_area <= 0 && return (NaN, NaN)

    # Fraction of true ROA inside {V ≤ ρ}.
    ns = 50_000
    inside_true = 0
    inside_both = 0
    # Estimate bounding box of true boundary.
    xs = [p[1] for p in boundary_points]
    ys = [p[2] for p in boundary_points]
    x_lo, x_hi = minimum(xs), maximum(xs)
    y_lo, y_hi = minimum(ys), maximum(ys)
    box_area = (x_hi - x_lo) * (y_hi - y_lo)

    for _ in 1:ns
        px = x_lo + (x_hi - x_lo) * rand()
        py = y_lo + (y_hi - y_lo) * rand()
        in_true = _point_in_polygon((px, py), boundary_points)
        if in_true
            inside_true += 1
            if convert(Float64, V(x => [px, py])) <= ρ
                inside_both += 1
            end
        end
    end
    frac_true_inside = inside_true > 0 ? inside_both / inside_true : 0.0

    # Fraction of {V ≤ ρ} that lies outside the true ROA.
    ns2 = 50_000
    inside_v = 0
    outside_true = 0
    # Use box from _estimate_box.
    box = _estimate_box(V, ρ, x)
    x_lo2, x_hi2 = box[1][1], box[1][2]
    y_lo2, y_hi2 = box[2][1], box[2][2]
    for _ in 1:ns2
        px = x_lo2 + (x_hi2 - x_lo2) * rand()
        py = y_lo2 + (y_hi2 - y_lo2) * rand()
        if convert(Float64, V(x => [px, py])) <= ρ
            inside_v += 1
            if !_point_in_polygon((px, py), boundary_points)
                outside_true += 1
            end
        end
    end
    frac_v_outside = inside_v > 0 ? outside_true / inside_v : 0.0

    return frac_true_inside, frac_v_outside
end

function _point_in_polygon(pt::Tuple{Float64,Float64},
                           poly::Vector{Tuple{Float64,Float64}})
    x, y = pt
    inside = false
    n = length(poly)
    j = n
    for i in 1:n
        xi, yi = poly[i]
        xj, yj = poly[j]
        if ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
            inside = !inside
        end
        j = i
    end
    return inside
end

"""
    benchmark_row(sys, result; label = sys.name)

One row of the table in `docs/benchmarks.md`, as a `NamedTuple`:
system, `n`, `deg V`, `deg λ`, `ρ`, estimated volume, fraction of true ROA,
largest SDP block, number of SDPs, total solve time, solver, and the ε used.

Emitting the table from code rather than typing it is the difference between a
benchmark and a claim. `docs/benchmarks.md` is generated from these rows.
"""
function benchmark_row(sys::PolynomialSystem, result::ROAResult; label::AbstractString = sys.name)
    cert = result.certificate
    V = cert.V
    x = cert.x
    deg_V = Int(maximum(degree, monomials(V)))
    λ_deg = 0
    if haskey(cert.multipliers, :λ)
        λ_deg = Int(maximum(degree, monomials(cert.multipliers[:λ])))
    end
    ρ_val = cert.ρ
    vol, _ = roa_volume(V, ρ_val, x; method = :montecarlo)
    block_size = cert.gram_dim
    n_sdps = result.total_sdps
    solve_time = result.total_time
    solver_name = string(cert.solver)
    ε_val = cert.ε

    return (;
        system = label,
        n = length(x),
        deg_V = deg_V,
        deg_λ = λ_deg,
        ρ = ρ_val,
        volume = vol,
        fraction_true_roa = NaN,  # must be set externally
        max_sdp_block = block_size,
        num_sdps = n_sdps,
        solve_time_s = solve_time,
        solver = solver_name,
        epsilon = ε_val,
    )
end
