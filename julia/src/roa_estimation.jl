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
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.4")
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
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.4")
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
function v_s_iteration(sys::PolynomialSystem; opts::SOSOptions = SOSOptions(), V0 = nothing, shape = nothing)
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.4")
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
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.4")
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
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.4")
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
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.4")
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
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.4")
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
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.4")
end
