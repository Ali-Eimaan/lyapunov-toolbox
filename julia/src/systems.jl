# systems.jl — the benchmark vector fields, defined exactly once.
#
# This file is not in the original repository specification. It was added
# because the same four systems appear in `tests/`, in all four `examples/`, in
# `matlab/*.m`, and in `docs/benchmarks.md`; without a single definition those
# five copies drift, and a cross-language agreement test between Julia and
# MATLAB that compares two *different* Van der Pol oscillators is worse than no
# test at all. The MATLAB scripts restate the dynamics by necessity — the guide
# requires a comment there pointing back at this file, and a test that compares
# the two numerically. Keep it.

"""
    PolynomialSystem

A polynomial vector field together with everything a certificate needs to be
interpretable.

# Fields
- `x::Vector` — state variables, in physical order. Never reorder these.
- `f::Vector` — `f[i]` is `ẋ[i]`, same length and order as `x`.
- `name::String` — used in plot titles, log lines, and benchmark table rows.
- `equilibrium::Vector{Float64}` — the point being certified; must satisfy
  `f(equilibrium) = 0`, checked at construction.
- `constraints::Vector` — semialgebraic set `{gⱼ(x) ≥ 0}` the certificate must
  stay inside (obstacles, actuator limits, the region where a Taylor recast is
  trusted). Empty means unconstrained.
- `invariants::Vector` — algebraic equalities `{hₖ(x) = 0}` from
  [`trig_recast`](@ref). These enter the SOS program with free-sign multipliers,
  not SOS ones. Getting this wrong makes the program infeasible in a way that
  looks like a solver problem.
- `state_labels::Vector{String}` — axis labels with units, for the plots.
- `notes::String` — one paragraph: where the model comes from, what was
  approximated, and what the certificate therefore does and does not mean.

# Constructor contract
The inner constructor validates: `length(f) == length(x)`, that every `f[i]` is
polynomial in `x` only, and that `‖f(equilibrium)‖∞ < 1e-12`. Raise `ArgumentError`
with the offending residual in the message.
"""
struct PolynomialSystem
    x::Vector
    f::Vector
    name::String
    equilibrium::Vector{Float64}
    constraints::Vector
    invariants::Vector
    state_labels::Vector{String}
    notes::String

    function PolynomialSystem(x, f, name, equilibrium, constraints, invariants, state_labels, notes)
        error("not implemented — IMPLEMENTATION_GUIDE.md §5.2")
    end
end

"""
    van_der_pol(; reversed = true, μ = 1.0)

The Van der Pol oscillator, by default **time-reversed**:

```
ẋ₁ = -x₂
ẋ₂ =  x₁ + μ (x₁² - 1) x₂
```

The reversed system is the canonical ROA benchmark because its origin is
asymptotically stable and the boundary of its region of attraction is exactly
the unstable limit cycle — which is the forward system's *stable* limit cycle,
and therefore computable to high accuracy by simply integrating the forward
system. That gives ground truth to compare an SOS inner approximation against,
which almost no other nonlinear benchmark offers.

With `reversed = false` you get the forward oscillator, whose origin is
*unstable*. It is provided for one purpose: a negative test. Every synthesis
routine must fail to certify it, and must fail by reporting infeasibility rather
than by returning a certificate with a tiny `ρ`.

This reproduces the setup in Tedrake's MIT 6.832 notes and in Packard/Topcu's
ROA papers; `matlab/van_der_pol_roa.m` is the independent YALMIP implementation
of the same problem and the two must agree.
"""
function van_der_pol(; reversed::Bool = true, μ::Real = 1.0)
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.2")
end

"""
    inverted_pendulum(; m = 1.0, l = 1.0, b = 0.1, g = 9.81, K = nothing, recast = :taylor, taylor_order = 3)

Damped pendulum about the **upright** equilibrium under a stabilising state
feedback `u = -K·[θ, θ̇]`, with `θ = 0` upright:

```
θ̈ = (g/l)·sin θ - (b/(m l²))·θ̇ + u/(m l²)
```

Two recastings, selected by `recast`:

- `:taylor` — `sin θ ≈ θ - θ³/6` (`taylor_order = 3`) or higher. Cheap, and the
  standard textbook treatment, but the resulting certificate is for the
  *truncated* system. The `notes` field must say so, and the example must show
  where the truncation and the true field diverge.
- `:trig` — exact, via [`trig_recast`](@ref): states become `[s, c, θ̇]` with
  `s² + c² = 1` carried as an invariant. Three states instead of two and a
  bigger SDP, but the certificate is about the real pendulum.

If `K === nothing`, compute the LQR gain for the linearisation with `Q = I`,
`R = 1` — but do it with an explicit Riccati solve written here, not by adding a
control-systems dependency.

**Ground-truth anchor:** the true ROA of the closed loop is bounded by the
stable manifold of the downward-swinging saddle, so any certified sublevel set
that contains `θ = ±π` is wrong. Assert this in the tests; it is the cheapest
possible check that the whole pipeline is sound.
"""
function inverted_pendulum(; m::Real = 1.0, l::Real = 1.0, b::Real = 0.1, g::Real = 9.81,
                           K = nothing, recast::Symbol = :taylor, taylor_order::Integer = 3)
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.2")
end

"""
    unicycle_with_obstacle(; obstacle_centre = (1.0, 0.0), obstacle_radius = 0.4, k_ρ = 1.0, k_α = 3.0, gains...)

Unicycle `(ẋ, ẏ, θ̇) = (v cos θ, v sin θ, ω)` under a smooth polynomial
stabilising feedback, with a circular obstacle expressed as the semialgebraic
constraint

```
g(x, y) = (x - x_o)² + (y - y_o)² - r² ≥ 0.
```

This is the example that exercises the constrained Positivstellensatz path: the
certificate is only required to hold on `{g ≥ 0}`, so a multiplier `σ(x)·g(x)`
enters the S-procedure. Without the constraint the problem is a plain ROA
computation and the obstacle is decoration — the point of the example is that
the certified set is provably obstacle-free, which is exactly the statement a
CBF-based controller needs as its terminal ingredient.

`sin θ`/`cos θ` are handled by [`trig_recast`](@ref) (exact), giving states
`[x, y, s, c]` with the invariant `s² + c² = 1`. Do not Taylor-expand the
heading here: the interesting states are precisely the ones with large `|θ|`.

The polar-tracking gains `k_ρ`, `k_α` follow the standard unicycle controller;
whichever form is used, the closed loop must be polynomial in the recast states
and the linearisation at the goal must be Hurwitz — verify that before running
any SOS program, because an unstable closed loop makes every downstream
infeasibility uninformative.
"""
function unicycle_with_obstacle(; obstacle_centre = (1.0, 0.0), obstacle_radius::Real = 0.4,
                                k_ρ::Real = 1.0, k_α::Real = 3.0)
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.2")
end

"""
    consensus_network(; n = 3, laplacian = nothing, coupling = :cubic, κ = 1.0)

`n`-agent consensus with nonlinear coupling,

```
ẋᵢ = -Σⱼ aᵢⱼ · φ(xᵢ - xⱼ),      φ(z) = z + κ z³   (`coupling = :cubic`)
```

on a connected graph with Laplacian `L`. Default: the 3-agent path graph.

**The trap, and the reason this example is here.** The consensus manifold
`span(𝟙)` is an equilibrium *set*, not a point: `L𝟙 = 0`. A Lyapunov function
that is positive definite on all of `ℝⁿ` cannot exist, and an SOS program that
asks for one returns infeasible with no hint as to why. The fix is to certify
the **disagreement dynamics**: pick an orthonormal basis `U ∈ ℝ^{n×(n-1)}` of
`𝟙^⊥`, change coordinates `z = U'x`, and run the whole pipeline on `ż`. The
origin in `z` is the consensus manifold in `x`, and it is an isolated
equilibrium there.

Return the system in `z` coordinates and store `U` in `notes` (and expose it via
the example) so the certified set can be mapped back to `x` for plotting.

**Citation hook.** A quadratic `V(z) = z'Pz` certified for the closed-loop
network is the single-mode special case of the composition argument used for
RRC-CBF in the thesis; the example is written to make that connection explicit
in prose, not just implied.
"""
function consensus_network(; n::Integer = 3, laplacian = nothing, coupling::Symbol = :cubic, κ::Real = 1.0)
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.2")
end

"""
    simulate(sys, x0, tspan; dt = 1e-3, max_norm = 1e6)

Fixed-step RK4 integration of `sys.f` from `x0`.

Deliberately hand-written rather than pulled from `DifferentialEquations.jl`:
the package must stay installable in a bare CI container, and every use here is
either a sanity check or a ground-truth boundary trace where a fixed step with a
stated `dt` is easier to defend than an adaptive solver's tolerances.

Terminates early and returns `(t, X, :diverged)` once `‖x‖∞ > max_norm`, which
is how the ROA ground-truth sweep detects escape without waiting out the full
`tspan`.
"""
function simulate(sys::PolynomialSystem, x0::AbstractVector, tspan::Tuple{<:Real,<:Real};
                  dt::Real = 1e-3, max_norm::Real = 1e6)
    error("not implemented — IMPLEMENTATION_GUIDE.md §5.2")
end
