"""
    LyapunovToolbox

Sum-of-Squares Lyapunov synthesis and region-of-attraction (ROA) estimation for
polynomial and polynomial-recastable vector fields.

# What this package certifies

Given `ẋ = f(x)` with `f(0) = 0` and `f` polynomial, the package searches for a
polynomial `V` and a level `ρ > 0` such that

```
V(x) - ε‖x‖²             is SOS          (V positive definite)
-(V̇(x) + ε‖x‖²) + λ(x)·(V(x) - ρ)   is SOS,  λ SOS
```

The second line is the S-procedure statement of "`V̇ < 0` everywhere on
`{x : V(x) ≤ ρ}`". Together with `V ≻ 0` and boundedness of the sublevel set,
this certifies `Ω_ρ = {x : V(x) ≤ ρ}` is an invariant inner approximation of the
region of attraction of the origin. The mathematics is written out in
`docs/derivations/sos_decomposition.tex` and `docs/derivations/positivstellensatz.tex`.

**`Ω_ρ` must be bounded.** SOS says nothing about unboundedness, and an
unbounded sublevel set makes the invariance conclusion false. Every entry point
here checks boundedness of the sublevel set and refuses to return a certificate
without it.

# Module layout

| file | contents |
|---|---|
| `polynomial_utils.jl` | monomial bases, Lie derivatives, Taylor/trigonometric recasting, numerical helpers |
| `systems.jl` | the four benchmark systems, defined once and shared by tests, examples, and MATLAB cross-checks |
| `sos_lyapunov.jl` | the SOS programs themselves: positivity, S-procedure, Putinar multipliers, certificate objects |
| `roa_estimation.jl` | outer loops: ρ-bisection, V–s iteration, volume estimation, ground-truth comparison |

# Quick start

```julia
using LyapunovToolbox

sys  = van_der_pol()                       # reversed Van der Pol
V0   = lyapunov_from_linearization(sys)    # quadratic seed from the Jacobian
cert = estimate_roa(sys, V0)               # bisection on ρ
verify_certificate(cert)                   # independent numerical re-check
```

See `examples/` for the four worked demos and `docs/theory_notes.md` for the
theory in the same notation as the code.
"""
module LyapunovToolbox

using LinearAlgebra
using Printf
using Random

using DynamicPolynomials
using MultivariatePolynomials
using SemialgebraicSets
using SumOfSquares
using JuMP

# Order matters: polynomial_utils defines the primitives the rest use, systems
# defines the plants, sos_lyapunov the programs, roa_estimation the outer loops.
include("polynomial_utils.jl")
include("systems.jl")
include("sos_lyapunov.jl")
include("roa_estimation.jl")

export
    # polynomial_utils
    monomial_basis,
    sos_basis,
    gram_dimension,
    lie_derivative,
    jacobian_at,
    taylor_recast,
    trig_recast,
    prune_terms,
    evaluate_on_grid,
    gram_residual,
    is_sublevel_set_bounded,
    # systems
    PolynomialSystem,
    van_der_pol,
    inverted_pendulum,
    unicycle_with_obstacle,
    consensus_network,
    simulate,
    # sos_lyapunov
    SOSOptions,
    LyapunovCertificate,
    default_solver,
    is_sos,
    sos_decomposition,
    lyapunov_from_linearization,
    verify_sublevel_set,
    multiplier_step,
    v_step,
    beta_step,
    synthesize_lyapunov,
    verify_certificate,
    # roa_estimation
    ROAResult,
    estimate_roa,
    bisect_level_set,
    v_s_iteration,
    roa_volume,
    roa_boundary_2d,
    true_roa_boundary_2d,
    containment_fraction,
    benchmark_row

end # module
