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
using LyapunovToolbox

function main()
    # ---------------------------------------------------------------- 1. system
    # TODO(impl): sys = van_der_pol(μ = 1.0)   (reversed by default)
    # Print sys.name, the state labels, and the Jacobian at the origin. Assert
    # the Jacobian is Hurwitz before anything else — every later infeasibility is
    # uninterpretable if it is not.

    # ------------------------------------------------- 2. quadratic V, fixed-V ρ
    # TODO(impl): V0 = lyapunov_from_linearization(sys)
    #             cert2 = estimate_roa(sys, V0)
    # This is the baseline. Print ρ, the SDP size, and the solve time.

    # ------------------------------------------------------ 3. ground truth ROA
    # TODO(impl): integrate the FORWARD oscillator from a point outside the limit
    # cycle, drop the transient, keep one period. That closed curve is the exact
    # ROA boundary of the reversed system.
    #
    # Prefer this to true_roa_boundary_2d's convergence sweep here, and say why
    # in the printed output: the limit cycle is exact up to integration error,
    # whereas the sweep is an empirical estimate at a finite horizon.

    # ----------------------------------------------------- 4. containment check
    # TODO(impl): assert every point of the certified boundary lies strictly
    # inside the limit cycle. A failure here means a false certificate — throw,
    # do not warn. Report the minimum clearance as a number, since a certificate
    # that only just fits is worth knowing about.

    # -------------------------------------------------- 5. degree sweep + V–s it.
    # TODO(impl): for V_degree in (2, 4, 6): run v_s_iteration and record
    # (ρ, certified area, area/true_area, n_sdps, solve time).
    # Assert the certified area is strictly increasing in the degree, and that
    # every certified set is contained in the true ROA. Print the table in the
    # exact markdown format docs/benchmarks.md uses, so it can be pasted.

    # ----------------------------------------------------------- 6. verification
    # TODO(impl): verify_certificate on every certificate produced above; assert
    # all three checks pass and print the worst sampled V̇ for each.

    # ------------------------------------------------------------- 7. the figure
    # TODO(impl): call the Python plotting layer, or write the boundary curves to
    # docs/figures/data/van_der_pol.csv for the notebook to plot. Pick one and be
    # consistent across all four examples — the guide picks the CSV route (§9.1)
    # so that the .jl scripts have no plotting dependency and stay CI-cheap.

    return 0
end

exit(main())
