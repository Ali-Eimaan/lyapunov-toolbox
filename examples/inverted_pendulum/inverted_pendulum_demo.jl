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
using LyapunovToolbox

function main()
    # ------------------------------------------------------- 1. two formulations
    # TODO(impl):
    #   sys_taylor = inverted_pendulum(recast = :taylor, taylor_order = 3)
    #   sys_trig   = inverted_pendulum(recast = :trig)
    # Print both state dimensions and the resulting Gram dimensions from
    # gram_dimension — the cost difference is part of the story.

    # ------------------------------------------------------- 2. certify each one
    # TODO(impl): synthesize_lyapunov on both with V_degree = 4. Record ρ, SDP
    # size, solve time, and the number of SDPs.

    # ------------------------------------------ 3. map both back to (θ, θ̇) and compare
    # TODO(impl): the trig certificate lives in (s, c, θ̇) on the manifold
    # s² + c² = 1. Parameterise it by θ to project onto (θ, θ̇), then compare
    # areas with the Taylor certificate. Both numbers get printed; do not report
    # only the larger one.

    # -------------------------------------------------- 4. the honest cross-check
    # TODO(impl): sample points inside the TAYLOR certified set and integrate the
    # TRUE trigonometric dynamics from each. Report:
    #   - the fraction that converge upright;
    #   - the smallest |θ| among those that do not.
    # If the fraction is 1, say so and note the largest |θ| that was actually
    # tested — "no counterexample below |θ| = 0.9 rad" is a claim; "the Taylor
    # certificate is valid" is not.
    #
    # Then do the same for the TRIG certified set, where the fraction must be 1
    # by construction. A failure there is a bug in the toolbox, not a modelling
    # gap — assert it rather than reporting it.

    # ---------------------------------------------------- 5. the saddle sanity check
    # TODO(impl): assert neither certified set contains (θ, θ̇) = (±π, 0). The
    # downward equilibrium is a saddle of the closed loop; a certified set
    # containing it is wrong. Cheapest possible end-to-end check — throw on failure.

    # ------------------------------------------------------- 6. Taylor order sweep
    # TODO(impl): taylor_order ∈ (3, 5, 7): report ρ, projected area, and the
    # converged fraction from step 4 for each. Higher order is a better model, so
    # the converged fraction should improve. Whether the area does is not obvious;
    # report it either way and comment on what you see rather than what you expected.

    # ----------------------------------------------------------- 7. verify + export
    # TODO(impl): verify_certificate on both; write boundary curves to
    # docs/figures/data/inverted_pendulum.csv for the notebook.

    return 0
end

exit(main())
