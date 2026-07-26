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
using LyapunovToolbox

function main()
    # ------------------------------------------------------------- 1. the system
    # TODO(impl): sys = unicycle_with_obstacle(obstacle_centre = (1.0, 0.0),
    #                                          obstacle_radius = 0.4)
    # Verify the closed-loop linearisation at the goal is Hurwitz first. Print
    # the recast state list and the invariant polynomial.

    # ------------------------------------------- 2. unconstrained vs. constrained
    # TODO(impl): certify twice — once ignoring the obstacle, once with it as a
    # Putinar constraint. Assert ρ_constrained ≥ ρ_unconstrained and report both
    # plus the ratio. Equality to within solver noise means the constraint is not
    # binding at this level, which is a legitimate outcome but must be stated
    # rather than presented as a gain.

    # ------------------------------------------------- 3. free vs. SOS invariant
    # TODO(impl): run the trig invariant both ways (free-sign μ, and σ SOS).
    # Report both levels. The free-sign version should be at least as large; if
    # the SOS version wins, the free-sign multiplier is not being wired in.

    # ----------------------------------------------------- 4. obstacle-free check
    # TODO(impl): rejection-sample 50_000 points from the certified set —
    # ON THE MANIFOLD s² + c² = 1, not from the ambient 4-D box; a sample off the
    # manifold is not a state of the unicycle and testing anything there proves
    # nothing. Assert min g(x) ≥ -1e-8 and report the minimum clearance in metres.

    # -------------------------------------------------------- 5. shape-set sweep
    # TODO(impl): compare p(x) = x'x against an anisotropic shape matched to the
    # geometry (the certified set is long in the goal-approach direction and thin
    # across it). Report both β values and both volumes. State which shape was
    # used for the figure — a volume comparison across different shape sets is
    # not a comparison.

    # ---------------------------------------------------------- 6. trajectories
    # TODO(impl): simulate from 20 sampled interior points. Assert every
    # trajectory reaches the goal (‖x‖ < 1e-3) and never enters the obstacle.
    # This is the closed-loop statement the certificate is supposed to imply; if
    # simulation contradicts it, the certificate is wrong and the script must fail.

    # ------------------------------------------------------- 7. verify + export
    # TODO(impl): verify_certificate; write the (x, y) slice at θ = 0 and at
    # θ = π/4 to docs/figures/data/unicycle.csv. Two slices, because a single
    # slice of a 4-D set invites exactly the overstatement this repository is
    # trying to avoid.

    return 0
end

exit(main())
