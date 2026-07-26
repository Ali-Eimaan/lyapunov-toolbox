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
using LyapunovToolbox

function main()
    # ------------------------------------------------------ 1. the trap, shown once
    # TODO(impl): build the network in x coordinates and show that the Jacobian at
    # the consensus point has a zero eigenvalue with eigenvector 𝟙. Then attempt
    # the SOS program in x and record that it is infeasible.
    #
    # Keep this section. It costs one SDP and it is the difference between a
    # reader believing the change of coordinates is a convention and understanding
    # that it is necessary.

    # ------------------------------------------- 2. disagreement coordinates
    # TODO(impl): sys = consensus_network(n = 3, coupling = :cubic, κ = 1.0)
    # (returns the system already in z). Print U, verify Uᵀ𝟙 = 0 and UᵀU = I,
    # and assert the Jacobian of ż at the origin is Hurwitz.

    # -------------------------------------------------------- 3. certify
    # TODO(impl): quadratic V first (the case that connects to standard consensus
    # theory, where V = ½ zᵀz gives global exponential consensus for the LINEAR
    # protocol), then degree 4. Report both levels and both volumes.
    #
    # Comment on why the certified set is bounded at all: for κ > 0 the cubic term
    # is stabilising and the linear analysis suggests global convergence, so a
    # bounded certified set here reflects the LIMITS OF THE RELAXATION rather than
    # a real boundary. Say that explicitly — presenting a relaxation artefact as a
    # property of the system is the kind of error a reviewer will find.

    # -------------------------------------------------------- 4. κ sweep
    # TODO(impl): κ ∈ (0, 0.5, 1, 2, 5). Report ρ and volume for each. For κ < 0
    # (destabilising coupling) the certified set must shrink and eventually the
    # linearisation stops being Hurwitz; include one negative κ to show the
    # failure is reported as a raised ArgumentError rather than a tiny ρ.

    # ------------------------------------------------- 5. map back to x and simulate
    # TODO(impl): sample interior points of {V(z) ≤ ρ}, lift to x via x = Uz + α𝟙
    # for several α (the certificate is invariant to α — verify that numerically,
    # it is the whole justification for the coordinate change), and simulate.
    # Assert max pairwise |xᵢ - xⱼ| → 0 for every one.

    # --------------------------------------- 6. switching: the negative result
    # TODO(impl): certify two graphs separately — the path graph and the graph with
    # one edge removed (both connected). Then simulate a trajectory that switches
    # between them at times chosen adversarially (switch whenever the state is
    # near the boundary of the current mode's certified set) and show that the
    # switched trajectory can leave the union of the two certified sets.
    #
    # Report it as a counterexample with the switching times and the escaping
    # state. If you cannot find one at these parameters, say that too — "no
    # counterexample found in 10^4 randomized switching sequences" is an honest
    # and useful sentence. What is not acceptable is quietly dropping the section
    # because it did not produce a clean result.
    #
    # This is the citation hook: the transition-viable machinery in the thesis is
    # what closes exactly this gap.

    # ------------------------------------------------------- 7. verify + export
    # TODO(impl): verify_certificate on every certificate; write the 2-D
    # disagreement slice to docs/figures/data/consensus.csv.

    return 0
end

exit(main())
