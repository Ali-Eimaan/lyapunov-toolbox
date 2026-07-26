<!--
STUB. Rewrite this file LAST, once the numbers exist — see
IMPLEMENTATION_GUIDE.md §11. Every heading below is a placeholder, and every
number that appears in the finished README must have been produced by code in
this repository. Delete this comment when you rewrite.
-->

# lyapunov-toolbox

Sum-of-Squares Lyapunov function synthesis and region-of-attraction estimation in
Julia (SumOfSquares.jl), Python (bridge + plotting), and MATLAB (YALMIP+Mosek).

<!-- TODO(impl): CI badge, coverage badge, licence badge, Julia version badge. -->

<!-- TODO(impl): the four hero figures from docs/figures/, as a 2x2 table:
     roa_van_der_pol.png, roa_inverted_pendulum.png,
     roa_unicycle_obstacle.png, roa_consensus.png.
     Each caption states the certified ρ and the fraction of the true ROA. -->

## What this does

<!-- TODO(impl): four or five sentences. What an SOS certificate is, what the
     toolbox produces (a Lyapunov function, a level, AND the Gram matrices that
     let a third party re-verify it without re-solving), and who it is for. -->

## What is implemented

<!-- TODO(impl): short table.
     | | Julia | MATLAB | Python |
     rows: fixed-V ρ bisection, V–s iteration, shape sets, Putinar constraints,
     trigonometric recasting, certificate verification, plotting. -->

## Install

<!-- TODO(impl): the Julia instantiate command, the pip install, and the one
     environment variable the bridge reads. Three commands, no prose. -->

## Quickstart

<!-- TODO(impl): a ten-line snippet that certifies the Van der Pol ROA and prints
     ρ. It must be copy-pasteable and must actually run — test it in a clean
     checkout before committing. -->

## Results

<!-- TODO(impl): the headline table from docs/benchmarks.md §1, abridged to the
     four systems. Quote the fraction of the true ROA you actually measured. If
     it is 78%, write 78%. -->

## Examples

<!-- TODO(impl): one line per example directory, linking both the .jl script and
     the notebook. -->

## Theory

<!-- TODO(impl): link docs/theory_notes.md and the two compiled derivations, with
     one line each on what they contain. -->

## What this does not do

<!-- TODO(impl): keep this section and make it specific. At minimum:
     - inner approximations only — nothing here bounds the ROA from outside;
     - no rational/exact certificates, these are floating-point;
     - no robust or common-Lyapunov formulation, so no disturbance guarantees;
     - the V–s iteration is a local heuristic with no optimality guarantee;
     - degree/dimension limits, with the measured wall where it stops being
       practical;
     - discrete-time is not implemented.
     Naming the limits is what makes the claims above credible. -->

## Citing

<!-- TODO(impl): the citation hook. This library is what produces the inner
     approximations of transition-viable funnels in `transition-viable-swarm`;
     give the exact sentence and the commit-pinned form to use in a paper. -->

## References

<!-- TODO(impl): matching docs/theory_notes.md §10. -->

## Licence

MIT — see [LICENSE](LICENSE).
