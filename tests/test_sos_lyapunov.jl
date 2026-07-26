# Tests for julia/src/sos_lyapunov.jl — the SDP formulations.
#
# These tests call a conic solver, so they are the slow part of the suite.
# Budget: the unmarked tests must finish in under two minutes on a laptop.
# Anything slower goes behind `if get(ENV, "LYAPUNOV_SLOW_TESTS", "0") == "1"`.
#
# The governing principle here: never test an SOS routine by asking the SOS
# routine whether it succeeded. Every feasibility claim is cross-checked either
# against a hand-verifiable instance, against pointwise sampling, or against
# forward simulation of the dynamics.

using Test
using LinearAlgebra
using Random
using DynamicPolynomials
using LyapunovToolbox

@testset "sos_lyapunov" begin

    @testset "is_sos" begin
        @testset "accepts an obvious sum of squares" begin
            # (x1 + x2)^2 + x2^4. Trivially SOS by construction.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "rejects a polynomial that is negative somewhere" begin
            # p = x1^2 - 1 takes the value -1 at the origin.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "rejects the Motzkin polynomial" begin
            # M = x^2*y^4 + x^4*y^2 - 3*x^2*y^2*z^2 + z^6 is nonnegative
            # everywhere (AM–GM) but is NOT a sum of squares — the canonical
            # separating example. is_sos must return false, and must do so by
            # detecting infeasibility rather than by throwing or hanging.
            #
            # Also assert nonnegativity numerically at 10_000 random points, so
            # that the test documents what it is actually distinguishing:
            # "not SOS" is strictly stronger than "not nonnegative".
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "rejects an odd-degree polynomial without solving" begin
            # No odd-degree polynomial is SOS. Should short-circuit.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end
    end

    @testset "sos_decomposition" begin
        @testset "reconstructs the polynomial from its Gram matrix" begin
            # gram_residual(p, basis, Q) < 1e-7 and eigmin(Q) > -1e-8.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "Cholesky factor gives an explicit sum of squares" begin
            # Q = L*L' → p = Σ (L'z)_i^2. Verify at random points to 1e-7.
            # This is the object the .tex derivation walks through by hand.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "reports failure rather than clipping a negative eigenvalue" begin
            # Hand a non-SOS polynomial in; status must not be a success value.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end
    end

    @testset "lyapunov_from_linearization" begin
        @testset "solves the Lyapunov equation to residual 1e-9" begin
            # ‖A'P + PA + Q‖_F < 1e-9 * max(1, ‖P‖_F) for the reversed Van der Pol.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "returns a positive definite P" begin
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "V̇ is negative in a neighbourhood of the origin" begin
            # Sample points with ‖x‖ ≤ 1e-2 and assert V̇ < 0 at all of them.
            # This is the property that makes the ρ-bisection well-posed.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "throws for a non-Hurwitz linearisation" begin
            # The FORWARD Van der Pol: origin unstable. Must raise ArgumentError
            # naming the offending eigenvalues, not return a P.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end
    end

    @testset "verify_sublevel_set" begin
        @testset "feasible for a small level on the reversed Van der Pol" begin
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "infeasible for an absurdly large level" begin
            # ρ = 1e6 encloses the limit cycle; must report infeasible, not throw
            # and not time out.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "feasibility is monotone in rho" begin
            # THE test that justifies bisection. Sweep ρ over a geometric grid;
            # once feasibility flips to false it must never flip back. Report the
            # whole boolean sequence in the failure message.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "certified set really has V̇ < 0 throughout" begin
            # The anti-tautology test. Take a feasible (V, ρ), rejection-sample
            # 20_000 points from {V ≤ ρ}, and assert V̇ < 0 at every one with
            # ‖x‖ > 1e-6. Uses no solver. If the SDP is right and this fails, the
            # Lie derivative or the sign convention is wrong.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "does not mutate the options object" begin
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "throws on an unrecognised solver status" begin
            # Inject a stub optimizer that reports OTHER_ERROR. Silent conversion
            # of an unknown status to `infeasible` would shrink every reported ROA
            # with no visible symptom.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end
    end

    @testset "globally stable linear system" begin
        @testset "no finite level is ever infeasible" begin
            # ẋ = A x with A Hurwitz has ROA = R^n. With V from the Lyapunov
            # equation, verify_sublevel_set must be feasible at ρ = 1e8, and
            # bisect_level_set must therefore hit its doubling cap and THROW
            # rather than report a large finite ρ as if it were the boundary.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end
    end

    @testset "constrained certificates (Putinar)" begin
        @testset "obstacle constraint enlarges the feasible level set" begin
            # Same system, certified on {g ≥ 0} vs on all of R^n: the constrained
            # program is a relaxation, so its ρ must be ≥ the unconstrained one.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "certified set contains no point violating the constraint" begin
            # Sample {V ≤ ρ} and assert g(x) ≥ -1e-8 throughout. This is the
            # statement the unicycle example advertises; test it directly.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "equality invariants use free-sign multipliers" begin
            # For the trig-recast pendulum, assert the s^2 + c^2 - 1 multiplier is
            # NOT constrained to be SOS: build both programs and check the free
            # one certifies a level at least as large. An SOS multiplier on an
            # equality is a strictly weaker (though not incorrect) relaxation.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end
    end

    @testset "verify_certificate" begin
        @testset "passes on a genuine certificate" begin
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "fails when rho is inflated after the fact" begin
            # Take a valid certificate, multiply ρ by 5 without re-solving, and
            # assert verification fails on the sampling check. A verifier that
            # only re-reads the stored Gram matrices would pass this and is
            # therefore not a verifier.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "fails when V is replaced by an unrelated polynomial" begin
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "calls no solver" begin
            # Run with SOSOptions(solver = a factory that throws on construction).
            # Verification must complete.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end

        @testset "reports the three checks separately" begin
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end
    end

    @testset "determinism" begin
        @testset "two runs with the same seed agree to 1e-6" begin
            # Same options, same seed, same ρ. Conic solvers are deterministic
            # given identical input; if this fails, something is reading a global
            # RNG or iterating over an unordered collection.
            error("not implemented — IMPLEMENTATION_GUIDE.md §8")
        end
    end
end
