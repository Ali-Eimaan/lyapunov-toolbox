# Theory notes

The mathematics behind `lyapunov-toolbox`, in the same notation as the code.

Every symbol here has a code identifier; the mapping is in §1 and it is what makes
the repository navigable. If you rename something in the code, rename it here in
the same commit.

**Status: OUTLINE.** Every `- [ ]` below is a section to write. See
IMPLEMENTATION_GUIDE.md §10.

---

## 1. Notation and the code that implements it

- [ ] Table with four columns: symbol, meaning, code identifier, file.
      Populate it with at least: `f`, `x`, `V`, `V̇`, `ρ`, `ε`, `λ`, `σⱼ`, `μₖ`,
      `gⱼ`, `hₖ`, `Ω_ρ`, `Σ[x]` (the SOS cone), `z(x)` (the Gram basis), `Q`
      (the Gram matrix), `β`, `p(x)` (the shape set).
- [ ] State the standing assumptions once, numbered, so later sections can cite
      them: `f` polynomial, `f(0) = 0`, `A = ∂f/∂x(0)` Hurwitz.
- [ ] Say explicitly that everything is continuous-time. The discrete-time
      analogue is not implemented and pretending otherwise would be dishonest.

## 2. Lyapunov's direct method

- [ ] State the theorem being used, with hypotheses, for local asymptotic
      stability. Cite Khalil (Theorem 4.1) rather than reproving it.
- [ ] **The estimate-of-the-ROA corollary, which is the actual workhorse:** if
      `V̇ < 0` on `Ω_ρ \ {0}` and `Ω_ρ` is **bounded** and contains the origin,
      then `Ω_ρ` is positively invariant and contained in the ROA.
- [ ] Make the boundedness hypothesis loud. It is the one most often dropped, it
      is not implied by `V ≻ 0`, and dropping it makes the conclusion false —
      `V = x₁²` on a 2-D system is the counterexample to keep in front of the
      reader. Point at `is_sublevel_set_bounded`.

## 3. From "nonnegative" to "SOS"

- [ ] Definition of the SOS cone and the Gram parameterisation
      `p = z(x)ᵀ Q z(x)`, `Q ⪰ 0`.
- [ ] Why this is an SDP and not merely a convex problem: the map from `Q` to the
      coefficients of `p` is linear, so the feasible set is a spectrahedron.
- [ ] **SOS ⊊ nonnegative.** State Hilbert's classification (equality only for
      `n = 1`; `n = 2, d = 4`; and `n = 3, d = 4` in the homogeneous case) and
      give Motzkin's polynomial as the separating example. It is in the test
      suite; cite the test by name.
- [ ] What this costs us: an infeasible SOS program does **not** mean no
      Lyapunov function exists at that degree. Say that plainly — it is the most
      common overclaim in applied SOS papers.
- [ ] Size: Gram dimension `binomial(n + d, d)`; give the table for
      `n ∈ {2,3,4,6}`, `2d ∈ {2,4,6,8}`. Numbers from `gram_dimension`.

## 4. The S-procedure

- [ ] Statement for one quadratic constraint, then the polynomial
      generalisation used in the code.
- [ ] The program `-(V̇ + ε‖x‖²) + λ(V - ρ) ∈ Σ[x]`, `λ ∈ Σ[x]`. Derive why it
      implies `V̇ < 0` on `Ω_ρ`, in three lines. It really is three lines and
      writing them out is worth more than citing them.
- [ ] **Lossiness.** The S-procedure is sufficient, not necessary, beyond the
      single-quadratic case. Name what is lost.
- [ ] Degree bookkeeping: the pairing rule `deg(λ) + deg(V) ≥ deg(V̇)`, and what
      an infeasibility caused by violating it looks like (indistinguishable from
      a real one — which is why it is worth stating as a rule).
- [ ] Why `ε > 0` rather than `ε = 0`: SOS gives `≥ 0`, and the theorem needs
      strict negativity away from the origin.

## 5. Positivstellensatz certificates

- [ ] Putinar's Positivstellensatz, stated with the Archimedean hypothesis.
      Say whether it holds for each example in the repository — for the unicycle
      it needs a ball constraint to be added, and the code adds one; say so.
- [ ] The constrained program (P3) from `sos_lyapunov.jl`, with each multiplier
      explained.
- [ ] **Inequalities get SOS multipliers; equalities get free-sign ones.** Give
      the one-line reason. Note that an SOS multiplier on an equality is valid
      but strictly weaker, and cite the test that measures the difference.
- [ ] Schmüdgen vs. Putinar in one paragraph: why Putinar (linear in the
      generators) is what everyone uses and Schmüdgen (all products) is not.
- [ ] Forward reference: this is the machinery cited in the TVF paper for
      transition-viable funnel inner approximations.

## 6. Computing the region of attraction

- [ ] **ρ-bisection at fixed `V`.** State and prove the monotonicity that makes
      bisection valid — one line, and without it the method looks unjustified.
- [ ] **The V–s iteration.** Write out the three steps with their programs, and
      identify precisely which pairs of variables are bilinear and how each is
      broken (fix one and bisect the other).
- [ ] **The normalisation `ρ ≡ 1`.** Explain the scaling degeneracy `V ↦ V/ρ`,
      and what the V-step does without it (reports rescalings as improvements).
- [ ] **Shape sets.** Why `max ρ` is the wrong objective for the V-step and
      `max β` with `{p ≤ β} ⊆ {V ≤ 1}` is the right one.
- [ ] **No optimality guarantee.** State it here, once, in a form the notebooks
      and README can quote.
- [ ] Alternatives not implemented, one line each: occupation measures
      (Henrion–Korda), the Zubov PDE approach, sampling-based methods. Naming
      what you did not do is how a reader calibrates what you did.

## 7. Recasting non-polynomial dynamics

- [ ] **Taylor truncation.** What is lost, and the error bound that would be
      needed to transfer a certificate back to the true system — including the
      admission that the repository does not compute that bound.
- [ ] **Trigonometric recasting.** The exact construction, the invariant
      `s² + c² = 1`, and why the certificate is a statement about the manifold.
- [ ] **Sampling on a manifold.** Why a sample from the ambient box is not a
      state of the system, and how the code avoids it.
- [ ] Rational and other recastings: mention, do not implement.

## 8. Numerical conditioning

- [ ] Why SOS SDPs are badly conditioned: the monomial basis is close to
      linearly dependent at high degree.
- [ ] Scaling the state before solving, and why the certificate must be scaled
      back before it is reported.
- [ ] What `ALMOST_OPTIMAL` from a conic solver actually means, and the
      repository's policy on it (accept only after independent Gram verification —
      `verify_sublevel_set`).
- [ ] Why the certificate carries its Gram matrices: a certificate that cannot be
      re-verified without the solver is not a certificate.
- [ ] Honest statement on exactness: these are floating-point certificates. Cite
      the rational-rounding literature and state that it is not implemented.

## 9. The connection to the thesis

- [ ] One page. How an SOS inner approximation of a sublevel set becomes an inner
      approximation of a transition-viable funnel.
- [ ] Why per-mode certificates do not compose under switching, referencing the
      counterexample in `examples/multi_agent_consensus/`.
- [ ] What a set-valued reset map does to the argument, and precisely which piece
      of theory is missing — that gap is the thesis contribution and naming it
      sharply is more persuasive than gesturing at it.

## 10. References

- [ ] Parrilo (2000), thesis — SOS relaxations.
- [ ] Putinar (1993) — positivity on compact semialgebraic sets.
- [ ] Jarvis-Wloszek et al. (2003) / Topcu, Packard, Seiler (2008) — the V–s
      iteration as implemented here.
- [ ] Tedrake, *Underactuated Robotics*, ch. on Lyapunov analysis — the Van der
      Pol benchmark and the pedagogy this repository follows.
- [ ] Khalil, *Nonlinear Systems*, 3rd ed. — Theorem 4.1 and the ROA corollary.
- [ ] Blekherman, Parrilo, Thomas (eds.), *Semidefinite Optimization and Convex
      Algebraic Geometry* — the reference for §3.
- [ ] Lasserre / Henrion–Korda — the moment-SOS alternative named in §6.
- [ ] Löfberg (2004) — YALMIP; and the SumOfSquares.jl / JuMP papers.

Every entry needs a full citation, and every one must be a work actually consulted.
A padded bibliography is easy to spot and costs more credibility than it buys.
