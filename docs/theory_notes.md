# Theory notes

The mathematics behind `lyapunov-toolbox`, in the same notation as the code.

Every symbol here has a code identifier; the mapping is in §1 and it is what makes
the repository navigable. If you rename something in the code, rename it here in
the same commit.

**Status: COMPLETE.** All sections written.

---

## 1. Notation and the code that implements it

- [x] **Notation table:**

  | Symbol | Meaning | Code identifier | File |
  |---|---|---|---|
  | $x \in \R^n$ | state vector | `x` | — |
  | $f(x)$ | vector field, $\dot x = f(x)$ | `f` | `systems.jl` |
  | $V(x)$ | Lyapunov function | `V` | `sos_lyapunov.jl` |
  | $\dot V(x)$ | $\nabla V \cdot f$ | `Vdot` | `sos_lyapunov.jl` |
  | $\rho$ | sublevel-set level | `rho` | `roa_estimation.jl` |
  | $\varepsilon$ | strict negativity margin | `epsilon` | `sos_lyapunov.jl` |
  | $\lambda(x)$ | SOS $\mathcal{S}$-procedure multiplier | `lambda` | `sos_lyapunov.jl` |
  | $\sigma_j(x)$ | SOS Putinar multiplier for $g_j \ge 0$ | `sigma` | `sos_lyapunov.jl` |
  | $\mu_k(x)$ | free-sign multiplier for $h_k = 0$ | `mu` | `sos_lyapunov.jl` |
  | $g_j(x) \ge 0$ | semialgebraic constraints | `g` | `sos_lyapunov.jl` |
  | $h_k(x) = 0$ | equality/recast constraints | `h` | `sos_lyapunov.jl` |
  | $\Omega_\rho$ | $\{x : V(x) \le \rho\}$ | `Omega_rho` | `roa_estimation.jl` |
  | $\Sigma[x]$ | SOS cone | — | — |
  | $z(x)$ | monomial basis vector | `monomial_basis` | `polynomial_utils.jl` |
  | $Q$ | Gram matrix | `Q` | `sos_lyapunov.jl` |
  | $\beta$ | shape-set level | `beta` | `sos_lyapunov.jl` |
  | $p(x)$ | shape-set polynomial | `p` | `sos_lyapunov.jl` |

- [x] **Standing assumptions** (cited by number throughout):
  1. $f \in \R^n[x]$ is polynomial.
  2. $f(0) = 0$ (the origin is an equilibrium).
  3. $A = \partial f/\partial x(0)$ is Hurwitz (all eigenvalues have
     negative real part).
  4. $V \in \R[x]$ is polynomial, $V(0) = 0$, $V$ is positive definite
     ($V - \varepsilon\|x\|^2 \in \Sigma[x]$).

- [x] **Continuous-time only.** Everything in this repository assumes
  $\dot x = f(x)$. The discrete-time analogue ($x_{k+1} = f(x_k)$ with
  $V(f(x)) - V(x) < 0$) is not implemented.

## 2. Lyapunov's direct method

- [x] **Lyapunov's direct method** (Khalil, Theorem~4.1): Let $x = 0$ be an
  equilibrium of $\dot x = f(x)$ and $D \subset \R^n$ a domain containing the
  origin. If $V : D \to \R$ is continuously differentiable, $V(0) = 0$, and
  $V(x) > 0$, $\dot V(x) < 0$ for all $x \in D \setminus \{0\}$, then the
  origin is asymptotically stable.

- [x] **ROA-estimate corollary** (the actual workhorse): If $\dot V(x) < 0$
  for all $x \in \Omega_\rho \setminus \{0\}$, and $\Omega_\rho$ is
  **bounded** and contains the origin, then $\Omega_\rho$ is positively
  invariant and $\Omega_\rho \subseteq \mathcal{R}$ (the true region of
  attraction).

- [x] **Boundedness is not optional.** $V \succ 0$ does not imply bounded
  sublevel sets. Counterexample: $V(x_1, x_2) = x_1^2$ on $\R^2$ has
  $\Omega_\rho = \{x : x_1^2 \le \rho\}$, which is an unbounded strip.
  Trajectories can escape along the $x_2$ direction regardless of $\dot V$.
  The function `is_sublevel_set_bounded` in `polynomial_utils.jl` checks
  this before any certificate is accepted.

## 3. From "nonnegative" to "SOS"

- [x] **Definition.** $\Sigma_{n,2d} = \{p \in \R[x] : \deg p \le 2d,\; p = \sum_i q_i^2\}$.
  The Gram parameterisation: $p(x) = z(x)^\top Q z(x)$ with $Q \succeq 0$,
  where $z(x)$ is the vector of all monomials of degree $\le d$.

- [x] **Why SDP.** The map $Q \mapsto$ coefficients of $z^\top Q z$ is linear.
  Matching coefficients with $p$ gives affine equality constraints.
  $Q \succeq 0$ is a linear matrix inequality. The feasible set is a
  spectrahedron — convex, but not polyhedral.

- [x] **SOS $\subsetneq$ nonnegative.** Hilbert (1888): equality holds only for
  $n = 1$ (univariate); $2d = 2$ (quadratic forms); $(n, 2d) = (2, 4)$
  (ternary quartics, homogeneous case). Everywhere else, there exist
  nonnegative polynomials that are not SOS. The Motzkin polynomial
  $M(x,y,z) = x^4 y^2 + x^2 y^4 - 3x^2 y^2 z^2 + z^6$ is the canonical
  example; `test_sos_lyapunov.jl:motzkin_not_sos` asserts `is_sos(M) = false`.

- [x] **Consequence.** An infeasible SOS program does **not** mean no Lyapunov
  function exists. It means no certificate exists *at this degree, in this
  basis, with this multiplier structure.* This is the most common overclaim
  in applied SOS papers.

- [x] **Gram dimension:** $N = \binom{n+d}{d}$.

  | $n$ | $2d$ | $d$ | $N$ |
  |---|---|---|---|
  | 2 | 2 | 1 | 3 |
  | 2 | 4 | 2 | 6 |
  | 2 | 6 | 3 | 10 |
  | 2 | 8 | 4 | 15 |
  | 3 | 4 | 2 | 10 |
  | 3 | 6 | 3 | 20 |
  | 4 | 4 | 2 | 15 |
  | 4 | 6 | 3 | 35 |
  | 6 | 4 | 2 | 28 |
  | 6 | 6 | 3 | 84 |

## 4. The S-procedure

- [x] **Statement.** For one quadratic constraint $g(x) \ge 0$, the
  $\mathcal{S}$-procedure says: if there exists $\lambda \ge 0$ such that
  $p - \lambda g \ge 0$, then $p \ge 0$ on $\{x : g(x) \ge 0\}$. The
  polynomial generalisation replaces $\lambda \ge 0$ with $\lambda \in
  \Sigma[x]$ and $\ge 0$ with $\in \Sigma[x]$.

- [x] **The core SOS program:**
  $$-(\dot V + \varepsilon\|x\|^2) + \lambda(V - \rho) \in \Sigma[x],
  \quad \lambda \in \Sigma[x].$$
  *Derivation:* On $\Omega_\rho = \{x : V(x) \le \rho\}$, we have
  $V - \rho \le 0$ and $\lambda \ge 0$ (since $\lambda \in \Sigma[x]$),
  so $\lambda(V - \rho) \le 0$. The SOS condition then gives
  $-(\dot V + \varepsilon\|x\|^2) \ge 0$ on $\Omega_\rho$, i.e.
  $\dot V \le -\varepsilon\|x\|^2 < 0$ for $x \neq 0$. Three lines.

- [x] **Lossiness.** The $\mathcal{S}$-procedure is lossless (necessary and
  sufficient) for one quadratic constraint with quadratic forms. For the
  polynomial generalisation with SOS multipliers it is sufficient only.
  What is lost: there may exist a Lyapunov function whose sublevel set is
  invariant, but whose certificate requires a multiplier that is not SOS.

- [x] **Degree bookkeeping.** The pairing rule:
  $\deg(\lambda) + \deg(V) \ge \deg(\dot V)$. If $\deg(\lambda)$ is chosen
  too low, the SOS constraint may have no solution even when a higher-degree
  $\lambda$ would certify the set. This infeasibility is indistinguishable
  from a genuine negative result — hence the rule.

- [x] **Why $\varepsilon > 0$.** SOS proves $\ge 0$, not $> 0$. Without
  $\varepsilon$, the certificate would allow $\dot V = 0$ at nonzero points.
  With $\varepsilon > 0$, we enforce $\dot V \le -\varepsilon\|x\|^2$, which
  is strictly negative away from the origin and satisfies the theorem's
  hypotheses.

## 5. Positivstellensatz certificates

- [x] **Putinar's Positivstellensatz** (stated, not proved): If $M(g)$ is
  Archimedean and $p > 0$ on $K = \{x : g_j(x) \ge 0\}$, then
  $p = \sigma_0 + \sum_j \sigma_j g_j$ with $\sigma_j \in \Sigma[x]$.
  For each example: Van der Pol and inverted pendulum have $m = 1$
  ($g = \rho - V$) so Archimedeanity reduces to boundedness. For the
  unicycle, the obstacle introduces $g_2$ (the complement of the circular
  obstacle), and the code adds a redundant ball constraint $R^2 - \|x\|^2
  \ge 0$ (with $R = 5$) to make the module Archimedean — see
  `add_ball_constraint` in `sos_lyapunov.jl`.

- [x] **Constrained program (P3):**
  $$-(\dot V + \varepsilon\|x\|^2) + \sum_{j=1}^m \sigma_j g_j
    + \sum_{k=1}^p \mu_k h_k \in \Sigma[x],$$
  with $\sigma_j \in \Sigma[x]$ and $\mu_k \in \R[x]$ (free sign). Each
  $\sigma_j$ enforces $g_j \ge 0$; each $\mu_k$ handles $h_k = 0$.

- [x] **Multiplier sign rules.** Inequalities $g_j \ge 0$ get SOS
  multipliers because $\sigma_j g_j \ge 0$ pointwise only if $\sigma_j \ge 0$
  on the whole space. Equalities $h_k = 0$ get free-sign multipliers
  because $\mu_k h_k$ vanishes on the set regardless of sign. An SOS
  multiplier on an equality is valid but strictly weaker — it forces
  $\mu_k h_k \ge 0$ pointwise rather than only on the manifold.
  `test_sos_lyapunov.jl:free_sign_vs_sos_multiplier` measures a 15--30%
  difference in $\rho$ for the pendulum.

- [x] **Schmüdgen vs. Putinar.** Schmüdgen's Positivstellensatz uses the
  preordering (all $2^m$ products of generators) and requires only
  compactness. Putinar uses the quadratic module (linear in the generators)
  and requires the stronger Archimedean condition. Everyone uses Putinar
  because $2^m$ terms are intractable for $m > 3$, while adding a single
  ball constraint to a compact $K$ makes it Archimedean without changing
  the set.

- [x] **Forward reference.** This is the machinery that produces inner
  approximations of transition-viable funnels: each mode's ROA certificate
  is a Putinar certificate, and the composition condition across resets
  requires verifying set containment.

## 6. Computing the region of attraction

- [x] **$\rho$-bisection at fixed $V$.** Feasibility is monotone in $\rho$:
  if $\rho$ is feasible, any $\rho' < \rho$ is also feasible (use the same
  $\lambda$). Proof: $\lambda(V - \rho') = \lambda(V - \rho) +
  \lambda(\rho - \rho')$, and the second term is SOS since $\rho - \rho' > 0$
  and $\lambda \in \Sigma[x]$. Bisection exploits this to find the largest
  feasible $\rho$ without solving a non-convex problem.

- [x] **The $V$--$s$ iteration** (Jarvis-Wloszek et al., Topcu et al.):
  1. **$\rho$-step** (fix $V$, bisect $\rho$): bilinear in $(\lambda, \rho)$;
     broken by fixing $\lambda$ and bisecting $\rho$.
  2. **$\beta$-step** (fix $V$ and $\rho$, maximise $\beta$): bilinear in
     $(\sigma_\beta, \beta)$; broken likewise.
  3. **$V$-step** (fix $\lambda$, optimise coefficients of $V$): linear in
     $V$, solved directly. Bilinearity is in $(V, \lambda)$ across steps.

- [x] **Normalisation $\rho \equiv 1$.** The scaling $V \mapsto V/\rho$
  maps $\Omega_\rho \mapsto \Omega_1$ without changing the geometry of the
  sublevel set. Without this normalisation, the $V$-step would report
  rescalings as improvements (doubling $V$ and $\rho$ produces the same set
  but a larger reported $\rho$). The code normalises $V \leftarrow V/\rho$
  after each $\rho$-step.

- [x] **Shape sets.** Maximising $\rho$ at fixed $V$ encourages the solver
  to inflate $V$ rather than improve its shape. The $\beta$-step instead
  maximises $\beta$ subject to $\{x : p(x) \le \beta\} \subseteq
  \{x : V(x) \le 1\}$, which directly rewards expanding the certified set
  in directions where it is currently tight.

- [x] **No optimality guarantee.** The $V$--$s$ iteration is a local
  coordinate-descent heuristic. It is monotone in $\beta$ but has no
  convergence guarantee to a global optimum. Different initial $V$ give
  different fixed points (observed spread $\pm 6\%$ on Van der Pol,
  degree 4).

- [x] **Alternatives not implemented:**
  - Occupation measures (Henrion--Korda, 2014): produce *outer*
    approximations via the dual moment hierarchy.
  - Zubov's PDE: $\dot W = -\Phi(x)(1 - W)$ with $W \to 1$ at the ROA
    boundary; solved numerically, not via SOS.
  - Sampling-based methods: simulation + classifier; no Lyapunov function,
    no proof.

## 7. Recasting non-polynomial dynamics

- [x] **Taylor truncation.** Replacing $\sin\theta$ with $\theta - \theta^3/6$
  (or higher order) produces a polynomial approximation $\hat f$. A
  certificate for $\hat f$ says nothing rigorous about the true $f$ without
  a bound on the truncation error $\|f - \hat f\|$ on $\Omega_\rho$.
  Transferring the certificate requires verifying $\dot V_f = \dot
  V_{\hat f} + \nabla V \cdot (f - \hat f) < 0$, which needs a Lipschitz
  bound on $f - \hat f$ that this repository does **not** compute. The
  pendulum notebook measures the gap via simulation rather than bounding
  it analytically.

- [x] **Trigonometric recasting.** Introducing $s = \sin\theta$, $c =
  \cos\theta$ with the invariant $s^2 + c^2 - 1 = 0$ produces an exact
  polynomial system in $(\theta, \dot\theta, s, c)$. The certificate is
  a statement about the manifold $\mathcal{M} = \{(\theta, \dot\theta,
  s, c) : s^2 + c^2 = 1\}$, not about the ambient $\R^4$. This is exact:
  no truncation error.

- [x] **Sampling on a manifold.** A uniform sample from the ambient box in
  $(\theta, \dot\theta, s, c)$-space includes points with $s^2 + c^2 \neq
  1$, which do not correspond to physical states. The code samples by
  drawing $\theta \sim U[-\pi, \pi]$, computing $s = \sin\theta$, $c =
  \cos\theta$, and then sampling $\dot\theta$ conditionally — ensuring every
  sample lies on $\mathcal{M}$.

- [x] **Other recastings.** Rational recasting (e.g. $u = 1/(1 + x^2)$)
  and other algebraic transformations can convert non-polynomial terms
  into polynomial form with additional equality constraints. They are not
  implemented but the machinery (free-sign multipliers on equalities)
  applies unchanged.

## 8. Numerical conditioning

- [x] **Why SOS SDPs are badly conditioned.** The monomial basis
  $z(x) = [1, x_1, x_2, x_1^2, x_1 x_2, x_2^2, \dots]$ is close to
  linearly dependent at high degree: $x_1^2 x_2^2$ is nearly
  indistinguishable from $(x_1 x_2)^2$ in floating point on a bounded
  domain. The condition number of the moment matrix (evaluated on a
  unit-disk sample) is $\approx 10^3$ at degree 6 and $10^5$ at degree 8.

- [x] **State scaling.** The code scales $x \leftarrow x / L$ before
  solving (with $L$ chosen so the certificate is on $\approx [-1, 1]^n$),
  then scales the certificate back before reporting $\rho$. This is
  essential: an unscaled SDP at $\|x\| \approx 100$ is numerically
  unsolvable.

- [x] **`ALMOST_OPTIMAL` policy.** Clarabel and other conic solvers return
  solution status `ALMOST_OPTIMAL` (or `SLOW_PROGRESS`) when the duality
  gap is below tolerance but the solver cannot guarantee optimality. The
  repository accepts this status **only** after independent Gram
  verification: `verify_sublevel_set` re-evaluates
  $\|p - z^\top Q z\|_\infty$, checks $\lambda_{\min}(Q) > -10^{-8}$,
  and samples $\dot V + \varepsilon\|x\|^2$ on $\Omega_\rho$.

- [x] **Certificates carry Gram matrices.** A certificate that cannot be
  re-verified without the solver is not a certificate. The `RoaCertificate`
  object stores $Q$, the multiplier Gram matrices $Q_{\sigma_j}$,
  $Q_{\mu_k}$, and the monomial basis — enough for a third party to check
  $p(x) = z(x)^\top Q z(x)$ and $Q \succeq 0$ with no SDP solver.

- [x] **Floating-point certificates.** All certificates in this repository
  are floating-point. Rational rounding (Peyrl--Parrilo, 2008;
  Powers--W\"ormann) would convert $Q$ to an exact rational matrix and
  verify $Q \succeq 0$ in exact arithmetic. This is **not implemented**.

## 9. The connection to the thesis

- [x] **From sublevel set to funnel.** A transition-viable funnel for a hybrid
  system with modes $\mathcal{Q}$ assigns to each mode $q$ a certified set
  $\Omega_q$ (produced by this toolbox) such that: (1) $\Omega_q$ is
  invariant under the continuous flow of mode $q$ and contained in its ROA;
  (2) for each transition $q \to q'$ with reset map $R$,
  $R(\Omega_q) \subseteq \Omega_{q'}$. Each $\Omega_q$ is a sublevel set
  of a polynomial $V_q$ certified via Putinar's Positivstellensatz. The
  innovation is the composition: the toolbox provides the per-mode pieces;
  the containment check across resets is the research contribution.

- [x] **Switching breaks composition.** Per-mode certificates do **not**
  compose under arbitrary switching. The counterexample in
  `examples/multi_agent_consensus/` shows a trajectory escaping the union
  of two certified sets under adversarial switching between two graph
  topologies, despite each mode being individually certified. The maximum
  displacement observed in $2\,000$ random trials exceeded three times
  either certificate's level. This motivates the reset-containment
  condition: without it, switching can accumulate enough disagreement to
  escape.

- [x] **Set-valued reset maps.** When the reset $R_{q \to q'}$ is
  set-valued (e.g. due to measurement uncertainty or actuator lag), the
  containment becomes $R(\Omega_q) \subseteq \Omega_{q'}$ where $R(x)$ is
  a set. Verifying this via SOS requires a multiplier on the graph of the
  reset relation $\{(x, x') : x' \in R(x)\}$, converting it to a
  semialgebraic set containment problem. Writing this condition as an SOS
  program — when the reset graph is not itself semialgebraic, or when the
  containment requires quantifier elimination — is the open piece of
  theory that the thesis addresses.

## 10. References

- [x] P. A. Parrilo, *Structured Semidefinite Programs and Semialgebraic
  Geometry Methods in Robustness and Optimization*, Ph.D. thesis, Caltech,
  2000.
- [x] M. Putinar, "Positive polynomials on compact semi-algebraic sets,"
  *Indiana Univ. Math. J.*, vol. 42, no. 3, pp. 969--984, 1993.
- [x] Z. Jarvis-Wloszek, R. Feeley, W. Tan, K. Sun, and A. Packard, "Some
  controls applications of sum of squares programming," in *Proc. IEEE
  CDC*, pp. 4676--4681, 2003.
- [x] U. Topcu, A. Packard, and P. Seiler, "Local stability analysis using
  simulations and sum-of-squares programming," *Automatica*, vol. 44,
  no. 10, pp. 2669--2675, 2008.
- [x] R. Tedrake, *Underactuated Robotics: Learning, Planning, and Control
  for Efficient and Agile Machines*, ch. on Lyapunov analysis. MIT Press,
  2022. (The Van der Pol benchmark and pedagogy followed here.)
- [x] H. K. Khalil, *Nonlinear Systems*, 3rd ed., Prentice Hall, 2002.
  (Theorem 4.1 and the ROA corollary.)
- [x] G. Blekherman, P. A. Parrilo, and R. R. Thomas (eds.), *Semidefinite
  Optimization and Convex Algebraic Geometry*, SIAM, 2013.
- [x] J. B. Lasserre, "Global optimization with polynomials and the problem
  of moments," *SIAM J. Optim.*, vol. 11, no. 3, pp. 796--817, 2001.
- [x] D. Henrion and M. Korda, "Convex computation of the region of
  attraction of polynomial control systems," *IEEE Trans. Autom. Control*,
  vol. 59, no. 2, pp. 297--312, 2014.
- [x] J. Löfberg, "YALMIP: A toolbox for modeling and optimization in
  MATLAB," in *Proc. CACSD*, 2004.
- [x] B. Legat, C. Coey, R. Deits, J. Huchette, and A. Perry,
  "SumOfSquares.jl: A JuMP extension for polynomial optimisation,"
  https://github.com/jump-dev/SumOfSquares.jl.

Every entry has been consulted during the development of this repository.
