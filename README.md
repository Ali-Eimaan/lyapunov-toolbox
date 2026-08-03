# lyapunov-toolbox

Sum-of-Squares Lyapunov function synthesis and region-of-attraction estimation in
Julia (SumOfSquares.jl), Python (bridge + plotting), and MATLAB (YALMIP+Mosek).

[![CI](https://github.com/ali-eiman/lyapunov-toolbox/actions/workflows/ci.yml/badge.svg)](https://github.com/ali-eiman/lyapunov-toolbox/actions/workflows/ci.yml)
[![Julia](https://img.shields.io/badge/julia-%3E%3D1.9-blue)](https://julialang.org)
[![Licence](https://img.shields.io/badge/licence-MIT-green)](LICENSE)

<table>
<tr>
  <td align="center"><strong>Van der Pol</strong><br>
    <em>ρ = 2.97, 78% of true ROA</em></td>
  <td align="center"><strong>Inverted Pendulum</strong><br>
    <em>ρ = 0.48, certified to |θ| ≤ 0.71 rad</em></td>
</tr><tr>
  <td align="center"><strong>Unicycle + Obstacle</strong><br>
    <em>ρ = 1.42, obstacle excluded via Putinar constraint</em></td>
  <td align="center"><strong>Multi-Agent Consensus</strong><br>
    <em>ρ = 2.14 in disagreement coordinates, N = 5 agents</em></td>
</tr>
</table>

## What this does

This toolbox synthesises **Lyapunov functions** for polynomial (and
polynomially-recast) nonlinear systems using sum-of-squares programming.
For a given polynomial vector field $\dot x = f(x)$ with a locally
asymptotically stable equilibrium at the origin, it produces:

1. a polynomial Lyapunov function $V(x)$,
2. a certified sublevel set $\Omega_\rho = \{x : V(x) \le \rho\}$ that is
   positively invariant and contained in the true region of attraction,
3. the Gram matrices and SOS multiplier certificates that let a third party
   **re-verify** the result without re-solving the SDP.

It is written for researchers who need inner approximations of ROAs as
building blocks for larger analysis problems (compositional stability of
hybrid systems, transition-viable funnels, safety verification).

## What is implemented

| Capability | Julia | MATLAB | Python |
|---|---|---|---|
| Fixed-$V$ $\rho$-bisection | ✅ | ✅ | via bridge |
| $V$--$s$ iteration | ✅ | ✅ | via bridge |
| Shape-set expansion ($\beta$-step) | ✅ | — | via bridge |
| Putinar constraints (obstacles, barriers) | ✅ | — | via bridge |
| Trigonometric recasting (pendulum) | ✅ | — | via bridge |
| Certificate verification (Gram check) | ✅ | ✅ | ✅ |
| Plotting (ROA, phase portrait, $\dot V$) | — | — | ✅ |
| Save/load certificates | ✅ | ✅ | ✅ |

## Install

```bash
# Julia dependencies
cd julia && julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Python bridge + plotting
cd python && pip install -e .

# Tell the Python bridge where Julia lives
export JULIA_PROJECT=$(pwd)/julia
```

## Quickstart

```python
from lyapunov_toolbox import get_session, certify_roa, synthesize

session = get_session()

# Van der Pol, μ = 1, degree-4 Lyapunov candidate
V, rho, certificate = synthesize(
    system="van_der_pol", mu=1.0, degree=4, epsilon=1e-4
)

print(f"Certified ρ = {rho:.3f}")
print(f"Gram residual = {certificate.gram_residual():.1e}")
certificate.save("vdp_cert.json")
```

## Results

| System | $n$ | Degree | $\rho$ | ROA fraction | Time (s) |
|---|---|---|---|---|---|
| Van der Pol ($\mu = 1$) | 2 | 4 | 2.97 | 78% of limit cycle | 0.8 |
| Van der Pol ($\mu = 1$) | 2 | 6 | 3.51 | 92% of limit cycle | 3.1 |
| Inverted pendulum (trig) | 2 | 4 | 0.48 | $|\theta| \le 0.71$ rad | 1.2 |
| Unicycle + obstacle | 3 | 4 | 1.42 | obstacle excluded | 4.7 |
| Consensus ($N = 5$, disc.) | 4 | 4 | 2.14 | $\kappa \ge -0.2$ | 2.3 |

The ROA fraction is the ratio of the certified area (Shoelace formula on
the $V = \rho$ boundary) to the true limit-cycle area for Van der Pol, and
the fraction of sampled interior points that converge under the true
$\sin\theta$ dynamics for the pendulum.

## Examples

| Example | Julia script | Python notebook |
|---|---|---|
| Van der Pol oscillator | `examples/van_der_pol/van_der_pol_demo.jl` | `examples/van_der_pol/van_der_pol_demo.ipynb` |
| Inverted pendulum (trig + Taylor) | `examples/inverted_pendulum/inverted_pendulum_demo.jl` | `examples/inverted_pendulum/inverted_pendulum_demo.ipynb` |
| Unicycle with circular obstacle | `examples/unicycle_with_obstacle/unicycle_demo.jl` | `examples/unicycle_with_obstacle/unicycle_demo.ipynb` |
| Multi-agent consensus (adversarial) | `examples/multi_agent_consensus/consensus_demo.jl` | `examples/multi_agent_consensus/consensus_demo.ipynb` |

## Theory

- **[`docs/theory_notes.md`](docs/theory_notes.md)** — notation, Lyapunov's
  direct method, SOS parameterisation, the $\mathcal{S}$-procedure, Putinar
  certificates, ROA computation, recasting, and numerical conditioning.
- **[`docs/derivations/sos_decomposition.pdf`](docs/derivations/sos_decomposition.pdf)** —
  Gram matrix parameterisation with a fully worked hand example, Motzkin's
  non-SOS polynomial, and numerical reality.
- **[`docs/derivations/positivstellensatz.pdf`](docs/derivations/positivstellensatz.pdf)** —
  Putinar's Positivstellensatz, the $\mathcal{S}$-procedure specialisation,
  the $V$--$s$ iteration, and the forward connection to transition-viable funnels.

## What this does not do

- **Inner approximations only.** Nothing here provides an outer bound on
  the ROA; the true ROA may be much larger than the certified set.
- **No exact certificates.** All results are floating-point; rational
  rounding (Peyrl--Parrilo, Powers--W\"ormann) is not implemented.
- **No robustness guarantees.** No robust or common-Lyapunov formulation
  exists in this codebase; additive disturbances, parameter mismatch, and
  model error invalidate the certificates.
- **The $V$--$s$ iteration is a local heuristic.** It has no convergence
  guarantee to a global optimum; different seeds give different fixed
  points (observed spread $\pm 6\%$ on Van der Pol, degree 4).
- **Degree and dimension limits.** The monomial basis condition number
  grows rapidly; beyond degree 8 in two variables or degree 4 in six
  variables the SDP becomes numerically unreliable.
- **Continuous-time only.** Discrete-time Lyapunov analysis is not
  implemented.

Naming the limits is what makes the claims above credible.

## Citing

This library produces the per-mode inner approximations that compose into
transition-viable funnels. To cite it in a paper:

> The SOS certificates and sublevel-set estimates were computed with
> `lyapunov-toolbox` (commit `abcdef`, https://github.com/ali-eiman/lyapunov-toolbox),
> which wraps SumOfSquares.jl with Clarabel.

Replace `abcdef` with the actual commit hash used.

## References

1. P. A. Parrilo, *Structured Semidefinite Programs and Semialgebraic Geometry
   Methods in Robustness and Optimization*, Ph.D. thesis, Caltech, 2000.
2. M. Putinar, "Positive polynomials on compact semi-algebraic sets," *Indiana
   Univ. Math. J.*, vol. 42, no. 3, pp. 969--984, 1993.
3. Z. Jarvis-Wloszek, R. Feeley, W. Tan, K. Sun, and A. Packard, "Some controls
   applications of sum of squares programming," in *Proc. IEEE CDC*, 2003.
4. U. Topcu, A. Packard, and P. Seiler, "Local stability analysis using
   simulations and sum-of-squares programming," *Automatica*, vol. 44, no. 10,
   2008.
5. H. K. Khalil, *Nonlinear Systems*, 3rd ed., Prentice Hall, 2002.
6. G. Blekherman, P. A. Parrilo, and R. R. Thomas (eds.), *Semidefinite
   Optimization and Convex Algebraic Geometry*, SIAM, 2013.
7. D. Henrion and M. Korda, "Convex computation of the region of attraction of
   polynomial control systems," *IEEE Trans. Autom. Control*, vol. 59, no. 2, 2014.
8. J. Löfberg, "YALMIP: A toolbox for modeling and optimization in MATLAB," in
   *Proc. CACSD*, 2004.

## Licence

MIT — see [LICENSE](LICENSE).
