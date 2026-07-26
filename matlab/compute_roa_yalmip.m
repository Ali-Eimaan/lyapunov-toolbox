function [rho, V, info] = compute_roa_yalmip(f, x, V0, opts)
%COMPUTE_ROA_YALMIP  SOS region-of-attraction estimate via YALMIP.
%
%   [RHO, V, INFO] = COMPUTE_ROA_YALMIP(F, X, V0, OPTS) returns the largest
%   level RHO for which the sublevel set {x : V(x) <= RHO} is certified by a
%   sum-of-squares argument to lie inside the region of attraction of the
%   origin of XDOT = F(X).
%
%   This is an INDEPENDENT implementation of the same mathematics as
%   julia/src/roa_estimation.jl, written against a different modelling layer
%   (YALMIP) and a different solver (Mosek). That independence is the entire
%   point: the two are cross-validated against each other in
%   tests/test_roa_estimation.jl, and an agreement between two implementations
%   sharing no code is evidence of correctness in a way that a repeated run of
%   the same code is not.
%
%   Do NOT refactor the two toward a shared formulation. If they disagree,
%   find out which is wrong; do not make them agree by construction.
%
%   INPUTS
%     F     - n-by-1 sdpvar polynomial vector field, F(0) = 0.
%     X     - n-by-1 sdpvar decision variables, in the same order as F.
%     V0    - candidate Lyapunov function (sdpvar polynomial). Typically
%             x'*P*x with P from LYAP on the linearisation.
%     OPTS  - struct, all fields optional:
%               .epsilon        positive-definiteness margin      (1e-6)
%               .lambdaDegree   degree of the SOS multiplier      (2)
%               .constraints    cell array of g_j(x) >= 0         ({})
%               .sigmaDegree    degree of the constraint multipliers (2)
%               .rhoMax         upper bracket for the bisection   (auto)
%               .bisectionTol   relative stopping width           (1e-4)
%               .maxIter        bisection iteration cap           (40)
%               .solver         'mosek' | 'sedumi' | 'sdpt3'      ('mosek')
%               .verbose        0 or 1                            (0)
%
%   OUTPUTS
%     RHO   - certified level (double). Always a VERIFIED feasible value, never
%             the midpoint of the final bracket.
%     V     - the Lyapunov function (unchanged from V0 in this routine; the
%             V-s iteration is Julia-side only, deliberately - see the guide).
%     INFO  - struct with .iterations, .solverTime, .solverStatus, .lambda,
%             .sigma, .gramV, .gramVdot, .config. INFO.config records every
%             OPTS field actually used, so a result can be reproduced from the
%             INFO struct alone.
%
%   THE PROGRAM SOLVED, for fixed rho:
%
%       lambda(x)                                   is SOS
%       -(Vdot(x) + epsilon*x'*x) + lambda(x)*(V(x) - rho)
%         - sum_j sigma_j(x)*g_j(x)                 is SOS,  sigma_j SOS
%
%   On {V <= rho} the term lambda*(V - rho) is nonpositive, hence
%   Vdot <= -epsilon*x'*x < 0 away from the origin. See
%   docs/derivations/positivstellensatz.tex.
%
%   REQUIREMENTS
%     YALMIP (tested with R20230622) and one SDP solver. Mosek needs a licence;
%     SeDuMi is the free fallback and is what the SOLVER field should be set to
%     on a machine without one. The routine must ERROR with an actionable
%     message when neither is on the path -- not fall through to a default that
%     produces a wrong answer quietly.
%
%   See also: VAN_DER_POL_ROA, PENDULUM_ROA.

% TODO(impl): IMPLEMENTATION_GUIDE.md §7.
%
% Implementation outline:
%
%   1. Validate inputs. length(f) == length(x); V0 vanishes at the origin
%      (replace(V0, x, zeros(n,1)) below 1e-12); every entry of f is polynomial
%      in x. Error with the offending value in the message.
%
%   2. Fill OPTS defaults into a local struct, then copy it into info.config.
%      Do not read OPTS again after this point -- one place where defaults are
%      resolved, so the recorded config cannot drift from what was used.
%
%   3. Vdot = jacobian(V0, x) * f.   One line, and the single most important
%      line in the file. Note YALMIP's jacobian returns a row vector, so this is
%      a row-times-column product; a transpose error here yields a scalar of the
%      wrong sign and a certificate that looks plausible.
%
%   4. Check that the sublevel set is bounded (top-degree form positive
%      definite). Refuse to return a certificate otherwise, matching the Julia
%      side's behaviour.
%
%   5. Bracket rho: if opts.rhoMax is empty, double from 1e-3 until the program
%      is infeasible, capped at 30 doublings; error if no infeasible level is
%      found (globally stable system -- see the Julia docstring for why that
%      must not be reported as a large finite rho).
%
%   6. Bisect. Inner feasibility solve, per candidate rho:
%
%        [lambda, cl] = polynomial(x, opts.lambdaDegree, 0);
%        F = [sos(lambda)];
%        expr = -(Vdot + opts.epsilon*(x'*x)) + lambda*(V0 - rho);
%        for each constraint g_j:
%            [sigma_j, cs_j] = polynomial(x, opts.sigmaDegree, 0);
%            F = [F, sos(sigma_j)];
%            expr = expr - sigma_j*g_j;
%        F = [F, sos(expr)];
%        sol = solvesos(F, [], sdpsettings('solver', opts.solver, ...
%                                          'verbose', opts.verbose), ...
%                       [cl; cs...]);
%
%      Feasible iff sol.problem == 0. Treat sol.problem == 4 (numerical
%      problems) as INFEASIBLE and record it in info.solverStatus -- accepting
%      it would inflate rho, which is the one direction of error this file must
%      not make. Any other nonzero sol.problem that is not 1 (infeasible) is an
%      error: rethrow it with yalmiperror(sol.problem).
%
%   7. Return the largest verified feasible rho, re-solving once at that value
%      to recover the multipliers and Gram matrices for INFO.
%
%   8. Extract Gram matrices with sosd/getbase and store them, so the MATLAB
%      result can be re-verified without MATLAB. The cross-validation fixture
%      written by van_der_pol_roa.m depends on this.

error('lyapunovToolbox:notImplemented', ...
      'compute_roa_yalmip is not implemented - see IMPLEMENTATION_GUIDE.md section 7');

end
