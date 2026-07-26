function results = pendulum_roa(varargin)
%PENDULUM_ROA  ROA of the upright equilibrium of a damped, LQR-stabilised pendulum.
%
%   RESULTS = PENDULUM_ROA() computes an SOS inner approximation of the region
%   of attraction of theta = 0 (upright) for
%
%       thetaddot = (g/l)*sin(theta) - (b/(m*l^2))*thetadot + u/(m*l^2),
%       u = -K*[theta; thetadot],
%
%   under a cubic Taylor recasting of sin(theta), and compares the certified
%   set against the true ROA obtained by simulation.
%
%   WHAT THIS EXAMPLE IS FOR, beyond a second data point: it is the one where
%   the model used for certification is NOT the model of the physical system.
%   sin(theta) ~ theta - theta^3/6 is a different vector field, so a certificate
%   proved for it is a statement about the truncated system. The script makes
%   that gap visible instead of glossing it:
%
%     - it certifies on the cubic model;
%     - it then simulates the TRUE trigonometric dynamics from points inside the
%       certified set;
%     - and it reports the fraction that actually converge upright.
%
%   That fraction should be 1 for a small certified set and will start to drop
%   as the set grows toward |theta| where the truncation error becomes
%   comparable to the field itself. Report where it drops. Reporting the
%   breakdown point is the difference between a demonstration and a claim.
%
%   RESULTS = PENDULUM_ROA('Name', Value, ...) accepts:
%       'M','L','B','G'    physical parameters      (1.0, 1.0, 0.1, 9.81)
%       'TaylorOrder'      3, 5, or 7               (3)
%       'Epsilon'          margin                   (1e-6)
%       'LambdaDegree'     multiplier degree        (2)
%       'Solver'           SDP solver               ('mosek')
%       'Plot'             draw the figure          (true)
%
%   RESULTS is a struct with RHO, V, P, K, TRUEFRACTION, BREAKDOWNTHETA,
%   SOLVERTIME, CONFIG.
%
%   See also: COMPUTE_ROA_YALMIP, VAN_DER_POL_ROA.

% TODO(impl): IMPLEMENTATION_GUIDE.md §7.
%
% Steps:
%
%   1. Parse options; record them in results.config.
%
%   2. Linearise about upright: A = [0 1; g/l, -b/(m*l^2)],
%      B = [0; 1/(m*l^2)]. Note A(2,1) is +g/l, not -g/l: upright is UNSTABLE,
%      the eigenvalues of A straddle the imaginary axis. If you get a Hurwitz A
%      before applying feedback, the sign is wrong.
%
%   3. K = lqr(A, B, diag([10 1]), 1). Verify eig(A - B*K) is Hurwitz before
%      going further; an SOS infeasibility caused by a non-stabilising gain is
%      indistinguishable from one caused by a modelling error.
%
%   4. Build the closed-loop polynomial field with the Taylor recast:
%
%        sdpvar th thd
%        x = [th; thd];
%        sinApprox = taylor_sin(th, TaylorOrder);   % local helper, below
%        u = -K*x;
%        f = [thd;
%             (g/l)*sinApprox - (b/(m*l^2))*thd + u/(m*l^2)];
%
%   5. V0 from lyap on (A - B*K) with Q = I; rho from compute_roa_yalmip.
%
%   6. TRUE ROA by simulation of the exact dynamics (sin, not the truncation):
%      grid theta in [-pi, pi] and thetadot in [-10, 10], integrate each with
%      ode45 for 20 s, mark convergence as final |x| < 1e-3. This is an
%      empirical estimate at a stated grid and horizon -- label it as such
%      everywhere it appears.
%
%   7. ASSERT the certified set excludes theta = +-pi. The downward equilibrium
%      is a saddle of the closed loop and cannot be in the ROA of the upright
%      one; a certified set containing it is wrong whatever the SDP reported.
%      This is the cheapest end-to-end soundness check available and it must
%      ERROR, not warn.
%
%   8. Simulate the TRUE dynamics from points sampled inside the certified set;
%      report the converged fraction and the smallest |theta| at which a
%      certified point fails. Both numbers go in the printed summary and in
%      docs/benchmarks.md.
%
%   9. Plot: certified sublevel set, empirical ROA boundary, the two separatrix
%      branches through (+-pi, 0), and the sampled trajectories. Label axes
%      'theta [rad]' and 'thetadot [rad/s]'.
%
%  10. Sweep TaylorOrder over [3 5 7] and report rho and the converged fraction
%      for each. A higher-order recast is a better model, so the converged
%      fraction should improve; whether rho does is not obvious and is worth
%      reporting either way.

% -------------------------------------------------------------------------
% function s = taylor_sin(th, order)
% TODO(impl): local helper. Odd orders only -- assert mod(order, 2) == 1 and
% error otherwise, since an even-order truncation of an odd function silently
% adds nothing and wastes a degree.
% -------------------------------------------------------------------------

error('lyapunovToolbox:notImplemented', ...
      'pendulum_roa is not implemented - see IMPLEMENTATION_GUIDE.md section 7');

end
