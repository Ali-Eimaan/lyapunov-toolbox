function results = van_der_pol_roa(varargin)
%VAN_DER_POL_ROA  Canonical Van der Pol ROA benchmark (Tedrake, MIT 6.832).
%
%   RESULTS = VAN_DER_POL_ROA() computes an SOS inner approximation of the
%   region of attraction of the origin of the TIME-REVERSED Van der Pol
%   oscillator
%
%       x1dot = -x2
%       x2dot =  x1 + mu*(x1^2 - 1)*x2
%
%   and compares it against the exact ROA boundary, which for this system is
%   the forward oscillator's stable limit cycle.
%
%   THIS IS THE SANITY CHECK FOR THE WHOLE TOOLBOX. It reproduces a result that
%   is published, widely replicated, and independently checkable, using a solver
%   stack (YALMIP + Mosek) that shares no code with the Julia implementation. If
%   this script and julia/ disagree, one of them is wrong, and the disagreement
%   is worth more than either number alone.
%
%   RESULTS = VAN_DER_POL_ROA('Name', Value, ...) accepts:
%       'Mu'            oscillator parameter                       (1.0)
%       'Epsilon'       positive-definiteness margin               (1e-6)
%       'LambdaDegree'  SOS multiplier degree                      (2)
%       'Solver'        'mosek' | 'sedumi' | 'sdpt3'               ('mosek')
%       'Plot'          draw the comparison figure                 (true)
%       'WriteFixture'  write the cross-validation CSV             (false)
%
%   RESULTS is a struct with fields RHO, V, P, AREA, TRUEAREA, AREARATIO,
%   SOLVERTIME, CONFIG.
%
%   CROSS-VALIDATION FIXTURE
%   With 'WriteFixture' true the script writes tests/fixtures/matlab_vdp.csv,
%   which tests/test_roa_estimation.jl reads to assert Julia/MATLAB agreement to
%   1e-3 relative. Regenerating that fixture is a manual, release-time step:
%   CI has neither MATLAB nor Mosek, so the fixture is committed and the commit
%   message must record the MATLAB and Mosek versions that produced it.
%
%   See also: COMPUTE_ROA_YALMIP, PENDULUM_ROA.

% TODO(impl): IMPLEMENTATION_GUIDE.md §7.
%
% Steps:
%
%   1. Parse the name/value options with inputParser and record them in
%      results.config.
%
%   2. Declare the dynamics.
%
%        sdpvar x1 x2
%        x = [x1; x2];
%        f = [-x2; x1 + mu*(x1^2 - 1)*x2];
%
%      NOTE: this restates julia/src/systems.jl's van_der_pol. That duplication
%      is unavoidable across the language boundary and is exactly what the
%      cross-validation fixture guards. If you change one, change both in the
%      same commit and re-run this script.
%
%   3. Seed candidate from the linearisation. A = [0 -1; 1 -mu] at the origin.
%      Solve the continuous Lyapunov equation A'P + PA = -I with LYAP (which
%      uses that sign convention -- verify the residual rather than trusting it,
%      because the sign convention of lyap is a recurring source of error) and
%      set V0 = x'*P*x.
%
%   4. rho = compute_roa_yalmip(f, x, V0, opts).
%
%   5. Ground truth. Integrate the FORWARD oscillator
%        [-x2; x1 + mu*(x1^2-1)*x2]  with the sign of time flipped
%      i.e. ode45 on [x2; -x1 - mu*(x1^2-1)*x2], from a point outside the limit
%      cycle, discard the transient, and keep one period. That closed curve is
%      the exact ROA boundary of the reversed system. Using an exact boundary
%      rather than a sampled convergence test is what makes the area ratio here
%      trustworthy.
%
%   6. Areas by polyarea on both closed curves; report the ratio.
%
%   7. ASSERT CONTAINMENT: every point of the certified boundary must lie
%      strictly inside the limit cycle (inpolygon). If it does not, the
%      certificate is false -- ERROR, do not warn and do not plot. A figure
%      showing a certified set spilling outside the true ROA, published in a
%      README, is worse than no figure.
%
%   8. Plot: limit cycle dashed, certified sublevel set filled, a few
%      trajectories from inside and outside. Axis labels 'x_1', 'x_2', legend
%      with the certified rho and the area ratio to three significant figures.
%
%   9. Print a summary block with rho, both areas, the ratio, the solver, the
%      solve time, and the MATLAB/YALMIP/solver versions. That block is pasted
%      verbatim into docs/benchmarks.md -- it is the provenance for a number CI
%      cannot regenerate.

error('lyapunovToolbox:notImplemented', ...
      'van_der_pol_roa is not implemented - see IMPLEMENTATION_GUIDE.md section 7');

end
