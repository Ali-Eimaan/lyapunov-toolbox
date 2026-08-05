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
%   MINCLEARANCE, PERIOD, SOLVERTIME, INFO, CONFIG.
%
%   GROUND TRUTH, AND WHY IT IS EXACT HERE
%   The reversed system's ROA boundary IS the forward system's stable limit
%   cycle, so it is computable to integrator accuracy rather than estimated by
%   a convergence sweep. Extracting it correctly matters: the cycle must be
%   exactly ONE period. A trace spanning several periods is a self-overlapping
%   vertex list, and both POLYAREA (shoelace, whose signed loop contributions do
%   not simply add) and INPOLYGON (even-odd winding) then return nonsense --
%   quietly, and in either direction. This script locates the period with a
%   Poincare section and verifies the orbit has actually settled before using it.
%
%   CROSS-VALIDATION FIXTURE
%   With 'WriteFixture' true the script writes tests/fixtures/matlab_vdp.csv via
%   WRITE_ROA_FIXTURE, which tests/test_roa_estimation.jl reads to assert
%   Julia/MATLAB agreement. Regenerating it is a manual, release-time step: CI
%   has neither MATLAB nor Mosek, so the fixture is committed and the commit
%   message must record the MATLAB, YALMIP, and solver versions that produced it.
%
%   See also: COMPUTE_ROA_YALMIP, PENDULUM_ROA, WRITE_ROA_FIXTURE.

% --- Implementation ---
results = struct();

% 1. Parse options.
p = inputParser;
addParameter(p, 'Mu', 1.0, @(v) isscalar(v) && v > 0);
addParameter(p, 'Epsilon', 1e-6, @(v) isscalar(v) && v > 0);
addParameter(p, 'LambdaDegree', 2, @(v) isscalar(v) && v >= 0 && mod(v,2)==0);
addParameter(p, 'Solver', 'mosek', @(v) ismember(v, {'mosek','sedumi','sdpt3'}));
addParameter(p, 'Plot', true, @islogical);
addParameter(p, 'WriteFixture', false, @islogical);
parse(p, varargin{:});
config = p.Results;
mu = config.Mu;

% 1b. Check the environment BEFORE the first SDPVAR. Without this the script
%     dies on "Unrecognized function or variable 'sdpvar'", which tells a
%     first-time user nothing about what to install.
config.Solver = assert_yalmip_env(config.Solver);
results.config = config;

% 2. Declare dynamics (reversed Van der Pol).
%    NOTE: this restates van_der_pol() in julia/src/systems.jl. The duplication
%    is unavoidable across the language boundary and is exactly what the
%    cross-validation fixture guards. Change one, change both in the same
%    commit, and re-run this script.
% Explicit assignment, not command syntax -- see the note in pendulum_roa.m:
% command-syntax sdpvar creates names at run time, but MATLAB binds
% variable-vs-function at parse time, so any name shadowing a toolbox function
% is parsed as a call. x1/x2 happen to be safe today; this does not rely on it.
x1 = sdpvar(1, 1);
x2 = sdpvar(1, 1);
x = [x1; x2];
f = [-x2; x1 + mu*(x1^2 - 1)*x2];

% 3. Seed from the linearisation. A = [0 -1; 1 -mu] is Hurwitz for mu > 0.
%    MATLAB's LYAP(A, Q) solves A*X + X*A' + Q = 0, so LYAP(A', I) gives
%    A'P + PA = -I. Verify the residual rather than trusting the convention:
%    the sign of lyap is a recurring source of error.
A = [0 -1; 1 -mu];
assert(all(real(eig(A)) < 0), ...
       'A is not Hurwitz for mu = %g; the reversed origin should be stable.', mu);
P = lyap(A', eye(2));
residual = A'*P + P*A + eye(2);
assert(norm(residual) < 1e-9 * max(1, norm(P)), ...
       'Lyapunov equation residual too large: %e', norm(residual));
V0 = x'*P*x;

% 4. Compute the certified level.
opts = struct('epsilon', config.Epsilon, ...
              'lambdaDegree', config.LambdaDegree, ...
              'solver', config.Solver);
[rho, V, info] = compute_roa_yalmip(f, x, V0, opts);
results.rho        = rho;
results.V          = V;
results.P          = P;
results.info       = info;
results.solverTime = info.solverTime;

% 5. Ground truth: exactly ONE period of the forward system's limit cycle.
[cycle, period] = limit_cycle_one_period(mu);
results.period   = period;
results.trueArea = polyarea(cycle(:,1), cycle(:,2));

% 6. Certified boundary by ray-marching from the origin (the equilibrium).
n_rays   = 400;
angles   = (0:n_rays-1)' * (2*pi / n_rays);   % no duplicated endpoint
boundary = zeros(n_rays, 2);
for i = 1:n_rays
    d  = [cos(angles(i)); sin(angles(i))];
    lo = 0; hi = 100;
    for k = 1:60
        mid = (lo + hi) / 2;
        if double(replace(V, x, mid * d)) <= rho
            lo = mid;
        else
            hi = mid;
        end
    end
    boundary(i,:) = (lo * d)';
end
results.area      = polyarea(boundary(:,1), boundary(:,2));
results.areaRatio = results.area / results.trueArea;

% 7. Assert containment -- the soundness check this script exists for.
%    Test membership, not distance: the distance from a certified point to the
%    nearest cycle point is positive whether the point is inside or outside.
inside = inpolygon(boundary(:,1), boundary(:,2), cycle(:,1), cycle(:,2));
if ~all(inside)
    bad = find(~inside, 1);
    error('lyapunovToolbox:falseCertificate', ...
          ['Certified boundary point (%.4f, %.4f) lies OUTSIDE the true ROA. ' ...
           'This is a false certificate: %d of %d boundary points escape.'], ...
          boundary(bad,1), boundary(bad,2), sum(~inside), n_rays);
end

% Report how much room is left. A certificate that only just fits is worth
% knowing about, and the number is meaningless if it is never printed.
minClearance = inf;
for i = 1:size(boundary, 1)
    d2 = (cycle(:,1) - boundary(i,1)).^2 + (cycle(:,2) - boundary(i,2)).^2;
    minClearance = min(minClearance, sqrt(min(d2)));
end
results.minClearance = minClearance;

% 8. Cross-validation fixture.
if config.WriteFixture
    write_roa_fixture('matlab_vdp', V, x, rho, info, config.Epsilon);
end

% 9. Plot.
if config.Plot
    figure;
    plot([cycle(:,1); cycle(1,1)], [cycle(:,2); cycle(1,2)], ...
         'k--', 'LineWidth', 2); hold on;
    fill(boundary(:,1), boundary(:,2), [0 0.4470 0.7410], ...
         'FaceAlpha', 0.25, 'EdgeColor', [0 0.4470 0.7410], 'LineWidth', 1.5);
    plot(0, 0, 'k.', 'MarkerSize', 15);
    xlabel('x_1'); ylabel('x_2');
    title(sprintf('Van der Pol ROA (\\rho = %.4f, %.1f%% of true area)', ...
          rho, results.areaRatio * 100));
    legend('True ROA (limit cycle, exact)', 'Certified ROA (SOS)', ...
           'Equilibrium', 'Location', 'best');
    axis equal; grid on;
end

fprintf(['Van der Pol (mu=%g): rho=%.6g, area=%.6g, trueArea=%.6g, ' ...
         'ratio=%.2f%%, clearance=%.4g, period=%.4f, solver=%s, time=%.2fs\n'], ...
        mu, rho, results.area, results.trueArea, results.areaRatio*100, ...
        results.minClearance, period, info.solverResolved, info.solverTime);
fprintf('  bounded via %s; gram residual %.2e; min eig Q %.2e\n', ...
        info.boundednessEvidence, info.gramResidual, info.minEigQ);

end % van_der_pol_roa


% =========================================================================
function [cycle, period] = limit_cycle_one_period(mu)
%LIMIT_CYCLE_ONE_PERIOD  Exactly one period of the forward VdP limit cycle.
%
%   Returns an (N x 2) open vertex list (first point not repeated at the end)
%   traversing the cycle once, and the period.
%
%   Method: integrate the forward oscillator from outside the cycle, record
%   crossings of the Poincare section {x2 = 0, x2 increasing} -- one crossing
%   per period -- and re-integrate densely between the last two. The section
%   direction matters: {x2 = 0} alone is crossed twice per period.

fwd = @(~, y) [y(2); -y(1) - mu*(y(1)^2 - 1)*y(2)];

odeOpts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12, ...
                 'Events', @section_x2_up);
[~, ~, te, ye] = ode45(fwd, [0 80], [2; 0], odeOpts);

if numel(te) < 4
    error('lyapunovToolbox:noLimitCycle', ...
          ['Found only %d Poincare crossings in 80 s; the forward oscillator ' ...
           'did not settle onto its limit cycle.'], numel(te));
end

% Use the last two crossings, so the transient is long past.
t0 = te(end-1);
t1 = te(end);
period = t1 - t0;

tGrid = linspace(t0, t1, 2001);
[~, orbit] = ode45(fwd, tGrid, ye(end-1,:)', ...
                   odeset('RelTol', 1e-10, 'AbsTol', 1e-12));
cycle = orbit(1:end-1, :);   % drop the closing point: the list is open

% Confirm the orbit has actually converged, by comparing this period's area
% against the previous one. An unsettled trace would give a spiral, and every
% number downstream would be wrong without any visible symptom.
tGridPrev = linspace(te(end-2), te(end-1), 2001);
[~, orbitPrev] = ode45(fwd, tGridPrev, ye(end-2,:)', ...
                       odeset('RelTol', 1e-10, 'AbsTol', 1e-12));
areaNow  = polyarea(cycle(:,1), cycle(:,2));
areaPrev = polyarea(orbitPrev(1:end-1,1), orbitPrev(1:end-1,2));
relDiff  = abs(areaNow - areaPrev) / max(1, areaPrev);
if relDiff > 1e-3
    error('lyapunovToolbox:cycleNotConverged', ...
          ['Consecutive periods enclose areas differing by %.2e (>1e-3); the ' ...
           'orbit has not settled onto the limit cycle. Integrate longer.'], relDiff);
end

% The shoelace area and INPOLYGON both assume a simple (non-self-intersecting)
% polygon. The VdP limit cycle encircles the origin, so a sufficient check is
% that the polar angle advances monotonically.
ang   = unwrap(atan2(cycle(:,2), cycle(:,1)));
dAng  = diff(ang);
if ~(all(dAng > 0) || all(dAng < 0))
    error('lyapunovToolbox:cycleNotSimple', ...
          ['The extracted cycle is not star-shaped about the origin, so its ' ...
           'shoelace area and inpolygon test are not trustworthy.']);
end

end


% =========================================================================
function [value, isterminal, direction] = section_x2_up(~, y)
%SECTION_X2_UP  Poincare section {x2 = 0} crossed with x2 increasing.
value      = y(2);
isterminal = 0;
direction  = 1;
end
