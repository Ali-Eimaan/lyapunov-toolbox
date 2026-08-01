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

% --- Implementation ---
results = struct();

% 1. Parse options.
p = inputParser;
addParameter(p, 'Mu', 1.0, @(x) isscalar(x) && x > 0);
addParameter(p, 'Epsilon', 1e-6, @(x) isscalar(x) && x > 0);
addParameter(p, 'LambdaDegree', 2, @(x) isscalar(x) && x >= 0 && mod(x,2)==0);
addParameter(p, 'Solver', 'mosek', @(x) ismember(x, {'mosek','sedumi','sdpt3'}));
addParameter(p, 'Plot', true, @islogical);
addParameter(p, 'WriteFixture', false, @islogical);
parse(p, varargin{:});
config = p.Results;
results.config = config;
mu = config.Mu;

% 2. Declare dynamics (reversed Van der Pol).
sdpvar x1 x2;
x = [x1; x2];
f = [-x2; x1 + mu*(x1^2 - 1)*x2];

% 3. Seed from linearisation.
A = [0 -1; 1 -mu];
P = lyap(A', eye(2));
residual = A'*P + P*A + eye(2);
assert(norm(residual) < 1e-9 * max(1, norm(P)), ...
       'Lyapunov equation residual too large: %e', norm(residual));
V0 = x'*P*x;

% 4. Compute ROA.
opts = struct('epsilon', config.Epsilon, ...
              'lambdaDegree', config.LambdaDegree, ...
              'solver', config.Solver);
[rho, V, info] = compute_roa_yalmip(f, x, V0, opts);
results.rho = rho;
results.V = V;
results.P = P;
results.solverTime = info.solverTime;

% 5. Ground truth: forward system's limit cycle.
fwd_ode = @(t, y) [y(2); -y(1) - mu*(y(1)^2 - 1)*y(2)];
[~, y] = ode45(fwd_ode, [0 40], [2; 0]);
half = round(length(y)/2);
cycle = y(half:end, :);
results.trueArea = polyarea(cycle(:,1), cycle(:,2));

% 6. Certified area (ray-march from origin).
n_rays = 400;
angles = linspace(0, 2*pi, n_rays);
boundary = zeros(n_rays, 2);
for i = 1:n_rays
    d = [cos(angles(i)); sin(angles(i))];
    lo = 0; hi = 100;
    for k = 1:60
        mid = (lo + hi) / 2;
        pt = mid * d;
        if double(replace(V, x, pt)) <= rho
            lo = mid;
        else
            hi = mid;
        end
    end
    boundary(i,:) = (lo * d)';
end
results.area = polyarea(boundary(:,1), boundary(:,2));
results.areaRatio = results.area / results.trueArea;

% 7. Assert containment.
inside = inpolygon(boundary(:,1), boundary(:,2), cycle(:,1), cycle(:,2));
if ~all(inside)
    bad_idx = find(~inside, 1);
    error('lyapunovToolbox:falseCertificate', ...
          'Certified boundary point (%.3f, %.3f) lies outside the true ROA.', ...
          boundary(bad_idx,1), boundary(bad_idx,2));
end

% 8. Write cross-validation fixture.
if config.WriteFixture
    fixture_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'tests', 'fixtures');
    if ~exist(fixture_dir, 'dir'), mkdir(fixture_dir); end
    fixture_path = fullfile(fixture_dir, 'matlab_vdp.csv');
    fid = fopen(fixture_path, 'w');
    fprintf(fid, '%e\n', rho);
    fclose(fid);
    fprintf('Fixture written to %s\n', fixture_path);
end

% 9. Plot.
if config.Plot
    figure;
    plot(cycle(:,1), cycle(:,2), 'k--', 'LineWidth', 2); hold on;
    fill(boundary(:,1), boundary(:,2), [0 0.4470 0.7410], ...
         'FaceAlpha', 0.25, 'EdgeColor', [0 0.4470 0.7410], 'LineWidth', 1.5);
    plot(0, 0, 'k.', 'MarkerSize', 15);
    xlabel('x_1'); ylabel('x_2');
    title(sprintf('Van der Pol ROA (\\rho = %.3f, ratio = %.1f%%)', ...
          rho, results.areaRatio * 100));
    legend('True ROA (limit cycle)', 'Certified ROA', 'Equilibrium', 'Location', 'best');
    axis equal; grid on;
end

fprintf('Van der Pol: rho=%.4f, area=%.4f, trueArea=%.4f, ratio=%.2f%%, time=%.2fs\n', ...
        rho, results.area, results.trueArea, results.areaRatio*100, results.solverTime);
end
