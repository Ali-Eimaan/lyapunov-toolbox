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
%   When nothing fails, RESULTS.MAXTHETATESTED is the number that matters:
%   "no counterexample below |theta| = 0.9 rad" is a usable claim, whereas
%   "the Taylor certificate is valid" is not one this script can support.
%
%   RESULTS = PENDULUM_ROA('Name', Value, ...) accepts:
%       'M','L','B','G'    physical parameters      (1.0, 1.0, 0.1, 9.81)
%       'TaylorOrder'      3, 5, or 7               (3)
%       'Epsilon'          margin                   (1e-6)
%       'LambdaDegree'     multiplier degree        (2)
%       'Solver'           SDP solver               ('mosek')
%       'Plot'             draw the figure          (true)
%       'WriteFixture'     write the cross-validation CSV   (false)
%       'Seed'             RNG seed for the sampling study   (1)
%
%   RESULTS is a struct with RHO, V, P, K, TRUEFRACTION, BREAKDOWNTHETA,
%   MAXTHETATESTED, MINCLEARANCE, SOLVERTIME, INFO, CONFIG.
%
%   See also: COMPUTE_ROA_YALMIP, VAN_DER_POL_ROA, WRITE_ROA_FIXTURE.

% --- Implementation ---
results = struct();

% 1. Parse options.
p = inputParser;
addParameter(p, 'M', 1.0, @(v) isscalar(v) && v > 0);
addParameter(p, 'L', 1.0, @(v) isscalar(v) && v > 0);
addParameter(p, 'B', 0.1, @(v) isscalar(v) && v > 0);
addParameter(p, 'G', 9.81, @(v) isscalar(v) && v > 0);
addParameter(p, 'TaylorOrder', 3, @(v) isscalar(v) && ismember(v, [3,5,7]));
addParameter(p, 'Epsilon', 1e-6, @(v) isscalar(v) && v > 0);
addParameter(p, 'LambdaDegree', 2, @(v) isscalar(v) && v >= 0 && mod(v,2)==0);
addParameter(p, 'Solver', 'mosek', @(v) ismember(v, {'mosek','sedumi','sdpt3'}));
addParameter(p, 'Plot', true, @islogical);
addParameter(p, 'WriteFixture', false, @islogical);
addParameter(p, 'Seed', 1, @(v) isscalar(v) && v >= 0);
parse(p, varargin{:});
config = p.Results;

% 1b. Check the environment BEFORE the first SDPVAR. Without this the script
%     dies on "Unrecognized function or variable 'sdpvar'", which tells a
%     first-time user nothing about what to install.
config.Solver = assert_yalmip_env(config.Solver);
results.config = config;

m = config.M;  l = config.L;  b = config.B;  g = config.G;
tOrder = config.TaylorOrder;

% 2. Linearise about upright. The (2,1) entry is +g/l: upright is UNSTABLE
%    open-loop, which is the whole reason a controller is present. A Hurwitz A
%    here means the sign convention is wrong.
A = [0 1; g/l, -b/(m*l^2)];
B = [0; 1/(m*l^2)];
if all(real(eig(A)) < 0)
    error('lyapunovToolbox:wrongLinearisation', ...
          'A is Hurwitz before feedback; upright should be unstable.');
end

% 3. LQR design. Verify the closed loop is Hurwitz before any SOS program: an
%    infeasibility caused by a non-stabilising gain is indistinguishable from a
%    real one, and it costs an afternoon.
Q = diag([10 1]); R = 1;
K = lqr(A, B, Q, R);
Acl = A - B*K;
if any(real(eig(Acl)) >= 0)
    error('lyapunovToolbox:nonStabilisingGain', ...
          'eig(A-B*K) has non-negative real part: LQR did not stabilise.');
end
results.K = K;

% 4. Closed-loop polynomial field with the Taylor recast.
sdpvar th thd;
x = [th; thd];
sinApprox = taylor_sin(th, tOrder);
u = -K*x;
f = [thd; (g/l)*sinApprox - (b/(m*l^2))*thd + u/(m*l^2)];

% 5. V0 from the Lyapunov equation, then the certified level.
P = lyap(Acl', eye(2));
residual = Acl'*P + P*Acl + eye(2);
assert(norm(residual) < 1e-9 * max(1, norm(P)), ...
       'Lyapunov equation residual too large: %e', norm(residual));
V0 = x'*P*x;
results.P = P;

opts = struct('epsilon', config.Epsilon, ...
              'lambdaDegree', config.LambdaDegree, ...
              'solver', config.Solver);
[rho, V, info] = compute_roa_yalmip(f, x, V0, opts);
results.rho        = rho;
results.V          = V;
results.info       = info;
results.solverTime = info.solverTime;

% 6. Assert the certified set excludes theta = +-pi.
%    The downward equilibrium is a saddle of the closed loop and cannot lie in
%    the ROA of the upright one. Cheapest end-to-end soundness check available,
%    and it must ERROR rather than warn.
Vpi   = double(replace(V, x, [ pi; 0]));
Vmipi = double(replace(V, x, [-pi; 0]));
if Vpi <= rho || Vmipi <= rho
    error('lyapunovToolbox:falseCertificate', ...
          ['Certified set includes theta = +-pi (V(pi,0)=%.4g, V(-pi,0)=%.4g, ' ...
           'rho=%.4g). The certificate is FALSE.'], Vpi, Vmipi, rho);
end
results.minClearance = min(Vpi, Vmipi) - rho;

% 7. The honest cross-check: sample INSIDE the certified set, then integrate the
%    TRUE trigonometric dynamics. The certificate was proved for the truncated
%    field; this measures how far that transfers.
exact_ode = @(~, y) [y(2);
    (g/l)*sin(y(1)) - (b/(m*l^2))*y(2) + (-K*y)/(m*l^2)];
odeTol = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

n_test  = 200;
rng(config.Seed, 'twister');
test_pts   = zeros(n_test, 2);
count      = 0;
attempts   = 0;
maxAttempts = 200 * n_test;
while count < n_test && attempts < maxAttempts
    attempts = attempts + 1;
    cand = [pi*(2*rand-1); 10*(2*rand-1)];
    if double(replace(V, x, cand)) <= rho
        count = count + 1;
        test_pts(count,:) = cand';
    end
end
if count < n_test
    % Do not silently report a fraction computed from three points.
    error('lyapunovToolbox:samplingFailed', ...
          ['Only %d of %d interior samples found in %d attempts; the certified ' ...
           'set is a vanishing fraction of the sampling box (rho=%.4g).'], ...
          count, n_test, attempts, rho);
end

test_converged = false(n_test, 1);
for i = 1:n_test
    [~, yy] = ode45(exact_ode, [0 20], test_pts(i,:)', odeTol);
    test_converged(i) = norm(yy(end,:)) < 1e-3;
end
results.trueFraction   = mean(test_converged);
results.maxThetaTested = max(abs(test_pts(:,1)));

failed = test_pts(~test_converged, :);
if isempty(failed)
    results.breakdownTheta = Inf;   % read together with maxThetaTested
else
    results.breakdownTheta = min(abs(failed(:,1)));
end

% 8. Certified boundary by ray-marching from the origin (the equilibrium).
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
results.area = polyarea(boundary(:,1), boundary(:,2));

% 9. Cross-validation fixture.
if config.WriteFixture
    write_roa_fixture('matlab_pendulum', V, x, rho, info, config.Epsilon);
end

% 10. Plot. The empirical ROA grid costs 1600 ODE solves, so it is computed
%     only when it is going to be drawn -- the sweep in step 11 calls this
%     function with 'Plot', false and would otherwise pay for it three times.
if config.Plot
    theta_grid = linspace(-pi, pi, 40);
    thd_grid   = linspace(-10, 10, 40);
    [TT, DD]   = meshgrid(theta_grid, thd_grid);
    converged  = false(size(TT));
    for i = 1:numel(TT)
        [~, yy] = ode45(exact_ode, [0 20], [TT(i); DD(i)], odeTol);
        converged(i) = norm(yy(end,:)) < 1e-3;
    end

    figure;
    plot(TT(converged), DD(converged), 'g.', 'MarkerSize', 2); hold on;
    fill(boundary(:,1), boundary(:,2), [0 0.4470 0.7410], ...
         'FaceAlpha', 0.25, 'EdgeColor', [0 0.4470 0.7410], 'LineWidth', 1.5);

    % Separatrix branches through (+-pi, 0), integrated backwards from just
    % either side of the saddle.
    [~, sep_a] = ode45(exact_ode, [0 5], [ pi-1e-6; 0], odeTol);
    [~, sep_b] = ode45(exact_ode, [0 5], [-pi+1e-6; 0], odeTol);
    plot(sep_a(:,1), sep_a(:,2), 'r-', 'LineWidth', 1.5);
    plot(sep_b(:,1), sep_b(:,2), 'r-', 'LineWidth', 1.5);

    if ~isempty(failed)
        plot(failed(:,1), failed(:,2), 'rx', 'MarkerSize', 8, 'LineWidth', 1.5);
    end

    plot(0, 0, 'k.', 'MarkerSize', 15);
    xlabel('\theta [rad]'); ylabel('$\dot\theta$ [rad/s]', 'Interpreter', 'latex');
    title(sprintf('Pendulum ROA (Taylor order %d, \\rho = %.4f, converged %.1f%%)', ...
          tOrder, rho, results.trueFraction*100));
    legend('True ROA (grid, empirical)', 'Certified ROA (SOS, truncated model)', ...
           'Separatrix', 'Location', 'best');
    axis equal; grid on;
end

fprintf(['Pendulum (order %d): rho=%.6g, area=%.6g, converged=%.1f%% of %d, ' ...
         'maxThetaTested=%.4f, breakdownTheta=%.4f, solver=%s, time=%.2fs\n'], ...
        tOrder, rho, results.area, results.trueFraction*100, n_test, ...
        results.maxThetaTested, results.breakdownTheta, ...
        info.solverResolved, info.solverTime);
fprintf('  bounded via %s; gram residual %.2e; min eig Q %.2e\n', ...
        info.boundednessEvidence, info.gramResidual, info.minEigQ);

% 11. Taylor-order sweep, when called as a script (no output requested).
%     A higher-order truncation is a better model, so the converged fraction
%     should improve. Whether rho does is not obvious -- report it either way.
if nargout == 0
    orders = [3 5 7];
    sweep  = struct('order', {}, 'rho', {}, 'trueFraction', {}, ...
                    'breakdownTheta', {}, 'maxThetaTested', {}, 'time', {});
    for oi = 1:numel(orders)
        if orders(oi) == tOrder
            sw = results;
        else
            sw = pendulum_roa('M', m, 'L', l, 'B', b, 'G', g, ...
                              'TaylorOrder', orders(oi), ...
                              'Epsilon', config.Epsilon, ...
                              'LambdaDegree', config.LambdaDegree, ...
                              'Solver', config.Solver, ...
                              'Seed', config.Seed, ...
                              'Plot', false);
        end
        sweep(end+1) = struct('order', orders(oi), 'rho', sw.rho, ...
                              'trueFraction', sw.trueFraction, ...
                              'breakdownTheta', sw.breakdownTheta, ...
                              'maxThetaTested', sw.maxThetaTested, ...
                              'time', sw.solverTime); %#ok<AGROW>
    end
    disp(struct2table(sweep));
end

end % pendulum_roa


% =========================================================================
function s = taylor_sin(th, order)
%TAYLOR_SIN  Odd Taylor truncation of sin, to the given odd order.
%
%   An even-order truncation of an odd function adds nothing and wastes a
%   degree, so even orders are rejected rather than silently accepted.

assert(mod(order, 2) == 1, ...
       'TaylorOrder must be odd, got %d (an even order wastes a degree).', order);
s   = 0;
sgn = 1;                      % not `sign`: that shadows the builtin
for k = 1:2:order
    s   = s + sgn * th^k / factorial(k);
    sgn = -sgn;
end

end
