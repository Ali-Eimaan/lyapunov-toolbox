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

% --- Implementation ---
results = struct();

% 1. Parse options.
p = inputParser;
addParameter(p, 'M', 1.0, @(x) isscalar(x) && x > 0);
addParameter(p, 'L', 1.0, @(x) isscalar(x) && x > 0);
addParameter(p, 'B', 0.1, @(x) isscalar(x) && x > 0);
addParameter(p, 'G', 9.81, @(x) isscalar(x) && x > 0);
addParameter(p, 'TaylorOrder', 3, @(x) isscalar(x) && ismember(x, [3,5,7]));
addParameter(p, 'Epsilon', 1e-6, @(x) isscalar(x) && x > 0);
addParameter(p, 'LambdaDegree', 2, @(x) isscalar(x) && x >= 0 && mod(x,2)==0);
addParameter(p, 'Solver', 'mosek', @(x) ismember(x, {'mosek','sedumi','sdpt3'}));
addParameter(p, 'Plot', true, @islogical);
parse(p, varargin{:});
config = p.Results;
results.config = config;

m = config.M;  l = config.L;  b = config.B;  g = config.G;
tOrder = config.TaylorOrder;

% 2. Linearise about upright (unstable).
A = [0 1; g/l, -b/(m*l^2)];
B = [0; 1/(m*l^2)];
check_eig = eig(A);
if all(real(check_eig) < 0)
    error('lyapunovToolbox:wrongLinearisation', ...
          'A is Hurwitz before feedback; upright should be unstable.');
end

% 3. LQR design.
Q = diag([10 1]); R = 1;
K = lqr(A, B, Q, R);
Acl = A - B*K;
cl_eig = eig(Acl);
if any(real(cl_eig) >= 0)
    error('lyapunovToolbox:nonStabilisingGain', ...
          'eig(A-B*K) has non-negative real part: LQR did not stabilise.');
end
results.K = K;

% 4. Build closed-loop polynomial field with Taylor recast.
sdpvar th thd;
x = [th; thd];
sinApprox = taylor_sin(th, tOrder);
u = -K*x;
f = [thd; (g/l)*sinApprox - (b/(m*l^2))*thd + u/(m*l^2)];

% 5. V0 from Lyapunov equation, then compute ROA.
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
results.rho = rho;
results.V = V;
results.solverTime = info.solverTime;

% 6. True ROA by simulation of exact (sin) dynamics.
exact_ode = @(t, y) [y(2);
    (g/l)*sin(y(1)) - (b/(m*l^2))*y(2) + (-K*y)/(m*l^2)];
theta_grid = linspace(-pi, pi, 40);
thd_grid = linspace(-10, 10, 40);
[TT, DD] = meshgrid(theta_grid, thd_grid);
converged = false(size(TT));
for i = 1:numel(TT)
    [~, yy] = ode45(exact_ode, [0 20], [TT(i); DD(i)]);
    converged(i) = norm(yy(end,:)) < 1e-3;
end

% 7. Assert certified set excludes theta = ±pi.
Vpi_val = double(replace(V, x, [pi; 0]));
Vmipi_val = double(replace(V, x, [-pi; 0]));
if Vpi_val <= rho || Vmipi_val <= rho
    error('lyapunovToolbox:falseCertificate', ...
          'Certified set includes theta = ±pi. Certificate is FALSE.');
end

% 8. Simulate true dynamics from points inside certified set.
n_test = 200;
rng(1);  % deterministic
test_pts = zeros(n_test, 2);
count = 0;
while count < n_test
    cand = [pi*(2*rand-1); 10*(2*rand-1)];
    if double(replace(V, x, cand)) <= rho
        count = count + 1;
        test_pts(count,:) = cand';
    end
end
test_converged = false(n_test, 1);
for i = 1:n_test
    [~, yy] = ode45(exact_ode, [0 20], test_pts(i,:)');
    test_converged(i) = norm(yy(end,:)) < 1e-3;
end
results.trueFraction = mean(test_converged);
failed = test_pts(~test_converged, :);
if isempty(failed)
    results.breakdownTheta = Inf;
else
    results.breakdownTheta = min(abs(failed(:,1)));
end

% 9. Plot.
if config.Plot
    figure;
    % True ROA grid.
    [TT, DD] = deal(TT, DD);
    theta_conv = TT(converged(:));
    thd_conv = DD(converged(:));
    plot(theta_conv, thd_conv, 'g.', 'MarkerSize', 2); hold on;

    % Certified boundary.
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
    fill(boundary(:,1), boundary(:,2), [0 0.4470 0.7410], ...
         'FaceAlpha', 0.25, 'EdgeColor', [0 0.4470 0.7410], 'LineWidth', 1.5);

    % Separatrices through (±pi, 0).
    sep_fwd = []; sep_bwd = [];
    opts_ode = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [~, y_sep] = ode45(exact_ode, [0 5], [pi-1e-6; 0], opts_ode);
    sep_fwd = y_sep;
    [~, y_sep2] = ode45(exact_ode, [0 5], [pi+1e-6; 0], opts_ode);
    sep_bwd = y_sep2;
    plot(sep_fwd(:,1), sep_fwd(:,2), 'r-', 'LineWidth', 1.5);
    plot(sep_bwd(:,1), sep_bwd(:,2), 'r-', 'LineWidth', 1.5);

    % Failing test points.
    if ~isempty(failed)
        plot(failed(:,1), failed(:,2), 'rx', 'MarkerSize', 8, 'LineWidth', 1.5);
    end

    plot(0, 0, 'k.', 'MarkerSize', 15);
    xlabel('theta [rad]'); ylabel('thetadot [rad/s]');
    title(sprintf('Pendulum ROA (order=%d, \\rho=%.3f, frac=%.2f)', ...
          tOrder, rho, results.trueFraction));
    axis equal; grid on;
    legend('True ROA (grid)', 'Certified ROA', 'Separatrix', ...
           'Location', 'best');
end

% 10. Sweep Taylor orders.
if nargout == 0
    fprintf('Pendulum: order=%d, rho=%.4f, trueFraction=%.2f%%, breakdownTheta=%.3f, time=%.2fs\n', ...
            tOrder, rho, results.trueFraction*100, results.breakdownTheta, results.solverTime);
    orders = [3 5 7];
    sweep = struct('order', {}, 'rho', {}, 'trueFraction', {}, 'breakdownTheta', {}, 'time', {});
    for oi = 1:numel(orders)
        if orders(oi) == tOrder, continue; end
        sw = pendulum_roa('M',m,'L',l,'B',b,'G',g, ...
                          'TaylorOrder', orders(oi), ...
                          'Epsilon', config.Epsilon, ...
                          'LambdaDegree', config.LambdaDegree, ...
                          'Solver', config.Solver, ...
                          'Plot', false);
        sweep(end+1) = struct('order', orders(oi), 'rho', sw.rho, ...
                              'trueFraction', sw.trueFraction, ...
                              'breakdownTheta', sw.breakdownTheta, ...
                              'time', sw.solverTime); %#ok<AGROW>
    end
    tbl = struct2table(sweep);
    disp(tbl);
end
end

% -------------------------------------------------------------------------
function s = taylor_sin(th, order)
    assert(mod(order, 2) == 1, ...
           'TaylorOrder must be odd, got %d (even order wastes a degree).', order);
    s = 0;
    sign = 1;
    for k = 1:2:order
        s = s + sign * th^k / factorial(k);
        sign = -sign;
    end
end
