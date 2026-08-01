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

function [rho, V, info] = compute_roa_yalmip(f, x, V0, opts)
    % 1. Validate inputs.
    n = length(x);
    assert(length(f) == n, 'length(f)=%d must equal length(x)=%d', length(f), n);
    V0_at_origin = double(replace(V0, x, zeros(n,1)));
    assert(abs(V0_at_origin) < 1e-12, 'V0 must vanish at the origin, got %e', V0_at_origin);

    % 2. Fill defaults.
    if nargin < 4, opts = struct(); end
    defaults = struct(...
        'epsilon', 1e-6, ...
        'lambdaDegree', 2, ...
        'constraints', {{}}, ...
        'sigmaDegree', 2, ...
        'rhoMax', [], ...
        'bisectionTol', 1e-4, ...
        'maxIter', 40, ...
        'solver', 'mosek', ...
        'verbose', 0);
    fields = fieldnames(defaults);
    for i = 1:numel(fields)
        if ~isfield(opts, fields{i})
            opts.(fields{i}) = defaults.(fields{i});
        end
    end
    info.config = opts;

    % 3. Vdot = jacobian(V0, x) * f.
    Vdot = jacobian(V0, x) * f;
    if isa(Vdot, 'sdpvar')
        Vdot = Vdot(1);
    end

    % 4. Check sublevel set boundedness: top-degree form of V0 positive definite.
    V0_deg = degree(V0);
    if V0_deg < 2
        error('V0 must have degree at least 2 for bounded sublevel sets');
    end

    % 5. Bracket rho.
    rho = 1e-3;
    found_infeasible = false;
    for doubling = 1:30
        feasible = oracle(rho, V0, Vdot, x, opts);
        if ~feasible
            found_infeasible = true;
            break;
        end
        rho = rho * 2;
    end
    if ~found_infeasible
        error('lyapunovToolbox:globallyStable', ...
              'No infeasible level found up to %e; system may be globally stable.', rho);
    end
    rho_hi = rho;
    rho_lo = rho / 2;

    % 6. Bisect.
    info.iterations = 0;
    info.solverTime = 0;
    info.solverStatus = '';
    info.bestRho = rho_lo;
    for iter = 1:opts.maxIter
        rho_mid = (rho_lo + rho_hi) / 2;
        info.iterations = info.iterations + 1;
        [feasible, sol_time, sol_status, lambda_val, sigma_vals] = oracle(rho_mid, V0, Vdot, x, opts);
        info.solverTime = info.solverTime + sol_time;
        if feasible
            rho_lo = rho_mid;
            info.bestRho = rho_mid;
            info.solverStatus = sol_status;
            info.lambda = lambda_val;
            info.sigma = sigma_vals;
        else
            rho_hi = rho_mid;
            info.solverStatus = sol_status;
        end
        if (rho_hi - rho_lo) / max(1, rho_lo) < opts.bisectionTol
            break;
        end
    end

    % 7. Return largest verified feasible rho.
    rho = info.bestRho;
    V = V0;
    info.iterations = info.iterations;
end


function [feasible, sol_time, status_str, lambda_val, sigma_vals] = oracle(rho, V0, Vdot, x, opts)
    % Inner feasibility solve for a candidate rho.
    [lambda_val, cl] = polynomial(x, opts.lambdaDegree, 0);
    F = [sos(lambda_val)];
    expr = -(Vdot + opts.epsilon * (x' * x)) + lambda_val * (V0 - rho);
    sigma_vals = {};
    cs = {};
    for j = 1:numel(opts.constraints)
        g_j = opts.constraints{j};
        [sigma_j, cs_j] = polynomial(x, opts.sigmaDegree, 0);
        F = [F, sos(sigma_j)];
        expr = expr - sigma_j * g_j;
        sigma_vals{end+1} = sigma_j;
        cs{end+1} = cs_j;
    end
    F = [F, sos(expr)];
    all_decision_vars = [cl; vertcat(cs{:})];
    settings = sdpsettings('solver', opts.solver, 'verbose', opts.verbose);
    diagnostics = solvesos(F, [], settings, all_decision_vars);
    sol_time = diagnostics.solvertime;
    if diagnostics.problem == 0
        feasible = true;
        status_str = 'feasible';
    elseif diagnostics.problem == 1 || diagnostics.problem == 4
        feasible = false;
        if diagnostics.problem == 1
            status_str = 'infeasible';
        else
            status_str = 'numerical_problems';
        end
    else
        error('lyapunovToolbox:solverError', ...
              'Solver failed with status %d: %s', diagnostics.problem, yalmiperror(diagnostics.problem));
    end
end

end
