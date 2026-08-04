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
%               .boundedRadius  sphere radius for the boundedness
%                               sampling fallback                 (50)
%               .boundedSamples number of sphere samples          (5000)
%
%   OUTPUTS
%     RHO   - certified level (double). Always a VERIFIED feasible value: it is
%             re-solved once at the end of the bisection, and the routine errors
%             rather than returning a level no solver ever accepted.
%     V     - the Lyapunov function (unchanged from V0 in this routine; the
%             V-s iteration is Julia-side only, deliberately - see README).
%     INFO  - struct with:
%               .iterations           bisection iterations run
%               .solverTime           summed solver time (s)
%               .solverStatus         status at the returned RHO
%               .lambda               multiplier at the returned RHO
%               .sigma                cell array of constraint multipliers
%               .gramV                Gram matrix of the V0 positivity
%                                     certificate (P1)
%               .gramVdot             Gram matrix of the derivative
%                                     constraint (P2/P3) at the returned RHO
%               .gramLambda           Gram matrix of the SOS multiplier
%               .gramResidual         worst SOS residual reported by SOLVESOS
%               .minEigQ              smallest eigenvalue over all Gram
%                                     matrices at the returned RHO
%               .bounded              logical
%               .boundednessEvidence  'top_form_pd' | 'sampling_only'
%               .solverResolved       solver actually used
%               .config               every OPTS field actually used, so a
%                                     result is reproducible from INFO alone
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
%   BOUNDEDNESS. The conclusion "Omega_rho is inside the ROA" needs Omega_rho to
%   be BOUNDED; without it, Vdot < 0 on the sublevel set implies neither
%   invariance nor convergence. This routine checks it and reports which kind of
%   evidence it obtained, mirroring is_sublevel_set_bounded on the Julia side:
%     'top_form_pd'   - conclusive (exact for a quadratic V0);
%     'sampling_only' - a necessary-condition screen, NOT a proof.
%   A witness of unboundedness is an error, not a warning.
%
%   REQUIREMENTS
%     YALMIP (tested with R20230622) and one SDP solver. Mosek needs a licence;
%     SeDuMi is the free fallback and is what the SOLVER field should be set to
%     on a machine without one. The routine ERRORS with an actionable message
%     when the requested solver is not on the path -- it does not silently
%     substitute another one, because two benchmark rows produced by different
%     solvers without saying so are worse than one missing row.
%
%   See also: VAN_DER_POL_ROA, PENDULUM_ROA, WRITE_ROA_FIXTURE.

% ---------------------------------------------------------------- validation
n = length(x);
assert(length(f) == n, 'length(f)=%d must equal length(x)=%d', length(f), n);

V0_at_origin = double(replace(V0, x, zeros(n,1)));
assert(abs(V0_at_origin) < 1e-12, ...
       'V0 must vanish at the origin, got %e', V0_at_origin);

f_at_origin = double(replace(f, x, zeros(n,1)));
assert(max(abs(f_at_origin)) < 1e-12, ...
       'f must vanish at the origin, got ||f(0)||_inf = %e', max(abs(f_at_origin)));

% ------------------------------------------------------------------ defaults
if nargin < 4, opts = struct(); end
defaults = struct( ...
    'epsilon',        1e-6, ...
    'lambdaDegree',   2, ...
    'constraints',    {{}}, ...
    'sigmaDegree',    2, ...
    'rhoMax',         [], ...
    'bisectionTol',   1e-4, ...
    'maxIter',        40, ...
    'solver',         'mosek', ...
    'verbose',        0, ...
    'boundedRadius',  50, ...
    'boundedSamples', 5000);
fields = fieldnames(defaults);
for i = 1:numel(fields)
    if ~isfield(opts, fields{i})
        opts.(fields{i}) = defaults.(fields{i});
    end
end

% Resolve the environment BEFORE any solve, so a missing YALMIP or licence is
% reported once, by name, rather than as a confusing failure inside the
% bisection. Entry-point scripts call this too, before their first SDPVAR --
% it is idempotent and cheap.
opts.solver = assert_yalmip_env(opts.solver);
info = struct();
info.config = opts;
info.solverResolved = opts.solver;
settings = sdpsettings('solver', opts.solver, 'verbose', opts.verbose);

% -------------------------------------------------------- Vdot and degrees
% YALMIP's jacobian returns a ROW vector, so this is row-times-column. A
% transpose slip here yields a scalar of the wrong sign and a certificate that
% looks entirely plausible.
Vdot = jacobian(V0, x) * f;
if numel(Vdot) ~= 1
    error('lyapunovToolbox:badVdot', ...
          'Vdot should be scalar, got %d entries; check the jacobian orientation.', ...
          numel(Vdot));
end

degV = degree(V0);
if mod(degV, 2) ~= 0 || degV < 2
    error('lyapunovToolbox:badCandidateDegree', ...
          'V0 must have even degree >= 2, got %d.', degV);
end

% ------------------------------------------- (P1) positive definiteness of V0
% V0 - epsilon*(x'x)^(degV/2) must be SOS. This is the hypothesis every later
% step assumes and it was previously never checked.
posExpr = V0 - opts.epsilon * (x' * x)^(degV/2);
[posDiag, ~, posQ, posRes] = solvesos(sos(posExpr), [], settings, []);
if posDiag.problem ~= 0
    error('lyapunovToolbox:candidateNotPositive', ...
          ['V0 - %g*(x''x)^%d is not SOS (solver status %d: %s). V0 is not a ' ...
           'positive definite candidate at this margin; reduce epsilon or fix ' ...
           'the Lyapunov solve.'], ...
          opts.epsilon, degV/2, posDiag.problem, yalmiperror(posDiag.problem));
end
info.gramV = posQ{1};
info.solverTime = posDiag.solvertime;

% --------------------------------------------------------------- boundedness
[isBounded, evidence] = check_boundedness(V0, x, degV, opts);
info.bounded = isBounded;
info.boundednessEvidence = evidence;
if ~isBounded
    error('lyapunovToolbox:unboundedSublevelSet', ...
          ['The sublevel sets of V0 are not bounded (evidence: %s). Vdot < 0 on ' ...
           'an unbounded set certifies nothing -- boundedness is a hypothesis of ' ...
           'the ROA theorem, not a technicality.'], evidence);
end

% ------------------------------------------------------------------ bracket
% Feasibility is MONOTONE in rho: a certificate for {V <= rho} is a certificate
% for every smaller sublevel set, because that set is a subset. That is the
% entire justification for bisecting a problem that is bilinear in (lambda, rho).
% The status at each probe is not recorded: the routine re-solves at the level
% it finally returns (see below), so info.solverStatus always describes THAT
% solve rather than whichever bisection step happened to run last.
bestRho       = [];       % assigned ONLY at a level a solver accepted
rhoLo         = 0;
info.iterations = 0;

if isempty(opts.rhoMax)
    rhoProbe       = 1e-3;
    foundInfeasible = false;
    for doubling = 1:30
        [ok, tsol] = oracle(rhoProbe, V0, Vdot, x, opts, settings);
        info.solverTime = info.solverTime + tsol;
        if ok
            bestRho    = rhoProbe;
            rhoLo      = rhoProbe;
            rhoProbe   = rhoProbe * 2;
        else
            foundInfeasible = true;
            rhoHi = rhoProbe;
            break;
        end
    end
    if ~foundInfeasible
        error('lyapunovToolbox:globallyStable', ...
              ['No infeasible level found up to rho = %e after 30 doublings. ' ...
               'For a globally stable system there is no finite boundary; pass ' ...
               'opts.rhoMax explicitly if you want a capped estimate. Returning ' ...
               'the cap as a certified level would be a silent over-estimate.'], ...
              rhoProbe);
    end
else
    rhoHi = opts.rhoMax;
end

% ------------------------------------------------------------------ bisect
for iter = 1:opts.maxIter
    if (rhoHi - rhoLo) <= opts.bisectionTol * max(1, rhoLo)
        break;
    end
    rhoMid = (rhoLo + rhoHi) / 2;
    info.iterations = info.iterations + 1;
    [ok, tsol] = oracle(rhoMid, V0, Vdot, x, opts, settings);
    info.solverTime = info.solverTime + tsol;
    if ok
        bestRho    = rhoMid;
        rhoLo      = rhoMid;
    else
        rhoHi = rhoMid;
    end
end

if isempty(bestRho)
    error('lyapunovToolbox:noFeasibleLevel', ...
          ['No feasible rho found in (0, %e]. V0 is probably not a Lyapunov ' ...
           'function for this field near the origin -- check the sign of Vdot ' ...
           'and that the linearisation is Hurwitz.'], rhoHi);
end

% -------------------------------- final re-solve at the returned level
% Re-solve at bestRho so the returned multipliers and Gram matrices belong to
% the level actually reported. Without this the witness would be whichever
% bisection step happened to run last, which is not necessarily bestRho.
[ok, tsol, status, lambdaVal, sigmaVals, gramQ, gramRes] = ...
    oracle(bestRho, V0, Vdot, x, opts, settings);
info.solverTime = info.solverTime + tsol;
if ~ok
    error('lyapunovToolbox:reverificationFailed', ...
          ['rho = %e was feasible during bisection but failed re-verification ' ...
           '(status: %s). This is a numerical-reproducibility failure, not a ' ...
           'tighter answer; do not use the level.'], bestRho, status);
end

rho = bestRho;
V   = V0;

info.solverStatus = status;
info.lambda       = lambdaVal;
info.sigma        = sigmaVals;
info.gramLambda   = gramQ{1};
info.gramVdot     = gramQ{end};
info.gramResidual = max(abs([gramRes(:); posRes(:)]));
minEig = min(eig((info.gramV + info.gramV')/2));
for kQ = 1:numel(gramQ)
    Qk = (gramQ{kQ} + gramQ{kQ}')/2;
    minEig = min(minEig, min(eig(Qk)));
end
info.minEigQ = minEig;

if opts.verbose
    fprintf(['compute_roa_yalmip: rho=%.6g, %d bisection iters, solver=%s, ' ...
             'time=%.2fs, gram residual=%.2e, min eig Q=%.2e, bounded via %s\n'], ...
            rho, info.iterations, info.solverResolved, info.solverTime, ...
            info.gramResidual, info.minEigQ, info.boundednessEvidence);
end

end % compute_roa_yalmip


% =========================================================================
function [feasible, solTime, statusStr, lambdaVal, sigmaVals, gramQ, gramRes] = ...
    oracle(rho, V0, Vdot, x, opts, settings)
%ORACLE  Inner feasibility solve for one candidate rho.

[lambdaVal, cl] = polynomial(x, opts.lambdaDegree, 0);
F    = sos(lambdaVal);
expr = -(Vdot + opts.epsilon * (x' * x)) + lambdaVal * (V0 - rho);

nCon       = numel(opts.constraints);
sigmaVals  = cell(1, nCon);
cs         = cell(1, nCon);
for j = 1:nCon
    [sigma_j, cs_j] = polynomial(x, opts.sigmaDegree, 0);
    F              = [F, sos(sigma_j)];      %#ok<AGROW>
    expr           = expr - sigma_j * opts.constraints{j};
    sigmaVals{j}   = sigma_j;
    cs{j}          = cs_j;
end
F = [F, sos(expr)];

decisionVars = [cl; vertcat(cs{:})];
[diagnostics, ~, gramQ, gramRes] = solvesos(F, [], settings, decisionVars);
solTime = diagnostics.solvertime;

switch diagnostics.problem
    case 0
        feasible  = true;
        statusStr = 'feasible';
    case 1
        feasible  = false;
        statusStr = 'infeasible';
    case 4
        % "Numerical problems". Treated as INFEASIBLE deliberately: accepting a
        % marginal solve inflates rho, and an over-estimate of a certified ROA
        % is a false certificate. This is the one direction of error this
        % routine must never make.
        feasible  = false;
        statusStr = 'numerical_problems';
    otherwise
        error('lyapunovToolbox:solverError', ...
              'Solver failed with status %d: %s', ...
              diagnostics.problem, yalmiperror(diagnostics.problem));
end

if ~feasible
    lambdaVal = [];
    sigmaVals = {};
    gramQ     = {};
    gramRes   = 0;
end

end


% =========================================================================
function [isBounded, evidence] = check_boundedness(V0, x, degV, opts)
%CHECK_BOUNDEDNESS  Is every sublevel set of V0 bounded?
%
%   Two tiers, matching is_sublevel_set_bounded in
%   julia/src/polynomial_utils.jl:
%
%     1. Structural, and CONCLUSIVE. For a quadratic V0 = x'Px the question is
%        exactly "P > 0", answerable from the Hessian with no solver.
%     2. Sampling fallback, NOT a proof. Points on a large sphere: if V0 <= 0
%        anywhere far out, the sublevel sets are certainly unbounded; if not,
%        we have found no counterexample and say only that.
%
%   The caller records which tier it got. Reporting tier 2 as proof is the
%   error this split exists to prevent.

n = length(x);

if degV == 2
    P = double(hessian(V0, x)) / 2;
    P = (P + P') / 2;
    if min(eig(P)) > 0
        isBounded = true;
        evidence  = 'top_form_pd';
        return;
    end
    % A quadratic with a non-positive-definite Hessian has unbounded (or empty)
    % sublevel sets. That is conclusive in the other direction.
    isBounded = false;
    evidence  = 'unbounded_witness';
    return;
end

% Higher degree: sample a range of radii, not a single sphere. One sphere can
% miss a narrow escape cone entirely, which would report an unbounded set as
% bounded -- the dangerous direction.
%
% The evidence must not vary run to run, so the draw is seeded -- but a library
% routine must not clobber the caller's stream, so the previous state is saved
% and restored. Callers downstream of this (pendulum_roa's sampling study) seed
% their own RNG and would otherwise inherit whatever was left here.
rngState = rng;
cleanupRng = onCleanup(@() rng(rngState));
rng(20260804, 'twister');
radii = opts.boundedRadius * [0.25, 0.5, 1.0];
perR  = max(1, floor(opts.boundedSamples / numel(radii)));
worst = Inf;
for r = radii
    for s = 1:perR
        u  = randn(n, 1);
        nu = norm(u);
        if nu == 0, continue; end
        pt = r * (u / nu);
        worst = min(worst, double(replace(V0, x, pt)));
    end
end

if worst <= 0
    isBounded = false;
    evidence  = 'unbounded_witness';
else
    isBounded = true;
    evidence  = 'sampling_only';
end

end
