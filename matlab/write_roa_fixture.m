function path = write_roa_fixture(name, V, x, rho, info, epsilon)
%WRITE_ROA_FIXTURE  Write a Julia/MATLAB cross-validation fixture.
%
%   PATH = WRITE_ROA_FIXTURE(NAME, V, X, RHO, INFO, EPSILON) writes
%   tests/fixtures/<NAME>.csv in the schema documented in
%   tests/fixtures/README.md and returns the path it wrote.
%
%   WHY THIS IS A COMMITTED ARTEFACT. CI has neither a MATLAB nor a Mosek
%   licence, so this file cannot be regenerated automatically. The commit that
%   adds or updates one MUST record, in its message, the MATLAB version, the
%   YALMIP release, the solver and its version, and the machine. That message
%   is the only provenance the fixture has.
%
%   LAYOUT (long format; one row per monomial of V, header row first):
%
%     rho,monomial,coeff,degree,epsilon,solver,solver_status,
%     gram_residual,min_eig_Q,timestamp
%
%   The scalar fields repeat on every row so that any single row is
%   self-describing and the file can be read with a plain CSV reader without
%   tracking state. `monomial` is a quoted, bracketed, comma-separated exponent
%   tuple in the variable order of X -- e.g. "[2,0]" is x1^2 for X = [x1; x2].
%
%   The exponent tuple, not a printed monomial string, is what makes the
%   fixture machine-comparable: the Julia reader reconstructs V from these rows
%   and compares coefficient by coefficient against its own result, which is a
%   stronger check than comparing rho alone.
%
%   ANY CHANGE TO THIS LAYOUT MUST UPDATE tests/fixtures/README.md AND THE
%   JULIA READER IN tests/test_roa_estimation.jl IN THE SAME COMMIT.
%
%   See also: COMPUTE_ROA_YALMIP, VAN_DER_POL_ROA, PENDULUM_ROA.

narginchk(6, 6);

fixtureDir = fullfile(fileparts(mfilename('fullpath')), '..', 'tests', 'fixtures');
if ~exist(fixtureDir, 'dir')
    mkdir(fixtureDir);
end
path = fullfile(fixtureDir, [name '.csv']);

[coeffs, monos] = coefficients(V, x);
E = monomial_exponents(monos, x);
nTerms = numel(coeffs);
assert(size(E, 1) == nTerms, ...
       'exponent matrix has %d rows but V has %d terms', size(E, 1), nTerms);

degV      = degree(V);
solverStr = field_or(info, 'solverResolved', 'unknown');
statusStr = field_or(info, 'solverStatus',  'unknown');
gramRes   = field_or(info, 'gramResidual',  NaN);
minEigQ   = field_or(info, 'minEigQ',       NaN);
stamp     = iso_timestamp();

fid = fopen(path, 'w');
if fid < 0
    error('lyapunovToolbox:fixtureWriteFailed', ...
          'Could not open %s for writing.', path);
end
cleaner = onCleanup(@() fclose(fid));

fprintf(fid, ['rho,monomial,coeff,degree,epsilon,solver,solver_status,' ...
              'gram_residual,min_eig_Q,timestamp\n']);
for i = 1:nTerms
    monoStr = ['"[' strjoin(arrayfun(@(e) sprintf('%d', e), E(i, :), ...
                                     'UniformOutput', false), ',') ']"'];
    fprintf(fid, '%.17g,%s,%.17g,%d,%.17g,%s,%s,%.17g,%.17g,%s\n', ...
            rho, monoStr, double(coeffs(i)), degV, epsilon, ...
            solverStr, statusStr, gramRes, minEigQ, stamp);
end

fprintf('Fixture written to %s (%d terms)\n', path, nTerms);

end


% =========================================================================
function E = monomial_exponents(monos, x)
%MONOMIAL_EXPONENTS  Exponent matrix of a YALMIP monomial vector.
%
%   Row i, column j is the exponent of x(j) in monos(i). Tries YALMIP's
%   EXPONENTS first and falls back to per-element DEGREE, because the exact
%   signature of EXPONENTS has varied across YALMIP releases and this file has
%   to keep working across them.

nVars  = numel(x);
nMonos = numel(monos);

try
    E = exponents(monos, x);
    if size(E, 1) == nMonos && size(E, 2) == nVars
        return;
    end
catch
    % fall through
end

E = zeros(nMonos, nVars);
for i = 1:nMonos
    for j = 1:nVars
        E(i, j) = degree(monos(i), x(j));
    end
end

end


% =========================================================================
function v = field_or(s, name, default)
if isfield(s, name)
    v = s.(name);
else
    v = default;
end
end


% =========================================================================
function s = iso_timestamp()
%ISO_TIMESTAMP  ISO-8601 local timestamp, without the deprecated DATESTR.
try
    s = char(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss'));
catch
    s = datestr(now, 'yyyy-mm-ddTHH:MM:SS'); %#ok<DATST,TNOW1>
end
end
