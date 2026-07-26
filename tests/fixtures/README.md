# Cross-validation fixtures

Recorded output from the MATLAB/YALMIP/Mosek implementation, used by
`tests/test_roa_estimation.jl` to check that the Julia and MATLAB paths — which
share no code — agree.

CI has neither a MATLAB nor a Mosek licence, so these files cannot be
regenerated automatically. They are committed artefacts, and the commit that
adds or updates one **must** record in its message: the MATLAB version, the
YALMIP release, the solver and its version, and the machine.

## Files

| file | produced by | contents |
|---|---|---|
| `matlab_vdp.csv` | `van_der_pol_roa('WriteFixture', true)` | `rho`, the coefficients of `V` with their monomial exponents, and the full config struct flattened to columns |
| `matlab_pendulum.csv` | `pendulum_roa('WriteFixture', true)` | same schema |

## Schema

TODO(impl): document the exact column layout once `compute_roa_yalmip` is
written, and keep the Julia reader in `test_roa_estimation.jl` in sync with it.
A fixture whose schema is only defined by the code that writes it is a fixture
nobody can regenerate in two years.

## When a fixture is absent

The Julia test emits `@info` naming the missing file and skips. It does **not**
pass silently — see the note in `tests/test_roa_estimation.jl`.
