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
| `matlab_vdp.csv` | `van_der_pol_roa('WriteFixture', true)` | `rho`, the coefficients of `V` with their monomial exponents, and the solver provenance |
| `matlab_pendulum.csv` | `pendulum_roa('WriteFixture', true)` | same schema |

Both are written by [`matlab/write_roa_fixture.m`](../../matlab/write_roa_fixture.m),
which is the single definition of the layout below.

## Schema

Each CSV fixture has the following columns, in order. The layout is **long
format**: a header row, then one row per monomial of $V$, with the scalar
fields repeated on every row so that any single row is self-describing and the
file reads with a plain CSV reader.

| Column | Type | Description |
|---|---|---|
| `rho` | Float64 | The certified sublevel-set level |
| `monomial` | String | Quoted exponent tuple in the variable order of `x`, e.g. `"[2,0]"` for $x_1^2$ |
| `coeff` | Float64 | Coefficient of that monomial in $V(x)$ |
| `degree` | Int | Degree of $V$ (redundant, for validation) |
| `epsilon` | Float64 | The $\\varepsilon$ margin used |
| `solver` | String | Solver actually used, as resolved by `compute_roa_yalmip`, e.g. `mosek` |
| `solver_status` | String | Terminal status at the returned level: `feasible` |
| `gram_residual` | Float64 | $\\|V - z^\\top Q z\\|_\\infty$ |
| `min_eig_Q` | Float64 | Minimum eigenvalue of the Gram matrix |
| `timestamp` | String | ISO-8601 timestamp of fixture generation |

The monomial/coefficient pairs are repeated for each term in $V$; the row
count therefore varies with the degree and number of variables. The Julia
reader in `test_roa_estimation.jl` parses these columns, reconstructs $V$
as a polynomial, and compares $\\rho$ and the coefficients to the Julia-side
result within a relative tolerance of $10^{-5}$.

**Any change to this schema must update both the MATLAB writer and the Julia
reader in the same commit.**

## When a fixture is absent

The Julia test emits `@info` naming the missing file and skips. It does **not**
pass silently — see the note in `tests/test_roa_estimation.jl`.
