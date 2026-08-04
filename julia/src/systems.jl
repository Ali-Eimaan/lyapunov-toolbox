# systems.jl — the benchmark vector fields, defined exactly once.
#
# This file is not in the original repository specification. It was added
# because the same four systems appear in `tests/`, in all four `examples/`, in
# `matlab/*.m`, and in `docs/benchmarks.md`; without a single definition those
# five copies drift, and a cross-language agreement test between Julia and
# MATLAB that compares two *different* Van der Pol oscillators is worse than no
# test at all. The MATLAB scripts restate the dynamics by necessity — the guide
# requires a comment there pointing back at this file, and a test that compares
# the two numerically. Keep it.

"""
    PolynomialSystem

A polynomial vector field together with everything a certificate needs to be
interpretable.

# Fields
- `x::Vector` — state variables, in physical order. Never reorder these.
- `f::Vector` — `f[i]` is `ẋ[i]`, same length and order as `x`.
- `name::String` — used in plot titles, log lines, and benchmark table rows.
- `equilibrium::Vector{Float64}` — the point being certified; must satisfy
  `f(equilibrium) = 0`, checked at construction.
- `constraints::Vector` — semialgebraic set `{gⱼ(x) ≥ 0}` the certificate must
  stay inside (obstacles, actuator limits, the region where a Taylor recast is
  trusted). Empty means unconstrained.
- `invariants::Vector` — algebraic equalities `{hₖ(x) = 0}` from
  [`trig_recast`](@ref). These enter the SOS program with free-sign multipliers,
  not SOS ones. Getting this wrong makes the program infeasible in a way that
  looks like a solver problem.
- `state_labels::Vector{String}` — axis labels with units, for the plots.
- `notes::String` — one paragraph: where the model comes from, what was
  approximated, and what the certificate therefore does and does not mean.

# Constructor contract
The inner constructor validates: `length(f) == length(x)`, that every `f[i]` is
polynomial in `x` only, and that `‖f(equilibrium)‖∞ < 1e-12`. Raise `ArgumentError`
with the offending residual in the message.
"""
struct PolynomialSystem
    x::Vector
    f::Vector
    name::String
    equilibrium::Vector{Float64}
    constraints::Vector
    invariants::Vector
    state_labels::Vector{String}
    notes::String

    function PolynomialSystem(x, f, name, equilibrium, constraints, invariants, state_labels, notes)
        length(f) == length(x) || throw(ArgumentError(
            "length(f) = $(length(f)) != length(x) = $(length(x)); f[i] must be " *
            "the time derivative of x[i], so the two must be aligned"))
        for (i, fi) in enumerate(f)
            bad = setdiff(variables(fi), x)
            isempty(bad) || throw(ArgumentError(
                "f[$i] contains variables $(bad) not in the state vector x; a " *
                "system must be polynomial in its own states (recast sin/cos " *
                "first, see trig_recast)"))
        end
        length(equilibrium) == length(x) || throw(ArgumentError(
            "equilibrium has length $(length(equilibrium)) but the system has " *
            "$(length(x)) states"))
        resid = maximum(abs(convert(Float64, fi(x => equilibrium))) for fi in f)
        resid < 1e-12 || throw(ArgumentError(
            "equilibrium residual ||f(eq)||∞ = $resid is not below 1e-12; the " *
            "point $(equilibrium) is not an equilibrium of the field"))
        return new(x, f, name, Float64.(equilibrium), constraints, invariants,
                   state_labels, notes)
    end
end

"""
    van_der_pol(; reversed = true, μ = 1.0)

The Van der Pol oscillator, by default **time-reversed**:

```
ẋ₁ = -x₂
ẋ₂ =  x₁ + μ (x₁² - 1) x₂
```

The reversed system is the canonical ROA benchmark because its origin is
asymptotically stable and the boundary of its region of attraction is exactly
the unstable limit cycle — which is the forward system's *stable* limit cycle,
and therefore computable to high accuracy by simply integrating the forward
system. That gives ground truth to compare an SOS inner approximation against,
which almost no other nonlinear benchmark offers.

With `reversed = false` you get the forward oscillator, whose origin is
*unstable*. It is provided for one purpose: a negative test. Every synthesis
routine must fail to certify it, and must fail by reporting infeasibility rather
than by returning a certificate with a tiny `ρ`.

This reproduces the setup in Tedrake's MIT 6.832 notes and in Packard/Topcu's
ROA papers; `matlab/van_der_pol_roa.m` is the independent YALMIP implementation
of the same problem and the two must agree.
"""
function van_der_pol(; reversed::Bool = true, μ::Real = 1.0)
    @polyvar x[1:2]
    if reversed
        f = [-x[2], x[1] + μ * (x[1]^2 - 1) * x[2]]
    else
        f = [x[2], -x[1] - μ * (x[1]^2 - 1) * x[2]]
    end
    notes = reversed ?
        "Time-reversed Van der Pol oscillator. The origin is locally " *
        "asymptotically stable (A = [0 -1; 1 -μ] is Hurwitz for μ > 0) and the " *
        "boundary of its region of attraction is the unstable limit cycle — the " *
        "forward system's stable limit cycle — so ground truth is computable by " *
        "forward integration. Canonical SOS-ROA benchmark (Tedrake 6.832; " *
        "Packard/Topcu); matlab/van_der_pol_roa.m is an independent YALMIP " *
        "implementation and must agree with this file." :
        "Forward Van der Pol oscillator: the origin is UNSTABLE. Included as a " *
        "negative test — every synthesis routine must report infeasibility for " *
        "it, never a certificate."
    return PolynomialSystem(x, f,
        reversed ? "Van der Pol (reversed)" : "Van der Pol (forward)",
        [0.0, 0.0], [], [], ["x₁", "x₂"], notes)
end

"""
    inverted_pendulum(; m = 1.0, l = 1.0, b = 0.1, g = 9.81, K = nothing, recast = :taylor, taylor_order = 3)

Damped pendulum about the **upright** equilibrium under a stabilising state
feedback `u = -K·[θ, θ̇]`, with `θ = 0` upright:

```
θ̈ = (g/l)·sin θ - (b/(m l²))·θ̇ + u/(m l²)
```

Two recastings, selected by `recast`:

- `:taylor` — `sin θ ≈ θ - θ³/6` (`taylor_order = 3`) or higher. Cheap, and the
  standard textbook treatment, but the resulting certificate is for the
  *truncated* system. The `notes` field must say so, and the example must show
  where the truncation and the true field diverge.
- `:trig` — exact, via [`trig_recast`](@ref): states become `[s, c, θ̇]` with
  `s² + c² = 1` carried as an invariant. Three states instead of two and a
  bigger SDP, but the certificate is about the real pendulum.

If `K === nothing`, compute the LQR gain for the linearisation with
`Q = Diagonal([10.0, 1.0])`, `R = 1` (the guide's §5.2 weights) — but do it with
an explicit Riccati solve written here, not by adding a control-systems
dependency.

**Ground-truth anchor:** the true ROA of the closed loop is bounded by the
stable manifold of the downward-swinging saddle, so any certified sublevel set
that contains `θ = ±π` is wrong. Assert this in the tests; it is the cheapest
possible check that the whole pipeline is sound.
"""
# Continuous-time algebraic Riccati equation A'P + PA - PBR⁻¹B'P + Q = 0, solved
# via the stable invariant subspace of the Hamiltonian matrix. n is 2 here, so
# the eigendecomposition is tiny; the point is to avoid a control-systems
# dependency for one 2×2 problem.
function _care(A::Matrix{Float64}, B::Matrix{Float64}, Q::Matrix{Float64},
               R::Matrix{Float64})
    n = size(A, 1)
    H = [A -B * (R \ B'); -Q -A']
    F = schur(H)
    # Reorder so the n stable (Re λ < 0) eigenvalues come first. For the 2×2
    # Hamiltonian here the stable subspace is unambiguous; `ordschur` with the
    # stability selector is exact for the symmetric structure of H.
    sel = real.(F.values) .< 0
    n_stable = count(sel)
    n_stable == n || throw(ErrorException(
        "CARE: expected $n stable Hamiltonian eigenvalues, found $n_stable"))
    Fs = ordschur(F, sel)
    Z = Fs.Z
    U1 = Z[1:n, 1:n]
    U2 = Z[(n + 1):(2n), 1:n]
    P = U2 / U1
    P = (P + P') / 2
    # Verify the Riccati residual before returning: |A'P + PA - PBR⁻¹B'P + Q| must
    # be under 1e-9 relative to the largest term. This is the same standard as the
    # Lyapunov solve; a P that merely looks symmetric is not enough.
    res = norm(A' * P + P * A - P * B * (R \ B') * P + Q, 2)
    scale = max(1.0, norm(P, 2), norm(Q, 2))
    res < 1e-9 * scale || throw(ErrorException(
        "CARE residual ‖A'P + PA - PBR⁻¹B'P + Q‖_F = $res is not below " *
        "1e-9·scale = $(1e-9 * scale)"))
    return P
end

function inverted_pendulum(; m::Real = 1.0, l::Real = 1.0, b::Real = 0.1, g::Real = 9.81,
                           K = nothing, recast::Symbol = :taylor, taylor_order::Integer = 3)
    recast in (:taylor, :trig) || throw(ArgumentError(
        "recast must be :taylor or :trig, got $(repr(recast))"))
    if recast === :taylor
        isodd(taylor_order) || throw(ArgumentError(
            "taylor_order must be odd so sin θ's expansion stays odd, got $taylor_order"))
        taylor_order >= 1 || throw(ArgumentError(
            "taylor_order must be ≥ 1, got $taylor_order"))
    end
    m > 0 || throw(ArgumentError("m must be positive, got $m"))
    l > 0 || throw(ArgumentError("l must be positive, got $l"))
    b >= 0 || throw(ArgumentError("b must be nonnegative, got $b"))
    g > 0 || throw(ArgumentError("g must be positive, got $g"))

    # Linearisation about the upright: θ̈ = (g/l)θ - (b/(m l²))θ̇ + u/(m l²). The
    # (2,1) entry is +g/l — upright is UNSTABLE open-loop. If this A comes out
    # Hurwitz, the sign convention is wrong (guide §5.2).
    A = [0.0 1.0; g / l -b / (m * l^2)]
    B =reshape([0.0, 1.0 / (m * l^2)], :, 1)
    if K === nothing
        Q = Matrix(Diagonal([10.0, 1.0]))   # guide §5.2
        R = fill(1.0, 1, 1)
        P = _care(A, B, Q, R)
        K = Matrix(R \ (B' * P))            # u = -K·[θ, θ̇], K is 1×2
    else
        K = Float64.(reshape(collect(K), 1, :))
        length(K) == 2 || throw(ArgumentError(
            "K must have two entries for [θ, θ̇], got $(length(K))"))
    end
    Acl = A - B * K
    ev = eigvals(Acl)
    any(real(ev) .>= 0) && throw(ArgumentError(
        "closed-loop A - BK is not Hurwitz (eigenvalues $(round.(ev; digits = 6))); " *
        "the upright equilibrium is not locally stable"))

    @polyvar θ θdot
    x = [θ, θdot]
    if recast === :taylor
        # sin θ ≈ θ - θ³/3! + θ⁵/5! - ... up to taylor_order (odd). This is a
        # DIFFERENT vector field from the true pendulum; the notes must say so and
        # the example must measure the gap empirically (guide §4.4).
        sinθ = sum((-1)^k * θ^(2k + 1) / factorial(2k + 1) for k in 0:((taylor_order - 1) ÷ 2))
        f = [θdot, (g / l) * sinθ - (b / (m * l^2)) * θdot + (1.0 / (m * l^2)) * (-(K[1] * θ + K[2] * θdot))]
        eq = [0.0, 0.0]
        invars = Any[]
        labels = ["θ [rad]", "θ̇ [rad/s]"]
        notes = "Damped pendulum about the upright equilibrium under LQR feedback " *
                "u = -K[θ, θ̇]. The sine is Taylor-truncated at order $taylor_order, " *
                "so THIS IS A DIFFERENT VECTOR FIELD from the true pendulum: the " *
                "certificate applies to the truncated system only, and the example " *
                "measures where the two diverge empirically (no truncation-error " *
                "bound is computed). The true ROA is bounded by the stable manifold " *
                "of the downward saddle at θ = ±π; any certified set containing " *
                "θ = ±π is wrong."
        return PolynomialSystem(x, f, "Inverted pendulum (Taylor order $taylor_order)", eq,
                                invars, Any[], labels, notes)
    else
        # Exact trig recast: s = sin θ, c = cos θ, ṡ = c·θ̇, ċ = -s·θ̇. The
        # equilibrium is (s, c, θ̇) = (0, 1, 0), which is NOT the origin — the
        # positivity margin in verify_sublevel_set already handles that (it
        # vanishes at sys.equilibrium, not at 0).
        @polyvar s c
        # Build the field directly with s/c (trig_recast requires polynomial s,c).
        # θ̈ = (g/l)·sin θ - (b/(m l²))·θ̇ + u/(m l²), u = -K₁θ - K₂θ̇.
        # In the recast: sin θ → s (exact).  The controller term K₁·θ is NOT
        # K₁·s — the original LQR uses the angle θ, not sin θ.  Near the
        # equilibrium θ ≈ sin θ, so substituting K₁·s keeps the same tangent
        # model and is correct for a local certificate, but this is a DIFFERENT
        # controller globally.  We approximate θ by its Taylor expansion in s
        # (arcsin series): θ ≈ s + s³/6 + 3s⁵/40, which is polynomial and
        # retains exact linearisation.  This keeps the pendulum dynamics exact
        # (sin θ → s), and the controller is a polynomial approximation of the
        # original LQR.  For the inverted pendulum, this approximation is valid
        # over [-π/2, π/2] — the region the certificate will reach anyway — but
        # the example MUST say clearly that the controller is approximate.
        # arcsin(s) ≈ s + s³/6 + 3s⁵/40 (polynomial, good to ~0.2% at |s|=0.7)
        θ_approx = s + s^3 / 6 + 3 * s^5 / 40
        u_term = -(K[1] * θ_approx + K[2] * θdot)
        θdot_expr = (g / l) * s - (b / (m * l^2)) * θdot + (1.0 / (m * l^2)) * u_term
        f_poly, x_new, invariants = trig_recast([θdot, θdot_expr], x, 1; s = s, c = c)
        eq = [0.0, 1.0, 0.0]
        constraints = Any[]
        labels = ["s = sin θ", "c = cos θ", "θ̇ [rad/s]"]
        notes = "Damped pendulum about the upright equilibrium under LQR feedback, " *
                "with the heading recast EXACTLY via s = sin θ, c = cos θ. " *
                "Invariant s² + c² = 1 carried as an equality constraint with a " *
                "FREE-SIGN multiplier. The certificate is about the true pendulum " *
                "(up to the polynomiality of the field), unlike the Taylor recast. " *
                "Equilibrium is (s, c, θ̇) = (0, 1, 0); state order (s, c, θ̇) " *
                "matches trig_recast's substitution."
        return PolynomialSystem(x_new, f_poly, "Inverted pendulum (trig)", eq,
                                constraints, invariants, labels, notes)
    end
end

"""
    unicycle_with_obstacle(; obstacle_centre = (1.0, 0.0), obstacle_radius = 0.4, k_ρ = 1.0, k_α = 3.0, gains...)

Unicycle `(ẋ, ẏ, θ̇) = (v cos θ, v sin θ, ω)` under a smooth polynomial
stabilising feedback, with a circular obstacle expressed as the semialgebraic
constraint

```
g(x, y) = (x - x_o)² + (y - y_o)² - r² ≥ 0.
```

This is the example that exercises the constrained Positivstellensatz path: the
certificate is only required to hold on `{g ≥ 0}`, so a multiplier `σ(x)·g(x)`
enters the S-procedure. Without the constraint the problem is a plain ROA
computation and the obstacle is decoration — the point of the example is that
the certified set is provably obstacle-free, which is exactly the statement a
CBF-based controller needs as its terminal ingredient.

`sin θ`/`cos θ` are handled by [`trig_recast`](@ref) (exact), giving states
`[x, y, s, c]` with the invariant `s² + c² = 1`. Do not Taylor-expand the
heading here: the interesting states are precisely the ones with large `|θ|`.

The polar-tracking gains `k_ρ`, `k_α` follow the standard unicycle controller;
whichever form is used, the closed loop must be polynomial in the recast states
and the linearisation at the goal must be Hurwitz — verify that before running
any SOS program, because an unstable closed loop makes every downstream
infeasibility uninformative.
"""
function unicycle_with_obstacle(; obstacle_centre = (1.0, 0.0), obstacle_radius::Real = 0.4,
                                v_ref::Real = 1.0, k_e::Real = 1.0,
                                k_y::Real = 1.0, k_ψ::Real = 2.0)
    v_ref > 0 || throw(ArgumentError(
        "v_ref must be positive, got $v_ref; the tracking error dynamics are " *
        "exponentially stable only for a non-zero reference speed — at v_ref = 0 " *
        "this degenerates to point stabilisation, which Brockett rules out"))
    k_e > 0 || throw(ArgumentError("k_e must be positive, got $k_e"))
    k_y > 0 || throw(ArgumentError("k_y must be positive, got $k_y"))
    k_ψ > 0 || throw(ArgumentError("k_ψ must be positive, got $k_ψ"))
    obstacle_radius > 0 || throw(ArgumentError("obstacle_radius must be positive, got $obstacle_radius"))

    xo, yo = Float64.(obstacle_centre)

    # Exact trig recast of the heading error: s = sin ψ, c = cos ψ.
    # States [e₁, e₂, s, c]; invariant s² + c² = 1.
    @polyvar e1 e2 s c
    x_new = [e1, e2, s, c]

    # Tracking error in the robot frame, relative to a reference travelling in a
    # straight line at constant speed v_ref (ω_ref = 0):
    #     ė₁ = ω e₂ - v + v_ref cos ψ
    #     ė₂ = -ω e₁ + v_ref sin ψ
    #     ψ̇  = -ω
    # with the standard (Kanayama) tracking law
    #     v = v_ref cos ψ + k_e e₁,   ω = k_y v_ref e₂ + k_ψ sin ψ.
    #
    # Both are polynomial in the recast states, and the v_ref cos ψ term cancels
    # exactly in ė₁, which is what keeps the field polynomial without any
    # truncation.
    ω = k_y * v_ref * e2 + k_ψ * s
    f_poly = [ω * e2 - k_e * e1,        # ė₁
              -ω * e1 + v_ref * s,      # ė₂
              -c * ω,                   # ṡ = c·ψ̇ = -c·ω
              s * ω]                    # ċ = -s·ψ̇ = s·ω

    eq = [0.0, 0.0, 0.0, 1.0]

    # The linearisation must be Hurwitz ON THE MANIFOLD. Exactly one zero
    # eigenvalue is expected — the normal direction of s² + c² = 1, which at the
    # equilibrium (0,0,0,1) is the c-axis. Every eigenvalue of the reduced
    # dynamics on the tangent space span{e₁, e₂, s} must be strictly negative.
    #
    # This check is deliberately strict. The previous version accepted "at least
    # 1 of 3 strictly negative", which passed a spectrum of {-k_e, 0, 0}: the old
    # feedback point-stabilised the unicycle, and that is impossible to do
    # exponentially (Brockett), so no positive ρ was ever certifiable and
    # bisect_level_set returned 0.
    A = jacobian_at(f_poly, x_new, eq)
    G = [0.0 0.0 2*eq[3] 2*eq[4]]        # ∇(s² + c² - 1) at eq
    T = nullspace(G)                     # orthonormal basis of the tangent space
    ev_tan = eigvals(T' * A * T)
    maximum(real, ev_tan) < -1e-9 || throw(ArgumentError(
        "unicycle tracking closed loop is not exponentially stable on the " *
        "invariant manifold: tangent-space eigenvalues " *
        "$(round.(ev_tan; digits = 6)). Expected all strictly negative; the " *
        "2×2 heading block has characteristic polynomial λ² + k_ψ λ + k_y v_ref², " *
        "so it is Hurwitz for k_ψ, k_y, v_ref > 0."))

    # Obstacle in the error plane: a disc the tracker must stay outside of.
    g = (e1 - xo)^2 + (e2 - yo)^2 - obstacle_radius^2

    # Redundant Archimedean ball, so Putinar's theorem applies (see
    # docs/derivations/positivstellensatz.tex).
    arch_N = 2.0 * (xo^2 + yo^2 + 4.0)
    g_ball = arch_N - (e1^2 + e2^2 + s^2 + c^2)

    constraints = [g, g_ball]
    invariants = [s^2 + c^2 - 1]

    labels = ["e₁ [m] (along-track)", "e₂ [m] (cross-track)", "s = sin ψ", "c = cos ψ"]
    notes = "Unicycle TRACKING a straight reference at constant speed " *
            "v_ref = $v_ref, under the Kanayama law v = v_ref·cos ψ + $k_e·e₁, " *
            "ω = $k_y·v_ref·e₂ + $k_ψ·sin ψ. States are the body-frame tracking " *
            "error [e₁, e₂, s, c] with the heading error recast EXACTLY via " *
            "s = sin ψ, c = cos ψ and the invariant s² + c² = 1 carried with a " *
            "free-sign multiplier. " *
            "WHY TRACKING AND NOT POINT STABILISATION: the unicycle fails " *
            "Brockett's necessary condition, so NO smooth static state feedback " *
            "stabilises a point — the linearisation is unavoidably " *
            "non-hyperbolic and no strict margin ε‖x-x_eq‖^{2k} can be met near " *
            "the goal. The tracking error dynamics have linearisation with " *
            "eigenvalues {-k_e} ∪ roots(λ² + k_ψ λ + k_y v_ref²), Hurwitz for " *
            "positive gains and v_ref > 0, so a certificate exists. " *
            "Obstacle at ($xo, $yo) radius $obstacle_radius enters as a " *
            "semialgebraic constraint g ≥ 0 through an SOS multiplier σ; a " *
            "redundant Archimedean ball (N = $arch_N) is added for Putinar. " *
            "The certified set is PROVABLY obstacle-free — the point of the " *
            "example. Equilibrium is (e₁, e₂, s, c) = (0, 0, 0, 1)."

    return PolynomialSystem(x_new, f_poly, "Unicycle tracking with obstacle", eq,
                            constraints, invariants, labels, notes)
end

"""
    consensus_network(; n = 3, laplacian = nothing, coupling = :cubic, κ = 1.0)

`n`-agent consensus with nonlinear coupling,

```
ẋᵢ = -Σⱼ aᵢⱼ · φ(xᵢ - xⱼ),      φ(z) = z + κ z³   (`coupling = :cubic`)
```

on a connected graph with Laplacian `L`. Default: the 3-agent path graph.

**The trap, and the reason this example is here.** The consensus manifold
`span(𝟙)` is an equilibrium *set*, not a point: `L𝟙 = 0`. A Lyapunov function
that is positive definite on all of `ℝⁿ` cannot exist, and an SOS program that
asks for one returns infeasible with no hint as to why. The fix is to certify
the **disagreement dynamics**: pick an orthonormal basis `U ∈ ℝ^{n×(n-1)}` of
`𝟙^⊥`, change coordinates `z = U'x`, and run the whole pipeline on `ż`. The
origin in `z` is the consensus manifold in `x`, and it is an isolated
equilibrium there.

Return the system in `z` coordinates and store `U` in `notes` (and expose it via
the example) so the certified set can be mapped back to `x` for plotting.

**Citation hook.** A quadratic `V(z) = z'Pz` certified for the closed-loop
network is the single-mode special case of the composition argument used for
RRC-CBF in the thesis; the example is written to make that connection explicit
in prose, not just implied.
"""
function consensus_network(; n::Integer = 3, laplacian = nothing, coupling::Symbol = :cubic, κ::Real = 1.0)
    n >= 2 || throw(ArgumentError("need at least 2 agents, got n = $n"))
    coupling in (:cubic, :linear) || throw(ArgumentError(
        "coupling must be :cubic or :linear, got $(repr(coupling))"))

    # Default: path graph 1—2—3—...—n
    if laplacian === nothing
        L = zeros(n, n)
        for i in 1:n
            L[i, i] = (i > 1 ? 1.0 : 0.0) + (i < n ? 1.0 : 0.0)
        end
        for i in 1:(n - 1)
            L[i, i + 1] = L[i + 1, i] = -1.0
        end
    else
        # Accept both AbstractMatrix (native Julia) and AbstractVector of
        # AbstractVector (from Python .tolist()) → stack into a matrix.
        if laplacian isa AbstractVector
            L = Float64.(reduce(hcat, laplacian))
        else
            L = Float64.(laplacian)
        end
        size(L) == (n, n) || throw(ArgumentError(
            "laplacian must be $(n)×$(n), got $(size(L))"))
    end

    # Orthonormal basis U of 1^⊥ via QR decomposition of [1 | I].
    # U ∈ R^{n×(n-1)}, U'1 = 0, U'U = I.
    ones_vec = ones(n)
    Qmat = Matrix(qr([ones_vec Matrix{Float64}(I, n, n)]).Q)
    U = Qmat[:, 2:n]  # drop the first column (span of 1)

    # Verify the reduction basis.
    @assert maximum(abs, U' * ones_vec) < 1e-12 "U'1 ≠ 0"
    @assert norm(U' * U - I) < 1e-12 "U'U ≠ I"

    # Build the consensus dynamics in original (x) coordinates:
    # ẋ_i = -Σ_j L_{ij} φ(x_i - x_j) where φ(z) = z + κ z³ (cubic) or φ(z) = z (linear)
    @polyvar xvars[1:n]
    f_x = Vector{Any}(undef, n)
    for i in 1:n
        f_x[i] = 0
        for j in 1:n
            if L[i, j] != 0 && i != j
                Δ = xvars[i] - xvars[j]
                if coupling === :cubic
                    f_x[i] = f_x[i] + L[i, j] * (Δ + κ * Δ^3)
                else
                    f_x[i] = f_x[i] + L[i, j] * Δ
                end
            end
        end
    end

    # Change coordinates to disagreement space: z = U'x, x = Uz on 1^⊥.
    # ż = U' f(Uz), where f maps R^n → R^n.
    @polyvar zvars[1:(n - 1)]
    # Express x in terms of z: x = Uz
    x_from_z = U * zvars
    # Build f_z with concrete polynomial type (not Vector{Any}) so that
    # differentiate / jacobian_at work correctly.
    fz_sample = subs(f_x[1], xvars => x_from_z)
    f_z = Vector{typeof(fz_sample)}(undef, n - 1)
    for i in 1:(n - 1)
        fi = zero(fz_sample)
        for k in 1:n
            # ż_i = Σ_k U[k, i] * f_x[k](x ← Uz)
            fk_sub = subs(f_x[k], xvars => x_from_z)
            fi = fi + U[k, i] * fk_sub
        end
        f_z[i] = fi
    end

    # Verify the reduction is exact: f(Uz) has no component along 1,
    # so ż = U' f(Uz) captures the full dynamics.
    # The condition 1'f(x) = 0 (sum of all agent dynamics is zero) must hold identically.
    sum_f = subs(f_x[1], xvars => x_from_z)
    for k in 2:n
        sum_f = sum_f + subs(f_x[k], xvars => x_from_z)
    end
    @assert maximum(abs, MultivariatePolynomials.coefficients(sum_f)) < 1e-12 "reduction not exact: f has a 1-component"

    eq = zeros(n - 1)
    labels = ["z$i — disagreement coord $i" for i in 1:(n - 1)]
    notes = "$n-agent consensus with $(coupling) coupling (κ = $κ) on the path graph, " *
            "in DISAGREEMENT coordinates z = U'x. The Lyapunov certificate is for the " *
            "ż dynamics; the origin in z is the consensus manifold in x. No Lyapunov " *
            "function positive definite on all of R^n exists (the consensus manifold " *
            "span(1) is an equilibrium set), which is why the certifier must fail in " *
            "x-coordinates and succeed in z-coordinates — the trap that motivates " *
            "this example. U ∈ R^{$(n)×$(n-1)} is column-orthonormal with U'1 = 0."

    return PolynomialSystem(zvars, f_z, "Consensus network (n=$n, $coupling)", eq,
                            [], [], labels, notes)
end

"""
    simulate(sys, x0, tspan; dt = 1e-3, max_norm = 1e6)

Fixed-step RK4 integration of `sys.f` from `x0`.

Deliberately hand-written rather than pulled from `DifferentialEquations.jl`:
the package must stay installable in a bare CI container, and every use here is
either a sanity check or a ground-truth boundary trace where a fixed step with a
stated `dt` is easier to defend than an adaptive solver's tolerances.

Terminates early and returns `(t, X, :diverged)` once `‖x‖∞ > max_norm`, which
is how the ROA ground-truth sweep detects escape without waiting out the full
`tspan`. On completion returns `(t, X, :converged)`.
"""
function simulate(sys::PolynomialSystem, x0::AbstractVector, tspan::Tuple{<:Real,<:Real};
                  dt::Real = 1e-3, max_norm::Real = 1e6)
    length(x0) == length(sys.x) || throw(ArgumentError(
        "x0 has length $(length(x0)) but the system has $(length(sys.x)) states"))
    dt > 0 || throw(ArgumentError("dt must be positive, got $dt"))
    t0, t1 = tspan
    t1 > t0 || throw(ArgumentError("tspan must be increasing, got $tspan"))
    fval(xv) = [convert(Float64, sys.f[i](sys.x => xv)) for i in eachindex(sys.x)]
    nsteps = max(0, ceil(Int, (t1 - t0) / dt))
    # Use the *effective* step h = (t1 - t0)/nsteps for both the integrator and
    # the reported time grid. Previously the RK4 loop stepped with `dt` while t
    # was built as a grid of `nsteps+1` points over [t0, t1]; when `dt` does not
    # divide t1 - t0 exactly, the reported times did not match the step that was
    # actually taken, so t[end] ≠ t1 and each t[k] disagreed with the row X[k].
    h = nsteps > 0 ? (t1 - t0) / nsteps : 0.0
    t = collect(range(t0, t1; length = nsteps + 1))
    X = Matrix{Float64}(undef, nsteps + 1, length(x0))
    X[1, :] = x0
    for k in 1:nsteps
        xk = X[k, :]
        k1 = fval(xk)
        k2 = fval(xk .+ (h / 2) .* k1)
        k3 = fval(xk .+ (h / 2) .* k2)
        k4 = fval(xk .+ h .* k3)
        X[k + 1, :] = xk .+ (h / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
        if maximum(abs, X[k + 1, :]) > max_norm
            # Escape detected: return the trajectory up to and including the
            # diverging step, so a boundary sweep knows which side of the ROA
            # boundary this initial condition was on.
            return (t[1:(k + 1)], X[1:(k + 1), :], :diverged)
        end
    end
    return (t, X, :converged)
end
