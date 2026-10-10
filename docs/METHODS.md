# Mathematical model and numerical method

## 1. Scope and normalization

We solve the normal-phase, elastic, nonrelativistic UU model in three velocity dimensions. The statistical sign is +1 for Bose, -1 for Fermi, and 0 for classical statistics. The angular kernel is isotropic. No condensate is represented.

Let `v0 = sqrt(k_B T / m)` and let `n` be the physical equilibrium number density. Define `v_hat = (v-u)/v0` and `f_hat = v0³ f_number/n`. Thus `integral f_hat dv_hat = 1`. The polynomial argument is the dimensionless energy `E=|v_hat|²/2`.

For `dσ/dΩ=d²/4`, use length unit `1/(n d²)` and time unit `1/(n d² v0)`. With diameter parameter 1 in the implementation:

```
mu_physical     = mu_code * sqrt(m k_B T) / d²
lambda_physical = lambda_code * k_B sqrt(k_B T/m) / d²
```

The CLI's optional `--diameter` multiplies the reference cross section; both coefficients must scale as its inverse square, unless `--strength` explicitly overrides the kernel prefactor. These restoration factors are for the specified effective model, not for every atomic-spin scattering convention.

Write `E=|v_hat|²/2` and `D(E)=1-epsilon*z*exp(-E)`. The dimensionless number-normalized equilibrium is

```
f0(v) = exp(-E) / (Z_norm D(E)),
Z_norm = integral exp(-E)/D(E) dv,
A_hat = z Z_norm,
r(v) = 1/D(E),
w(v) = f0(v) r(v).
```

This gives `A_hat f0 = 1/(z^-1 exp(E)-epsilon)`. The dimensionless degeneracy per internal state is `G_(3/2)(z)` and `A_hat=(2 pi)^(3/2) G_(3/2)(z)`.

Physical `h`, particle mass and internal degeneracy are not changed by the code. Choosing z selects a dimensionless degeneracy; converting a scan to a fixed real gas requires choosing a physical density/temperature path. In particular, matching Bose and Fermi points at equal z does not mean equal physical degeneracy.

## 2. Collision brackets

The classical hard-sphere convention is

```
B = |v-v1| d² / 4,  integral_S² B dOmega = pi d² |v-v1|.
```

With `W=f0 f01 r' r1'` and `Delta phi=phi+phi1-phi'-phi1'`,

```
[phi,psi] = (1/4) integral B W Delta(phi) Delta(psi) dv dv1 dOmega.
```

The thermal and shear components used for the solves are

```
Phi_n^T = P_n(E) v_x,
Phi_n^S = Q_n(E) v_x v_y.
```

The thermal polynomials begin with `L_1^(3/2)`. Subtract their weighted mean with measure `w |v|² dv` to enforce momentum matching. Shear polynomials begin with `L_0^(5/2)`. Their matching is automatic by angular symmetry.

Using `J_x=v_x(E-5p/2)` and `n=T=m=k_B=1`:

```
M^T_mn = [Phi_m^T, Phi_n^T],   d^T_m = integral w Phi_m^T J_x dv,
M^S_mn = [Phi_m^S, Phi_n^S],   d^S_m = integral w Phi_m^S v_x v_y dv.
M^T alpha=d^T,                 M^S beta=d^S.
lambda = (d^T)^T alpha,        mu = (d^S)^T beta.
```

Moment formulas provide an algebraically equivalent check. Neither transport coefficient is supplied as a fit parameter.

## 3. Reduction of the collision integral

Introduce the center velocity `V=(v+v1)/2` and half-relative velocity `q=(v-v1)/2`. For equal masses,

```
v=V+q, v1=V-q, v'=V+q', v1'=V-q', |q'|=|q|.
```

Energy and momentum conservation are exact at the level of this geometry. The equilibrium detailed-balance identity is verified separately in tests.

Global rotations can be integrated analytically: take V along the z axis, the incoming q in the xz plane, and integrate the two polar cosines and the relative outgoing azimuth. Together with `x=|V|²` and `y=|q|²`, this leaves five integration variables.

For the unit-density Maxwellian `phi(v)=(2pi)^(-3/2) exp(-|v|²/2)`, V and q are independent Gaussian vectors with component variance 1/2. Dividing W by `phi(v)phi(v1)` gives

```
R = 1 / [(G_(3/2)(z)/z)² D(E)D(E1)D(E')D(E1')].
```

The radial measures after including relative speed and the bracket prefactor are

```
2 d² * x^(1/2) y^1 exp(-x-y) dx dy
```

times the normalized angular measures `d cos(theta_in)/2`, `d cos(theta_out)/2`, `dphi/(2pi)` and the factor R.

Generalized Gauss--Laguerre integrates x and y, Gauss--Legendre integrates the polar cosines, and a periodic trapezoid rule integrates the azimuth. The implementation has no hard velocity-domain cutoff, but quadrature resolution still limits the accuracy of sharp degenerate distributions.

After global rotational averaging, a single vector component bracket equals one third of the sum over the three vector components. A single off-diagonal traceless-tensor component bracket equals one tenth of the sum over all tensor components. The latter follows from isotropy on the five-dimensional symmetric traceless representation, including the two equal off-diagonal entries in the Frobenius sum.

## 4. Independent quantum variational benchmark

`wu_reference` in `src/collisions.cpp` evaluates Eq. (50) of Wu arXiv:1810.03090v1 through a different, four-dimensional quadrature. Its reduced integrals I_A and I_B are evaluated directly, not by calling the collision-bracket solver.

The code divides I_A and I_B by z² analytically before integration. Their classical limits in this convention are `sqrt(pi)/2` and `sqrt(pi)`. The G_nu factors are likewise divided by z to avoid unnecessary underflow.

The comparison uses **ratios to the classical first-profile limits**. Constant spin/cross-section multiplicity factors in the original physical Bose and balanced Fermi cases cancel in these ratios. This check validates the one-profile quantum integrals. It does not certify higher-profile convergence, reproduce the full Fourier solver, or compare absolute viscosities of identical real gases.

## 5. Homogeneous dynamics

For a finite perturbation basis write `chi=sum c_n(t)Phi_n`. The homogeneous linearized equation is

```
H c_dot = -M c,
H_mn = integral w Phi_m Phi_n dv.
```

Use Cholesky whitening of H and diagonalize the resulting symmetric collision matrix. Positive generalized eigenvalues give decay and monotonicity of the quadratic perturbation entropy `c^T H c/2`.

With initial coefficients `c(0)=H^-1 d`, the response `d^T c(t)` has integral `d^T M^-1 d`, equal to the corresponding steady transport coefficient in this normalization. This is an internal response identity; it is not an independent validation and is not a derivation of the full quantum many-body Kubo formula.

Polynomial perturbations need not keep `f0+w chi` nonnegative at every velocity for finite amplitudes. The implementation is a tangent/linear-response calculation, and does not claim a positivity-preserving nonlinear evolution.

## 6. Spatial kinetic shear wave

Retain a conserved transverse momentum mode `v_x`, the matched thermal-vector profiles `P_n(E)v_x`, and shear profiles `Q_n(E)v_x v_y`. Allow their coefficients to vary along y. Projection of the linearized streaming equation gives

```
H c_t + S c_y = -M c,
S_mn = integral w Phi_m v_y Phi_n dv.
```

For a spatial Fourier mode, `c(y,t)=Re(c_k(t) exp(iky))`,

```
H c_k_dot = -(M+ikS)c_k.
```

The code solves this finite-dimensional evolution by matrix eigendecomposition. It does not time-step a Navier--Stokes equation with a prescribed viscosity. Viscosity is inferred from the slow eigenvalue s(k) as `mu_dynamic=-rho Re(s(k))/k²`. As k tends to zero it approaches the transport viscosity.

For these transverse modes, perturbations of mass and energy vanish by parity. The momentum basis is explicitly retained, unlike the matched homogeneous collision-only solves.

Important limit: only vector and shear angular sectors are retained. Streaming also excites higher angular sectors in the full kinetic equation. Therefore finite-k dispersion is a model truncation; small-k agreement is a consistency check of the projected hydrodynamic limit, not full finite-Knudsen validation.
