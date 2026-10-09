# From the collision integral to measured sound: the Hilker et al. comparison

This note explains how the Galerkin transport coefficients are turned into a prediction for first sound in a normal (non-condensed) 3D Bose gas, and how that prediction is compared with Fig. 4 of Hilker et al., PRL 128, 223601 (2022). Only the five observations with L/ℓ ≥ 3 are used: below that threshold the gas is not hydrodynamic and Navier–Stokes–Fourier modes do not apply. The condensate, second sound and the crossover to the collisionless regime are not modelled.

Sources: [paper and supplement](https://arxiv.org/abs/2112.14763), [original data, CC BY 4.0](https://doi.org/10.17863/CAM.83837). The published theory curve is used only as a separate check of the sound-mode algebra, not as data.

## 1. Power-law collision kernel

The solver accepts an isotropic kernel B(g) = C g^γ, with g = |v − v₁|, C > 0 and γ > −3. With V = (v+v₁)/2, q = (v−v₁)/2, x = V², y = q² and g = 2√y, the radial weight in y becomes y^((1+γ)/2) exp(−y) and the bracket prefactor is 4·2^γ·C. For γ = 1 and C = d²/4 this reduces to the hard-sphere weight y·exp(−y) and prefactor 2d².

Classical first-Sonine values for the same kernel are

```
S_γ = 4C · 2^(γ−1) · Γ((γ+7)/2) / 6,   η₁ = 5 / (16 √π S_γ),   κ₁ = (15/4) η₁
```

in code units n = m = k_B T = 1. The second Sonine approximation has closed forms (Chapman–Cowling): 1 + 3/202 and 1 + 1/44 for γ = 1, and 1 + 3/190 and 1 + 1/40 for γ = −1. Both are unit tests.

## 2. Physical kernel for identical bosons

For identical bosons in the same internal state, at low energy and with s-wave scattering, σ = 8πa². Integrating over the full sphere this corresponds to dσ/dΩ = 2a², i.e. B = 2a² g: the γ = 1 case with d² = 8a². The scattering length a is taken from the dataset and is never adjusted. From a result computed with C = 1/4:

```
η = η_code · √(m k_B T) / (8a²),     κ = κ_code · k_B √(k_B T / m) / (8a²)
```

The classical check is η₁ = 5√(m k_B T / π) / (128 a²). Bose enhancement is already inside the collision integral through the (1 + f) factors; the cross section is not multiplied by an extra stimulation factor. The mass is that of ³⁹K. The outputs record n a³, g n/(k_B T) and a√(m k_B T)/ħ to show the physical regime; a small value of the last parameter only roughly justifies an energy-independent cross section.

## 3. Reconstructing the state of each measurement

The dataset gives L, a, T and I = L/ℓ. With the paper's definition of the mean free path,

```
n = I / (8π a² L),     λ_T = √(2πħ² / (m k_B T)),     n λ_T³ = G_{3/2}(z)
```

and the last relation is inverted for z. The Bose functions G_s(z) = Σ_j z^j / j^s are summed directly with a tail-error check, independently of the solver's equilibrium quadrature. The density is therefore derived from the published abscissa, not an independent measurement.

z = exp[(μ − 2g n)/(k_B T)] with g = 4πħ²a/m is the effective Hartree–Fock fugacity.

## 4. Kinetic and total pressure

The velocity moment gives the kinetic pressure p_kin = n k_B T G_{5/2}/G_{3/2}. The mean field adds

```
p = p_kin + g n²,     E/V = (3/2) p_kin + g n²
```

so E/V = 3p/2 does not hold for the total pressure. The thermodynamic derivatives are

```
p_n = k_B T G_{3/2}/G_{1/2} + 2 g n
p_T = n k_B ( (5/2) G_{5/2}/G_{3/2} − (3/2) G_{3/2}/G_{1/2} )
C_V = (3/2) p_T,     C_P = C_V + T p_T² / (n p_n)
c_ad² = (C_P / C_V) p_n / m = 5 p_kin / (3 m n) + 2 g n / m
```

with C_V and C_P per unit volume. The tests compare the derivatives of p and E and the two expressions for c_ad.

## 5. Why the Galerkin profiles can be reused

In the bulk the mean field acts as a force X = −2g∇n/m. Euler's equation gives

```
D_t u = −∇(p_kin + g n²) / (m n),     so     D_t u − X = −∇p_kin / (m n)
```

The mean-field force cancels the g n² contribution in the Chapman–Enskog source, and the internal kinetic energy obeys the usual adiabatic equation (D_t ln z = 0, D_t ln T = −(2/3)∇·u at first order). The thermal source stays U_i (c² − 5 Z_kin / 2), with Z_kin = p_kin / (n k_B T), and the viscous source stays the traceless symmetric tensor. The local energy shift 2gn is velocity-independent and cancels in each collision's energy balance. The UU matrices and radial profiles can therefore be used at the effective fugacity. This holds for the normal Hartree–Fock kinetic model with the vacuum cross section; it adds no critical or in-medium corrections.

## 6. From transport coefficients to measured sound

The linearized Navier–Stokes–Fourier equations are

```
∂_t δn + n ∂_x u = 0
m n ∂_t u = −∂_x (p_n δn + p_T δT) + (4/3) η ∂_x² u
C_V ∂_t δT = −T p_T ∂_x u + κ ∂_x² δT
```

With k = π/L and exp[i(kx − ωt)], define Q_η = 4ηk / (3 m n c_ad), Q_κ = κk / (C_V c_ad) and R = C_V / C_P. The dispersion relation is

```
C³ + i(Q_η + Q_κ) C² − (1 + Q_η Q_κ) C − i Q_κ R = 0,     C = ω / (k c_ad)
```

(Eq. S52 of the supplement). The code checks that its roots match the eigenvalues of the 3×3 system built directly from the balance equations, that it reduces to Euler sound without dissipation, and that it gives the Stokes–Kirchhoff damping γ ≈ (k²/2)[(4/3)η/ρ + κ(1/C_V − 1/C_P)] at small k.

The dataset's c_I is |ω|/k, not the phase velocity Re(ω)/k; the two differ when damping is large. The damping is γ_I = −Im(ω). The experimental damping column is in units of 2π·Hz and is multiplied by 2π to obtain s⁻¹.

## 7. Three separate comparisons

1. **Published curve.** z = 0.75, C_V/C_P = 0.5, q_η = 0.73, Pr = 0.67 inserted into S52 reproduce the authors' theory file within 0.08% for L/ℓ ≥ 3 (differences come from the rounded published parameters). This checks the algebra, not the physics.
2. **Fixed-state Galerkin curve.** The same curve with our N = 7 transport coefficients at z = 0.75: speed within ~2% and damping ~6.5% above the published curve for L/ℓ ≥ 3.
3. **Point-by-point prediction (main result).** z and Hartree–Fock thermodynamics reconstructed for each observation, N = 7 coefficients, refined quadrature, compared with the measured speed and damping. Two model variants are also shown: classical transport with the same Hartree–Fock thermodynamics, and Bose transport without the mean field.

## 8. Result, convergence and uncertainties

- **Speed:** RMS residual / experimental error = 1.03, i.e. agreement within errors. The speed is set mainly by thermodynamics, so this tests the state reconstruction rather than the collision integral.
- **Damping:** signed differences from −21% to +16%; RMS residual / error = 2.65. This is the part that tests the transport coefficients.
- **Quantum vs classical transport:** Bose statistics lowers η by 13–37% and the predicted damping by 9–25% on these points. Classical transport gives the same RMS (2.65) with the opposite average sign (+1.7σ instead of −1.3σ). With five points and this scatter, the data do not discriminate between the two.
- **Numerical convergence:** N = 5 → 7 changes η by at most 0.04% and κ by at most 0.17%; refining the quadrature from (24,16,20) to (32,20,24) changes them by < 0.001%.
- **Input sensitivity:** the published errors of I, a and T are varied one at a time by ±1σ, recomputing density, fugacity, transport and sound. These are sensitivity ranges, not confidence intervals: I already contains a and the input covariance is not published. The error on L is not tabulated per row and is not included.
- **Data consistency:** row 13 (zero-based; L = 70 μm, a = 600 a₀) has a damping ratio γ/(ck) consistent with L = 72 μm instead. The original data are not changed; the row is outside the selected subset. The other fifteen rows reproduce the tabulated ratios.

Model errors are not removed by numerical convergence: finite Knudsen number, box walls, density inhomogeneity, energy-dependent cross section and approximate thermodynamics all remain. A disagreement is reported, never removed by changing a.

The standing-wave figure evolves the NSF equations from δn/n = 0.01 cos(πx/L), u = 0, δT = 0, with zero normal velocity at the walls. It illustrates the predicted mode; it does not reproduce the experiment's shaking protocol.

## 9. Reproduction

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --parallel 2
build/qtransport hilker --root .
build/qtransport figures --root .
```

Outputs: `results/hilker/comparison.json`, `comparison.svg`, `standing_wave.json`, `standing_wave.svg`. The input CSV files are unmodified copies from the attributed archive in `references/hilker2022/`.
