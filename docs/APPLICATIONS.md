# Applications of the collision matrices

The transport solver and homogeneous relaxation form the core project. The two applications below add spatial dynamics under explicit extra assumptions.

## Transverse kinetic mode

For variation along y, retain transverse momentum, thermal-vector and shear-tensor velocity modes. Galerkin projection of streaming and collisions gives

```text
H dc_k/dt = -(M + i k S) c_k,
S_mn = integral w Phi_m v_y Phi_n dv.
```

At small wave number the slow eigenvalue recovers shear diffusion, and its inferred viscosity approaches the steady transport coefficient. This is an internal consistency check: the two calculations share collision brackets. At finite k, streaming would excite angular sectors outside the retained basis; the present dispersion is not a converged full kinetic prediction.

```sh
build/qtransport transport --statistics bose --fugacity 0.5 --basis 4 --output runs/bose.json
build/qtransport wave --input runs/bose.json --k 0.1 --output runs/bose-wave.json
```

The full reproduction command already generates representative modes and scans.

## Normal-Bose sound and Hilker et al.

The optional application uses the published dataset for Hilker et al., *First and Second Sound in a Compressible 3D Bose Fluid*, [PRL 128, 223601 (2022)](https://doi.org/10.1103/PhysRevLett.128.223601). Data are available from the [Cambridge repository](https://doi.org/10.17863/CAM.83837) under CC BY 4.0; full attribution is retained in `references/hilker2022/ATTRIBUTION.json`.

The calculation maps the low-energy bosonic cross section to `B=2 a^2 |v-v1|`, computes transport with the Galerkin solver, adds normal-phase Hartree--Fock thermodynamics, and solves longitudinal Navier--Stokes--Fourier modes. The [technical note](technical_note/main.pdf) derives this chain and explains its normalization.

Five observations are selected by the published inverse Knudsen number L/ell>=3. Density is reconstructed from that quantity, the tabulated scattering length and box length. No model parameter is fitted. The selection, conversions and unselected row 13's inconsistent damping-ratio entry are recorded explicitly in the output manifest.

```sh
build/qtransport hilker --root .
build/qtransport figures --root .
```

The saved study has two separate comparisons:

- The published theory curve evaluated with its rounded parameters is recovered within about 0.08% in sound speed over the selected hydrodynamic range. This checks the fluid/thermodynamic calculation.
- The project's five predictions have experimental damping discrepancies of approximately -21% to +16%. This tests the adequacy of the additional physical assumptions, not just numerical integration.

Input variations are one-at-a-time sensitivity studies, not confidence intervals; input covariances are unavailable. The illustrative standing wave starts with a specified density perturbation and does not reproduce the experimental driving protocol. These limitations remain even if the collision integrals and polynomial basis are numerically converged.

## Classical-transport control

To isolate the effect of transport statistics, `classical_transport_same_hf` in each selected `comparison.json` record replaces the quantum viscosity and conductivity by their classical counterparts **at the same basis size, physical state and cross-section convention**, while retaining the same Hartree--Fock equation of state and thermodynamic derivatives. It is not a fully classical gas model. The independent `no_mean_field` calculation answers a different question and is not this control.

For observable y, we report R = sqrt(mean(((y_model - y_data)/sigma_data)^2)), using each point's published measurement uncertainty:

| Statistic, five selected observations | Quantum transport + HF | Classical transport + HF |
|---|---:|---:|
| Sound-speed normalized RMS residual | 1.0308 | 1.0404 |
| Damping normalized RMS residual | 2.6544 | 2.6475 |
| Largest absolute damping discrepancy | 20.85% | 34.18% |

The RMS agreement is almost unchanged, although individual predictions and the largest relative discrepancy differ. Thus this selected comparison does not demonstrate that quantum transport improves the overall agreement. It also does not establish that quantum corrections are absent: the transport coefficients differ, but the sound observables also depend on the shared thermodynamics and the fluid closure.

These are descriptive standardized residuals, not a formal significance test or a reduced chi-square. Uncertainty and covariance in reconstructed inputs are not propagated into the denominator; input variations remain separate sensitivity calculations. Some individual damping residuals exceed four reported measurement uncertainties. The five points should not be described as all lying within 1.03 uncertainties.

The plot uses inverse Knudsen number L/ell on the horizontal axis. Circles and triangles distinguish 50 and 70 micrometre boxes; colors distinguish measured, quantum-transport and classical-transport values. No lines connect observations with different physical inputs. Vertical data bars are reported measurement uncertainties; quantum-model bars are one-input-at-a-time sensitivity ranges. Classical-control markers show central values only.

`qtransport_report_export` computes the metrics directly from the saved C++ predictions and writes [hilker_metrics.json](technical_note/data/hilker_metrics.json), including per-observation residuals, and a LaTeX table. The figures and metrics can therefore be rebuilt without a new collision calculation.
