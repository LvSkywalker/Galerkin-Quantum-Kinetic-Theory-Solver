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
