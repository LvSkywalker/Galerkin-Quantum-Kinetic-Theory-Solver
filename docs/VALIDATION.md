# Verification, convergence and model assessment

These answer three different questions: whether the implementation solves the intended finite problem, how its result changes under refinement, and how the physical model compares with an observation.

## Evidence behind the numbers

| Claim | Evidence | Interpretation |
|---|---|---|
| 26 test groups | `results/cpp/test_report.json`; `tests/tests.cpp` | Automated physics and regression checks at their recorded settings |
| 36 transport calculations | `results/cpp/transport.json` and `transport_same_order.csv` | 18 states, each computed with N=1 and N=4 |
| Maximum N=1 reference discrepancy <0.0031% | `results/cpp/benchmark.json` | Same one-profile approximation evaluated by two integration routes |
| Largest displayed N=4 to N=5 change: 0.3134% | `results/cpp/convergence.json` | Fermi conductivity at z=10, normalized to N=4 |
| Separate quadrature refinement | `kind=quadrature` records in `convergence.json` | Change at fixed N=4 under integration refinement |
| Linear relaxation and entropy decay | `results/cpp/dynamics.json` and tests | Finite-dimensional dynamics and internal response consistency |
| Optional five-point sound comparison | `results/cpp/hilker/comparison.json` | Assessment of the additional Hartree--Fock/fluid model |

The basis figure divides by N=5, so the same largest change appears there as approximately -0.3124%. All curves reach zero at N=5 by construction. Neither number certifies the remaining truncation error.

## One-profile verification

For classical hard spheres in the documented units, the code recovers

```text
mu_1 = 5/(16 sqrt(pi))
lambda_1 = 75/(64 sqrt(pi)).
```

These are benchmarks for the first-profile approximation, not exact coefficients of the full hard-sphere collision equation.

The quantum check uses classical, Bose z=0.5 and 0.9, and Fermi z=1 and 10. `brackets()` forms the five-variable Galerkin collision integral. Separately, `wu_reference()` evaluates the four-variable integrals I_A and I_B in Wu's Eq. (50), after matching normalization and taking ratios to the corresponding classical result. Both target the same N=1 trial space. They share equilibrium helper functions; this is a distinct integration route, not a wholly independent implementation of every input.

This check does not implement Wu's Fourier collision algorithm, compare N=4 against a continuum solution, or measure agreement with an experiment.

## Test coverage

The C++ suite checks Gaussian-rule moments, Laguerre recurrence, equilibrium normalization, quantum matching, detailed balance, Pauli occupation bounds, invalid inputs, elastic collision geometry, classical first- and second-profile benchmarks, power-law prefactors and kernel-strength scaling. It also checks quantum N=1 integrals, nested Galerkin estimates, generalized relaxation, entropy and response identities, and the long-wavelength shear limit.

The optional application has checks for Bose thermodynamics, unit conversion, sound dispersion and published data conventions. The standalone Fourier diagnostic is tested only as a representation of known functions. Frozen prototype outputs additionally test matrices, coefficients and transient curves; their provenance is retained in `tests/reference/README.md`.

The CTest configuration adds two CLI smoke checks. Local build evidence is recorded in [BUILD_VERIFICATION.md](BUILD_VERIFICATION.md). GitHub Actions outcomes must be read from an actual remote run.

## Refinement and limits

N increases at fixed quadrature, while quadrature is separately refined at fixed N. For exact nested Galerkin spaces the variational estimates are nondecreasing; inaccurate integration or ill conditioning can spoil that property. The saved sequence is useful evidence but does not establish a continuum error bound.

The production scan uses Bose z<=0.9 and Fermi z<=10. Accepting a larger range in the API is not a claim of verification near condensation or in arbitrarily strong degeneracy. Positive finite-dimensional decay rates do not prove a spectral gap for the full collision operator. Polynomial linear perturbations need not preserve distribution positivity at finite amplitude.

The projected shear model retains only selected angular sectors. The sound application adds thermodynamic and hydrodynamic assumptions; its experimental discrepancies are much larger than the one-profile benchmark discrepancy. These are distinct numerical and physical questions.
