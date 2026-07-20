# Deferred items (manual parity plan)

Items from `2026-07-19-viennaps-manual-parity-implementation-plan.md` tracked here.

## Completed (framework + smoke verification)

| Item | Implementation | Smoke marker |
|------|----------------|--------------|
| MC BCA implant (channeling + cascade) | `models/psMCBcaImplant.hpp` | `bca-check` |
| MFEM vector elasticity + traction + stress tensor | `fields/MfemElasticityKernel.hpp` | `mfem-elastic-check` |
| SUNDIALS field/MFEM dofs + amgcl linear solves | `SundialsTimeIntegrator` + `AmgclSolver.hpp` | `cvode-field-check`, `amgcl-check` |
| Domain material-map → mesh marking | `GeometryFieldCoupler::markFromMaterialMap/LayerStack` | `material-map-check` |
| LOCOS bird’s-beak + doping qualitative regression | `LocosDopingValidator.hpp` | `locos-check` |
| 3D mesh + band-limited solves | `PhysicsField::initMesh3D` + `BandLimitedSolver` | `band3d-check` |
| Phase 2 silicidation | `models/psSilicidation.hpp` | `silicide-check` |
| Phase 2 lithography (aerial/threshold mask) | `models/psLithography.hpp` | `litho-check` |
| Phase 2 advanced SPER | `fields/SPERKernel.hpp` | `sper-check` |
| Phase 2 parameter DB (inheritance + blend) | `fields/ParameterDatabase.hpp` | `paramdb-check` |
| Python multiphysics polish | `python/viennaps/multiphysics.pyi` stubs + C++ umbrella | (stubs) |
| Debug ↔ Release mfem CRT note | CMake `VIENNAPS_MFEM_REQUIRE_MATCHING_CRT` | `crt-check` |

## Remaining limits (not commercial calibration)

- ATHENA/SProcess **numerical** golden-profile matching (needs proprietary tables / lab data)
- Full crystal BCA with real silicon lattice potentials and full cascade recoils (engine is Kinchin–Pease + channeling model)
- Full 3D production performance (CUDA amgcl backend, adaptive band mesh from live LS)
- Complete pybind11 export of every multiphysics type (stubs document C++ API)

These limits do **not** block the deferred Phase 1–2 framework deliverables above.

## Earlier foundation (reference)

- Multiphysics foundation, cluster models, oxidation adapter OED, field-dof CVODE, pair/charged kernels, dose conservation fixes
