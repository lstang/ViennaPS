# Diffusion Engine — Implementation Progress Ledger

Tracking Phase 1 implementation via subagent-driven-development.
Source plan: `docs/superpowers/plans/2026-07-20-diffusion-phase1.md`.
Architecture decision: see `docs/superpowers/specs/adr-0001-mfem-architecture.md`.

Build dir for this work: `build_phase1/`.
Configure command (rerun if cache dropped):
```
cmake -B build_phase1 -G "Visual Studio 17 2022" -A x64 \
  -DVIENNAPS_BUILD_TESTS=ON \
  -DCMAKE_TOOLCHAIN_FILE=F:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake \
  -DVCPKG_TARGET_TRIPLET=x64-windows \
  -DVCPKG_INSTALLED_DIR=F:/dev/vcpkg/installed \
  -DVCPKG_MANIFEST_INSTALL=OFF
```
Build a single test target:
```
cmake --build build_phase1 --config Release --target testDiffusion
```
Run a single test:
```
ctest -R testDiffusion --test-dir build_phase1 -C Release --output-on-failure
```

Verified at skill start: MFEM 4.9.1 at `f:/dev/mfem/build` detected, SUNDIALS detected, full solution builds clean, sundials_cvode.dll present.

## Task ledger

- Task 0 (ADR): complete (commit pending)
- Task 1 (MeshAttributes): pending
- Task 2 (DiffusionModel): pending
- Task 3 (ConstantDiffusion): pending
- Task 3.5 (DiffusionPhysics): pending
- Task 4 (LevelSetToMesh): pending
- Task 5 (DiffusionEngine): pending
- Task 6 (umbrella header): pending
- Task 7 (dose conservation test): pending
- Final review: pending
