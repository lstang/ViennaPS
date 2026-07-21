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

- Task 1 (MeshAttributes): complete (commits 8c907ef..909ff2f, review clean — spec ✅ + quality Approved)

- Task 2 (DiffusionModel): complete (commits 909ff2f..48dfaba, review clean — spec ✅ + quality Approved; fixed stale filename in plan)

- Task 3 (ConstantDiffusion): complete (commits 9003343..1cd035e, review clean — spec ✅ + quality Approved; attrs_ shadow accepted as intentional)

- Task 3.5 (DiffusionPhysics): complete (commits 1cd035e..0ffaf30, review clean — spec ✅ + quality Approved; Fermi/Cdd stubbed in test per controller decision)

- Task 4 (LevelSetToMesh): complete (commits 0ffaf30..71d6557, review found Critical unconditional <mfem.hpp> include → fixed in 71d6557 → re-reviewed clean — spec ✅ + quality Approved; getMinBounds fix for INFINITE_BOUNDARY verified against hrleGrid.hpp:122-133)

- Task 5 (DiffusionEngine): complete — see task-5-report.md. Two pre-existing bugs in ConstantDiffusion fixed along the way (dangling Coefficient& references); engine exercised through implicit-Euler path because the linked MFEM 4.9.1 was built without MFEM_USE_SUNDIALS and without MFEM_USE_MPI (so neither mfem::CVODESolver nor mfem::HypreBoomerAMG is available). SUNDIALS+HypreBoomerAMG code paths compile-clean under the appropriate #ifdef gates and will activate when MFEM is rebuilt with those options.
