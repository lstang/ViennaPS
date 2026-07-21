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

- Task 5 (DiffusionEngine): complete (commits 71d6557..1aced0f, review clean — spec ✅ + quality Approved; ConstantDiffusion dangling-coef bug from Task 3 fixed; SUNDIALS/HYPRE paths gated by MFEM_USE_SUNDIALS/MPI — production solver path needs MFEM rebuild as follow-up but fallback is in-spec)

- Task 6 (umbrella header): complete (commits 1aced0f..31b94a4, review clean — spec ✅ + quality Approved; 2 pre-existing failures in intermediate/removeStrayPoints are unrelated psAnalyticImplant.hpp include-path bug)

- Task 7 (dose conservation): complete (commits 31b94a4..78b2ed3, review clean — spec ✅ + quality Approved; ratio 1.28e-16 far below 1% tolerance)
- All Phase 1 implementation tasks complete. Final whole-branch review next.

## Final whole-branch review: APPROVED (no Phase-1 blockers)

Commits 777f645..78b2ed3 (11 commits, 78KB). All Phase 1 spec ✅ + quality Approved across every task. `testDiffusion` PASS end-to-end. Dose conservation validated to 1.28e-16 (FP precision, far below 1% spec).

### Important findings (Phase 2 pre-work — do NOT block Phase 1 merge)
1. **Engine never consumes `DiffusionPhysics` BCs** — `DiffusionPhysics.hpp:70-98` defines BC API but `DiffusionEngine::solveImplicitEuler` keeps `essTdof` empty. Correct for Neumann-zero (all Phase 1 tests), but the first Phase 2 Dirichlet/segregation/Robin BC test will silently behave as Neumann-zero. Fix: call `physics.allBoundaryConditions()` and apply per-type, or add a TODO warning if non-empty.
2. **`shouldCreateTimeDerivative` not re-entrant across `solve()` calls** — `DiffusionPhysics.hpp:53-56` claims the species into `timeDerivativeClaimed_` and never clears it. On a second `solve()` on the same physics object, every species gets denied the mass matrix → silent fallback to identity mass. Phase 1 tests use fresh physics per test so this is latent. Fix: clear `timeDerivativeClaimed_` at top of `assembleAllSpecies()`, or make call-scoped.
3. **No multi-species integration test** — both engine tests use exactly 1 species. Species-outer loop, `allSpecies_` map plumbing, packed-block CVODE layout, `modelTargetsSpecies` filter are all multi-species code paths exercised only by inspection. Add a 2-species smoke test in Phase 2 Task 1 before relying on this for Fermi+Cdd composition.

### Minor findings (style/doc)
- `DiffusionEngine.hpp:391-392` — CVODE abstol `1e-9 * 1e18` reads as magic number; add a comment.
- `DiffusionPhysics.hpp:62-63` — comment references "the previous flat std::vector<BCSpec>" as if refactoring prior code; this file is new. Doc drift.
- `DiffusionModel.hpp:7` — `<map>` over-include on non-MFEM builds.
- `tests/diffusion/testDiffusion.cpp:57-82` — `FermiDiffusionStub`/`CddDiffusionStub` are correctly test-only; add a TODO citing Phase 3 Task 10 as the replacement site.

### Known follow-up (from Task 5)
- Rebuild MFEM at `f:/dev/mfem/build` with `MFEM_USE_MPI=ON MFEM_USE_SUNDIALS=ON MFEM_USE_HYPRE=ON` to exercise the production CVODE+HypreBoomerAMG path. The implicit-Euler+DSmoother fallback path is in-spec for Phase 1's linear ConstantDiffusion; production path is compile-verified but runtime-untested. HYPRE + SUNDIALS libs are already in vcpkg.

### Pre-existing failures (NOT from this branch)
- `intermediate` and `removeStrayPoints` tests fail on a `psAnalyticImplant.hpp:9` include-path bug (`psProcessModel.hpp` lives at `include/viennaps/process/psProcessModel.hpp`). Predates this branch.

## Status: Phase 1 implementation complete. Ready for Phase 2.

## Post-review fixes (Phase 1 polish)

Four issues raised after the final review, all resolved:

1. **Picard strategy documented on DiffusionModel.hpp** — added a
   `\section jacobian-strategy` doc block enumerating the three strategies
   (Picard default, manual Jacobian, MFEM NL kernel) per Phase 1 Task 5
   spec. Models opting into strategy (b) must declare so in their header.

2. **Dirichlet/Neumann BCs now applied by the engine** — previously the
   engine consumed `DiffusionPhysics::BCSpec` only as documentation and
   always ran with natural zero-flux. Now:
   - Dirichlet: `BilinearForm::EliminateEssentialBC(ess_attr_marker,
     presc_values, rhs, DIAG_ONE)` per step.
   - Neumann (non-zero): `BoundaryLFIntegrator(flux_coef)` added to RHS.
   - Boundary spec `"all"` (all attrs) or numeric string `"1"`, `"2"`, ...
   - segregation/robin: warning + treated as natural (Phase 2 Task 5
     territory).
   New tests: `TestDirichletBC` (interior relaxes to clamped 1e18),
   `TestNeumannBC` (dose increases at expected flux*perimeter*t rate,
   4% spatial-discretization error). All 34 tests pass.

3. **System-matrix caching** — DROPPED in favor of per-step rebuild.
   Rationale: the per-step `BilinearForm` rebuild enables clean Dirichlet
   elimination via MFEM's idiomatic API. The cache was the source of a
   subtle bug (`EliminateRow(DIAG_ONE)` requires the diagonal to exist
   in sparsity pattern, fragile). For Phase 1's small meshes the cost is
   negligible; for Phase 2 nonlinear D the K matrix changes per step
   anyway so caching wouldn't help. Trade-off documented inline.

4. **"Conforming cut-cell deferred" tracked as follow-up F6** — see
   `docs/superpowers/specs/diffusion-phase1-followups.md`. Lists the
   three deferred pieces (conforming cut mesh, transition layer,
   MFEMCutTransitionSubMesh integration) with triggers for revisit.
