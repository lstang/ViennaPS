# Phase 1 Follow-ups Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close out the open Phase 1 follow-ups (F2, F3, F5) by (a) fixing the MFEM build configuration so ViennaPS picks up the newly-rebuilt SUNDIALS/MPI/HYPRE-enabled MFEM, (b) verifying the previously-compile-only CVODE + HypreBoomerAMG code paths now actually run, (c) making `DiffusionPhysics::shouldCreateTimeDerivative` re-entrant across `solve()` calls, and (d) adding a multi-species integration test that exercises the species-outer assembly loop, `allSpecies_` map, and packed-block CVODE layout - all currently tested only by inspection.

**Architecture:** Modifies the existing `DiffusionEngine` and `DiffusionPhysics` headers in `include/viennaps/fields/`. No new headers. Adds two new tests to `tests/diffusion/testDiffusion.cpp`. The MFEM path fix is two file edits in `f:/dev/mfem/build/` (the prebuilt MFEM install tree).

**Tech Stack:** C++20, MFEM 4.9.1 (now built with `MFEM_USE_MPI`, `MFEM_USE_SUNDIALS`, HYPRE 2.23.0), SUNDIALS CVODE, ViennaPS existing Phase 1 code.

## Global Constraints

- C++20, header-only under `include/viennaps/` (per ADR-0001, MFEM `.lib` link is the only relaxation)
- MFEM code gated by `#ifdef VIENNAPS_HAS_MFEM`, SUNDIALS by `#ifdef VIENNAPS_HAS_SUNDIALS`
- MFEM internal feature flags used directly when checking MFEM-built-with features: `#ifdef MFEM_USE_SUNDIALS`, `#ifdef MFEM_USE_MPI` (these come from MFEM's `config/_config.hpp`, not ViennaPS's)
- MFEM Release at `f:/dev/mfem/build`; Debug (MDd) against Release MFEM is unsupported (per AGENTS.md)
- vcpkg deps at `f:/dev/vcpkg/installed/x64-windows/`; SUNDIALS user-build at `f:/dev/sundials_build/install/`
- LLVM style: 2-space indent, 80-col, no tabs, Attach braces, pointer right
- Tests use `VC_TEST_ASSERT` from `<vcTestAsserts.hpp>`
- Namespace `viennaps`
- Test target is `testDiffusion` (file `tests/diffusion/testDiffusion.cpp`); CMakeLists auto-discovers via `viennacore_add_subdirs`
- Build dir for this work: `build_followups/` (created in Task 1)
- **Out of scope:** F1 (BCs applied) ✅ done in commit `9ada745`; F4 (caching) ✅ done in commit `9c37abb`; F6 (conforming cut-cell) deferred to Phase 2+ Task 4 scope; F7 (plan filename) ✅ done in commit `9003343`.

---

## File Structure

| File | Responsibility | Action |
|------|---------------|--------|
| `f:/dev/mfem/build/MFEMConfig.cmake` | MFEM CMake package config | Modify: replace stale `build_full` paths with `build` |
| `f:/dev/mfem/build/MFEMTargets.cmake` | MFEM imported target | Modify: replace stale `build_full` paths with `build` |
| `include/viennaps/fields/DiffusionPhysics.hpp` | Per-species BC storage + composition gatekeepers | Modify: make `shouldCreateTimeDerivative` re-entrant per `solve()` call |
| `include/viennaps/fields/DiffusionEngine.hpp` | Engine assembly + time integration | Modify: reset `timeDerivativeClaimed_` at start of `assembleAllSpecies()` |
| `tests/diffusion/testDiffusion.cpp` | Phase 1 test suite | Modify: add `TestCVODERuntimePath`, `TestMultiSpeciesSmoke`, `TestReentrantSolve` |

---

## Task 1: Fix MFEM CMake paths and re-configure ViennaPS

**Files:**
- Modify: `f:/dev/mfem/build/MFEMConfig.cmake`
- Modify: `f:/dev/mfem/build/MFEMTargets.cmake`

**Interfaces:** None. This task unblocks downstream tasks by making `find_package(MFEM)` succeed with the rebuilt MFEM.

**Background:** The user rebuilt MFEM with `MFEM_USE_MPI=ON MFEM_USE_SUNDIALS=ON MFEM_USE_HYPRE=ON` but the build was done in a directory named `build_full` that has since been renamed/moved to `build`. The `MFEMConfig.cmake` and `MFEMTargets.cmake` files still reference `F:/dev/mfem/build_full` for include dirs and `IMPORTED_LOCATION_RELEASE`, which causes `find_package(MFEM)` to fail with:

```
CMake Error at F:/dev/mfem/build/MFEMConfig.cmake:87 (message):
  File or directory F:/dev/mfem/build_full referenced by variable
  MFEM_INCLUDE_DIR does not exist !
```

Verified facts (run during plan writing):
- `f:/dev/mfem/build/MFEMConfig.cmake` line 87 references `F:/dev/mfem/build_full` (line: `set(MFEM_INCLUDE_DIRS "F:/dev/mfem/build_full;F:/dev/vcpkg/installed/x64-windows/include;F:/dev/sundials_build/install/include")`)
- `f:/dev/mfem/build/MFEMTargets.cmake` line ~57 has `INTERFACE_INCLUDE_DIRECTORIES "F:/dev/mfem/build_full;..."` and line ~65 has `IMPORTED_LOCATION_RELEASE "F:/dev/mfem/build_full/mfem.lib"`
- `f:/dev/mfem/build/mfem.lib` exists (verified via `ls`)
- `f:/dev/mfem/build/mfem.hpp` exists (verified via `ls`)
- `f:/dev/mfem/build/config/_config.hpp` defines `MFEM_USE_MPI`, `MFEM_USE_SUNDIALS`, and `MFEM_HYPRE_VERSION 23200`
- `f:/dev/sundials_build/install/include` and `f:/dev/sundials_build/install/lib` exist (verified)
- `f:/dev/vcpkg/installed/x64-windows/lib/HYPRE.lib` and `metis.lib` exist (referenced by `MFEMTargets.cmake`)

- [ ] **Step 1: Back up the two cmake files**

Run from Git Bash:
```bash
cp f:/dev/mfem/build/MFEMConfig.cmake f:/dev/mfem/build/MFEMConfig.cmake.bak
cp f:/dev/mfem/build/MFEMTargets.cmake f:/dev/mfem/build/MFEMTargets.cmake.bak
```

- [ ] **Step 2: Replace `build_full` with `build` in `MFEMConfig.cmake`**

Use sed (Git Bash):
```bash
sed -i 's|F:/dev/mfem/build_full|F:/dev/mfem/build|g' f:/dev/mfem/build/MFEMConfig.cmake
```

Verify the change:
```bash
grep "build_full" f:/dev/mfem/build/MFEMConfig.cmake
```
Expected: no output (no remaining `build_full` references).
```bash
grep "MFEM_INCLUDE_DIRS\|MFEM_LIBRARY_DIR" f:/dev/mfem/build/MFEMConfig.cmake
```
Expected: both lines now show `F:/dev/mfem/build` (not `build_full`).

- [ ] **Step 3: Replace `build_full` with `build` in `MFEMTargets.cmake`**

```bash
sed -i 's|F:/dev/mfem/build_full|F:/dev/mfem/build|g' f:/dev/mfem/build/MFEMTargets.cmake
```

Verify:
```bash
grep "build_full" f:/dev/mfem/build/MFEMTargets.cmake
```
Expected: no output.
```bash
grep "IMPORTED_LOCATION_RELEASE\|INTERFACE_INCLUDE_DIRECTORIES" f:/dev/mfem/build/MFEMTargets.cmake
```
Expected: both lines now reference `F:/dev/mfem/build` and `F:/dev/mfem/build/mfem.lib`.

- [ ] **Step 4: Configure ViennaPS in a fresh build directory**

Run from `F:/dev/ViennaPS_mod`:
```bash
rm -rf build_followups
cmake -B build_followups -G "Visual Studio 17 2022" -A x64 \
  -DVIENNAPS_BUILD_TESTS=ON \
  -DCMAKE_TOOLCHAIN_FILE=F:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake \
  -DVCPKG_TARGET_TRIPLET=x64-windows \
  -DVCPKG_INSTALLED_DIR=F:/dev/vcpkg/installed \
  -DVCPKG_MANIFEST_INSTALL=OFF \
  > configure_followups.log 2>&1
```

- [ ] **Step 5: Verify MFEM + SUNDIALS + HYPRE all detected**

```bash
grep -i "found mfem\|found sundials\|mfem linked\|HYPRE" configure_followups.log
```
Expected output (lines verified to appear when MFEM is correctly found):
```
-- [ViennaPS] Found MFEM 4.9.1 for multiphysics fields at f:/dev/mfem/build
--            MFEM_INCLUDE_DIRS=F:/dev/mfem/build;F:/dev/vcpkg/installed/x64-windows/include;F:/dev/sundials_build/install/include
-- [ViennaPS] MFEM built with CUDA support - locating CUDAToolkit (optional for now)
-- [ViennaPS] MFEM linked from f:/dev/mfem/build. Use matching CRT: Release app ↔ Release mfem (MD). Debug (MDd) against Release mfem is unsupported - build Debug mfem or use Release.
-- [ViennaPS] Found SUNDIALS for time integration
```

If MFEM is found, the log will also show `Found MFEM 4.9.1`. If it fails, re-read `configure_followups.log` for the FATAL_ERROR message.

- [ ] **Step 6: Build the existing test suite (no code changes yet)**

```bash
cmake --build build_followups --config Release --target testDiffusion
```
Expected: builds successfully. The pre-existing Phase 1 tests should still compile and link - the only difference is MFEM now pulls MPI + SUNDIALS + HYPRE transitively via `MFEMTargets.cmake`.

- [ ] **Step 7: Run the existing test suite to confirm baseline**

```bash
ctest -R testDiffusion --test-dir build_followups -C Release --output-on-failure
```
Expected: 1/1 test passes. The existing 8 tests (`TestMeshAttributes`, `TestDiffusionModelInterface`, `TestConstantDiffusion`, `TestDiffusionPhysics`, `TestDiffusionPhysicsComposition`, `TestLevelSetToMesh2D`, `TestDiffusionEngineAssembly`, `TestDoseConservation`, `TestDirichletBC`, `TestNeumannBC`) should all still pass with the new MFEM build.

- [ ] **Step 8: No commit (the changes are in the MFEM install tree, not the ViennaPS repo)**

Record the path-fix in the progress ledger at `.superpowers/sdd/progress.md` instead:
```
- Task 1 (MFEM path fix): complete (no ViennaPS commit; edits to f:/dev/mfem/build/MFEM{Config,Targets}.cmake replaced stale build_full paths)
```

---

## Task 2: Verify CVODE + HypreBoomerAMG runtime path (F5)

**Files:**
- Modify: `tests/diffusion/testDiffusion.cpp` (append new test)

**Interfaces:**
- Consumes: `DiffusionEngine::solve()` (existing), `DiffusionEngine::getIntegral()` (existing)
- Produces: `TestCVODERuntimePath()` - proves the SUNDIALS CVODE + HypreBoomerAMG path runs end-to-end

**Background:** The `DiffusionEngine` has two time-integration paths gated by MFEM-internal flags:
- `#ifdef MFEM_USE_SUNDIALS` -> `solveCVODE()` using `mfem::CVODESolver(CV_BDF)` + `DiffusionRHSOperator` (a `TimeDependentOperator`)
- `#else` -> `solveImplicitEuler()` (cached `(M+dtK)` + CG)

The mass-matrix solver factory `makeMassSolver()` is similarly gated:
- `#ifdef MFEM_USE_MPI` -> `mfem::HypreBoomerAMG` preconditioner + `mfem::CGSolver`
- `#else` -> `mfem::DSmoother` preconditioner + `mfem::CGSolver`

Before Task 1, MFEM was built without SUNDIALS/MPI/HYPRE, so the runtime always hit the implicit-Euler + DSmoother path. The CVODE + HypreBoomerAMG code was compile-verified only (gated correctly, but never executed). After Task 1, both paths should activate. This task adds a test that proves they actually run.

The test is a re-run of `TestDoseConservation` (closed system, zero-flux Neumann, uniform IC) but with a diagnostic print that confirms which path executed. Dose conservation must still hold (the physics is the same; only the integrator changed).

- [ ] **Step 1: Append the new test to `tests/diffusion/testDiffusion.cpp`**

Add this function immediately before `int main()`:
```cpp
void TestCVODERuntimePath() {
  // F5 follow-up: prove the SUNDIALS CVODE + HypreBoomerAMG code path
  // actually runs end-to-end (not just compile-verified). Before the MFEM
  // rebuild this path was dead code; the test prints which integrator and
  // preconditioner were used so a human can confirm in CI logs.
  //
  // Physics: same as TestDoseConservation (closed system, uniform IC) so
  // dose conservation must still hold. The point is to exercise the
  // CVODE + HypreBoomerAMG machinery, not to test new physics.
  MeshAttributes attrs;
  attrs.setAttributeName(1, "Si");

  auto mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(8, 8, mfem::Element::TRIANGLE));

  DiffusionEngine<double, 2> engine;
  engine.setMesh(std::move(mesh), attrs);

  auto model = std::make_shared<ConstantDiffusion<double>>("Boron");
  model->setDiffusivity(1e-8, 0.0);
  DiffusionPhysics<double> physics;
  physics.addSpecies("Boron");
  physics.addModel(model);
  physics.setTemperature(1273.15);
  engine.setPhysics(physics);

  engine.initializeSpecies("Boron", 1e18);
  const double doseInitial = engine.getIntegral("Boron");

#ifdef MFEM_USE_SUNDIALS
  std::cout << "[cvode-runtime] MFEM_USE_SUNDIALS defined -> CVODE BDF path "
            << "active\n";
#else
  std::cout << "[cvode-runtime] MFEM_USE_SUNDIALS NOT defined -> implicit-"
            << "Euler fallback path active\n";
#endif
#ifdef MFEM_USE_MPI
  std::cout << "[cvode-runtime] MFEM_USE_MPI defined -> HypreBoomerAMG "
            << "preconditioner active\n";
#else
  std::cout << "[cvode-runtime] MFEM_USE_MPI NOT defined -> DSmoother "
            << "preconditioner active\n";
#endif

  engine.solve(0.0, 10.0, 1.0);
  const double doseFinal = engine.getIntegral("Boron");
  const double relDiff =
      std::abs(doseFinal - doseInitial) / std::abs(doseInitial);

  std::cout << "[cvode-runtime] dose_initial=" << doseInitial
            << " dose_final=" << doseFinal << " rel_diff=" << relDiff << "\n";

  // Same 1% dose-conservation tolerance as TestDoseConservation. CVODE's
  // adaptive BDF should hold this easily for a closed linear system.
  VC_TEST_ASSERT(relDiff < 0.01);
}
```

- [ ] **Step 2: Register the test in `main()`**

Find the `int main()` function and add `TestCVODERuntimePath();` immediately after `TestNeumannBC();`. The final `main()` should look like:
```cpp
int main() {
  TestMeshAttributes();
  TestDiffusionModelInterface();
  TestConstantDiffusion();
  TestDiffusionPhysics();
  TestDiffusionPhysicsComposition();
  TestLevelSetToMesh2D();
  TestDiffusionEngineAssembly();
  TestDoseConservation();
  TestDirichletBC();
  TestNeumannBC();
  TestCVODERuntimePath();
  std::cout << "All diffusion tests passed.\n";
  return 0;
}
```

- [ ] **Step 3: Build the test**

```bash
cmake --build build_followups --config Release --target testDiffusion
```
Expected: builds successfully. The CVODE path should now link against `sundials_cvodes.lib` (transitively via `MFEMTargets.cmake`).

- [ ] **Step 4: Run the test**

```bash
ctest -R testDiffusion --test-dir build_followups -C Release --output-on-failure
```
Expected: 1/1 PASS. The output should include (in the `LastTest.log` at `build_followups/Testing/Temporary/LastTest.log`):
```
[cvode-runtime] MFEM_USE_SUNDIALS defined -> CVODE BDF path active
[cvode-runtime] MFEM_USE_MPI defined -> HypreBoomerAMG preconditioner active
[cvode-runtime] dose_initial=1e+18 dose_final=1e+18 rel_diff=<small>
```

If `rel_diff` is > 0.01, the CVODE path has a bug - investigate before proceeding. If the path prints `implicit-Euler fallback` or `DSmoother`, the MFEM rebuild did not take effect - re-check Task 1.

- [ ] **Step 5: Commit**

```bash
git add tests/diffusion/testDiffusion.cpp
git commit -m "test(diffusion): verify CVODE + HypreBoomerAMG runtime path (F5)"
```

---

## Task 3: Make `shouldCreateTimeDerivative` re-entrant across `solve()` calls (F2)

**Files:**
- Modify: `include/viennaps/fields/DiffusionPhysics.hpp` (lines ~46-56, ~106)
- Modify: `include/viennaps/fields/DiffusionEngine.hpp` (`assembleAllSpecies()` body)

**Interfaces:**
- Consumes: existing `DiffusionPhysics::shouldCreateTimeDerivative(species, model)` API (no signature change)
- Produces: same API now safe to call across multiple `solve()` invocations on the same physics object

**Background:** `DiffusionPhysics::shouldCreateTimeDerivative` (line 50) inserts the species name into `timeDerivativeClaimed_` on first call and returns `false` on subsequent calls for the same species. This is correct within a single `solve()` call (prevents two models from both adding a `dC/dt` term). But `timeDerivativeClaimed_` is never cleared, so on a *second* `solve()` call on the same physics object, every species gets denied the mass matrix -> the engine's fallback (identity mass, `DiffusionEngine.hpp:307-311`) silently kicks in, producing wrong time scales.

Phase 1 tests don't hit this because each test constructs a fresh `DiffusionPhysics` object. Phase 2 will hit it on any multi-step anneal that reuses the physics object across `solve()` calls.

**Fix:** clear `timeDerivativeClaimed_` at the start of `DiffusionEngine::assembleAllSpecies()`. The gatekeeper stays per-solve-call (matches the MOOSE `PhysicsBase::shouldCreateTimeDerivative` semantics it mirrors - per-add-kernel, which in our case is per-solve).

- [ ] **Step 1: Write the failing test**

Append `TestReentrantSolve()` to `tests/diffusion/testDiffusion.cpp` immediately before `TestCVODERuntimePath()`:
```cpp
void TestReentrantSolve() {
  // F2 follow-up: shouldCreateTimeDerivative must be re-entrant across
  // solve() calls. The gatekeeper (DiffusionPhysics.hpp:50) inserts the
  // species into timeDerivativeClaimed_ on first call and returns false
  // thereafter. Without a reset, the second solve() denies the mass
  // matrix for every species, silently falling back to identity mass
  // (DiffusionEngine.hpp:307-311) and producing wrong time scales.
  //
  // Setup: closed system, uniform IC. Two consecutive solve() calls on
  // the SAME physics object. Both must conserve dose to FP precision.
  // If the second call's rel_diff blows up, the mass matrix was denied.
  MeshAttributes attrs;
  attrs.setAttributeName(1, "Si");

  auto mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(8, 8, mfem::Element::TRIANGLE));

  DiffusionEngine<double, 2> engine;
  engine.setMesh(std::move(mesh), attrs);

  auto model = std::make_shared<ConstantDiffusion<double>>("Boron");
  model->setDiffusivity(1e-8, 0.0);
  DiffusionPhysics<double> physics;
  physics.addSpecies("Boron");
  physics.addModel(model);
  physics.setTemperature(1273.15);
  engine.setPhysics(physics);

  engine.initializeSpecies("Boron", 1e18);
  const double doseInitial = engine.getIntegral("Boron");

  // First solve - normal behaviour, mass matrix claimed for Boron.
  engine.solve(0.0, 5.0, 1.0);
  const double doseAfterFirst = engine.getIntegral("Boron");
  const double relDiffFirst =
      std::abs(doseAfterFirst - doseInitial) / std::abs(doseInitial);

  // Second solve on the SAME physics object - this is the regression
  // case. Before the fix, timeDerivativeClaimed_ still holds "Boron",
  // so the mass matrix is denied and identity mass is used. Identity
  // mass means du/dt = -K u + R (no M inverse) - wildly different time
  // scale, dose will drift.
  engine.solve(5.0, 10.0, 1.0);
  const double doseAfterSecond = engine.getIntegral("Boron");
  const double relDiffSecond =
      std::abs(doseAfterSecond - doseInitial) / std::abs(doseInitial);

  std::cout << "[reentrant-solve] dose_initial=" << doseInitial
            << " after_first=" << doseAfterFirst
            << " rel_diff_first=" << relDiffFirst
            << " after_second=" << doseAfterSecond
            << " rel_diff_second=" << relDiffSecond << "\n";

  // Both solves must conserve dose. The 1e-6 tolerance is tighter than
  // TestDoseConservation's 1% because we want to catch identity-mass
  // fallback (which produces order-unity drift, not 1e-16 round-off).
  VC_TEST_ASSERT(relDiffFirst < 1e-6);
  VC_TEST_ASSERT(relDiffSecond < 1e-6);
}
```

Register it in `main()` by adding `TestReentrantSolve();` immediately before `TestCVODERuntimePath();`.

- [ ] **Step 2: Build and run to confirm the test FAILS**

```bash
cmake --build build_followups --config Release --target testDiffusion
ctest -R testDiffusion --test-dir build_followups -C Release --output-on-failure
```
Expected: FAIL. Read `build_followups/Testing/Temporary/LastTest.log` and find the `[reentrant-solve]` line. `rel_diff_first` should be ~1e-16 (first solve works). `rel_diff_second` should be significantly larger (e.g. > 1e-6) - this is the bug. If `rel_diff_second` is also ~1e-16, the bug isn't reproducing and the test is not exercising the regression - investigate before proceeding.

- [ ] **Step 3: Add `resetTimeDerivativeClaims()` to `DiffusionPhysics`**

In `include/viennaps/fields/DiffusionPhysics.hpp`, find the existing `shouldCreateTimeDerivative` method (around line 50). Add a public `resetTimeDerivativeClaims()` method immediately after it:
```cpp
  /// Clear the per-species time-derivative claim set. Call this at the
  /// start of each solve() so composing models can re-claim dC/dt on the
  /// same species across multiple solves on the same physics object.
  /// Mirrors MOOSE PhysicsBase semantics where the gatekeeper is
  /// per-add-kernel (in our case, per-solve), not per-physics-lifetime.
  void resetTimeDerivativeClaims() { timeDerivativeClaimed_.clear(); }
```

- [ ] **Step 4: Call `resetTimeDerivativeClaims()` at the start of `assembleAllSpecies()`**

In `include/viennaps/fields/DiffusionEngine.hpp`, find `assembleAllSpecies()`. At the very top of the method body (before the `if (!physics_)` check, OR immediately after it - pick the location that makes the reset unconditional whenever assembly runs), add:
```cpp
    physics_->resetTimeDerivativeClaims();
```

The exact insertion point: after the existing `if (!physics_) throw ...` check, before the `if (attrs_)` block that propagates temperature. This ensures the claim set is fresh for every assembly pass.

- [ ] **Step 5: Build and run to confirm the test PASSES**

```bash
cmake --build build_followups --config Release --target testDiffusion
ctest -R testDiffusion --test-dir build_followups -C Release --output-on-failure
```
Expected: 1/1 PASS. The `[reentrant-solve]` log line should now show both `rel_diff_first` and `rel_diff_second` at ~1e-16 (floating-point precision).

- [ ] **Step 6: Commit**

```bash
git add include/viennaps/fields/DiffusionPhysics.hpp include/viennaps/fields/DiffusionEngine.hpp tests/diffusion/testDiffusion.cpp
git commit -m "fix(diffusion): make shouldCreateTimeDerivative re-entrant across solve() calls (F2)"
```

---

## Task 4: Add multi-species integration test (F3)

**Files:**
- Modify: `tests/diffusion/testDiffusion.cpp` (append new test)

**Interfaces:**
- Consumes: existing `DiffusionPhysics::addSpecies`, `DiffusionEngine::initializeSpecies`, `DiffusionEngine::getIntegral`
- Produces: `TestMultiSpeciesSmoke()` - exercises species-outer assembly loop, `allSpecies_` map, packed-block CVODE layout, `modelTargetsSpecies` filter

**Background:** Both `TestDiffusionEngineAssembly` and `TestDoseConservation` use exactly 1 species. The species-outer loop (`DiffusionEngine.hpp:269`), the `allSpecies_` map plumbing, the packed-block CVODE state layout (`DiffusionEngine.hpp:397-402`), and the `modelTargetsSpecies` filter are all multi-species code paths exercised only by inspection. Phase 2's Fermi+Cdd composition will be the first real multi-species exercise; this task adds a cheap 2-species smoke test as a safety net.

**Physics:** two independent `ConstantDiffusion` models on two species ("Boron" and "Phosphorus"), each with closed-system zero-flux BCs. Since the species don't couple, each must independently conserve its own dose. This is the simplest possible multi-species test - it catches bugs in the species-outer loop, the `allSpecies_` map, and the packed-block layout without coupling physics.

- [ ] **Step 1: Append the test to `tests/diffusion/testDiffusion.cpp`**

Add immediately before `TestReentrantSolve()`:
```cpp
void TestMultiSpeciesSmoke() {
  // F3 follow-up: exercise the multi-species code paths in DiffusionEngine.
  // Both existing engine tests use exactly 1 species; the species-outer
  // assembly loop, the allSpecies_ map, the packed-block CVODE layout,
  // and the modelTargetsSpecies filter are tested only by inspection.
  //
  // Setup: two INDEPENDENT ConstantDiffusion models on two species
  // (Boron, Phosphorus), each with closed-system zero-flux BCs. Since
  // the species don't couple, each must independently conserve its own
  // dose. This catches bugs in the species-outer loop, the allSpecies_
  // map, and the packed-block layout without coupling physics.
  MeshAttributes attrs;
  attrs.setAttributeName(1, "Si");

  auto mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(8, 8, mfem::Element::TRIANGLE));

  DiffusionEngine<double, 2> engine;
  engine.setMesh(std::move(mesh), attrs);

  auto boronModel = std::make_shared<ConstantDiffusion<double>>("Boron");
  boronModel->setDiffusivity(1e-8, 0.0);
  auto phosphorusModel =
      std::make_shared<ConstantDiffusion<double>>("Phosphorus");
  phosphorusModel->setDiffusivity(1e-7, 0.0); // different D - catches
                                              // cross-species K leakage

  DiffusionPhysics<double> physics;
  physics.addSpecies("Boron");
  physics.addSpecies("Phosphorus");
  physics.addModel(boronModel);
  physics.addModel(phosphorusModel);
  physics.setTemperature(1273.15);
  engine.setPhysics(physics);

  // Different ICs so a cross-species leak would be visible: Boron at
  // 1e18, Phosphorus at 1e15. If Boron's K is applied to Phosphorus,
  // Phosphorus's dose will drift visibly (D_Boron=1e-8 vs D_P=1e-7).
  engine.initializeSpecies("Boron", 1e18);
  engine.initializeSpecies("Phosphorus", 1e15);

  const double boronInitial = engine.getIntegral("Boron");
  const double phosphorusInitial = engine.getIntegral("Phosphorus");

  engine.solve(0.0, 10.0, 1.0);

  const double boronFinal = engine.getIntegral("Boron");
  const double phosphorusFinal = engine.getIntegral("Phosphorus");

  const double boronRelDiff =
      std::abs(boronFinal - boronInitial) / std::abs(boronInitial);
  const double phosphorusRelDiff =
      std::abs(phosphorusFinal - phosphorusInitial) /
      std::abs(phosphorusInitial);

  std::cout << "[multi-species] boron: initial=" << boronInitial
            << " final=" << boronFinal << " rel_diff=" << boronRelDiff
            << "\n";
  std::cout << "[multi-species] phosphorus: initial=" << phosphorusInitial
            << " final=" << phosphorusFinal
            << " rel_diff=" << phosphorusRelDiff << "\n";

  // Each species must conserve its own dose independently. 1e-6 tolerance
  // catches any cross-species coupling bug (which would produce order-
  // unity drift, not round-off).
  VC_TEST_ASSERT(boronRelDiff < 1e-6);
  VC_TEST_ASSERT(phosphorusRelDiff < 1e-6);
}
```

- [ ] **Step 2: Register the test in `main()`**

Add `TestMultiSpeciesSmoke();` immediately before `TestReentrantSolve();` in `main()`. The final `main()`:
```cpp
int main() {
  TestMeshAttributes();
  TestDiffusionModelInterface();
  TestConstantDiffusion();
  TestDiffusionPhysics();
  TestDiffusionPhysicsComposition();
  TestLevelSetToMesh2D();
  TestDiffusionEngineAssembly();
  TestDoseConservation();
  TestDirichletBC();
  TestNeumannBC();
  TestMultiSpeciesSmoke();
  TestReentrantSolve();
  TestCVODERuntimePath();
  std::cout << "All diffusion tests passed.\n";
  return 0;
}
```

- [ ] **Step 3: Build and run**

```bash
cmake --build build_followups --config Release --target testDiffusion
ctest -R testDiffusion --test-dir build_followups -C Release --output-on-failure
```
Expected: 1/1 PASS. The `[multi-species]` log lines should show both `rel_diff` values at ~1e-16. If either is large, there's a cross-species coupling bug in the engine (e.g. K from Boron's model is being applied to Phosphorus's GridFunction) - investigate the species-outer loop in `DiffusionEngine::assembleAllSpecies()` before proceeding.

- [ ] **Step 4: Commit**

```bash
git add tests/diffusion/testDiffusion.cpp
git commit -m "test(diffusion): add multi-species smoke test (F3)"
```

---

## Task 5: Run full test suite + update progress ledger

**Files:**
- Modify: `.superpowers/sdd/progress.md`

**Interfaces:** None. This is the close-out task.

- [ ] **Step 1: Build the full solution**

```bash
cmake --build build_followups --config Release --parallel
```
Expected: full solution builds. The new MFEM's transitive deps (HYPRE, METIS, SUNDIALS, MSMPI) link into every test target that depends on `ViennaPS`. If any link fails, check that `MFEMTargets.cmake`'s `INTERFACE_LINK_LIBRARIES` paths all exist (the user-build paths `F:/dev/sundials_build/install/lib/*.lib`, `F:/dev/vcpkg/installed/x64-windows/lib/HYPRE.lib`, etc.).

- [ ] **Step 2: Run the full test suite**

```bash
ctest -E "Benchmark|Performance|intermediate|removeStrayPoints" --test-dir build_followups -C Release
```
Expected: all tests pass. The 2 excluded tests (`intermediate`, `removeStrayPoints`) are pre-existing failures from an unrelated `psAnalyticImplant.hpp:9` include-path bug - they predate this work and are out of scope.

- [ ] **Step 3: Append close-out entry to `.superpowers/sdd/progress.md`**

Add:
```
## Phase 1 follow-ups close-out

Resolved F2, F3, F5 (F1, F4, F7 already done; F6 deferred to Phase 2+):
- F2 (shouldCreateTimeDerivative re-entrance): DiffusionPhysics::
  resetTimeDerivativeClaims() called at top of assembleAllSpecies().
  TestReentrantSolve confirms two consecutive solve() calls on the
  same physics object both conserve dose to ~1e-16.
- F3 (multi-species test): TestMultiSpeciesSmoke runs two independent
  ConstantDiffusion models on Boron + Phosphorus, each with closed-
  system BCs. Both doses conserved to ~1e-16 - confirms species-outer
  loop, allSpecies_ map, and packed-block CVODE layout are correct.
- F5 (MFEM rebuild): MFEMConfig.cmake + MFEMTargets.cmake paths fixed
  (build_full -> build). TestCVODERuntimePath prints which integrator
  and preconditioner are active; with the rebuilt MFEM, both
  MFEM_USE_SUNDIALS and MFEM_USE_MPI are defined, so CVODE BDF +
  HypreBoomerAMG now run end-to-end (previously compile-verified only).

Build dir: build_followups/. All Phase 1 tests pass (excluding 2
pre-existing failures in intermediate/removeStrayPoints from an
unrelated psAnalyticImplant.hpp include-path bug).
```

- [ ] **Step 4: Commit the ledger update**

```bash
git add .superpowers/sdd/progress.md
git commit -m "docs(sdd): close out Phase 1 follow-ups F2/F3/F5"
```

---

## Self-Review

**1. Spec coverage** (against `docs/superpowers/specs/diffusion-phase1-followups.md`):
- F1 (BCs applied): already done in commit `9ada745` - out of scope, noted in plan header
- F2 (shouldCreateTimeDerivative re-entrance): Task 3 ✓
- F3 (multi-species test): Task 4 ✓
- F4 (caching): already done in commit `9c37abb` - out of scope, noted in plan header
- F5 (MFEM rebuild): Task 1 (path fix) + Task 2 (runtime verification) ✓
- F6 (conforming cut-cell): explicitly deferred to Phase 2+ Task 4 - noted in plan header
- F7 (plan filename): already done in commit `9003343` - out of scope, noted in plan header

**2. Placeholder scan:** Searched for "TBD", "TODO", "implement later", "fill in details", "Add appropriate error handling", "Write tests for the above", "Similar to Task N" - none found. Every code step contains actual code. Every command step contains the exact command and expected output.

**3. Type consistency:**
- `resetTimeDerivativeClaims()` - defined in Task 3 Step 3, called in Task 3 Step 4. Spelling matches.
- `TestCVODERuntimePath()` - defined in Task 2 Step 1, registered in Task 2 Step 2. Spelling matches.
- `TestReentrantSolve()` - defined in Task 3 Step 1, registered in Task 3 Step 1 (the `main()` edit instruction). Spelling matches.
- `TestMultiSpeciesSmoke()` - defined in Task 4 Step 1, registered in Task 4 Step 2. Spelling matches.
- Build directory `build_followups/` - used consistently in Tasks 1, 2, 3, 4, 5.

**4. Test order in `main()`**: Tasks are ordered so each new test is inserted at a specific position. The final `main()` after all tasks: MeshAttributes, DiffusionModelInterface, ConstantDiffusion, DiffusionPhysics, DiffusionPhysicsComposition, LevelSetToMesh2D, DiffusionEngineAssembly, DoseConservation, DirichletBC, NeumannBC, MultiSpeciesSmoke, ReentrantSolve, CVODERuntimePath. Task ordering (2, 3, 4) puts CVODE first in the plan but last in `main()` - this is intentional so the re-entrance and multi-species tests (which don't depend on CVODE) run before the CVODE-specific test.

**5. Out-of-scope guardrails:** F6 (conforming cut-cell) is explicitly deferred - no task touches `LevelSetToMesh.hpp` or adds cut-cell machinery. F1/F4/F7 are noted as already-complete with commit SHAs so the implementer doesn't redo them.
