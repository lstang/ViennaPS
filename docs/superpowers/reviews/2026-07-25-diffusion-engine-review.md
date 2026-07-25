# Code Review: Diffusion Engine (origin/master..b23d3a4)

**Date:** 2026-07-25
**Scope:** 37 commits, 103 files, ~16,276 lines. Multiphysics field layer for ViennaPS implementing FEM-based diffusion engine with MFEM + SUNDIALS CVODE.
**Base:** 2956ed5 (origin/master)
**Head:** b23d3a4 (zcode)

---

## Strengths

- **Excellent documentation**: `DiffusionEngine.hpp` thoroughly documents design decisions, MOOSE references, caching strategy, and the F9 Heisenbug root cause. Every model cites its MOOSE counterpart per ADR-0001.
- **Correct MFEM coefficient lifetime management**: `ConstantDiffusion.hpp:63-67` caches coefficients in `mutable unique_ptr` to prevent dangling references when `BilinearForm::Assemble()` runs after `assembleStiffness()` returns. This was a Phase 1 bug, now correctly fixed.
- **F9 Heisenbug fix**: `DiffusionEngine.hpp:736` correctly sets `iterative_mode = false` on inner CGSolver, preventing stale `tmp_` buffer from causing NaN.
- **F2 re-entrance fix**: `DiffusionPhysics.hpp:68` `resetTimeDerivativeClaims()` called at `DiffusionEngine.hpp:413` at the top of `assembleAllSpecies()`.
- **System-matrix caching**: The `EliminateVDofs`/`EliminateVDofsInRHS` pattern (`DiffusionEngine.hpp:951-986`) is the canonical MFEM pattern.
- **Robin BC**: `DiffusionEngine.hpp:470-473` correctly adds `BoundaryMassIntegrator` implementing the weak form of `-D dC/dn = hC`.
- **Species-outer assembly loop**: Mirrors MOOSE `MultiSpeciesDiffusionCG::addFEKernels()`.
- **Test quality**: Dose conservation (1.28e-16), Dirichlet/Neumann/Robin BC verification, re-entrance regression test, multi-species K-leakage test - all verify real physics.
- **Physics correctness**: `IntrinsicCarrier` ni formula, Fermi-Dirac activity coefficient, Arrhenius diffusivity, segregation two-compartment model - all mathematically correct.

---

## Issues

### Critical (Must Fix)

**1. Broken include paths break existing tests**
- **File:** `include/viennaps/models/psAnalyticImplant.hpp:9`, `include/viennaps/models/psBasicDiffusion.hpp:11`
- **Issue:** `#include "psProcessModel.hpp"` doesn't resolve - the file is at `include/viennaps/process/psProcessModel.hpp`. All other model files use `"../process/psProcessModel.hpp"`. Since `viennaps.hpp:45,67` includes these unconditionally, any TU including `viennaps.hpp` fails to compile, breaking `intermediate` and `removeStrayPoints` tests.
- **Fix:** Change to `#include "../process/psProcessModel.hpp"` in both files.

**2. Non-standard `#include <viennals.hpp>` umbrella**
- **File:** `include/viennaps/models/psAnalyticImplant.hpp:15`
- **Issue:** Only file in the codebase using the `<viennals.hpp>` umbrella. ViennaLS via CPM may not expose it, and `psDomain.hpp` (line 10) already pulls in all needed ViennaLS headers.
- **Fix:** Remove line 15.

**3. Hardcoded `VIENNAPS_HAS_SUNDIALS` in test source**
- **File:** `tests/domain/domain.cpp:1`
- **Issue:** `#define VIENNAPS_HAS_SUNDIALS 1` before includes forces SUNDIALS code paths regardless of CMake detection. Non-portable.
- **Fix:** Remove the define. `CMakeLists.txt:319` already propagates it via the INTERFACE target.

**4. Unconditional SUNDIALS define in test CMake**
- **File:** `tests/domain/CMakeLists.txt:5`
- **Issue:** `target_compile_definitions(... PRIVATE VIENNAPS_HAS_SUNDIALS)` is unconditional, even when SUNDIALS isn't found.
- **Fix:** Remove or gate with `if(SUNDIALS_FOUND)`.

**5. CddDiffusion: use-after-free with multiple bimolecular reactions**
- **File:** `include/viennaps/fields/models/CddDiffusion.hpp:149-151`
- **Issue:** In `assembleReaction()`, the loop over `reactions_` overwrites `prodCoef_` (a `mutable unique_ptr<ProductCoef>`) on each iteration. Each iteration registers `DomainLFIntegrator(*prodCoef_)` with R, which stores a `Coefficient&` reference. When the next iteration creates a new `prodCoef_`, the old `ProductCoef` is destroyed, but the `DomainLFIntegrator` from the previous iteration still holds a dangling reference. When `R.Assemble()` is called later (`DiffusionEngine.hpp:478`), it evaluates all integrators including the one with the destroyed coefficient - use-after-free.
- **Trigger:** Any `CddDiffusion` instance with 2+ bimolecular reactions in `reactions_` that match the current species. The API explicitly supports this via `addReaction()`.
- **Fix:** Use a `std::vector<std::unique_ptr<ProductCoef>>` to own all coefficients simultaneously, instead of a single `prodCoef_` that gets overwritten.

### Important (Should Fix)

**6. ConstantDiffusion returns D0 at T<=0 instead of 0**
- **File:** `include/viennaps/fields/models/ConstantDiffusion.hpp:40-42`
- **Issue:** `if (this->T_ <= 0) return D0_;` returns the pre-exponential factor when T<=0. Physically, D->0 as T->0 (frozen). The guard prevents division by zero but gives a wrong result that could mask bugs.
- **Fix:** Return `NumericType(0)` for T<=0.

**7. CddDiffusion tedTimeline C_I_eq guard uses 1**
- **File:** `include/viennaps/fields/models/CddDiffusion.hpp:178`
- **Issue:** `D_pair * CI / std::max(C_I_eq, NumericType(1))` - if C_I_eq is 0, using 1 as denominator gives unphysically large D. C_I_eq is typically ~1e10-1e15 cm^-3.
- **Fix:** Return 0 if C_I_eq <= 0.

**8. TestDirichletBC tolerance too loose**
- **File:** `tests/diffusion/testDiffusion.cpp:1286`
- **Issue:** `VC_TEST_ASSERT(interiorMean > 0.5e18)` allows 50% of clamped value. Observed baseline is 7.6e17 (76%), so the threshold wouldn't catch a regression to 55%.
- **Fix:** Tighten to `> 0.7e18`.

**9. Monolithic domain.cpp test with shared mutable state**
- **File:** `tests/domain/domain.cpp:40-509`
- **Issue:** ~470-line `RunTest()` function with ~15 sub-tests sharing a `field` object. If one block fails, subsequent blocks may produce misleading failures.
- **Fix:** Extract sub-tests into separate functions with their own field instances.

**10. Progress.md incorrectly claims include-path bug "predates this branch"**
- **File:** `.superpowers/sdd/progress.md:79`
- **Issue:** States the `psAnalyticImplant.hpp` include-path bug "Predates this branch." However, `psAnalyticImplant.hpp` was introduced by commit `c53d776` which is part of this branch.
- **Fix:** Correct the documentation.

**11. Duplicate MFEM path detection in CMakeLists.txt**
- **File:** `CMakeLists.txt:55-71` and `CMakeLists.txt:239-246`
- **Issue:** MFEM path is detected in two places with the same `EXISTS` check. The second block's fallback is dead code.
- **Fix:** Consolidate into a single block.

**12. ReactDiffusion: same mutable coefficient pattern as CddDiffusion**
- **File:** `include/viennaps/fields/models/ReactDiffusion.hpp:81-83`
- **Issue:** `recombCoef_` is a mutable `unique_ptr` overwritten in `assembleReaction`. Currently safe because there's only one recombination reaction, but the pattern is the same latent bug as CddDiffusion #5. If `assembleReaction` is called for both I and V species, the first call's `recombCoef_` is destroyed before `R.Assemble()` is called for that species.
- **Wait** - actually the engine calls `R.Assemble()` within the same species iteration (`DiffusionEngine.hpp:478`), so the coefficient is alive during Assemble. This is safe for the current single-reaction case. But it's a fragile pattern.
- **Fix:** Use a vector of owned coefficients for robustness, matching the fix for #5.

### Minor (Nice to Have)

**13. CVODE abstol magic number**
- **File:** `DiffusionEngine.hpp:613`
- **Issue:** `cvode.SetSStolerances(1e-6, 1e5)` - abstol=1e5 is appropriate for dopant concentrations (1e10-1e20) but too large for defect species (1e5-1e15). Add a comment or make it species-dependent.

**14. DiffusionModel.hpp over-include `<map>`**
- **File:** `include/viennaps/fields/DiffusionModel.hpp:7`
- **Issue:** `<map>` only needed when `VIENNAPS_HAS_MFEM` is defined (for `std::map<std::string, mfem::GridFunction*>`). Move inside the `#ifdef`.

**15. CMakeLists.txt redundant case check on Windows**
- **File:** `CMakeLists.txt:68`
- **Issue:** `if(EXISTS "f:/dev/mfem/build/MFEMConfig.cmake" OR EXISTS "F:/dev/mfem/build/MFEMConfig.cmake")` - Windows paths are case-insensitive.

**16. KMC atomize/deatomize round-trip test only checks size**
- **File:** `tests/diffusion/testDiffusion.cpp:391`
- **Issue:** `VC_TEST_ASSERT(back.size() == lat2.size())` verifies size but not that concentrations are preserved.

**17. `const_cast` in Mult and Segregation**
- **Files:** `DiffusionEngine.hpp:719`, `Segregation.hpp:88`
- **Issue:** `const_cast` to strip const from MFEM Vector/FiniteElementSpace. Common MFEM workaround, but should be commented as intentional.

**18. Per-step solver allocation in implicit Euler**
- **File:** `DiffusionEngine.hpp:990-994`
- **Issue:** `makeMassSolver()` creates a new CGSolver+DSmoother per step per species. Could cache the solver like the CVODE path does. Negligible cost for Phase 1 mesh sizes.

---

## Recommendations

1. **Fix Critical #1-5 immediately** before any merge. Issues #1-4 are one-line fixes. Issue #5 requires changing `prodCoef_` to a vector.
2. **Apply the same coefficient ownership pattern from #5 to ReactDiffusion** for consistency and future-proofing.
3. **Consider a clang-tidy run** with `bugprone-dangling-handle` and `cppcoreguidelines-pro-bounds-constant-array-index` to catch similar patterns.
4. **Run the format target** (`cmake --build build --target format`) before merge - some files may exceed 80-col limit.

---

## Assessment

**Ready to merge? No - with fixes.**

**Reasoning:** The core engine architecture is sound with correct MFEM API usage, well-documented design decisions, and physics-meaningful tests. However, 5 Critical issues must be fixed first: broken include paths (#1-2) that prevent compilation, hardcoded SUNDIALS defines (#3-4) that break portability, and a use-after-free in CddDiffusion (#5) that will crash when multiple reactions are configured. All are straightforward fixes - #1-4 are one-line changes, #5 requires changing a single `unique_ptr` to a `vector<unique_ptr>`.
