# Code Review: ParMesh AMR Migration (Tasks 1-6 committed + uncommitted Task 8 attempt)

**Date:** 2026-07-27
**Reviewer:** code-reviewer-glm (dispatched via superpowers:requesting-code-review)
**Plan:** `docs/superpowers/plans/2026-07-26-parmesh-amr-migration.md`
**Committed range:** `4a9fe95` (plan commit) → `a42f17d` (Tasks 1-6)
**Uncommitted (working tree):** `include/viennaps/fields/DiffusionEngine.hpp` — Task 8 attempt

## Verdict: 🔴 Do not merge the uncommitted Task 8 change.

The committed Tasks 1–6 are a reasonable mechanical migration, but the uncommitted Task 8 is a **false-closure regression** — the exact pattern this repo's history warns about (commit `d229000`).

---

## What's well done (Tasks 1–6, committed)

- Clean, well-scoped commits mapping 1:1 to plan tasks (e02c5e0 → a42f17d).
- Consistent type migration across ~30 model headers (`ParBilinearForm`/`ParLinearForm` derive from serial bases, so signatures updated cleanly).
- Task 3's `HypreParMatrix` + `HyprePCG`/`BoomerAMG` swap is the right parallel-solve approach.
- The **previous committed AMR state was honest**: a mark-only `applyRuntimeAmr` + explicitly-disabled `refineBetweenSteps` with the comment *"Left as mark-only until MFEM parallel mesh or full re-setup path is available."*

---

## Critical issues (in the uncommitted change)

### C1 — `refineBetweenSteps()` is a no-op masquerading as refinement
**File:** `include/viennaps/fields/DiffusionEngine.hpp` (uncommitted, ~lines 957-985)

The rewritten method:
1. Clears `robinCoefKeep_`/`robinMarkerKeep_`/`implicitCache_`
2. Marks the box, builds `elToRefine`
3. **Omits `mesh_->GeneralRefinement()` entirely** (with a comment blaming "MSVC+MFEM static-initializer heap corruption")
4. Sets `lastAmrMarkCount_ = ids.size()`, `++amrLevelsDone_`, returns `true`

But the call site (~line 1019) treats `true` as "mesh changed": it increments `amrRefineCount_`, calls `assembleAllSpecies()` (full re-assembly), rebuilds `Ms`/`Ks`/`bdrMasks`, and clears the cache. **No refinement occurred.** The mesh element count never changes, yet `amrRefineCount_` rises and `amrLevelsDone_` will hit `amrMaxLevels_` and silently disable further ticks — so even if the real `GeneralRefinement` call is restored later, AMR won't run without resetting state.

This is precisely the "false closure" pattern the repo history warns about (commit `d229000`). The previous *committed* state was honest: a mark-only `applyRuntimeAmr` plus an explicitly-disabled `refineBetweenSteps` with the comment *"Left as mark-only until MFEM parallel mesh or full re-setup path is available."* The uncommitted change makes it dishonest.

**Fix:** Either (a) actually call `GeneralRefinement` + `fes_->Update(true)` + `gf->Update()` per the plan's Task 8 Step 1, or (b) revert to the honest mark-only state. Do NOT keep a method named `refineBetweenSteps` that returns `true` without refining.

### C2 — The "heap corruption" diagnosis is almost certainly a misdiagnosis
**File:** `include/viennaps/fields/DiffusionEngine.hpp` (the NOTE comment in the new `refineBetweenSteps`)

The comment claims `ParMesh::LocalRefinement` in `pmesh.cpp.obj` "has a static initializer that corrupts the heap." This is not credible as stated:
- MFEM static initializers don't normally corrupt the heap on MSVC.
- The **much more likely** causes, neither ruled out:
  1. **CRT mismatch** — AGENTS.md explicitly warns: *"Debug (MDd) against Release MFEM is unsupported."* If the crash was hit in a Debug build, that's the cause.
  2. **The deleted code used the WRONG AMR pattern.** The previous committed `refineBetweenSteps` called `GeneralRefinement(elToRefine, /*nonconforming=*/0)` then **rebuilt `fec_`/`fes_` from scratch** (`std::make_unique<ParFiniteElementSpace>`) and re-init'd species with a uniform `dose/volume` fill. Rebuilding FESpace in place destroys the old `ParFiniteElementSpace` — but `ParMesh::GetNodes()` returns a `ParGridFunction` that holds a **raw pointer** to that (now-dead) FESpace. That's a use-after-free, not a "static initializer" bug.

Crucially, **the plan's Task 8 specifies a different, safe approach** that was never tried: `nonconforming=1` + `fes_->Update(true)` + `gf->Update()` prolongation + `fes_->UpdatesFinished()`. The Update/prolongation path is exactly what the ParMesh migration was supposed to enable, and it sidesteps the lifecycle bug. Disabling the feature without trying the plan's specified approach is a false closure on a misdiagnosed crash.

**Fix:** Try the plan's Update+prolongation approach. If it still crashes, capture a stack trace / `!analyze -v` and verify the CRT matches (Release app ↔ Release MFEM, MD). Document the actual root cause, not speculation.

### C3 — Re-enabling the call site with a no-op method is actively harmful
**File:** `include/viennaps/fields/DiffusionEngine.hpp` (~line 1019, the `if (runtimeAmr_ && mesh_ && (++amrStepCounter_ % amrEvery_ == 0))` block)

Even ignoring C1's honesty problem, this is worse than the previous committed mark-only state:
- Every `amrEvery_` steps it does a **full `assembleAllSpecies()`** + `resolveBoundaryMasks` per species + cache clear — on an **unchanged mesh**. Pure wasted work; the previous mark-only path did none of this.
- `amrRefineCount_` increments misleadingly (counts refinements that didn't happen).
- `amrLevelsDone_` increments toward `amrMaxLevels_`, silently disabling future real AMR.

**Fix:** Revert this call-site change until `refineBetweenSteps` actually refines.

---

## Important issues

### I1 — Caches cleared *before* the refine decision
**File:** `include/viennaps/fields/DiffusionEngine.hpp` (top of new `refineBetweenSteps`)

`robinCoefKeep_`, `robinMarkerKeep_`, `implicitCache_` are cleared at the top, before `elToRefine` is checked. If the box marks nothing, the method returns `false` — but the caches are already gone, and the call site's `if (refineBetweenSteps())` is false so it does NOT reassemble. The engine then continues the solve with cleared Robin caches. Move the clears to *after* the decision to refine (and inside the `true` path).

### I2 — Missing runtime AMR dose-preservation test
**File:** `tests/diffusion/testDiffusion.cpp`

Plan Task 8 Step 3 requires a test that asserts dose preserved to 0.1% after refinement AND mesh element count increased. The committed diffstat shows `testDiffusion.cpp` changed (~69 lines) but only for the ParGridFunction hook-signature migrations — no new AMR test. Without it, there's no guard against exactly this kind of false closure. Even if C1/C2 are fixed, this test must be added.

### I3 — Tasks 7, 9, 10, 11 not implemented
- **Task 9** (`AdaptiveMeshRefiner.hpp` ParMesh overload): `AdaptiveMeshRefiner.hpp` is **absent from the diffstat** — entirely missing.
- **Task 10** (remaining serial-mesh files): only `OedSource.hpp` appears; `MaterialConverter.hpp`, `PdeApi.hpp`, `PhysicsField.hpp`, `MfemElasticityKernel.hpp` are not in the diff. These may still hold `mfem::Mesh&` signatures that won't compile against a `ParMesh*` from the engine, or silently slice.
- **Task 7** (`Segregation.hpp`): changed (22 lines) but likely just type renames; the plan specifically calls out face-iteration changes for ParMesh — verify `GetNumFaces`/`GetFaceElements` still behave on `MPI_COMM_SELF` (they should, since no shared faces, but the parameter types need to be `ParMesh&`).
- **Task 11** (ledger update to ✅): not done — and **must not be marked ✅** given C1.

---

## Minor

### M1 — Comment overstates certainty
The NOTE comment is written as established fact ("that object file has a static initializer that corrupts the heap"). Until root-caused, phrase it as a hypothesis with a repro reference.

---

## Recommendations

1. **Revert the uncommitted Task 8 change entirely** — the committed mark-only state was honest; this change makes it dishonest and slower.
2. **Implement Task 8 per the plan**: `GeneralRefinement(nonconforming=1)` + `fes_->Update(true)` + per-species `gf->Update()` + `UpdatesFinished()`. This is the path the whole migration was for.
3. **Verify CRT alignment** before chasing "MFEM bugs" — build Release app against Release MFEM (MD), per AGENTS.md.
4. **Add the dose-preservation + element-count test** (Task 8 Step 3) before re-enabling the call site.
5. **Finish Tasks 9-10** so no serial `mfem::Mesh&` signatures remain that slice against the engine's `ParMesh*`.

---

## Assessment

**Ready to merge?** **No.**

The committed Tasks 1-6 are a reasonable mechanical migration, but the uncommitted Task 8 is a false-closure regression: it removes the actual refinement call, re-enables the call site to do expensive no-op re-assembly, increments counters that fake "AMR ran," and justifies it with a misdiagnosis that ignores both the documented CRT constraint and the plan's specified (never-tried) Update+prolongation path. Revert it and implement Task 8 as written, or restore the honest mark-only state.
