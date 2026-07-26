# Implementation Plan: Close Audit-Trail Carry-Forward Items

**Scope:** All 7 carry-forward items from the gap-analysis audit trail (`docs/superpowers/reviews/2026-07-25-spec-gap-analysis.md`). Organized into 5 dependency-ordered waves; each wave is independently committable and testable. Tests run via `cmake --build build_phase2 --config Release --target testDiffusion && ./build_phase2/tests/testDiffusion.exe`.

**Decisions made:**
- **III-V:** Add GaAs/InP material data + implement eq. 3-239/3-240
- **KMC epitaxy:** Real BKL Arrhenius events integrated into `step()` (replaces deterministic sweeps)
- **KMC lattice:** Keep 🟡 honest note, do NOT restructure to 2×FCC (existing recomb/cluster tests rely on cubic-stencil face-adjacency; a 2×FCC restructure would break them - documented limit, not a correctness bug)
- **Runtime AMR:** Full in-place refine + prolongation via `FESpace::Update(want_transform=true)` + `GridFunction::Update()`

---

## Wave 1 — Test gaps (turns 🟡 rows into verified ✅, no new physics) ✅ COMPLETED

### 1.1 Engine-solve tests for 5 untested FEM models (I2)
For each, register the model in a `DiffusionEngine`, run a short solve, assert a physics-sensitive outcome:

| Model | Test assertion |
|---|---|
| `SiGeCDiffusion` (FEM) | Register C + I + CI species; assert `CI > 0` AND `ΔC_s + ΔC_sI ≈ 0` to 0.1% |
| `GeBPairingModel` | Register Ge + B + GeB; assert mobile B decays, GeB grows, `C_GeB/(C_Ge·C_B) → k_f/k_b` |
| `StrainDiffusionModel` | Register with strain field; assert `D_eff(strained) ≠ D_eff(unstrained)` |
| `MeltDiffusion` (FEM `assembleStiffness`) | Register Boron + MeltFraction=0 and =1 variants; assert dose evolves faster with φ=1 |
| `PolysiliconDiffusion` (FEM GB path) | Call `enableGbSegregationSpecies` + `setAttributes`; assert GI/GB coupling fires |

### 1.2 `TestLevelSetToMesh3D` (I3)
Mirror `TestLevelSetToMesh2D` with a 3D domain: `Domain<double,3>::New()`, `MakePlane<double,3>`, `LevelSetToMeshConverter<double,3>`, assert `GetNV()>0`, `GetNE()>0` with tet elements, non-empty attribute set.

### 1.3 IDW deatomize dose conservation (I4)
Extend the IDW test: after `deatomizeIDW`, assert `Σ(out) ≈ original atom count` within ~30%. Add a `deatomizeToGridFunction` helper that scatters atoms onto an MFEM GridFunction (via `ProjectCoefficient` with a coefficient sampling atom densities from `KmcSite::x/y/z`), then assert `SolutionTransfer::integrate(gf) ≈ atomCount/volumePerSite` to 0.1%.

**Commit:** `"test(diffusion): close Wave 1 test gaps (5 FEM closures + 3D LS2Mesh + IDW dose)"` ✅

---

## Wave 2 — Cheap physics fixes (reuse existing patterns)

### 2.1 Poly GB face residual (C2) — reuses `Segregation::assembleInterfaceResidual` ✅ COMPLETED
**Files:** `include/viennaps/fields/models/PolysiliconDiffusion.hpp`
- Replace `assembleReaction`'s `DomainLFIntegrator` with an override of the **`finalizeReaction`** hook (called post-`Assemble` at `DiffusionEngine.hpp:514`).
- In `finalizeReaction`, call `SegregationCondition::assembleInterfaceResidual(R, C_gi_gf, C_gb_gf, mesh, interiorAttr_, boundaryAttr_, side)` for each side (1 and 2). This walks `Poly_GI:Poly_GB` faces and writes the mass-conserving residual directly into `R[id]`.
- Remove the old domain-integrated `SegCoef` path; keep `setFemSegregationRate`/`enableGbSegregationSpecies` API stable.
- **Test:** Wave 1.1 Poly GB test asserts the face-residual fires (GI loses mass, GB gains, sum conserved).
- **Ledger:** Poly GB row → ✅.

### 2.2 PDE `applyTo` auto-applyICs overload (M1)
**Files:** `include/viennaps/fields/PdeApi.hpp`
- Add an overload `applyTo(DiffusionPhysics<NumericType> &physics, Engine *engine = nullptr)` that, when `engine != nullptr`, also calls `applyICs(*engine)`. Keep the existing 1-arg overload (for backward compat with the 2 existing call sites at testDiffusion.cpp:719, 833).
- Update `TestPhase10PdeApi` (line 719) to use the 2-arg form and remove the manual `initializeSpecies` at line 722.
- **Ledger:** note M1 resolved.

### 2.3 HeatTransfer latent heat (M2)
**Files:** `include/viennaps/fields/models/FlashLaserAnneal.hpp`
- `HeatTransfer::assembleReaction` currently ignores `allSpecies` (line 55). Change it to look up `"MeltFraction"` from `allSpecies` (the same convention `MeltDiffusion` uses at line 144).
- Add a model-side `previous_phi_` GridFunction pointer + `setPreviousPhiStep()` API. Compute `∂φ/∂t ≈ (φ_current - φ_previous)/dt` per QP; add `ρ·L·∂φ/∂t` to the heat source reaction term.
- **Test:** with φ transitioning 0→1 (melting), assert T rise is less than the no-latent-heat case (energy goes into phase change).
- **Ledger:** Flash heat row note updated to include latent-heat coupling.

**Commit:** `"feat(diffusion): poly GB face residual (finalizeReaction), PDE auto-applyICs, HeatTransfer latent heat"`

---

## Wave 3 — Bigger physics implementations

### 3.1 III-V material data + eq. 3-239/3-240 (C1)
**Files:** `include/viennaps/fields/MaterialPropertySystem.hpp`, `include/viennaps/fields/models/IIIVDiffusion.hpp`
- **MaterialPropertySystem** (after line 150): add `GaAs` band data: `Nc=4.7e17`, `Nv=7.0e18`, `Eg=1.424` (300K); `InP`: `Nc=5.7e17`, `Nv=1.1e19`, `Eg=1.344`. Add `default` fallback for unknown III-V.
- **IIIVDiffusion::getDiffusivity**: implement:
  - Donor path (Si, Se, Ge): `D = D_AV·(n/ni) + D_AV²·(n/ni)²` (eq. 3-239)
  - Acceptor path (Be, Mg, Zn, C): `D = D_AI·(p/ni) + D_AI²·(p/ni)²` (eq. 3-240)
  - Use `IntrinsicCarrier::electronConcentration(C, T, material)` / `holeConcentration(...)` for n, p; `ni(T, material)` for ni.
  - Add `D_AV`, `D_AV²`, `D_AI`, `D_AI²` members set per (material, species) in the constructor; default values for GaAs Si/Se/Be from ATHENA Table 3-2.
  - Add a `setDonorSpecies(bool)` flag (or species-name lookup) to pick donor vs acceptor formula.
- **QP coefficient:** replace `ConstantCoefficient` at line 70 with a `IIIVDCoef` (mirror `ChargedFermiDCoef`) that reads the dopant GF + n/p per QP. Reuse the `IntrinsicCarrier` lookup.
- **Test:** Si in GaAs → D concentration-independent (first term only); Zn in GaAs → D ∝ (p/ni)² (assert `D(C=1e20) >> D(C=1e15)` for acceptor).
- **Ledger:** III-V row → ✅.

### 3.2 KMC epitaxy BKL events (C3)
**Files:** `include/viennaps/fields/kmc/KmcEvent.hpp`, `include/viennaps/fields/kmc/KmcAtomisticEngine.hpp`, `include/viennaps/fields/kmc/KmcEpitaxy.hpp`
- **KmcParameters** (`KmcEvent.hpp:25-50`): add `attachPreFactor`/`attachBarrier`, `desorbPreFactor`/`desorbBarrier`, `twinPreFactor`/`twinBarrier` + matching `attachRate()`/`desorbRate()`/`twinRate()` methods (Arrhenius form, same pattern as existing `hopRate()` etc.).
- **KmcAtomisticEngine::step()**: after the existing event loop, add a surface scan — for each column, find the topmost occupied site + the empty site above it; emit `Deposit` (attach) events at rate `attachRate()` for the empty surface site, `Desorb` events for occupied surface sites at `desorbRate()`, and `Twin` events at `twinRate()` with a per-site random draw.
- **apply()** (`KmcAtomisticEngine.hpp:200-234`): handle the already-declared `Deposit`/`Desorb` enum values (currently unused). `Deposit` fills the empty site with the growth species; `Desorb` clears it. Add a new `Twin` case that flips the stacking flag on a {111} site (use the existing `species` field with a new `KmcTwin` code).
- **Keep** `KmcEpitaxyModel` as-is but mark it deprecated/🟡 (deterministic helper). Add a `KmcEpitaxyEngine` thin wrapper that configures `KmcAtomisticEngine` with attach/desorb/twin rates for epitaxial growth.
- **Test:** run `KmcEpitaxyEngine` for N steps at high T; assert `depositCount > 0`, growth occurs (column heights increase), `twinCount` scales with `twinRate`.
- **Ledger:** KMC epitaxy row → ✅ (BKL path); deterministic helper stays 🟡.

**Commit:** `"feat: III-V eq. 3-239/3-240 with GaAs material data; KMC BKL epitaxy events"`

---

## Wave 4 — Engine infrastructure

### 4.1 Runtime AMR with prolongation (I1)
**Files:** `include/viennaps/fields/DiffusionEngine.hpp`, `include/viennaps/fields/AdaptiveMeshRefiner.hpp`
- **`applyRuntimeAmr`** (line 927): replace the mark-only stub with the full sequence:
  1. `AdaptiveMeshRefiner::markBox` (existing) → `ids`
  2. `mesh_->GeneralRefinement(el_to_refine)` (existing offline path)
  3. `fes_->Update(/*want_transform=*/true)` — computes the prolongation transform
  4. For each species GridFunction: `gf->Update()` — applies the transform in place
  4. `fes_->UpdatesFinished()` — releases the stored transform
  5. Invalidate the cached `M`/`K` matrices (`implicitCache_`) so they rebuild on the next step
- **`AdaptiveMeshRefiner::refineMarked`** (line 99): add a sibling `refineMarkedWithProlongation(mesh, fes, gridFunctions, ids)` that does steps 2-5 for the offline path too.
- **Test:** `TestProductionAmr` (line 1364): after refine, assert the pre-refine dose is preserved on the post-refine GridFunction to 0.1% (currently the pre-refine GridFunction is discarded). Add `runtimeAmrRefineCount() > 0` AND a dose-preservation assertion.
- **Ledger:** runtime AMR row → ✅.

### 4.2 KMC incremental event rebuild (audit item 4)
**Files:** `include/viennaps/fields/kmc/KmcAtomisticEngine.hpp`
- Replace the flat `events_` vector + per-step full rebuild with:
  - Per-site event lists: `std::vector<std::vector<KmcEvent>> eventsBySite_` indexed by site linear index `(k*ny+j)*nx+i`.
  - A Fenwick tree (binary indexed tree) `fenwick_` over per-site total rates for O(log N) selection + O(log N) update.
- After `apply(event)`, recompute event lists only for the 2 touched sites + their neighbors (≤14 sites under cubic stencil), updating the Fenwick tree incrementally.
- Keep the public API (`step()`, `run()`, counters) stable.
- **Test:** double the lattice size (8³ → 16³); assert per-step cost grows ~linearly with affected-site count, not with total sites (timing or operation-count assertion). Existing recomb/cluster/dissoc tests still pass.
- **Ledger:** KMC event tree row → ✅.

**Commit:** `"feat: runtime AMR with GridFunction prolongation; KMC incremental Fenwick event tree"` ✅ WAVE 4 DONE

---

## Wave 5 — Template migration

### 5.1 MobileImpurity SpeciesTag template (audit item 5)
**Files:** `include/viennaps/fields/models/MobileImpurity.hpp`, `include/viennaps/fields/models/CopperDiffusion.hpp`, NEW `include/viennaps/fields/models/MobileImpurityTags.hpp`
- Add a second template parameter with a default: `template <class NumericType, class SpeciesTag = GenericImpurityTag> class MobileImpurity`. The tag struct carries compile-time constants: `name`, `D0`, `Ea`, `charge`. `GenericImpurityTag` defaults to current runtime-string behavior (so existing 4 call sites at testDiffusion.cpp:1470, 1474, 1532, 1565, 1636 compile unchanged).
- `CopperDiffusion<NumericType>` becomes `using CopperDiffusion = MobileImpurity<NumericType, CopperTag>`. The `MobileImpurity` constructor initializes `D0_` and `z_` from `SpeciesTag::D0` / `SpeciesTag::charge`.
- Add `SodiumTag`, `IronTag` in `MobileImpurityTags.hpp` wired to Phase 10 DB values.
- **Test:** `MobileImpurity<double, SodiumTag>` assertion on D0 matches tag value; existing Copper tests use `setIonPairing(1.0)` explicitly since the tag does not carry beta.
- **Ledger:** audit item 5 resolved.

**Commit:** `"refactor: MobileImpurity<NumericType, SpeciesTag> template with Cu/Na/Fe tags"` ✅ WAVE 5 DONE

---

## Final wave — Ledger reconciliation

### 6.1 Update gap-analysis ledger + audit trail
**Files:** `docs/superpowers/reviews/2026-07-25-spec-gap-analysis.md`
- After each wave's commit, update the corresponding ledger rows to ✅ with the file:line evidence.
- After Wave 5, add a Pass-5 audit-trail entry summarizing what closed.
- KMC lattice row stays 🟡 (honest note, not restructured).

**Commit:** `"docs(gap-analysis): reconcile ledger after Pass-5 implementation cycle"`

---

## Per-item ledger row closures (summary)

| Wave | Item | Ledger row(s) closed | Status change |
|---|---|---|---|
| 1.1 | 5 FEM model tests | SiGeC, GeB, Strain, MeltDiffusion, Poly-GB (test side) | 🟡→✅ (verified) |
| 1.2 | TestLevelSetToMesh3D | 3D LevelSetToMesh | 🟡→✅ |
| 1.3 | IDW dose conservation | KMC IDW deatomize | 🟡→✅ |
| 2.1 | Poly GB face residual | Poly GB segregation | 🟡→✅ (physics) |
| 2.2 | PDE auto-applyICs | M1 | resolved |
| 2.3 | HeatTransfer latent heat | Flash heat | note updated |
| 3.1 | III-V eq. 3-239/3-240 | III-V FEM + I/V eq | 🟡→✅ |
| 3.2 | KMC BKL epitaxy | KMC epitaxy | 🟡→✅ (BKL) |
| 4.1 | Runtime AMR prolongation | Runtime AMR | 🟡→✅ |
| 4.2 | KMC Fenwick tree | KMC event tree | 🟡→✅ |
| 5.1 | SpeciesTag template | audit item 5 | resolved |
| — | KMC 2×FCC lattice | diamond lattice | stays 🟡 (honest) |

---

## Build/test discipline
- After each wave: `cmake --build build_phase2 --config Release --target testDiffusion && ./build_phase2/tests/testDiffusion.exe`. All tests must pass before committing.
- After each wave: `git add -p` the relevant files + the ledger update, commit with the message above.
- If a wave's tests fail, fix before moving on (don't accumulate regressions - the prior review cycles showed that pattern compounds).

## Estimated scope
- ~12 files modified, ~2-3 new files, ~600-900 lines of code + tests across all 5 waves.
- Each wave is 1 commit; total 6-7 commits including the final ledger reconciliation.