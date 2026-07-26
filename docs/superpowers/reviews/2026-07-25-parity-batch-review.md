# Code Review: 3-Commit Parity Batch (`0b9ca6d` + `8d414eb` + `367629e`)

**Date:** 2026-07-25
**Scope:** 3 commits, 21 files, +1677/−449. Closes gap-analysis rows across mesh, PDE, KMC, results, polysilicon, SiGeC, flash/laser, III-V, SPER. Includes the `367629e` drift-sign + SiGe `&&` Critical bug fixes from the prior review pass.
**Base:** `ac18975d` (last reviewed commit)
**Head:** `367629e` (current HEAD)
**Build/test:** `testDiffusion` Release - "All diffusion tests passed." (verified externally; judged by reading).

**Prior reviews in this cycle:**
- `2026-07-25-phases3-11-fulldepth-review.md` (review of `c475863`)
- Audit trail section in `2026-07-25-spec-gap-analysis.md`

**Recurring failure pattern (carried from prior reviews):** gap-ledger ✅ marks repeatedly overstate the code - FEM scaffolding without the cited physics, or tests that bypass the physics they claim to verify. This batch improves on `ac18975` (the bug fixes are real and well-tested) but the pattern persists in 3 rows.

---

## Strengths

- **`367629e` drift-sign fix is correct and well-tested.** `MobileImpurity.hpp:144` sets `a(0) = Drep * bNP * Ex_` (positive), with a clear IBP-identity derivation in the comment (lines 128-135). The E2 test (`testDiffusion.cpp:1557-1628`) is exactly the kind of direction-sensitive check the prior review demanded: non-uniform IC peaked at x_min, E along +x, z=+1, asserts `xCentroid1 > xCentroid0` AND `delta > 1e-8`. Would catch a sign inversion. Genuine closure.
- **`367629e` SiGe `&&` fix is correct and fail-loud.** `SiGeDiffusion.hpp:135` requires `useDefectMediated_ && CI && CV`; lines 140-147 add `MFEM_VERIFY(false, ...)` when only one defect field is registered. `DefectDCoef::Eval` (lines 101-108) also guards `if (CI_ && CV_)`. Genuine fix.
- **`LinearReactionDiffusion.hpp` is a real FEM reaction model.** `ScaledGF` coefficient (lines 47-58) reads the species GridFunction at QPs and scales by `-k`, fed to `DomainLFIntegrator`. `PdeApi.hpp:120-124` wires `ReactionPdeTerm` to it; `applyICs` (lines 175-179) calls `engine.initializeSpecies`. The PDE-reaction test (lines 820-842) asserts `d1 < d0` with a strong decay rate. Genuine closure for PDE-reaction row.
- **Jacobian strategy (b) is genuinely wired and tested.** `FermiDiffusion.hpp:91-99` implements `assembleStiffnessJacobian` with real analytic `dD/dC` via `MixedGradGradIntegrator`. Engine calls it at `DiffusionEngine.hpp:478` when `enableJacobian_` is true; test asserts `jacobianAssemblyCalls() > 0`. Genuine closure.
- **`SiGeCDiffusion` correctly inherits the carbon-suppression FEM residual.** Inherits from `CarbonDiffusion`, which assembles the `C_s + I ⇌ C_sI` residual with QP-local `TrapCoef`. Physics cited in Phase 6 plan is genuinely present (code-correct; test is partial - see I2).
- **`HeatTransfer` FEM is genuinely assembled and tested.** `FlashLaserAnneal.hpp:45-66` assembles `α·∇T·∇v` stiffness + `G` source. Test asserts `t1 > t0`. Genuine closure for the FEM-heat sub-claim.
- **Sub-cycling is genuinely in the solve loop.** `DiffusionEngine.hpp:990-998` divides `dtOuter` by `subCycles_` and runs the micro-step loop. Tested at line 971. Genuine closure.
- **L2 solution transfer is the real path and tested to 1e-16.** `SolutionTransfer::transferL2` (lines 91-141) builds `M` + `b` and runs CG; engine's `transferIntegralPreserving` defaults to it (line 147). Test asserts `relativeDoseError <= 1e-3`. Genuine closure.
- **KMC diamond-lattice row is honestly marked.** Ledger line 55 says "Body-diagonal opposite-sublattice (not full 2×FCC coords)" rather than claiming a full 2×FCC rewrite.

---

## Gap-Row Closure Verification

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| 1 | 3D LevelSetToMesh Cartesian tets | 🟡 Partial | `LevelSetToMesh.hpp:178-270` genuinely generates 3D tets via `MakeCartesian3D(..., TETRAHEDRON)`. Never instantiated in any test/example/engine call. "Code right, never run." |
| 2 | MovingMesh remesh trigger + Laplacian | 🟡 Partial | `laplacianSmooth` (line 129) and `needsRemesh` (line 207) exist and are unit-tested, BUT `DiffusionEngine::solve()` never calls them - standalone helpers only. |
| 3 | SolutionTransfer L2 Mx=b | ✅ Confirmed | `transferL2` is the path the engine calls. Tested to 1e-3. |
| 4 | Runtime AMR in solve | ❌ **False closure** | `applyRuntimeAmr` (`:927-955`) marks + counts but NEVER refines. Comment at 943-945 admits "Full GeneralRefinement + FESpace::Update can corrupt serial H1 spaces." Test asserts `runtimeAmrRefineCount() > 0` - passes because hook runs, no actual AMR. |
| 5 | Poly GB FEM segregation residual | ❌ **False closure** | `PolysiliconDiffusion.hpp:138` uses `DomainLFIntegrator` - a DOMAIN integrator, not a two-sided interior-face residual at `Poly_GI:Poly_GB` faces (Phase 5 plan lines 183-196 require face integrator mirroring MOOSE `InterfaceReaction`). `enableGbSegregationSpecies` never called in any test. |
| 6 | SiGeC carbon suppression FEM | 🟡 Partial | FEM residual genuinely inherited from `CarbonDiffusion` (correct physics). But the only SiGeC test (lines 494-499) calls static `tedFactor`, never instantiates `SiGeCDiffusion` in an engine solve. |
| 7 | GeB pairing FEM | 🟡 Partial | `GeBPairingModel` (`SiGeCDiffusion.hpp:43-115`) is a genuine 3-species FEM reaction with `kf·B·G − kr·P` (matches eq. 332). But the only GeB test (lines 501-505) uses the OLD 1D `GeBPairing::applyStep`, not the FEM model. |
| 8 | Strain FEM | 🟡 Partial | `StrainDiffusionModel` (`SiGeCDiffusion.hpp:118-165`) is genuine FEM with `D0·exp(-alpha·eps/kT)`. But the only strain test (lines 531-533) uses OLD `StrainDiffusionModifier` coefficient helper. |
| 9 | III-V FEM + I/V equilibrium | ❌ **False closure** | `IIIVDiffusion::getDiffusivity` (`:40-45`) is plain `D0·exp(-Ea/kT)` - NO `(n/ni)` or `(p/ni)` dependence, no Ga/As sublattice distinction. Plan requires eq. 3-239/3-240. `C_I_eq`/`C_V_eq` are higher-Ea Arrhenius proxies, not 4-sublattice equilibrium. Plan line 383 defers I/V eq to 🟡 but ledger line 46 marks ✅. |
| 10 | SPER orientation | ✅ Code / 🟡 test | `SPERKernel.hpp:28-36` implements orientation factors (100)=1, (110)=0.7, (111)=0.5; line 55 applies to velocity. Test only asserts the setter, not orientation-dependent velocity in a solve. |
| 11 | Flash FEM heat + MeltDiffusion FEM | 🟡 Partial | Heat FEM genuine and tested. `MeltDiffusion::assembleStiffness` (`:181-191`) uses `MeltDCoef` reading phi from GridFunction - φ-dependent D genuinely assembled. But the only MeltDiffusion test (lines 677-678) asserts `getDiffusivity(1.0) > getDiffusivity(0.0)` - tests the math formula, never the FEM path. |
| 12 | Diamond A-B neighbor stencil | 🟡 Partial (honestly marked) | `KmcAtomisticEngine.hpp:177-183` uses 4 body diagonals on cubic grid; `KmcLattice.hpp:29-31` places sites at `i*a0` (simple cubic), not 2×FCC offsets. Ledger honestly notes "not full 2×FCC coords" but marks ✅. |
| 13 | KMC amorphous pocket | 🟡 Partial | `KmcAmorphousPocket::implant` (`:326-345`) marks a cubic region as amorphous code 7. Static helper, NOT a KMC event with Arrhenius rate. |
| 14 | KMC IDW deatomize | 🟡 Partial | `deatomizeIDW` (`:292-320`) does IDW smoothing onto a 1D depth-bin array, NOT an L2 projection onto an MFEM GridFunction. Plan cited `MultiAppProjectionTransfer::assembleL2`. Test only asserts `idw.size() == 8`; plan requires dose conservation to 0.1%. |
| 15 | KMC epitaxy (planar/coord/twin/surf-seg) | ❌ **False closure** | `KmcEpitaxy.hpp:25-117`: `planarGrow`/`coordinationGrow`/`surfaceSegregate`/`formTwin` are deterministic per-pass sweeps; `formTwin` uses `(i+j+k)%7==0` - a modulo, not an Arrhenius rate. Phase 8 plan requires stochastic KMC events. |
| 16 | PDE Reaction term + ICs | ✅ Confirmed | `buildModels` handles `ReactionPdeTerm` via `LinearReactionDiffusion` (`PdeApi.hpp:120-124`); `applyICs` (lines 175-179) calls `engine.initializeSpecies`. Test asserts `d1 < d0` with strong decay. |
| 17 | Jacobian strategy (b) | ✅ Confirmed | See Strengths. |
| 18 | Sub-cycling | ✅ Confirmed | See Strengths. |
| 19 | 3D LevelSetToMesh Cartesian | 🟡 | Same as #1. |
| 20 | PDE reaction + IC | ✅ Confirmed | Same as #16. |
| 21 | Poly GB segregation FEM residual | ❌ | Same as #5. |

**Score: 5 confirmed (3, 16, 17, 18, 20), 9 partial, 4 false closures (4, 5, 9, 15) plus 10/21 with code-correct-but-untested-FEM-path.**

---

## Issues

### Critical (Must Fix)

**C1. III-V `IIIVDiffusion` ships a stub Arrhenius formula but the ledger claims "III-V FEM + I/V eq ✅".**
- **File:** `include/viennaps/fields/models/IIIVDiffusion.hpp:40-45`
- **Issue:** `getDiffusivity` returns `D0·exp(-Ea/kT)` - no `(n/ni)` or `(p/ni)` carrier dependence, no Ga-vs-As sublattice distinction. Phase 6 plan (lines 295-330) requires eq. 3-239/3-240. `C_I_eq`/`C_V_eq` (lines 50-58) are higher-Ea Arrhenius proxies, not a 4-sublattice equilibrium. Plan line 383 explicitly defers I/V eq to 🟡, yet ledger line 46 marks ✅.
- **Why it matters:** "Stub formula behind real-looking assembly" - silent physics absence. A user reading "III-V ✅" expects ATHENA eq. 3-239/3-240; they get a generic Arrhenius that cannot reproduce concentration-dependent III-V diffusion.
- **Fix:** Either (a) implement eq. 3-239/3-240 with `D_AV`/`D_AI`/`D_AV²`/`D_AI²` prefactors and `IntrinsicCarrier` for n/ni, p/ni; or (b) downgrade the ledger row to 🟡 with note "Arrhenius D-form only; eq. 3-239/3-240 carrier-dependent D and 4-sublattice I/V eq deferred (per plan line 383)".

**C2. Poly GB segregation uses a domain integrator, not the cited two-sided interior-face residual.**
- **File:** `include/viennaps/fields/models/PolysiliconDiffusion.hpp:115-139`
- **Issue:** `assembleReaction` adds `DomainLFIntegrator(*segCoefs_.back())` (line 138), integrating `kf·C_gi − kb·C_gb` over the WHOLE domain. Phase 5 plan (lines 183-196) explicitly requires a two-sided interior-face residual at `Poly_GI:Poly_GB` faces, mirroring MOOSE `InterfaceReaction` and reusing Phase 2 `Segregation`'s face-integrator pattern. The domain integrator applies the exchange flux everywhere, not just at the GB interface - physically wrong (couples GI and GB in the grain interior too). `enableGbSegregationSpecies`/`setFemSegregationRate` never called in any test.
- **Why it matters:** Silently produces a different PDE than cited. The "GB FEM segregation residual ✅" claim is false.
- **Fix:** Either (a) replace with `FaceIntegrator` (MFEM `TraceDomainLFIntegrator` or custom interior-face integrator) scoped to `Poly_GI:Poly_GB` faces via attribute marker, following `Segregation.hpp`'s pattern; or (b) downgrade ledger to 🟡 "domain-integrated exchange proxy, not face-local residual; FEM GB path untested".

**C3. KMC epitaxy methods are deterministic sweeps, not KMC events with Arrhenius rates.**
- **File:** `include/viennaps/fields/kmc/KmcEpitaxy.hpp:25-117`
- **Issue:** `planarGrow`, `coordinationGrow`, `surfaceSegregate`, `formTwin` deposit/mark every qualifying site deterministically. `formTwin` uses `(i+j+k)%7==0` - a modulo, not a stochastic rate. Phase 8 plan requires KMC events with Arrhenius rates `ν₀·exp(-E_m/kT)`.
- **Why it matters:** Ledger marks "Epitaxy planar/coord/twin/surf-seg ✅". These are not KMC at all - deterministic geometry routines. "Letter not spirit" anti-pattern.
- **Fix:** Either (a) implement as stochastic KMC events integrated into `KmcAtomisticEngine`'s BKL loop; or (b) downgrade ledger to 🟡 "deterministic geometric deposition helpers, not Arrhenius KMC events".

### Important (Should Fix)

**I1. Runtime AMR hook is a no-op (mark-count only, no mesh refinement).**
- **File:** `include/viennaps/fields/DiffusionEngine.hpp:927-955`
- **Issue:** `applyRuntimeAmr` marks elements, increments counters, but never calls `refineMarked`/`GeneralRefinement`. Comment (lines 943-945) admits intentional due to MFEM serial-H1 corruption concerns. Test asserts `runtimeAmrRefineCount() > 0` - passes because hook runs, no refinement occurs.
- **Fix:** Downgrade ledger to 🟡: "Mark hook + offline refine; runtime GeneralRefinement deferred (serial H1 corruption risk per comment)".

**I2. Five newly-added FEM model classes never instantiated in any test.**
- **Files:** `SiGeCDiffusion.hpp` (FEM path), `SiGeCDiffusion.hpp:43-115` (`GeBPairingModel`), `SiGeCDiffusion.hpp:118-165` (`StrainDiffusionModel`), `FlashLaserAnneal.hpp:141-207` (`MeltDiffusion` FEM `assembleStiffness`), `PolysiliconDiffusion.hpp:115-139` (`assembleReaction` FEM path).
- **Issue:** Each has real `assembleStiffness`/`assembleReaction`, but tests use the OLD 1D `applyStep`/`applySegregationStep`/`modifyD`/`getDiffusivity(phi)` paths instead. `enableGbSegregationSpecies` is never called. A regression breaking any FEM path would not be caught.
- **Why it matters:** Repeats the exact "test bypasses the physics" pattern the prior reviews flagged. Ledger marks these rows ✅.
- **Fix:** For each, add an engine-solve test that registers the FEM model, runs a short solve, asserts a physics-sensitive outcome.

**I3. 3D LevelSetToMesh has no test.**
- **File:** `include/viennaps/fields/LevelSetToMesh.hpp:178-270`; no `TestLevelSetToMesh3D`.
- **Issue:** 3D specialization genuinely implemented, never exercised. A bug in `attributeAtPoint3` (3D index clamping at lines 259-261) would not be caught.
- **Fix:** Add `TestLevelSetToMesh3D` mirroring `TestLevelSetToMesh2D`.

**I4. KMC IDW deatomize lacks the L2 projection half and the plan's acceptance tests.**
- **File:** `include/viennaps/fields/kmc/KmcAtomisticEngine.hpp:292-320`
- **Issue:** `deatomizeIDW` does IDW smoothing onto a 1D depth-bin array. Plan (line 150) cites "IDW smoothing + L2 projection" (MOOSE `MultiAppGeometricInterpolationTransfer` + `MultiAppProjectionTransfer`). L2 Mx=b projection onto MFEM GridFunction is missing. Test (lines 948-949) asserts only `idw.size() == 8`; plan line 149 requires dose conservation (0.1%) and peak-integral checks.
- **Fix:** Either add L2 projection path onto MFEM GridFunction (reuse `SolutionTransfer::transferL2`), or downgrade to 🟡 "IDW smoothing only; L2 projection deferred". Add dose-conservation assertion.

**I5. KMC event rebuild still O(N) per step.**
- **File:** `include/viennaps/fields/kmc/KmcAtomisticEngine.hpp:48-162`
- **Issue:** `step()` rebuilds entire `events_` vector every call (O(N_sites · neighbors)), then prefix-sum + binary search for O(log N) selection. Ledger (line 56) honestly notes "Rebuild still O(N) - incremental heap stretch" - so honestly marked, not overclaimed.
- **Fix:** No action needed for this batch; the 🟡 mark is honest.

### Minor

- **M1.** `PdeEquation::applyTo` does not call `applyICs` (`PdeApi.hpp:166-172`). User must remember to call `eq.applyICs(engine)` separately. Consider having `applyTo` accept an optional engine pointer, or document the two-call contract prominently.
- **M2.** `HeatTransfer` omits the latent-heat coupling `ρ·L·∂φ/∂t` and `ρ·c_P` mass scaling from eq. 213 (`FlashLaserAnneal.hpp:45-66`). Basic `α∇T·∇v + G` form is correct but coupled phase-field feedback is absent. Worth a 🟡 note rather than implying full eq. 213.
- **M3.** `KmcAmorphousPocket::implant` is a static initializer, not a KMC event (`KmcAtomisticEngine.hpp:323-345`). Either integrate as a KMC event type or document as a static damage initializer.
- **M4.** Ledger row "AdaptiveMeshRefiner mark/refine/ZZ/threshold/derefine API ✅" notes "NC derefine no-op on serial Cartesian" (line 34). `tryDerefine` (line 200-206) returns `0 * ne0` - explicit no-op, honest. `zzIndicator` (line 159-179) is a gradient-magnitude proxy, not the true Zienkiewicz-Zhu flux-recovery estimator (no smoothed-gradient patch projection). Name slightly overclaimed.

---

## Recommendations

**Rework required (false closures - choose implement OR downgrade):**
- **C1 (III-V):** Implement eq. 3-239/3-240 with carrier-dependent D on Ga/As sublattices, OR downgrade ledger row to 🟡.
- **C2 (Poly GB):** Replace domain integrator with face-local two-sided residual at `Poly_GI:Poly_GB` faces (reuse `Segregation.hpp`'s pattern), OR downgrade to 🟡.
- **C3 (KMC epitaxy):** Implement stochastic Arrhenius KMC events, OR downgrade to 🟡 "deterministic geometric helpers".

**Add tests (half-closures):**
- **I2:** Add engine-solve tests for `GeBPairingModel`, `StrainDiffusionModel`, `MeltDiffusion` FEM, `SiGeCDiffusion` FEM, and `PolysiliconDiffusion` FEM GB path. Each should assert a physics-sensitive outcome, not just "runs".
- **I3:** Add `TestLevelSetToMesh3D`.
- **I4:** Add dose-conservation assertion to the IDW deatomize test.

**Adjust ledger wording (honest limits):**
- **I1:** Runtime AMR -> 🟡 "mark hook + offline refine; runtime refine deferred".
- **M2:** Flash heat -> 🟡 "basic α∇T·∇v + G; latent-heat φ-coupling deferred".

**Merge-ready as-is:**
- Rows 3 (L2 transfer), 16/20 (PDE reaction+IC), 17 (Jacobian), 18 (sub-cycling), and the `367629e` bug fixes (drift sign, SiGe `&&`) are genuine closure with adequate tests.
- Rows 6 (SiGeC code), 7 (GeB code), 8 (strain code), 11 (heat FEM) are code-correct but need the I2 tests added before they can be called verified.

---

## Assessment

**Do the 3 commits achieve the claimed parity?** **Partially.** The `367629e` bug fixes and several engine/PDE rows (3, 16, 17, 18, 20) are genuine closure. But three rows are false closures (III-V I/V eq with a stub Arrhenius, Poly GB segregation with a domain integrator instead of a face residual, KMC epitaxy as deterministic sweeps), one row is a false closure (runtime AMR no-op), and five code-correct FEM models ship without any test exercising their FEM path - repeating the exact "test bypasses the physics" pattern the prior reviews flagged.

**Ready to merge?** **With fixes.**

**Reasoning:** The bug-fix commit (`367629e`) and the genuinely-wired rows are merge-ready, but the ledger's ✅ marks on III-V (C1), Poly GB (C2), KMC epitaxy (C3), and runtime AMR (I1) overstate what the code does. Either implement the cited physics or downgrade those four rows to 🟡 with honest notes, and add engine-solve tests for the five untested FEM models (I2) before calling them verified. The codebase has improved materially since `ac18975`, but the recurring "✅ marks taken at face value" failure pattern persists in this batch.

---

## Relevant file paths (all absolute)

- `D:\dev\ViennaPS_mod\include\viennaps\fields\models\MobileImpurity.hpp` (drift sign fix ✅)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\models\SiGeDiffusion.hpp` (`&&` fix ✅)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\models\LinearReactionDiffusion.hpp` (PDE reaction ✅)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\models\FermiDiffusion.hpp` (Jacobian ✅)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\SolutionTransfer.hpp` (L2 ✅)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\DiffusionEngine.hpp` (sub-cycles ✅; runtime AMR ❌ I1)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\models\IIIVDiffusion.hpp` (**C1 false closure**)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\models\PolysiliconDiffusion.hpp` (**C2 false closure**)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\kmc\KmcEpitaxy.hpp` (**C3 false closure**)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\models\SiGeCDiffusion.hpp` (I2 untested FEM)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\models\FlashLaserAnneal.hpp` (heat ✅; MeltDiffusion I2)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\LevelSetToMesh.hpp` (3D code 🟡 I3)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\kmc\KmcAtomisticEngine.hpp` (IDW I4; rebuild I5)
- `D:\dev\ViennaPS_mod\tests\diffusion\testDiffusion.cpp` (I2 test gaps)
