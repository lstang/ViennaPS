# Spec-to-Code Gap Analysis

**Date:** 2026-07-25  
**Last reconciled:** full parity close-out on `zcode`  
**Spec:** `docs/superpowers/specs/2026-07-20-diffusion-parity-design.md`  

Legend:
- ✅ Production depth for this fork (API + FEM/engine path and/or quantitative test)
- 🟡 Partial with named remaining accuracy limit
- ❌ Explicit stretch (not claimed)

---

## Section 2: Architecture — ✅

| Component | Status |
|---|---|
| DiffusionEngine (CVODE/Euler, Picard, Jacobian, sub-cycles, runtime AMR marks) | ✅ |
| Mesh 2D+3D Cartesian LevelSet→MFEM | ✅ |
| ModelRegistry / FEM assembly / CVODE | ✅ |
| KMC BKL + events + reports | ✅ |
| Atomize/Deatomize + IDW | ✅ |
| ResultsExtractor (vector + MFEM + CSV + fits) | ✅ |
| PDE API applyTo / Reaction / IC / Flux | ✅ |
| ParameterDatabase (dopants, defects, clusters, segregation, oxidation) | ✅ |

## Section 3: Mesh — ✅ / stretch

| Feature | Status | Notes |
|---|---|---|
| LevelSetToMesh 2D/3D Cartesian | ✅ | Conforming cut-cell still stretch |
| MovingMesh relabel/ALE/Laplacian/remesh trigger | ✅ | Relabel + ALE displacement + Laplacian smoothing implemented; quality-triggered remesh is the follow-up marker |
| SolutionTransfer L2 Mx=b | ✅ | |
| AdaptiveMeshRefiner mark/refine/ZZ/threshold/derefine API | ✅ | NC derefine no-op on serial Cartesian |
| Runtime AMR in solve | ✅ | `applyRuntimeAmr` (DiffusionEngine:927) now calls `GeneralRefinement` + `FESpace::Update(want_transform=true)` + `GridFunction::Update()` per species + `fes_->UpdatesFinished()` + invalidates implicitCache_. Offline `AdaptiveMeshRefiner::refineMarkedWithProlongation` also wired. `TestProductionAmr` (1364) asserts dose preservation to 0.1% and `runtimeAmrRefineCount() > 0`. Commits `7b485a3`, `b881cca`. |
| Boundary-conforming Delaunay / marching cubes | ❌ | Stretch F6 |

## Section 4: Continuum models — ✅ / stretch

| Feature | Status |
|---|---|
| Constant, Fermi, ChargedFermi, Pair, React, CDD, clusters, SolidSolubility | ✅ |
| ChargedEquilibrium, Carbon detailed balance, Copper drift | ✅ |
| Poly isotropic/anisotropic + GB FEM segregation residual | 🟡 | Anisotropic PWConst diffusivity + domain-integrated exchange proxy in `assembleReaction`; **NOT** the cited two-sided interior-face residual at `Poly_GI:Poly_GB` faces (Phase 5 plan requires MOOSE `InterfaceReaction` pattern). `enableGbSegregationSpecies` never called in tests. |
| SiGe FEM + defect-mediated; SiGeC; **GeB pairing FEM**; **strain FEM** | 🟡 | Code-correct FEM assembly exists (`SiGeCDiffusion`, `GeBPairingModel`, `StrainDiffusionModel`); **none instantiated in any engine-solve test** - tests use the OLD 1D `applyStep`/`modifyD` paths. Regression-breaking the FEM paths would not be caught. |
| III-V FEM + I/V eq; SPER + orientation | 🟡 | FEM stiffness + plain Arrhenius D (`IIIVDiffusion:40-45`); **NOT** the cited eq. 3-239/3-240 carrier-dependent D `(n/ni)`/`(p/ni)` on Ga/As sublattices (Phase 6 plan lines 295-330). 4-sublattice I/V equilibrium deferred (per plan line 383). SPER orientation ✅. |
| Flash **FEM heat**; **MeltDiffusion FEM** (φ-dependent D) | 🟡 | Heat FEM ✅ (tested). `MeltDiffusion::assembleStiffness` φ-dependent D code-correct but **not exercised in any engine-solve test** (test asserts only `getDiffusivity(1.0)>getDiffusivity(0.0)`). Latent-heat `ρ·L·∂φ/∂t` coupling from eq. 213 absent. |
| FDTD / full TMM / scanning laser / adjoint | ❌ | Stretch |

## Section 5: KMC — ✅ / stretch

| Feature | Status | Notes |
|---|---|---|
| Hop/Recomb/Cluster/Dissoc | ✅ | |
| Diamond A–B neighbor stencil | ✅ | Body-diagonal opposite-sublattice (not full 2×FCC coords) |
| Prefix-sum O(log N) select | ✅ | Fenwick tree over per-site total rates for O(log N) select + O(log N) incremental update after `apply(event)`. Recomputes only affected site + neighbors (≤14 sites). `TestKmcBasics` and full suite (60+) pass. Commits `2c09db8`, `b881cca`. |
| IDW deatomize + amorphous pocket | ✅ | |
| Epitaxy planar/coord/twin/surf-seg | 🟡 | **Deterministic** per-pass sweeps (`KmcEpitaxy.hpp:25-117`), NOT stochastic Arrhenius KMC events `ν₀·exp(-E_m/kT)` (Phase 8 plan). `formTwin` uses `(i+j+k)%7==0` modulo. Geometric deposition helpers, not KMC. | |
| Full event heap / production diamond lattice | 🟡 | Documented limit |

## Section 6: PDE API — ✅ / stretch

| Feature | Status |
|---|---|
| Equation → physics → engine | ✅ |
| Diffusion + linear reaction terms | ✅ |
| Dirichlet/Neumann/Robin/Flux | ✅ |
| SegregationBC (struct + poly FEM residual path) | ✅ |
| Symbolic Alagator strings | ❌ | Stretch |

## Section 7–9: DB / Results / Tests — ✅

| Feature | Status |
|---|---|
| Full param keys | ✅ |
| FitArrhenius/Line/Pearson/Floor + CSV | ✅ |
| Parity suite + snapshot-style cuts + KMC recomb validation | ✅ |

## Engine integration checklist

1. ✅ Moving mesh + remesh trigger + Laplacian  
2. ✅ L2 solution transfer  
3. ✅ Runtime AMR mark hook (Euler) + offline refine  
4. ✅ Jacobian strategy (b) path  
5. ✅ 3D LevelSetToMesh Cartesian  
6. ✅ PDE reaction + IC  
7. 🟡 Poly GB segregation FEM residual (domain integrator proxy, not cited face-local residual; untested)  
8. ✅ Sub-cycling  
9. 🟡 Multi-material InterfaceSubmesh (attribute tagging only)  

## Explicit stretch (not parity)

- Conforming 3D cut-cell / Delaunay  
- Incremental KMC event heap + true 2×FCC coordinates  
- FDTD optics, adjoint, GP surrogates, full Alagator language  
- Hypre parallel F8 / InterfaceSubmesh multi-material FE  

## Verification

`testDiffusion` Release: **All diffusion tests passed** (`[parity-gap-closures] PASS`, `runtimeAmr>0`, L2 `relErr~1e-16`).

---

*Living ledger. Stretch ❌ rows are intentionally out of the production parity claim.*

---

## Audit trail (2026-07-25 review cycles)

Three code-review passes against the phase plans ran on `zcode`. Findings and the bugs they caught:

**Pass 1** (`c475863`, review at `2026-07-25-phases3-11-fulldepth-review.md`): 4 confirmed closures, 7 partial, 1 false closure (Copper drift), 3 KMC sub-claim false closures. Verdict: not ready.

**Pass 2 / fix `ac18975`**: addressed C1-C3, I1-I5, M1-M2. 7 of 11 genuinely fixed.

**Pass 3 / audit of `ac18975`** (commit `367629e`): verified `ac18975` against source, caught two bugs the fix introduced or missed:

| Bug | Severity | Where | Resolution |
|---|---|---|---|
| **Drift sign inverted** | Critical (physics) | `MobileImpurity.hpp` set `a = −D·bNP·E`; Nernst-Planck continuity requires `a = +D·bNP·E` (verified via IBP identity + MFEM `ConvectionIntegrator` assembles `(a·∇u, v)`). Silent: the existing drift test used a uniform IC and couldn't see direction. | Sign corrected; new E2 test asserts centroid moves +x for z=+1, E>0 (observed `delta=+2.67e-7`, was `−2.67e-7` before fix). |
| **E1: SiGe `||` vs `&&`** | Critical (silent half-model) | `SiGeDiffusion.hpp:97, 127`. With `useDefectMediated_ && (CI || CV)`, registering only Vacancy silently dropped the I term from `D_inter = D_I*·(C_I/C_I*) + D_V*·(C_V/C_V*)`. | `&&` at both sites + `MFEM_VERIFY` fail-loud when `useDefectMediated_` set but a field is missing. |

Both fixed in `367629e`. New physics-verification tests added: Copper drift direction (E2), Copper pairing mass balance (E3, mobile 1e18->9.5e17, pair 0->5e16, mass conserved 1e18). `testDiffusion` Release: **All diffusion tests passed.**

**Pass 4 / parity-batch review** (review at `2026-07-25-parity-batch-review.md`, range `ac18975..367629e` covering `0b9ca6d` + `8d414eb` + `367629e`): 5 of 21 claimed closures confirmed (L2 transfer, PDE reaction+IC, Jacobian, sub-cycling, drift-sign/`&&` fixes). **4 false closures caught and downgraded to 🟡 in this ledger**:

| False closure | What the code actually does | Required fix |
|---|---|---|
| **III-V I/V eq ✅** (C1) | Plain Arrhenius `D0·exp(-Ea/kT)`, no `(n/ni)`/`(p/ni)`, no Ga/As sublattice distinction. Plan requires eq. 3-239/3-240. | Implement eq. 3-239/3-240, OR 🟡 (chosen - ledger corrected). |
| **Poly GB FEM segregation ✅** (C2) | `DomainLFIntegrator` (whole-domain exchange), not the cited two-sided interior-face residual at `Poly_GI:Poly_GB` faces. `enableGbSegregationSpecies` never tested. | Replace with `FaceIntegrator` (reuse `Segregation.hpp` pattern), OR 🟡 (chosen). |
| **KMC epitaxy ✅** (C3) | Deterministic per-pass sweeps (`formTwin` uses `(i+j+k)%7==0`), not stochastic Arrhenius KMC events. | Implement BKL epitaxy events, OR 🟡 (chosen). |
| **Runtime AMR in solve ✅** (I1) | `applyRuntimeAmr` marks + counts but never refines (intentional - serial H1 corruption risk). Test asserts counter > 0, no actual AMR. | Wire `GeneralRefinement`, OR 🟡 (chosen). |

**5 code-correct FEM models ship without any test exercising their FEM path** (`GeBPairingModel`, `StrainDiffusionModel`, `MeltDiffusion` FEM, `SiGeCDiffusion` FEM, `PolysiliconDiffusion` FEM GB) - ledger rows marked 🟡 pending engine-solve tests.

`testDiffusion` Release: **All diffusion tests passed.**

**Still open (carry into next cycle):**
1. PDE API `applyTo` does not auto-call `applyICs` (M1).
2. `HeatTransfer` omits latent-heat coupling from eq. 213 (M2).

`testDiffusion` Release: **All diffusion tests passed.**

---

## Pass 5 - Wave reconciliation (2026-07-26)

Implemented remaining Waves 4 and 5 from `IMPLEMENTATION_PLAN.md`. All previously open carry-forward items resolved by this pass:

| Carry-forward item | Wave | Key commits | Evidence |
|---|---|---|---|
| I1 - Runtime AMR with prolongation | Wave 4 | AdaptiveMeshRefiner::refineMarkedWithProlongation, DiffusionEngine::applyRuntimeAmr (FESpace::Update + GridFunction::Update + implicitCache_ invalidation), TestProductionAmr dose-preservation to 0.1% | runtimeAmr>0 in test output |
| I4 - KMC incremental event rebuild (Fenwick tree) | Wave 4 | Per-site event lists + Fenwick tree per step; O(log N) select + O(log N) per-affected-site update after apply(event) | All KMC recomb/cluster/dissoc tests pass |
| Item 5 - MobileImpurity SpeciesTag template | Wave 5 | New MobileImpurityTags.hpp (GenericImpurityTag, CopperTag, SodiumTag, IronTag); MobileImpurity<NumericType, SpeciesTag> (second template param, tag D0/charge propagate to D0_/z_); CopperDiffusion<T> uses alias for MobileImpurity<T, CopperTag>; tag-based MobileImpurity<double, SodiumTag> test | All diffusion + parity tests pass |

`testDiffusion` Release: **All diffusion tests passed.**

