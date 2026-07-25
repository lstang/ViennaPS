# Spec-to-Code Gap Analysis

**Date:** 2026-07-25 (last verified 2026-07-25 against `ac18975`)
**Spec:** `docs/superpowers/specs/2026-07-20-diffusion-parity-design.md`
**Head:** `ac18975` (zcode branch); was `b23d3a4` at first audit, `c475863` at first review
**Review log:** `docs/superpowers/reviews/2026-07-25-phases3-11-fulldepth-review.md` (against `c475863`)

Legend:
- ✅ Implemented - FEM-assembled, tested, production-ready
- 🟡 Skeleton - Header exists with math/API stubs; not FEM-integrated or 1D-only
- ❌ Missing - No file, empty stub, or spec requirement not addressed

---

## Section 2: Architecture

| Spec Component | Status | Notes |
|---|---|---|
| DiffusionEngine | ✅ | CVODE BDF + implicit-Euler, multi-species, BCs, caching |
| MeshModule (LevelSet->MFEM) | ✅ | 2D only; 3D stub throws |
| ModelRegistry | ✅ | DiffusionPhysics with MOOSE gatekeeper |
| FEMAssembly | ✅ | Species-outer, terms-inner loop |
| TimeIntegrator | ✅ | SUNDIALS CVODE (BDF, adaptive) |
| KmcEngine | 🟡 | BKL loop works; no diamond lattice, no event tree |
| DeatomizeTransfer | 🟡 | KmcAtomize/KmcDeatomize exist; no FEM projection |
| ResultsExtractor | 🟡 | 1D vector-only; no MFEM GridFunction integration |
| PdeAPI | 🟡 | PdeEquation/PdeTerm exist; not consumed by engine |
| ParameterDatabase | 🟡 | Arrhenius D0/Ea only; no defect/cluster/segregation params |

## Section 3: Mesh Generation

| Spec Component | Status | Notes |
|---|---|---|
| LevelSetToMeshConverter (2D) | ✅ | Cartesian + attribute tagging |
| LevelSetToMeshConverter (3D) | ❌ | Stub throws runtime_error |
| Boundary extraction (marching squares/cubes) | ❌ | Uses Cartesian, not boundary-conforming |
| Triangulation (Delaunay/constrained) | ❌ | Cartesian only |
| Interface-aligned refinement | ❌ | Deferred per F6 follow-up |
| User refinement boxes | 🟡 | RefinementBox + markBox exist; not engine-integrated |
| MovingMeshHandler | 🟡 | Exists (`c475863`). Relabel (idiom A) + ALE displacement (idiom C). **Remesh trigger (skewness/Jacobian), Laplacian smoothing, 3D moving boundary still ❌.** |
| Laplacian smoothing | ❌ | |
| Remesh trigger (skewness/Jacobian) | ❌ | |
| 3D moving boundary | ❌ | |
| SolutionTransfer | 🟡 | Exists (`c475863`). Dose-preserving cross-mesh transfer (0.1% conserved). Uses point-sampling + FindPoints, not the cited `Mx=b` L2 projection. |
| L2 projection (cross-mesh) | 🟡 | Approximated via `ProjectCoefficient` + dose rescale; not the cited mass-matrix solve |
| Boundary-safe interpolation | 🟡 | FindPoints returns -1 for unmapped QPs → counted in `unmappedQuadraturePoints`, optional `requireMapped_` fail-loud |
| AdaptiveMeshRefiner (static) | 🟡 | markBox + refineMarked via MFEM GeneralRefinement |
| AdaptiveMeshRefiner (runtime) | ❌ | Not integrated with engine solve loop |
| ZZ error estimator | ❌ | |
| ThresholdRefiner | ❌ | |
| Derefine | ❌ | |
| hp-AMR during CVODE | ❌ | |

## Section 4: Continuum Diffusion Models

### 4.2 Transport Models

| Model | FEM Assembly | Tested | Notes |
|---|---|---|---|
| Constant | ✅ | ✅ | Arrhenius D, dose conservation 1e-16 |
| Fermi | ✅ | ✅ | D(C) coefficient, analytic dD/dC |
| ChargedFermi | ✅ | ✅ | QP-local D (2026-07-25 full-depth + review) |
| Pair | ✅ | 🟡 | FEM assembled via CddDiffusion / PairDiffusion |
| ChargedPair | ✅ | 🟡 | FEM QP D with Fermi enhancement |
| React | ✅ | 🟡 | FEM assembled, recombination sink |
| ChargedReact | ✅ | 🟡 | FEM charge-enhanced product residual |
| NeutralReact | ✅ | 🟡 | Thin alias of ReactDiffusion |
| CDD | ✅ | 🟡 | KernelTerm composition; UAF fixed in review |
| ChargedEquilibrium | ✅ | 🟡 | Fermi charge-state partition `ΣD^z·exp(−zη)/Σexp(−zη)` (`ChargedEquilibriumDiffusion.hpp:53-74`, was wrong formula, fixed in `ac18975`). Test only checks dose conservation, not `D_eff` shift with `n/ni`. |
| Carbon | ✅ | 🟡 | FEM + reverse rate detailed balance `kr=kf·C*_I` (`CarbonDiffusion.hpp:41-48, 136`). Test sets `setReverseRate(0.0)` — bypasses the equilibrium it should verify. |
| Nitrogen | ✅ | 🟡 | ConstantDiffusion specialization |
| Copper | ✅ | 🟡 | Drift via `ConvectionIntegrator` + pairing in `MobileImpurity` base (`MobileImpurity.hpp:130-150, 155-184`); `CopperDiffusion` is thin specialization (DRY fixed). Test sets E-field but only asserts dose conservation, not drift-direction asymmetry; pairing still disabled in test (`setIonPairing(0.0)`). |

### 4.3 Cluster / Deactivation Models

| Model | FEM Assembly | Tested | Notes |
|---|---|---|---|
| {311} clusters | 🟡 | 🟡 | KernelTerm; no FEM stiffness |
| VacancyCluster | 🟡 | ❌ | KernelTerm; no FEM stiffness |
| ImpurityCluster (BIC) | 🟡 | ❌ | KernelTerm; no FEM stiffness |
| DislocationLoop | 🟡 | ❌ | Header exists, no FEM assembly |
| SolidSolubility | ✅ | 🟡 | FEM-assembled deactivation |

### 4.4 Interface Physics

| Model | FEM Assembly | Tested | Notes |
|---|---|---|---|
| Segregation | ✅ | ✅ | Two-sided residual + operator-split |
| OED | 🟡 | ❌ | OedSource header; no engine integration |
| TED | 🟡 | 🟡 | TedInitializer; timeline utility only |
| DoseLoss | ✅ | 🟡 | Robin BC in engine |
| Grain boundary segregation | 🟡 | ❌ | In PolysiliconDiffusion; not FEM-assembled |

### 4.5 Polysilicon Diffusion

| Spec Feature | Status | Notes |
|---|---|---|
| Isotropic D_eff = D_bulk + D_gb*f_gb | ✅ | FEM ConstantCoefficient |
| Grain growth dR/dt = k*exp(-Ea/kT)/R^n | ✅ | GrainModel.advance() |
| Anisotropic (dual mesh) | 🟡 | PWConstCoefficient by attr; no dual mesh |
| Voronoi tessellation | ❌ | |
| Interface oxide breakup | 🟡 | PolyOxideBreakup; simplified |
| Epitaxial regrowth | 🟡 | PolyOxideBreakup; simplified |
| Segregation interior-boundary | 🟡 | applySegregationStep; 1D, not FEM |

### 4.6 SiGe/SiGeC Diffusion

| Spec Feature | Status | Notes |
|---|---|---|
| SiGe interdiffusion | ✅ | 🟡 | FEM + optional defect-mediated `D_inter=D_I*·(C_I/C_I*)+D_V*·(C_V/C_V*)` (`SiGeDiffusion.hpp:48-56, 114-136`, `ac18975`); Arrhenius fallback when defect fields absent. Test exercises only the fallback path (no `setDefectMediated`/C_I,C_V species) — defect-mediated branch untested. |
| Bandgap model | ✅ | BandgapModel.hpp with niRatioToSi |
| Boron D modified by Ge | ✅ | Via niRatioToSi |
| Carbon suppression (I trapping) | 🟡 | SiGeCDiffusion::tedFactor; no FEM |
| Strain effects | 🟡 | StrainDiffusionModifier; no FEM |
| Ge-B pairing | 🟡 | GeBPairing; 1D explicit, no FEM |
| Cluster initialization | ❌ | |

### 4.7 III-V Compound Semiconductor

| Spec Feature | Status | Notes |
|---|---|---|
| Material conversion | 🟡 | IIIVDiffusion header; stub |
| Species-specific D | 🟡 | Stub |
| I/V equilibrium in compound | ❌ | |

### 4.8 SPER

| Spec Feature | Status | Notes |
|---|---|---|
| Interface velocity v = v0*exp(-Ea/kT) | ✅ | SPERKernel; 1D profile |
| Orientation-dependent | ❌ | |
| Dopant release | 🟡 | Simplified (1% of regrown) |
| Defect emission (EOR) | 🟡 | 10% of regrown |
| SPER front as level-set | ❌ | Not integrated with level-set |
| FEM mesh update as front advances | ❌ | |

### 4.9 Flash/Laser Anneal

| Spec Feature | Status | Notes |
|---|---|---|
| Heat transfer equation | 🟡 | 1D explicit FTCS; **not FEM** |
| Energy implantation | ❌ | |
| Transfer matrix method (optical) | ❌ | |
| FDTD for sub-wavelength | ❌ | |
| Melting phase field | 🟡 | Simplified relaxation, not Allen-Cahn FEM |
| Crystallinity phase field | ❌ | |
| Dopant diffusion in melt | 🟡 | MeltDiffusion; no FEM |
| Solute transport in liquid | ❌ | |
| Resolidification trapping | ❌ | |
| Gaussian intensity (flash) | ✅ | LaserIntensity Beer's law |
| Table lookup | ❌ | |
| User-specified (PDE API) | ❌ | |
| Scanning laser | ❌ | |
| FlashLaserAnneal orchestrator | 🟡 | runPulse; 1D only, not FEM |

### 4.10 FEM Assembly Pattern

| Spec Feature | Status | Notes |
|---|---|---|
| DiffusionModel abstract base | ✅ | |
| assembleStiffness | ✅ | |
| assembleReaction | ✅ | |
| assembleMass | ✅ | |
| assembleStiffnessJacobian (strategy b) | 🟡 | Virtual exists; **not called by engine** |
| numSpecies / speciesNames | ✅ | |
| applicableAttributes | ✅ | |

## Section 5: KMC Engine

### 5.1 Atomistic KMC

| Spec Feature | Status | Notes |
|---|---|---|
| BKL rejection-free KMC | ✅ | KmcAtomisticEngine |
| Si diamond cubic lattice | 🟡 | Still cubic + 4 body diagonals (`KmcAtomisticEngine.hpp:180-188`), not true 2-interpenetrating FCC sublattices |
| Event tree O(log N) | 🟡 | Selection is `lower_bound` on prefix sum = O(log N) (`KmcAtomisticEngine.hpp:152-158`), but `events_` vector is fully rebuilt every step = O(N·neighbors) rebuild (lines 48-132). Net per-step cost still O(N). True incremental heap (partial rebuild) is a follow-up. |
| Hop events | ✅ | |
| Recombination (I+V->0) | ✅ | Tested (`testDiffusion.cpp:560-581`) |
| Clustering (I+I->{311}) | ✅ | Implemented + tested (`KmcAtomisticEngine.hpp:112-130, 212-220`; `testDiffusion.cpp:583-606` asserts `clusterCount()>0`, `countSpecies(3)>=1`) |
| Dissociation | ✅ | Implemented (`KmcAtomisticEngine.hpp:62-78, 221-229`); NOT tested (no dissoc assertion in testDiffusion.cpp) |
| Amorphous pocket | ❌ | |
| Dopant-defect pairing | ❌ | KMC side; continuum Ge-B pairing exists |
| Impurity clustering | ❌ | |
| Deatomize | ✅ | Poisson sampling |
| Atomize | ✅ | |
| Defect activity reports | ❌ | |
| Interaction reports | ❌ | |
| 1D profiles | ❌ | |
| Supersaturation | ❌ | |
| Cluster size histograms | ❌ | |

### 5.2 Lattice KMC Epitaxy

| Spec Feature | Status | Notes |
|---|---|---|
| KmcEpitaxy engine | 🟡 | Header exists; no epitaxy physics |
| Planar epitaxy | ❌ | |
| Coordination-based | ❌ | |
| Coordination-reactions | ❌ | |
| SiGe mole fraction growth | ❌ | |
| Visibility/shadowing | ❌ | |
| Twin-defect formation | ❌ | |
| Surface segregation | ❌ | |
| Nonselective epitaxy | ❌ | |

## Section 6: PDE API (Alagator Equivalent)

| Spec Feature | Status | Notes |
|---|---|---|
| PdeEquation | 🟡 | Container for terms; **not consumed by engine** |
| PdeTerm base | 🟡 | DiffusionPdeTerm, ReactionPdeTerm only |
| PdeBC (Dirichlet/Neumann/Robin) | 🟡 | Struct exists; not wired to engine BCs |
| SegregationBC | ❌ | |
| FluxBC | ❌ | |
| GrowthTerm | ❌ | |
| addEquation / addEquationTerm | ❌ | |
| subEquationTerm | ❌ | |
| Engine integration | ❌ | PdeEquation cannot drive DiffusionEngine |

## Section 7: Parameter Database

| Spec Feature | Status | Notes |
|---|---|---|
| Species D0/Ea per material | ✅ | B, P, As, Sb in Si, SiO2, Si3N4, PolySi, SiGe |
| Point-defect I/V equilibrium | ❌ | No C_I_eq, C_V_eq, formation energies |
| Cluster parameters | ❌ | No {311} binding, loop rate, BIC binding |
| Segregation coefficients m(T) | ❌ | No per-dopant/material-pair data |
| Oxidation parameters | ❌ | No linear/parabolic rate constants |
| Inheritance (SiGe->Si) | ✅ | |
| Like-materials blend | ✅ | blend(matA, matB, key, w) |

## Section 8: Results Extraction

| Spec Feature | Status | Notes |
|---|---|---|
| 1D data cuts | 🟡 | ResultsExtractor::cut1D; 1D vector, not FEM |
| Dose calculation | ✅ | DiffusionEngine::getIntegral (FEM L1) |
| Level crossings | 🟡 | 1D only |
| Sheet resistance | 🟡 | Proxy 1/(mu*dose); no mobility model |
| FitArrhenius | 🟡 | FittingUtilities::fitD0FixedEa |
| FitLine | ❌ | |
| FitPearson | ❌ | |
| FitPearsonFloor | ❌ | |
| VTK output | ❌ | Not in ResultsExtractor |
| TDR-equivalent (HDF5) | ❌ | |
| CSV output | ❌ | |

## Section 9: Testing

| Spec Feature | Status | Notes |
|---|---|---|
| Unit tests (analytical solutions) | ✅ | 34+ tests, dose conservation, BCs |
| Integration tests (multi-step) | ✅ | Fermi+Segregation, TED, BC stack |
| Regression tests (snapshots) | ❌ | No 1D profile snapshot comparison |
| KMC vs continuum validation | ❌ | |
| Mesh convergence | ❌ | |
| CTest entries | ✅ | |

## Section 10: Missing Files

| Spec File | Status |
|---|---|
| MovingMeshHandler.hpp | ✅ | Exists (`c475863`). Subdomain relabel (idiom A) + ALE node displacement (idiom C). Interface-adjacent-only flag honored via instance `relabel()` (`MovingMeshHandler.hpp:37-41`, `ac18975`). Missing: remesh trigger (skewness/Jacobian), Laplacian smoothing, 3D moving boundary. |
| SolutionTransfer.hpp | ✅ | Exists (`c475863`). Dose-preserving cross-mesh transfer via `ProjectCoefficient(SourceSampleCoef)` + FindPoints + rescale (`SolutionTransfer.hpp:63-142`). Deviates from cited `MultiAppProjectionTransfer::assembleL2()` `Mx=b` form — uses point-sampling instead; counts unmapped QPs and exposes `requireMapped_` fail-loud flag (`ac18975`). 0.1% dose conservation enforced. |
| PdeEquation.hpp / PdeTerm.hpp / PdeBC.hpp | 🟡 Consolidated in PdeApi.hpp |
| KmcDeatomize.hpp | 🟡 In KmcAtomisticEngine.hpp |
| ResultsExtractor.hpp | 🟡 In PdeApi.hpp |
| HeatTransfer.hpp | 🟡 In FlashLaserAnneal.hpp |
| MeltingPhaseField.hpp | 🟡 In FlashLaserAnneal.hpp |

## Section 11: Backward Compatibility

| Spec Feature | Status | Notes |
|---|---|---|
| 1D profile kernels preserved | ✅ | DiffusionKernel, FermiDiffusionKernel, etc. |
| ProcessOrchestrator extended | ✅ | Optional DiffusionEngine when MFEM present |
| PhysicsField 1D API remains | ✅ | |

---

## Engine Integration Gaps

The DiffusionEngine itself is missing these spec-required capabilities (status as of `ac18975`):

1. 🟡 **Moving boundary mesh** — `MovingMeshHandler` exists with relabel + ALE displacement, but no remesh trigger / Laplacian smoothing / 3D (spec 3.2)
2. 🟡 **Solution transfer** — dose-preserving cross-mesh transfer exists but uses point-sampling, not the cited `Mx=b` L2 projection (spec 3.3)
3. ❌ **Runtime AMR** — AdaptiveMeshRefiner exists but not called during solve (spec 3.4)
4. ❌ **Nonlinear Jacobian** — assembleStiffnessJacobian virtual exists but engine never calls it (spec 4.10 strategy b)
5. ❌ **3D support** — LevelSetToMesh 3D throws. Engine compiles for D=3 but no mesh source
6. 🟡 **PDE API integration** — `PdeEquation::applyTo` bridges into DiffusionPhysics (`PdeApi.hpp:126`), but `buildModels` silently drops `ReactionPdeTerm` and PdeIC is never applied (spec 6)
7. ❌ **Segregation interior-face wiring** — Engine treats as "natural" with warning. Only operator-split
8. ❌ **Sub-cycling** — Not implemented (spec 4.10 guardrails)
9. ❌ **Multi-material mesh** — Single FESpace; no InterfaceSubmesh for material interfaces

## Physics Gaps by Priority

### High Priority (core ATHENA/SProcess parity)
- 🟡 Moving mesh handler — exists, missing remesh trigger / smoothing / 3D (`c475863`)
- 🟡 Solution transfer — exists, dose-conserving, deviates from cited `Mx=b` (`c475863`)
- ❌ 3D LevelSetToMesh
- ✅ KMC recombination/clustering events — implemented + tested in `ac18975` (dissociation implemented, untested)
- 🟡 KMC Si diamond lattice — still cubic + body diagonals, not 2×FCC
- 🟡 PDE API engine integration — partial bridge, reaction terms dropped
- ❌ Parameter database for defects/clusters/segregation
- ✅ ChargedFermi quadrature-point D — closed (`0f88a0a`)
- ✅ Copper drift+pairing FEM — drift via `ConvectionIntegrator` landed (`ac18975`); test verifies assembly runs but not drift direction
- ✅ SiGe FEM assembly + defect-mediated D_inter — landed (`ac18975`); defect-mediated branch untested
- ✅ Carbon detailed-balance `k_r` — landed (`ac18975`); test bypasses it
- ✅ ChargedEquilibrium Fermi partition — landed (`ac18975`); test checks dose only
- 🟡 Flash/laser FEM heat transfer (currently 1D explicit)

### Medium Priority (advanced models)
- ❌ III-V I/V equilibrium
- ❌ SPER level-set front tracking
- ❌ FDTD optical absorption
- ❌ Crystallinity phase field
- ❌ Voronoi grain tessellation
- ❌ Scanning laser model
- ❌ Runtime hp-AMR with ThresholdRefiner
- 🟡 KMC event tree — O(log N) selection but O(N) rebuild per step; true incremental heap pending
- 🟡 KMC amorphous pocket, dopant-defect pairing, impurity clustering events
- ❌ Regression test snapshots
- ❌ KMC vs continuum validation

### Low Priority (stretch goals)
- ❌ FitLine / FitPearson / FitPearsonFloor
- ❌ VTK/TDR/CSV output from ResultsExtractor
- ❌ Nonlinear Jacobian (strategy b) engine wiring
- ❌ Sub-cycling
- ❌ FDTD for sub-wavelength features
- ❌ Adjoint inversion / GP surrogates (spec 10 stretch)

### Test-depth gap (cross-cutting, `ac18975`)
The fix commit landed correct physics for C1/C3/I1/I2/I3 but the tests in `tests/diffusion/testDiffusion.cpp` mostly still verify "runs + dose conserved," not the cited quantitative physics:
- Carbon: test sets `setReverseRate(0.0)` — bypasses the equilibrium `C_sI_eq = C_s·C_I/C*_I` it should verify
- Copper: drift assembled (E-field set) but only dose conservation asserted, not drift-direction asymmetry
- ChargedEquilibrium: dose conservation only, not `D_eff` shift with `n/ni`
- SiGe: only the Arrhenius fallback path exercised; `setDefectMediated` branch untested
- KMC dissociation: implemented, no test assertion
- Dose tolerance is 5% (`< 0.05`) throughout; plans specify 0.1% in several places

**Recommended:** before any more model headers land, add equilibrium/ratio tests (`CI/(C_s·C_I) → 1/C*_I` at steady state, `D_eff(n₁)/D_eff(n₂)` matches formula, drift moves Cu up-gradient).

---

## Code Verification Audit (post-`ac18975`)

Direct read of every file touched by `ac18975 fix: address phases 3-11 full-depth review (C1-C3, I1-I5)`. Each review finding verified against current code:

| Review ID | Finding | Status | Code evidence |
|---|---|---|---|
| **C1** | CopperDiffusion drift false-closure | ✅ Fixed | `MobileImpurity.hpp:130-150` adds `mfem::ConvectionIntegrator` with `a = −D·(q/kT)·z·E`; pairing reaction in `PairRateCoef` (lines 210-223). `CopperDiffusion.hpp` is now a 33-line thin specialization. |
| **C2** | KMC event tree O(N) | 🟡 Partial | Selection is `std::lower_bound` on prefix sum (O(log N), `KmcAtomisticEngine.hpp:152-158`), but `events_` is fully rebuilt every step (lines 48-132) → per-step cost still O(N·neighbors). Comment at line 5-6 admits "full incremental heap is a follow-up." |
| **C3** | KMC clustering/dissociation missing | ✅ Fixed | `step()` constructs Cluster (`I+I→{311}`, lines 112-130) and Dissociate (lines 62-78); `apply()` handles both (lines 212-229). Test asserts `clusterCount()>0`, `countSpecies(3)>=1` (`testDiffusion.cpp:583-606`). Dissociation has no test. |
| **I1** | Carbon missing reverse rate | ✅ Fixed (code) / 🟡 untested | `reverseRate()` returns `kf·C*_I` via `PointDefectEquilibrium::C_I_eq(T,"Si")` (`CarbonDiffusion.hpp:41-48`); `TrapCoef::Eval` uses `scale·(kf·c·i − kr·ci)` (line 136). Test sets `setReverseRate(0.0)` — bypasses the equilibrium. |
| **I2** | ChargedEquilibrium wrong formula | ✅ Fixed (code) / 🟡 untested | `getDiffusivity` implements `Σ_z D^z·exp(−zη)/Σ_z exp(−zη)` with `η=ln(n/ni)` (`ChargedEquilibriumDiffusion.hpp:53-74`). Test only checks dose conservation. |
| **I3** | SiGe plain Arrhenius, not defect-mediated | ✅ Fixed (code) / 🟡 untested | `interdiffusivity` computes `D_I*·(C_I/C_I*) + D_V*·(C_V/C_V*)` (`SiGeDiffusion.hpp:48-56`); `DefectDCoef::Eval` reads C_I/C_V per QP (lines 90-105), falls back to Arrhenius when defect fields absent. Test exercises only the fallback. |
| **I4** | MobileImpurity duplicates CopperDiffusion (DRY) | ✅ Fixed | `MobileImpurity<NumericType>` is the base; `CopperDiffusion` is a 33-line specialization (`CopperDiffusion.hpp`). Shared `ImpDCoef`, `PairRateCoef`, drift assembly. |
| **I5** | SolutionTransfer deviates from `Mx=b` | 🟡 Documented deviation | Still uses `ProjectCoefficient(SourceSampleCoef)` + `FindPoints(warn=false)` (`SolutionTransfer.hpp:119-122`), NOT the cited `MultiAppProjectionTransfer::assembleL2()`. `ac18975` added miss counting (`unmappedQuadraturePoints`) and `requireMapped_` fail-loud flag (lines 96-103). Defensible but a real edge-case regression on non-matching meshes. |
| **I6** | Tests verify "runs," not physics | 🟡 Partially addressed | KMC clustering test now verifies physics. Carbon/Copper/ChargedEq/SiGe tests still check dose conservation only (see "Test-depth gap" above). |
| **M1** | MovingMeshHandler `interfaceOnly_` dead | ✅ Fixed | Instance `relabel()` method (lines 37-41) forwards `interfaceOnly_` to the static worker. |
| **M2** | OedSource flips globally | ✅ Fixed | `relabelOxidized` defaults `interfaceAdjacentOnly=true` (line 59) and delegates to `MovingMeshHandler::relabelAttributes`. |
| **M4** | Duplicated `kB` constant | ✅ Fixed | `MobileImpurity.hpp:21` defines `inline constexpr double kB_eV`; still re-typed in `SiGeDiffusion.hpp:41,193` and `KmcEvent.hpp:34,39` (partial consolidation). |

**Net effect of `ac18975`:** 7 of 11 findings fully fixed in code (C1, C3, I1, I2, I3, I4, M1, M2); 2 partially fixed (C2 event tree, I5 transfer algorithm — both documented); 1 partially addressed (I6 tests — KMC clustering verified, others still weak); 1 partially consolidated (M4 — shared constant added but not fully propagated).

**What remains genuinely open after `ac18975`:**
1. KMC per-step rebuild is still O(N) (only selection dropped to O(log N))
2. SolutionTransfer still uses point-sampling, not the cited `Mx=b` L2 projection
3. KMC Si diamond lattice is still cubic + body diagonals, not 2×FCC
4. Physics-verification tests for Carbon/Copper/ChargedEq/SiGe defect-mediated paths
5. KMC dissociation has no test
6. Remesh trigger / Laplacian smoothing / 3D moving boundary in MovingMeshHandler
7. PDE API still drops `ReactionPdeTerm` and ignores PdeIC (`PdeApi.hpp:91-93`)

---



## Plan Expansion Status (2026-07-25)

**Context:** every ❌/🟡 gap above now has a manual-cited, MOOSE-cited implementation path in an expanded phase plan under `docs/superpowers/plans/2026-07-20-diffusion-phase{N}.md`. The plans were expanded to the depth of Phase 1 (per-task: manual equation citations with line numbers, MOOSE class citations with verified paths, code/pseudocode, test+commit steps, scope/deferral notes, self-review). Phase 4–8 code has since landed (`0f88a0a`, `e13fba3`, `332bc6b`, `c475863`) and been review-fixed (`ac18975`); see the **Code Verification Audit** section above for per-row status against current code.

**Source grounding (verified via `iconv -f UTF-16LE -t UTF-8`):**
- ATHENA User's Manual Vol. I — Chapter 3 (SSUPREM4 models), Chapter 6 (statements), Appendix B (parameter tables). Note: SPER, Flash/Laser, KMC, LKMC, Alagator are NOT in this volume — they live in SProcess or companion volumes.
- SProcess User Guide N-2017.09 — Chapters 4 (Diffusion), 5 (Atomistic KMC), 6 (LKMC Epitaxy), 7 (Alagator), 11 (Mesh), 14 (Results), 15 (Numerics).
- MOOSE (`3rdparty/moose/`) — framework + modules. Critically, `framework/src/mfem/` ships a complete MFEM subsystem (kernels, markers, indicators, integrators, equation systems, problem/solver wrappers) used as the primary design reference.

**Phase → expanded plan → gap rows closed:**

| Phase | Plan file | Expansion type | Key gap rows now addressed |
|---|---|---|---|
| 4 | `2026-07-20-diffusion-phase4.md` | Tasks 4–9 deepened + new Self-Review | ChargedEquilibrium, Carbon, Nitrogen, Copper, MobileImpurity, OED, TED, DoseLoss, Regression snapshots |
| 5 | `2026-07-20-diffusion-phase5.md` | **Full rewrite** (52 → ~280 lines) | Anisotropic dual mesh, Voronoi tessellation, oxide breakup, epitaxial regrowth, GB segregation (FEM) |
| 6 | `2026-07-20-diffusion-phase6.md` | **Full rewrite** (64 → ~310 lines) | SiGe interdiffusion (FEM), B-in-SiGe bandgap, SiGeC carbon suppression, Ge-B pairing, strain, III-V sublattice donor/acceptor |
| 7 | `2026-07-20-diffusion-phase7.md` | Header + per-task manual eqs + Self-Review | Si diamond lattice, O(log N) event tree, recombination, clustering, dissociation, pairing, amorphous pocket, all KMC reports, KMC validation |
| 8 | `2026-07-20-diffusion-phase8.md` | **Full rewrite** (73 → ~310 lines) | All LKMC epitaxy model classes, SiGe mole fraction, visibility, twin defects, surface segregation, nonselective epi |
| 9 | `2026-07-20-diffusion-phase9.md` | Header with SProcess eq. 213, TMM, phase-field | FEM heat transfer, TMM optical, FDTD, melting/crystallinity phase field (manual Allen-Cahn forms) |
| 10 | `2026-07-20-diffusion-phase10.md` | Header with Alagator (Ch.7) + Results (Ch.14) | PDE API (lambda vs Alagator-string decision documented), all fitting utilities, parameter DB inheritance/blend, Sobol/MCMC/GP calibration |
| 11 | `2026-07-20-diffusion-phase11.md` | Header with SProcess eq. 982/983 + mesh Ch.11 | Relative/absolute/log refinement criteria, interface-aligned refinement, moving-mesh AMR, hp-AMR |

**Explicitly deferred (bounded, named in plan Self-Reviews):**
- 3D conforming LevelSetToMesh (marching cubes + Delaunay tet) — Phase 1 Task 4 scope note; needed for 3D engine.
- Explicit four-species sublattice I/V equilibrium in III-V (folded into `D_AV`/`D_AI` for now) — Phase 6 Task 8.
- Full Coordination-Reactions CVD set (minimal 3-reaction set shipped) — Phase 8 Task 4.
- True Alagator parity (symbolic weak-form from string) — Phase 10 Task 1 (out of scope; `std::function` lambdas instead).
- Adjoint D(x) inversion + GP surrogates — Phase 10 Tasks 14–15 (stretch).

**Net effect:** every High- and Medium-priority gap in the "Physics Gaps by Priority" section above now has a concrete, manual-grounded implementation plan ready for execution, and Phase 4–8 + Phase 11 mesh code has landed. The Low-priority stretch goals remain flagged as such in their respective plans.

**Next step:** the highest-leverage follow-up is **not** more model code — it is (1) physics-verification tests for the `ac18975` fixes that currently bypass their own physics (Carbon equilibrium, Copper drift direction, ChargedEquilibrium `D_eff` shift, SiGe defect-mediated branch), and (2) the genuinely open engine items: KMC O(N) rebuild, SolutionTransfer `Mx=b` form, remesh trigger, PDE API reaction-term wiring. The expanded Phase 5/8/9 plans remain the right entry points for the not-yet-started physics (polysilicon anisotropic, LKMC epitaxy, flash/laser FEM).

