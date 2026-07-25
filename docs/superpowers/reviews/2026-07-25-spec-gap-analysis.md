# Spec-to-Code Gap Analysis

**Date:** 2026-07-25
**Spec:** `docs/superpowers/specs/2026-07-20-diffusion-parity-design.md`
**Head:** b23d3a4 (zcode branch)

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
| MovingMeshHandler | ❌ | **No file.** No mesh deformation, no remesh trigger |
| Laplacian smoothing | ❌ | |
| Remesh trigger (skewness/Jacobian) | ❌ | |
| 3D moving boundary | ❌ | |
| SolutionTransfer | ❌ | **No file.** Naive copy in AdaptiveMeshRefiner only |
| L2 projection (cross-mesh) | ❌ | |
| Boundary-safe interpolation | ❌ | |
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
| ChargedFermi | 🟡 | 🟡 | Mean-concentration Picard, not quadrature-point D |
| Pair | ✅ | 🟡 | FEM assembled via CddDiffusion |
| ChargedPair | 🟡 | ❌ | Header exists, no FEM assembly |
| React | ✅ | 🟡 | FEM assembled, recombination sink |
| ChargedReact | 🟡 | ❌ | Header exists, no FEM assembly |
| NeutralReact | 🟡 | ❌ | Header exists, no FEM assembly |
| CDD | ✅ | 🟡 | KernelTerm composition + reactions; use-after-free bug |
| ChargedEquilibrium | 🟡 | ❌ | Header exists, no FEM assembly |
| Carbon | 🟡 | ❌ | Header exists, no FEM assembly |
| Nitrogen | 🟡 | ❌ | Header exists, no FEM assembly |
| Copper | 🟡 | ❌ | Header exists, no FEM assembly |

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
| SiGe interdiffusion | 🟡 | 1D explicit step; **no FEM assembly** |
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
| Si diamond cubic lattice | ❌ | Uses simple cubic (6-neighbor) |
| Event tree O(log N) | ❌ | Linear scan O(N) |
| Hop events | ✅ | |
| Recombination (I+V->0) | ❌ | |
| Clustering (I+I->{311}) | ❌ | |
| Dissociation | ❌ | |
| Amorphous pocket | ❌ | |
| Dopant-defect pairing | ❌ | |
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
| MovingMeshHandler.hpp | ❌ Missing |
| SolutionTransfer.hpp | ❌ Missing |
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

The DiffusionEngine itself is missing these spec-required capabilities:

1. **Moving boundary mesh** - No mesh deformation/remesh during solve (spec 3.2)
2. **Solution transfer** - No L2 projection between meshes (spec 3.3)
3. **Runtime AMR** - AdaptiveMeshRefiner exists but not called during solve (spec 3.4)
4. **Nonlinear Jacobian** - assembleStiffnessJacobian virtual exists but engine never calls it (spec 4.10 strategy b)
5. **3D support** - LevelSetToMesh 3D throws. Engine compiles for D=3 but no mesh source
6. **PDE API integration** - PdeEquation cannot drive the engine (spec 6)
7. **Segregation interior-face wiring** - Engine treats as "natural" with warning. Only operator-split
8. **Sub-cycling** - Not implemented (spec 4.10 guardrails)
9. **Multi-material mesh** - Single FESpace; no InterfaceSubmesh for material interfaces

## Physics Gaps by Priority

### High Priority (core ATHENA/SProcess parity)
- ❌ Moving mesh handler (required for oxidation-coupled diffusion)
- ❌ Solution transfer (required for mesh regeneration)
- ❌ 3D LevelSetToMesh
- ❌ KMC recombination/clustering/pairing events
- ❌ KMC Si diamond lattice
- ❌ PDE API engine integration
- ❌ Parameter database for defects/clusters/segregation
- 🟡 ChargedFermi quadrature-point D (currently mean-concentration Picard)
- 🟡 SiGe FEM assembly (currently 1D explicit)
- 🟡 Flash/laser FEM heat transfer (currently 1D explicit)

### Medium Priority (advanced models)
- ❌ III-V I/V equilibrium
- ❌ SPER level-set front tracking
- ❌ FDTD optical absorption
- ❌ Crystallinity phase field
- ❌ Voronoi grain tessellation
- ❌ Scanning laser model
- ❌ Runtime hp-AMR with ThresholdRefiner
- ❌ KMC event tree (O(log N))
- ❌ Regression test snapshots
- ❌ KMC vs continuum validation

### Low Priority (stretch goals)
- ❌ FitLine / FitPearson / FitPearsonFloor
- ❌ VTK/TDR/CSV output from ResultsExtractor
- ❌ Nonlinear Jacobian (strategy b) engine wiring
- ❌ Sub-cycling
- ❌ FDTD for sub-wavelength features
- ❌ Adjoint inversion / GP surrogates (spec 10 stretch)

---

## Plan Expansion Status (2026-07-25)

**Context:** every ❌/🟡 gap above now has a manual-cited, MOOSE-cited implementation path in an expanded phase plan under `docs/superpowers/plans/2026-07-20-diffusion-phase{N}.md`. The plans were expanded to the depth of Phase 1 (per-task: manual equation citations with line numbers, MOOSE class citations with verified paths, code/pseudocode, test+commit steps, scope/deferral notes, self-review). No code was written — this is the planning prerequisite for closing the gaps in implementation.

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

**Net effect:** every High- and Medium-priority gap in the "Physics Gaps by Priority" section above now has a concrete, manual-grounded implementation plan ready for execution. The Low-priority stretch goals remain flagged as such in their respective plans.

**Next step (per the brainstorming skill's terminal state):** invoke the writing-plans / executing-plans workflow against any one expanded phase to begin closing gaps in code. Recommended starting point: **Phase 4** (smallest scope, builds on the already-complete Phase 3; closes 8 gap rows including the High-priority ChargedFermi quadrature-point D fix path).

