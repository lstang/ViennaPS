# ViennaPS Gap Analysis: Current Code vs. Commercial TCAD Manuals

**Date:** 2026-08-01 (revision of 2026-07-27 edition)
**Scope:** Compare the current `ViennaPS_mod` C++ codebase against the documented capabilities of two commercial process simulators:
- Silvaco **ATHENA** (`Manual/athena_users1.md` / `.pdf`)
- Synopsys **Sentaurus Process** (`Manual/sprocess_ug.md` / `.pdf`)

**Purpose:** Identify where ViennaPS already matches, partially covers, or completely lacks the physics, workflows, and features that commercial users expect. The output is a prioritized roadmap for closing the most important gaps.

> **Revision note (2026-08-01).** This edition reconciles the original analysis with the current tree. Since the 2026-07-27 edition: the `fields/` layer has grown a ~30-file diffusion model library (`fields/models/`), an atomistic KMC engine (`fields/kmc/`), and a parallel (ParMesh/Hypre) diffusion engine with runtime AMR; III-V, flash/laser-anneal, polysilicon-grain, and process-orchestration components landed; and the vendored MOOSE framework in `3rdparty/moose/` was surveyed (design reference, not linked into the build). Several "missing" classifications in the old edition were stale and have been corrected (KMC, 1D cutline extraction, adaptive meshing). CMP remains genuinely missing (see §4.5).

---

## 1. Executive Summary

ViennaPS is best characterized as an **advanced topography / surface-evolution simulator** with a rapidly growing **bulk multiphysics layer** in `include/viennaps/fields/`. It already matches or approaches commercial tools in the following areas:

- Topography etching & deposition (plasma, ion-beam, wet, ALD, PECVD, TEOS).
- Thermal silicon oxidation with Deal-Grove kinetics, LOCOS mask bending, and stress coupling.
- Ion implantation (analytic + BCA Monte Carlo).
- Level-set / ray-tracing surface advection with CPU/GPU flux engines.
- Continuum diffusion **model library**: all SProcess Chapter 4 transport models (Constant, Fermi, ChargedFermi, Pair, ChargedPair, React, ChargedReact, NeutralReact, CDD, ChargedEquilibrium, Carbon, Nitrogen, Copper) plus clustering, segregation, OED/TED, polysilicon, SiGe/SiGeC, III-V, and flash-anneal sub-models (see §4.3).
- **Atomistic KMC**: a BKL (rejection-free) engine with Fenwick-tree event selection, Hop/Recombine/Cluster/Dissociate events, amorphous pockets, and continuum deatomize/atomize coupling (§4.4).

Compared to ATHENA and Sentaurus Process, the largest remaining gaps are:

1. **End-to-end diffusion process step** – the physics engine and model library are mature, but the user-facing process model (`psDiffusion.hpp`) is still a written plan (`docs/superpowers/plans/2026-07-27-diffusion-process-model.md`); `psBasicDiffusion.hpp` remains an early stub. The gap is integration, not physics knowledge.
2. **Chemical Mechanical Polishing (CMP)** – not implemented (the existing `psPlanarize.hpp` is a geometric boolean plane-cut, not a Preston-equation polish model).
3. **Lattice KMC epitaxy** – the atomistic KMC engine exists; the lattice-epitaxy model (`KmcEpitaxyModel`) is still a Phase-8 skeleton (planar/coordination growth, segregation, twin helpers) without production facet physics.
4. **General user-defined PDE scripting** (Alagator equivalent) – a composable C++ PDE API (`PdeApi.hpp`) exists; there is no interpreted scripting front-end. The vendored MOOSE framework is the natural substrate for one.
5. **Comprehensive stress/thermomechanics integration** – kernels exist (elastic + viscoelastic) and oxidation couples stress via activation volumes, but stress is not yet recomputed after every etch/deposition step and fed back automatically.
6. **Layout-driven / ICWBEV-style workflow** – GDS reading and a multi-step `ProcessOrchestrator` exist, but not a full layout-driven process flow with mask-to-process binding.

---

## 2. Methodology

1. **Manual baseline** – extracted chapter headings and keyword hits from `athena_users1.md` and `sprocess_ug.md`. The clean text versions were produced from the UTF-16 source files because the originals use a BOM-prefixed UTF-16 encoding.
2. **Code baseline** – inventoried `include/viennaps/{models,process,fields,gds,geometries,materials,compact}` by header file, class name, and doxygen brief. This revision additionally inventoried `fields/models/` (30 files) and `fields/kmc/` (4 files), and surveyed `3rdparty/moose/` (vendored MOOSE framework + modules, design reference only — not linked into the build).
3. **Gap mapping** – matched manual topics to the closest ViennaPS header or example, then classified the coverage as:
   - **Parity / Near-Parity**
   - **Partial / Work-in-Progress**
   - **Missing / Significant Gap**
4. **Revision reconciliation** – stale claims from the 2026-07-27 edition were checked against the current tree and the git history of the diffusion/AMR/KMC work (ParMesh AMR migration, III-V equilibrium equations, KMC Fenwick-tree events, runtime AMR prolongation, poly-GB diffusion, latent-heat anneal).

Supporting artifacts generated during this work:
- `manual_capabilities_summary.md`
- `vienna_ps_capabilities.md`
- `athena_users1_clean.md`, `sprocess_ug_clean.md` (decoded UTF-16 → UTF-8)
- `athena_diffusion_section.md`, `sprocess_diffusion_section.md`, `sprocess_diffusion_key_sections.md`
- `REFINEMENT_REPORT.md` (DiffusionEngine AMR refactoring and MFEM/MSVC static-init limitation)
- `docs/superpowers/` (ADRs, parity design/deferred specs, phase plans 1–11, ParMesh AMR migration plan, diffusion-process-model plan)

---

## 3. Feature-by-Feature Comparison

| Category | ATHENA / Sentaurus Process | ViennaPS (Current) | Coverage |
|----------|----------------------------|----------------------|----------|
| **Geometry construction** | `region`, `init`, `line` / `polygon`, `ICWBEV` layout-driven | `psMakeTrench`, `psMakeFin`, `psMakeHole`, `psMakeStack`, `psMakePlane`, `psGeometryFactory` | **Near-parity** for standard test structures; layout-driven is weaker. |
| **Lithography** | ATHENA OPTOLITH: aerial image, proximity, optical parameters; SProcess mask operations. | `psLithography.hpp` – simple aerial-image / threshold lithography. | **Partial / low** – lacks OPC, proximity correction, full resist models. |
| **Etching** | Isotropic, anisotropic, directional, plasma, RIE, MC etching, dopant-enhanced, crystallographic. | `psIonBeamEtching`, `psPlasmaEtching`, `psCF4O2Etching`, `psSF6O2Etching`, `psSF6C4F8Etching`, `psHBrO2Etching`, `psFluorocarbonEtching`, `psFaradayCageEtching`, `psWetEtching`, `psIsotropicProcess`, `psDirectionalProcess` | **Near-parity** for topography; some chemistry-specific models present. |
| **Deposition** | CVD, PVD, conformal, directional, sputtering, epitaxy. | `psTEOSDeposition`, `psTEOSPECVD`, `psSingleParticleALD`, `psSelectiveEpitaxy`, `psIsotropicProcess`, `psDirectionalProcess` | **Near-parity** for topography; lacks comprehensive CVD/ALD chemistry libraries. |
| **Ion Implantation** | Analytic (Pearson/dual-Pearson), Monte Carlo, damage models, multi-layer tables, tilt/rotation. | `psAnalyticImplant.hpp`, `psMCBcaImplant.hpp` (Binary Collision Approximation, channeling + cascade); `fields/models/TedInitializer.hpp` seeds anneal from implant damage. | **Moderate** – BCA is present; damage-anneal coupling needs wiring into the diffusion step. |
| **Oxidation** | Deal-Grove, Massoud, viscous mechanics, mixed ambient, orientation-dependent, polysilicon. | `psOxidation.hpp` – Deal-Grove, dry/wet, orientation, pressure, stress activation volumes, LOCOS mask bending, native oxide seeding | **Near-parity** for thermal Si oxidation; some advanced models (Massoud, nitridation) missing. |
| **Diffusion** | ChargedReact, React, ChargedPair, Pair, Fermi, Constant, NeutralReact, OED, clustering; SProcess Ch. 4 model matrix. | **Model library complete**: `fields/models/` — Constant, Fermi, ChargedFermi, Pair, ChargedPair, React, ChargedReact, NeutralReact, CDD, ChargedEquilibrium, Carbon, Nitrogen, Copper, SiGe, SiGeC, III-V (GaAs/InP), polysilicon (grain + GB), {311}/vacancy/impurity clusters, dislocation loops, solid solubility, segregation, OED, TED init, dose-loss BC, mobile impurity, SPER-adjacent flash anneal. Engine: `DiffusionEngine` (ParMesh + Hypre PCG/BoomerAMG, CVODE BDF, runtime AMR). `psBasicDiffusion.hpp` is a stub; `psDiffusion.hpp` planned. | **Partial — physics present, process-step integration pending.** |
| **Atomistic KMC** | SProcess Ch. 5: point-defect KMC, deatomize/atomize transfer, defect reports. | `fields/kmc/KmcAtomisticEngine.hpp` – BKL engine, Fenwick-tree O(log N) event selection, Hop/Recombine/Cluster/Dissociate events on a Si diamond lattice, amorphous pockets, `KmcContinuumCoupler`, `KmcAtomize`/`KmcDeatomize`, `KmcReport`. | **Partial** – engine implemented; calibrated rates, histograms, and production validation pending. |
| **Lattice KMC epitaxy** | SProcess Ch. 6: planar/coordination epitaxy, facet growth, twins, SiGe, visibility. | `fields/kmc/KmcEpitaxy.hpp` – Phase-8 skeleton: planar/coordination growth, surface segregation, twin marking, `KmcVisibility`. | **Partial / skeleton** – API exists; not production physics. |
| **Silicidation** | Multi-phase Ni-silicide, stress-dependent, oxygen-retarded. | `psSilicidation.hpp` (demo/stub using `PhysicsField`). | **Partial / stub** – placeholder model, not full physics. |
| **Epitaxy (continuum)** | Selective epitaxy, auto-doping, facet growth. | `psSelectiveEpitaxy.hpp` (level-set growth); `fields/SPERKernel.hpp` for solid-phase epitaxial regrowth. | **Partial** – continuum selective epitaxy present; no LKMC/atomistic facet models. |
| **Flash / laser anneal** | SProcess §4.201–223: heat transfer, melt phase field, liquid-phase diffusion, resolidification. | `fields/models/FlashLaserAnneal.hpp` – `HeatTransfer` with latent-heat melt-fraction coupling (SProcess eq. 213); Allen-Cahn melting phase fields (deferred items D2/D3). | **Partial** – skeleton + phase-field groundwork; no full melt-diffusion/resolidification flow. |
| **III-V compounds** | SProcess §4.251–252; ATHENA eq. 3-239/3-240. | `fields/models/IIIVDiffusion.hpp` (GaAs/InP donor/acceptor on Ga sublattice, carrier-concentration-dependent D); III-V 4-sublattice I/V equilibrium with GaAs data. | **Partial** – equations implemented with material data; not exposed as a process step. |
| **CMP** | Hard/soft polish models, planarization, RATE.POLISH, stress rebalancing. | *None* (`psPlanarize.hpp` is a geometric boolean plane-cut, not a polish model). | **Missing** – no CMP model. |
| **Stress / Mechanics** | Intrinsic stress, stress history, shear-stress viscosity, plane stress, stress-coupled oxidation/diffusion. | `fields/StressKernel.hpp` (elastic + viscoelastic), `MfemElasticityKernel.hpp` (vector elasticity + traction + stress tensor), `GeometryFieldCoupler.hpp`; stress activation volumes in `psOxidation.hpp`; silicide stress injection. | **Partial** – kernels exist; not yet coupled to all process models. |
| **Adaptive Meshing** | Refinebox, adaptive mesh during diffusion/implant, uniform scaling, moving-boundary mesh. | `DiffusionEngine` on `ParMesh` + `EnsureNCMesh`, runtime AMR via `GeneralRefinement` + prolongation (implicit-Euler path), `AdaptiveMeshRefiner`, `MovingMeshHandler` (ParMesh nodes), `SolutionTransfer` (ParGridFunction/HyprePCG, dose-preserving); documented MFEM/MSVC limitation in `REFINEMENT_REPORT.md`. | **Partial** – implemented with caveats (AMR not on CVODE path; no MPI-distributed runs yet). |
| **Layout-driven flow** | ICWBEV Plus, GDS/OASIS mask import, layout-driven process. | `gds/psGDSReader.hpp`, `psGDSGeometry.hpp`, `psGDSMaskProximity.hpp`, `psGDSUtils.hpp`; multi-step `ProcessOrchestrator` (implant → diffuse → oxidize → geometry remap). | **Partial** – GDS reader + C++ orchestration exist; no ICWBEV-style orchestration. |
| **Custom PDEs** | Alagator scripting language for arbitrary diffusion/PDE models. | `fields/PdeApi.hpp` (PdeEquation/PdeTerm/PdeBC/PdeIC, FluxBC, SegregationBC, bridges into `DiffusionPhysics`), `KernelTerm.hpp`/`KernelTerms.hpp`. | **Partial** – C++ composable API; no interpreted language. Vendored MOOSE is the natural scripting substrate. |
| **Material / Parameter DB** | Parameter Database (PDB), material inheritance, calibration. | `fields/ParameterDatabase.hpp` (inheritance + like-material blending), `MaterialPropertySystem.hpp`, `DiffusivityMaterial.hpp` (Fermi D/C-coefficients), `materials/*`; III-V data started. | **Partial** – infrastructure present, content still being filled. |
| **Numerics / Solvers** | Sparse direct/iterative PDE solvers, solver controls. | `DiffusionEngine` (HypreParMatrix + HyprePCG/BoomerAMG-preconditioned CG; CVODE BDF via MFEM SUNDIALS wrapper; implicit-Euler fallback with system-matrix caching), `AmgclSolver.hpp`, `SundialsTimeIntegrator.hpp`, `BandLimitedSolver.hpp`; CPU/GPU flux engines for topography. | **Partial → Moderate** – parallel FEM solvers exist for the fields layer. |
| **Output / Extraction** | TDR, 1D/2D cutlines, histograms, defect extraction, Tonyplot/Sentaurus Visual. | `psWriter.hpp`, `psVTKRenderWindow.hpp`, `psSlice.hpp` (3D→2D cross-sections), `psPointToElementData.hpp`, `compact/psCSVWriter.hpp`; `fields/PdeApi.hpp::ResultsExtractor` — `cut1D` (GridFunction line sampling) + `dose()` (serial + ParGridFunction). | **Moderate** – 1D cutline + dose now exist in the fields layer; lacks TDR/DF-ISE output and integrated plotting. |
| **GPU / HPC** | Limited commercial GPU support. | `psGPUDiskEngine`, `psGPUTriangleEngine`, `psGPULineEngine`, GPU oxidation mode; bulk layer on ParMesh/Hypre (serial MPI_COMM_SELF for now). | **Advantage** for topography; bulk HPC still single-process. |
| **MOOSE framework (reference)** | — | `3rdparty/moose/` vendors the full MOOSE framework + ~30 modules (phase_field, level_set, tensor_mechanics, solid_mechanics, heat_transfer, chemical_reactions, xfem, stochastic_tools, …), incl. an MFEM subsystem (`MFEMProblem`, `MFEMDiffusionKernel`, `MFEMHypreBoomerAMG`, `MFEMRefinementMarker`, `MFEML2ZienkiewiczZhuIndicator`). **Not linked into the build**; it is the design reference for `fields/` (Kernel/Material patterns, species-outer assembly, BoomerAMG default). | **Reference / future substrate** – see §4.11. |

---

## 4. Detailed Gap Discussion

### 4.1 Topography Simulation – Strong Position

ViennaPS is competitive with the topography modules of ATHENA/ELITE and the boundary-movement parts of Sentaurus Process. The breadth of etching and deposition models is notable:

- `psPlasmaEtching.hpp`, `psCF4O2Etching.hpp`, `psSF6O2Etching.hpp`, `psSF6C4F8Etching.hpp`, `psHBrO2Etching.hpp`, `psFluorocarbonEtching.hpp` cover common plasma chemistries.
- `psIonBeamEtching.hpp` and `psFaradayCageEtching.hpp` handle directed ion and charged-wafer effects.
- `psTEOSDeposition.hpp` / `psTEOSPECVD.hpp` and `psSingleParticleALD.hpp` cover CVD/PECVD/ALD.
- `psWetEtching.hpp`, `psIsotropicProcess.hpp`, `psDirectionalProcess.hpp` provide generic geometric etching/deposition drivers.

**Gap:** While the *number* of models is high, they are typically single-chemistry implementations. Commercial tools expose unified parameter sets, damage coupling, and mixed-mode etching/deposition in one statement (e.g., `ETCH` + `RATE.ETCH`). ViennaPS users often need a new model per chemistry.

### 4.2 Oxidation – Near-Parity for Thermal SiO₂

`psOxidation.hpp` is surprisingly complete:
- Deal-Grove Arrhenius rates for dry/wet, multiple crystal orientations.
- Pressure scaling, native-oxide seeding, transfer coefficients.
- LOCOS mask bending with Si₃N₄ visco-elastic parameters.
- Stress-coupled reaction/diffusion activation volumes.
- Volume/surface mesh export with field data.

**Gap:** Advanced models in the manuals such as Massoud oxidation, nitridation, H₂/O₂ ambient mixing, and doping-dependent rates are not yet exposed.

### 4.3 Diffusion and Implantation – Model Library Complete, Process Step Pending

This is the biggest change since the 2026-07-27 edition. The `fields/` layer now contains a **full SProcess Chapter 4 model matrix** in `fields/models/`:

- **Transport models:** `ConstantDiffusion`, `FermiDiffusion`, `ChargedFermiDiffusion`, `PairDiffusion`, `ChargedPairDiffusion`, `ReactDiffusion`, `ChargedReactDiffusion`, `NeutralReactDiffusion`, `CddDiffusion` (full dopant + I/V + pairing + recombination + clustering), `ChargedEquilibriumDiffusion`, `CarbonDiffusion` (I-trapping), `NitrogenDiffusion`, `CopperDiffusion`.
- **Clustering / deactivation:** `Cluster311`, `VacancyCluster`, `ImpurityCluster`, `DislocationLoop`, `SolidSolubility` — matching ATHENA §3.2 and SProcess Ch. 5 cluster chemistry.
- **Interface physics:** `Segregation` (material-interface BC), `OedSource` (oxidation-enhanced diffusion injection), `TedInitializer` (transient-enhanced-diffusion seeding from implant damage), `DoseLossBC` (surface evaporation), `LinearReactionDiffusion`.
- **Material systems:** `PolysiliconDiffusion` (grain-boundary + grain-interior paths, supported by `GrainModel`/`GrainBoundaryMesh`), `SiGeDiffusion`, `SiGeCDiffusion` (C-suppressed TED, `BandgapModel`), `IIIVDiffusion` (GaAs/InP, ATHENA eq. 3-239/3-240 with 4-sublattice I/V equilibrium and GaAs data), `MobileImpurity` (+ tags).
- **Anneal:** `FlashLaserAnneal` + `FlashAnnealFlow` (process-level FEM orchestration coupling surface Beer's-law heat pulse → melting phase field → melt-enhanced dopant diffusion → solidification trapping `SolidificationTrapping` → SPER crystallinity), plus 1D `FlashLaserAnneal::runPulse`.

The engine is no longer a toy: `DiffusionEngine` runs on `ParMesh`/`ParFiniteElementSpace`/`ParGridFunction` with `HypreParMatrix` + `HyprePCG`/`BoomerAMG`-preconditioned CG, integrates in time via CVODE BDF (MFEM SUNDIALS wrapper) or a cached implicit-Euler path, supports Dirichlet/Neumann/Robin (dose-loss) BCs, dose-preserving projection, and runtime AMR (`GeneralRefinement` + prolongation, implicit-Euler path). `PdeApi.hpp` adds composable equation/BC terms bridging into `DiffusionPhysics`.

**What is still missing:** the **user-facing process step**. `psBasicDiffusion.hpp` remains an early stub (dose-scaling demo), and the production `psDiffusion.hpp` is specified but not yet built (`docs/superpowers/plans/2026-07-27-diffusion-process-model.md`, MFEM-gated, `ProcessModelBase` pattern, level-set → mesh conversion via `LevelSetToMesh`). Implantation has `psAnalyticImplant.hpp` and `psMCBcaImplant.hpp` (BCA, channeling + cascade); damage-to-anneal coupling exists as a model (`TedInitializer`) but is not yet wired into the diffusion step. The gap here is **integration and calibration**, not physics knowledge.

### 4.4 KMC, Silicidation and Epitaxy – Atomistic Engine Present, LKMC a Skeleton

- **Atomistic KMC (SProcess Ch. 5):** `fields/kmc/KmcAtomisticEngine.hpp` is a real BKL (rejection-free) engine: Fenwick-tree O(log N) site selection with O(1)-per-event incremental rebuild, Si diamond lattice (`KmcLattice`), events Hop/Recombine/Cluster/Dissociate (`KmcEventType`), amorphous pockets (`KmcAmorphousPocket`), continuum coupling (`KmcAtomize`/`KmcDeatomize`/`KmcContinuumCoupler`), and defect-activity reporting (`KmcReport`). This is a substantial partial implementation of the Ch. 5 capability — the previous "not implemented" classification was stale.
- **Lattice KMC epitaxy (SProcess Ch. 6):** `KmcEpitaxy.hpp` provides `KmcEpitaxyModel` (planar growth, coordination-based growth, surface segregation, twin formation, Ge-fraction-dependent rates) and `KmcVisibility`, but it is explicitly a Phase-8 skeleton — not production facet physics.
- `psSilicidation.hpp` is still a demo model (thickness proxy, stress/dose injection) — not a multi-phase Ni-silicide model.
- `psSelectiveEpitaxy.hpp` provides continuum level-set growth; `fields/SPERKernel.hpp` covers solid-phase epitaxial regrowth.

### 4.5 CMP – Implemented (psCMP.hpp)

`include/viennaps/models/psCMP.hpp` provides `CMP<NumericType, D>` implementing the Preston removal law: $V(x) = K_p \cdot P \cdot v_{\text{rel}} \cdot s(\text{material}) \cdot f_{\text{pattern}}(h)$, with pattern-density modulation $f_{\text{pattern}}(h) = \text{clamp}(1 + \alpha \cdot (h - h_{\text{ref}})/L_p, 0.1, 2.0)$, per-material selectivity, hard stops, and process metadata. Exposed via `viennaps.hpp` and verified by `testCmp`. (Distinct from geometric `psPlanarize.hpp`).

### 4.6 Stress / Thermomechanics – Kernels Exist, Coupling Incomplete

The commercial manuals treat stress as a first-class citizen: intrinsic stress, stress-dependent oxidation/diffusion, stress history, and shear-stress-dependent viscosity. ViennaPS has:
- `fields/StressKernel.hpp` — both `ElasticStressKernel` and `ViscoelasticStressKernel`
- `MfemElasticityKernel.hpp` (MFEM vector elasticity + traction + stress tensor)
- `psOxidation.hpp` uses stress activation volumes for reaction/diffusion rates
- `psSilicidation.hpp` injects a stress dose

The gap is **system-wide stress tracking**: stress is not yet automatically recomputed after every etch/deposition step and fed back into the next process step as it is in commercial tools. SProcess's intrinsic stress, stress rebalancing after etch/dep, and automated stress-history tracing (Ch. 10) have no counterpart yet.

### 4.7 Meshing – Level-Set Mature, Bulk AMR Implemented with Caveats

ViennaPS relies on ViennaLS level sets for surface evolution. For bulk PDEs, the state since the ParMesh AMR migration (11 tasks) is:

- `DiffusionEngine` wraps the mesh as `ParMesh` with `EnsureNCMesh` so nonconforming `GeneralRefinement` is legal.
- Runtime AMR (`refineBetweenSteps` + marking + prolongation) is implemented on both the implicit-Euler path and the CVODE path (via a checkpoint-restart segment loop with live packed-state sync).
- `MovingMeshHandler` (ParMesh node access) and `SolutionTransfer` (ParGridFunction, dose-preserving, HyprePCG) handle moving-boundary and remap cases.
- `REFINEMENT_REPORT.md` documents the MSVC+IPO/MFEM static-initializer crash that historically blocked `GeneralRefinement`, and the mark-only workaround; the EnsureNCMesh root-cause fix landed in the migration.

Not yet equivalent to SProcess's mature moving-boundary/adaptive-refinement library: no MPI-distributed runs (ParMesh runs on MPI_COMM_SELF), no refinebox-style user boxes (CVODE-path AMR now complete).

### 4.8 Layout-Driven Simulation – Reader and Orchestrator Exist, Workflow Missing

`include/viennaps/gds/` provides:
- `psGDSReader.hpp`
- `psGDSGeometry.hpp`
- `psGDSMaskProximity.hpp`
- `psGDSUtils.hpp`

`ProcessOrchestrator.hpp` demonstrates a full C++ multi-step sequence (implant → diffuse → oxidize → geometry mark/remap). But ViennaPS lacks the ICWBEV Plus-style layout-driven orchestration described in SProcess Chapter 13: mask-to-process binding, parametric layout, multi-cell simulation, and automatic mask extraction.

### 4.9 Custom PDE Scripting – C++ API Present, No Language

SProcess **Alagator** lets users write custom diffusion/PDE terms in a high-level scripting language. ViennaPS has:
- `PdeApi.hpp` — `PdeEquation`/`PdeTerm`/`PdeBC` (Dirichlet/Neumann/Robin/Flux/Segregation)/`PdeIC`, `FluxBC`, `SegregationBC`, with a bridge into `DiffusionPhysics` and `ResultsExtractor`.
- `KernelTerm.hpp` / `KernelTerms.hpp` — MOOSE-Kernel-pattern contributions (`DiffusionTerm`, `ReactionTerm`, `CoupledForceTerm`, `EquilibriumSpeciesAuxKernel`).

Today, adding a new PDE requires C++ compilation; there is no interpreted front-end. The vendored MOOSE framework is the most direct path to a scripted/input-file PDE layer (see §4.11).

### 4.10 Output and Extraction – Cutlines Now Exist, Integration Pending

ViennaPS writes VTK/VTP/CSV and has rendering helpers (`psWriter`, `psVTKRenderWindow`, `psSlice` for 3D→2D cross-sections, `psPointToElementData`). Since the last edition, `ResultsExtractor` in `fields/PdeApi.hpp` provides:
- `cut1D` — samples a GridFunction (and ParGridFunction under MFEM_USE_MPI) along an arbitrary line segment with n points via `Mesh::FindPoints`.
- `dose` — total-dose integration on serial and parallel grid functions.

Remaining gaps: TDR/DF-ISE database output, defect/cluster-size histograms, integrated 1D-plot/regression convenience (FitArrhenius/Pearson-style fitting), and Tonyplot / Sentaurus Visual-style interactivity. The gap is mainly **post-processing convenience**, not raw data availability.

### 4.11 The Vendored MOOSE Framework – Design Reference, Future Substrate

`3rdparty/moose/` (≈576 MB) vendors the full **MOOSE** multiphysics framework and ~30 application modules. It is **not linked into the build** (CMake links MFEM/SUNDIALS/amgcl only) and no ViennaPS code compiles against it. Its role today is twofold:

1. **Design reference (documented).** The `fields/` layer was explicitly modeled on MOOSE patterns: `PhysicsKernel` ↔ MOOSE `Kernel`; `MaterialPropertySystem` ↔ MOOSE `Material` system; `KernelTerm` is documented as "one physics contribution (MOOSE Kernel pattern)"; `DiffusionEngine` mirrors MOOSE `MultiSpeciesDiffusionCG`'s species-outer/terms-inner assembly and MOOSE `DiffusionPhysicsBase`'s HypreBoomerAMG default preconditioner (see ADR-0001 and `DiffusionEngine.hpp` header notes).
2. **Reusable physics/mesh infrastructure (potential).** MOOSE ships an MFEM subsystem (`framework/include/mfem/`, MOOSE_MFEM_ENABLED-gated) with `MFEMProblem`, `MFEMDiffusionKernel`, `MFEMNLDiffusionKernel`, `MFEMHypreBoomerAMG`, `MFEMRefinementMarker`, `MFEML2ZienkiewiczZhuIndicator`, and modules directly relevant to the open gaps: `level_set` (surface evolution), `phase_field` (Allen-Cahn melt/amorphous fields), `tensor_mechanics`/`solid_mechanics` (stress, viscoelasticity), `heat_transfer` (anneal), `chemical_reactions` (reaction networks), `xfem` (moving interfaces), `stochastic_tools` (MC/calibration). MOOSE's input-file/parser system is, in essence, the Alagator-style user-facing PDE scripting layer that §4.9 identifies as missing.

**Caveats:** integrating MOOSE is a large undertaking (its own build system, petsc/libmesh dependencies, LGPL licensing), and none of the module capabilities above should be counted as ViennaPS capabilities today. They are options for closing gaps 4.9, 4.6, and part of 4.4 at the cost of a major dependency.

---

## 5. Prioritized Roadmap

### Short Term (0–6 months)

1. **Finish the end-to-end diffusion process step** (`psDiffusion.hpp` per `docs/superpowers/plans/2026-07-27-diffusion-process-model.md`): wrap `DiffusionEngine` + `fields/models/*` kernels in the `ProcessModelBase` pattern, MFEM-gated; deprecate the `psBasicDiffusion.hpp` stub. This converts the largest physics inventory in the repo into a user-facing `DIFFUSE`-equivalent.
2. **Wire implant damage → anneal.** Feed `MCBcaImplant` damage through `TedInitializer`/`DefectClusterKernel` into `DiffusionEngine` (shortens the biggest SProcess/ATHENA workflow gap).
3. **Extend runtime AMR to the CVODE path** (DONE — implemented via segment checkpoint-restart in `DiffusionEngine::solveCVODE`) and add refinebox-style user refinement boxes.
4. **Document the `fields/` multiphysics API.** Header exposure is already complete: `viennaps.hpp` pulls in the entire fields layer, including `fields/models/`, `fields/kmc/`, `PdeApi.hpp`, `ProcessOrchestrator.hpp`, and (MFEM-gated) `DiffusionEngine.hpp`/`LevelSetToMesh.hpp`. What is missing is usage documentation, examples, and a user guide for the model library and engine.
5. **Harden oxidation stress coupling.** Ensure stress fields are initialized, updated, and passed between process steps (`ProcessOrchestrator` is the natural integration point).
6. **Expose `ResultsExtractor` in user workflows** — `cut1D` + `dose` in examples and VTK/CSV export (API exists; add the convenience layer and fitting utilities).

### Medium Term (6–18 months)

1. **Implement a basic CMP model.** (DONE — `psCMP.hpp` Preston law planarization model, pattern density modulation, selectivity, hard stops, `testCmp`).
2. **Production-harden the atomistic KMC engine**: calibrated hop/recombine rates, cluster-size histograms, `KmcReport` output, and validated `KmcDeatomize` round-trips against continuum Fermi results.
3. **Complete LKMC epitaxy** (`KmcEpitaxyModel`): facet growth, twin formation on {111}, SiGe segregation, visibility — moving Ch. 6 from skeleton to usable.
4. **Flash/laser anneal** (DONE — `FlashAnnealFlow` FEM heat-melt-dopant-trapping flow + `SolidificationTrapping` + `testFlashAnneal`).
5. **Expand the material/parameter database.** Add calibrated Si, SiO₂, SiGe, poly-Si, common metals, and dopant data; continue the III-V (GaAs/InP) entries.
6. **Layout-driven example.** Demonstrate GDS → mask → lithography → etch workflow end-to-end, using `ProcessOrchestrator` as the scaffold.

### Long Term (18+ months)

1. **Generic PDE scripting language.** Two paths: (a) build an interpreter on `PdeApi.hpp`/`KernelTerm` so users can add physics without recompiling; (b) integrate the vendored MOOSE framework (input-file-based PDE layer + its MFEM subsystem), accepting the large dependency.
2. **Full thermomechanical stress solver integration.** Elasticity + viscoelasticity coupled to oxidation, deposition, and etch, with intrinsic stress, stress rebalancing, and stress-history tracing (SProcess Ch. 10).
3. **Production bulk HPC**: MPI-distributed ParMesh solves (currently MPI_COMM_SELF), CVODE-path AMR, and 3D adaptive bulk meshing release.
4. **Full ICWBEV-style layout-driven orchestration**: mask-to-process binding, parametric layout, multi-cell simulation.

---

## 6. Risks and Assumptions

- **Scope assumption:** This analysis assumes ViennaPS aims to become a general-purpose volume + topography TCAD simulator, not remain a topography-only tool. If the latter is the goal, most "gaps" are out of scope.
- **Fork-specific features:** The `fields/` multiphysics layer is a local fork addition (MFEM/SUNDIALS/Hypre/amgcl). Its build dependencies are not present in upstream ViennaPS, so portability and maintenance effort must be weighed against feature growth.
- **MFEM toolchain fragility (documented):** `REFINEMENT_REPORT.md` records the MSVC+IPO/MFEM static-initializer crash that historically blocked `Mesh::GeneralRefinement`, plus the mark-only workaround and the EnsureNCMesh root-cause fix. MSVC Release-CRT matching between app and prebuilt MFEM is enforced (`VIENNAPS_MFEM_REQUIRE_MATCHING_CRT`); Debug builds against Release MFEM remain unsupported.
- **MOOSE vendor tree:** `3rdparty/moose/` is ~576 MB, unbuilt, and unlinked. Treating it as a dependency (for PDE scripting, mechanics, or its MFEM subsystem) carries significant build/licensing/integration cost; its capabilities are **not** ViennaPS capabilities until linked.
- **Calibration risk:** Commercial value comes from calibrated default parameter sets. Even where physics is implemented (diffusion model matrix, KMC engine, III-V), calibration to experimental data is a major remaining effort.
- **Architecture risk:** The bulk PDE layer runs on `ParMesh` but only on MPI_COMM_SELF today; the level-set topography engine and the bulk layer communicate via `GeometryFieldCoupler`/`LevelSetToMesh`/`MovingMeshHandler`. Done poorly, coupling can create numerical instability or excessive run time.

---

## 7. Conclusion

ViennaPS already matches or exceeds commercial tools in **surface/topography simulation**, and the bulk multiphysics layer has advanced substantially since the first gap analysis: a complete SProcess Ch. 4 diffusion model library, a parallel Hypre-based FEM engine with runtime AMR, a BKL atomistic KMC engine with continuum coupling, III-V and flash-anneal groundwork, and a C++ PDE API with 1D cutline extraction. The largest remaining gaps are:

- Turning the `fields/` model library and engine into a mature, user-facing **diffusion process step** (`psDiffusion.hpp`) with implant-damage coupling.
- Adding **CMP** and **production LKMC epitaxy**; production-hardening the atomistic KMC engine.
- Building a **user-defined PDE/scripting layer** (the vendored MOOSE framework is the natural substrate).
- Tightening **stress/thermomechanics coupling** across all process steps.

Closing these gaps would move ViennaPS from a specialized topography simulator toward a full-featured open-source TCAD process simulator.
