# Spec-to-Code Gap Analysis

**Date:** 2026-07-25  
**Last reconciled:** parity wave after `ac18975` (this document rewrite)  
**Spec:** `docs/superpowers/specs/2026-07-20-diffusion-parity-design.md`  
**Head:** post-parity commit on `zcode`

Legend:
- ✅ Implemented — FEM-assembled and/or tested at production depth for this fork
- 🟡 Partial — real API + tests, remaining accuracy/stretch gaps noted
- ❌ Out of scope / deferred stretch (named)

---

## Section 2: Architecture

| Spec Component | Status | Notes |
|---|---|---|
| DiffusionEngine | ✅ | CVODE + implicit Euler, multi-species, BCs, Picard, Jacobian hook, runtime AMR counter |
| MeshModule (LevelSet→MFEM) | ✅ | 2D + **3D Cartesian** tetrahedral tagging (`LevelSetToMeshConverter<3>`) |
| ModelRegistry | ✅ | DiffusionPhysics gatekeeper |
| FEMAssembly | ✅ | Species-outer assembly |
| TimeIntegrator | ✅ | SUNDIALS CVODE when built |
| KmcEngine | ✅ | Hop/Recombine/Cluster/Dissociate + reports |
| DeatomizeTransfer | ✅ | Atomize/Deatomize + continuum coupler |
| ResultsExtractor | ✅ | Vector + **MFEM cut1D/dose**, CSV |
| PdeAPI | ✅ | `applyTo` + ReactionPdeTerm + ICs + Flux BC alias |
| ParameterDatabase | ✅ | Defects, clusters, segregation, oxidation Deal–Grove keys |

## Section 3: Mesh Generation

| Spec Component | Status | Notes |
|---|---|---|
| LevelSetToMesh 2D | ✅ | Cartesian triangles + attributes |
| LevelSetToMesh 3D | ✅ | Cartesian tets + attributes (not marching-cubes conforming) |
| MovingMeshHandler | ✅ | Relabel, ALE lift, **Laplacian smooth**, **aspect remesh trigger** |
| SolutionTransfer | ✅ | **L2 `M x = b`** + dose rescale (`transferL2`) |
| AdaptiveMeshRefiner static | ✅ | markBox / refineMarked |
| Runtime AMR in engine | 🟡 | Hook counts refine opportunities; FES rebuild opt-in |
| Boundary-conforming Delaunay | ❌ | Stretch (F6) |
| ZZ / hp-AMR / derefine | 🟡 | Static box path; full ThresholdRefiner stretch |

## Section 4: Continuum models

| Model / feature | Status | Notes |
|---|---|---|
| Constant / Fermi / ChargedFermi | ✅ | QP D where required |
| Pair / ChargedPair / React / ChargedReact / CDD | ✅ | FEM |
| ChargedEquilibrium | ✅ | Fermi charge-state partition; ratio-tested |
| Carbon | ✅ | Detailed balance `kr=kf·C*_I`; eq residual tested |
| Copper / MobileImpurity | ✅ | Drift + pairing base |
| Clusters {311}/VC/BIC/Loop | ✅ | FEM residuals (not bulk stiffness of immobile clusters) |
| SolidSolubility / Segregation / DoseLoss / OED / TED | ✅ | Engine or FEM paths |
| Polysilicon isotropic + grain growth | ✅ | Dual-attr anisotropic partial |
| SiGe FEM + defect-mediated D | ✅ | Ratio-tested |
| SiGeC / strain / GeB | 🟡 | Host 1D + factors; FEM Ge path primary |
| III-V | ✅ | FEM D + **C_I_eq / C_V_eq** proxies |
| SPER | ✅ | Arrhenius v + **orientation factor** |
| Flash heat | ✅ | **FEM HeatTransfer** model + 1D FTCS |
| Melting / melt D | 🟡 | Phase-field host; FEM heat closed |
| FDTD / TMM / scanning laser | ❌ | Stretch |

## Section 5: KMC

| Feature | Status | Notes |
|---|---|---|
| BKL | ✅ | |
| Event select O(log N) | 🟡 | Prefix-sum selection; rebuild O(N) |
| Hop / Recomb / Cluster / Dissoc | ✅ | Dissoc tested |
| Diamond topology | 🟡 | Cubic + diagonals + A/B sublattice parity |
| Reports (profile, supersat, histogram) | ✅ | `KmcReport` |
| Epitaxy planar / coord / twin / surf-seg | ✅ | `KmcEpitaxyModel` |
| Amorphous pocket / full IDW deatomize | 🟡 | Partial |

## Section 6: PDE API

| Feature | Status | Notes |
|---|---|---|
| PdeEquation → engine | ✅ | buildModels + applyTo + applyICs |
| Diffusion + Reaction terms | ✅ | `LinearReactionDiffusion` |
| Dirichlet/Neumann/Robin/Flux | ✅ | |
| SegregationBC struct | 🟡 | Struct + m; operator-split path for FEM |
| Symbolic Alagator strings | ❌ | Stretch (`std::function`/terms instead) |

## Section 7–9: DB / Results / Tests

| Feature | Status | Notes |
|---|---|---|
| Defect/cluster/segregation/oxidation params | ✅ | ParameterDatabase builtins |
| FitArrhenius / FitLine / FitPearson / Floor | ✅ | |
| CSV profile write | ✅ | |
| Unit + integration + parity suite | ✅ | `TestParityGapClosures` |
| Snapshot regression / KMC-vs-continuum gold | 🟡 | Manual-style checks; no gold files |
| Mesh convergence CTest | 🟡 | AMR + L2 tests |

## Engine integration (former gaps)

1. ✅ Moving mesh + remesh trigger + Laplacian  
2. ✅ Solution transfer L2  
3. 🟡 Runtime AMR hook (count); full rebind deferred  
4. ✅ Jacobian strategy (b) call path  
5. ✅ 3D LevelSetToMesh Cartesian  
6. ✅ PDE API reaction + IC  
7. 🟡 Segregation still primarily operator-split  
8. ❌ Sub-cycling (stretch)  
9. 🟡 Multi-material: attribute tagging; InterfaceSubmesh stretch  

## Explicit stretch (not claimed as parity)

- Conforming 3D cut-cell / Delaunay remesh  
- Full incremental KMC heap + true 2×FCC diamond geometry  
- FDTD optics, adjoint inversion, GP surrogates  
- Full Alagator string language  
- Production TCAD-calibrated parameter tables  

## Verification

`testDiffusion` Release: **All diffusion tests passed** including `[parity-gap-closures] PASS`.

---

*This document is the living parity ledger for the multiphysics fork. Stretch ❌ items remain explicitly out of the production claim.*
