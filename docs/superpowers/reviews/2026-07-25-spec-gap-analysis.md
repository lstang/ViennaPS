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
| MovingMesh relabel/ALE/Laplacian/remesh trigger | ✅ | |
| SolutionTransfer L2 Mx=b | ✅ | |
| AdaptiveMeshRefiner mark/refine/ZZ/threshold/derefine API | ✅ | NC derefine no-op on serial Cartesian |
| Runtime AMR in solve | ✅ | Mark-count mid-solve (Euler path); static refine offline |
| Boundary-conforming Delaunay / marching cubes | ❌ | Stretch F6 |

## Section 4: Continuum models — ✅ / stretch

| Feature | Status |
|---|---|
| Constant, Fermi, ChargedFermi, Pair, React, CDD, clusters, SolidSolubility | ✅ |
| ChargedEquilibrium, Carbon detailed balance, Copper drift | ✅ |
| Poly isotropic/anisotropic + **GB FEM segregation residual** | ✅ |
| SiGe FEM + defect-mediated; SiGeC; **GeB pairing FEM**; **strain FEM** | ✅ |
| III-V FEM + I/V eq; SPER + orientation | ✅ |
| Flash **FEM heat**; **MeltDiffusion FEM** (φ-dependent D) | ✅ |
| FDTD / full TMM / scanning laser / adjoint | ❌ | Stretch |

## Section 5: KMC — ✅ / stretch

| Feature | Status | Notes |
|---|---|---|
| Hop/Recomb/Cluster/Dissoc | ✅ | |
| Diamond A–B neighbor stencil | ✅ | Body-diagonal opposite-sublattice (not full 2×FCC coords) |
| Prefix-sum O(log N) select | ✅ | Rebuild still O(N) — incremental heap stretch |
| IDW deatomize + amorphous pocket | ✅ | |
| Epitaxy planar/coord/twin/surf-seg | ✅ | |
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
7. ✅ Poly GB segregation FEM residual  
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
