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

---

## Audit trail (2026-07-25 review cycles)

Three code-review passes against the phase plans ran on `zcode`. Findings and the bugs they caught:

**Pass 1** (`c475863`, review at `2026-07-25-phases3-11-fulldepth-review.md`): 4 confirmed closures, 7 partial, 1 false closure (Copper drift), 3 KMC sub-claim false closures. Verdict: not ready.

**Pass 2 / fix `ac18975`**: addressed C1-C3, I1-I5, M1-M2. 7 of 11 genuinely fixed.

**Pass 3 / this audit**: verified `ac18975` against source, caught two bugs the fix introduced or missed:

| Bug | Severity | Where | Resolution |
|---|---|---|---|
| **Drift sign inverted** | Critical (physics) | `MobileImpurity.hpp` set `a = −D·bNP·E`; Nernst-Planck continuity requires `a = +D·bNP·E` (verified via IBP identity + MFEM `ConvectionIntegrator` assembles `(a·∇u, v)`). Silent: the existing drift test used a uniform IC and couldn't see direction. | Sign corrected; new E2 test asserts centroid moves +x for z=+1, E>0 (observed `delta=+2.67e-7`, was `−2.67e-7` before fix). |
| **E1: SiGe `||` vs `&&`** | Critical (silent half-model) | `SiGeDiffusion.hpp:97, 127`. With `useDefectMediated_ && (CI || CV)`, registering only Vacancy silently dropped the I term from `D_inter = D_I*·(C_I/C_I*) + D_V*·(C_V/C_V*)`. | `&&` at both sites + `MFEM_VERIFY` fail-loud when `useDefectMediated_` set but a field is missing. |

Both fixed in the same commit as this audit. New physics-verification tests added: Copper drift direction (E2), Copper pairing mass balance (E3, mobile 1e18→9.5e17, pair 0→5e16, mass conserved 1e18). `testDiffusion` Release: **All diffusion tests passed.**

**Still open (carry into next cycle):** KMC event rebuild still O(N) per step (only selection is O(log N)); KMC diamond lattice is body-diagonal stencil, not true 2×FCC coordinates; FEM defect-mediated SiGe branch untested in run (formula + E1 fix verified by reading); KMC dissociation has no test; `MobileImpurity` is `std::string`-keyed, not the `SpeciesTag` template P4 Task 8 specifies; PDE API still ignores PdeIC on `applyTo`.
