# PETSc Feasibility Assessment for ViennaPS Diffusion Engine

**Date:** 2026-07-22
**Author:** Research Agent
**Scope:** Evaluate whether adding PETSc to ViennaPS_mod2's multiphysics stack (MFEM + SUNDIALS + AMGCL) is justified by the 11-phase diffusion engine plan.

---

## 1. Executive Summary

**Recommendation: Partial adopt.** Add PETSc as an optional backend gated by `VIENNAPS_HAS_PETSC`, using MFEM's existing PETSc wrappers (`mfem::PetscParMatrix`, `PetscNonlinearSolver`, `PetscODESolver`). PETSc directly solves the most critical bottleneck - the dense O(n^3) SUNDIALS cap at 256 dofs (`SundialsTimeIntegrator.hpp:167`) - by providing sparse, matrix-free, and MPI-parallel solvers through KSP/SNES/TS. PETSc's SNES replaces hand-rolled Picard with Newton-Krylov for Phases 2-3, TS provides native DAE support for Phase 3's reaction networks, and TAO provides the optimization framework for Phase 10's adjoint-based D(x) inversion. The marginal build cost is low: PETSc is already in vcpkg, and MFEM already ships PETSc wrappers (currently disabled: `/* #undef MFEM_USE_PETSC */` in `_config.hpp`). The main cost is rebuilding MFEM with `MFEM_USE_PETSC=ON` and accepting MPI as a runtime dependency. PETSc does NOT replace SUNDIALS (retained for serial CVODE), AMGCL (retained for header-only fallback), or MFEM's AMR (superior to PETSc DMForest for p-refinement).

---

## 2. Current Stack Assessment

### 2.1 What Works

| Component | Status | Evidence |
|-----------|--------|----------|
| MFEM FEM assembly | Functional | `DiffusionEngine.hpp:54-107` assembles M+K, solves `(M+dt*K)u = M*u_prev` |
| Per-species GridFunction storage | Functional | `PhysicsField.hpp:416-420`: `std::map<string, unique_ptr<GridFunction>>` |
| AMGCL sparse solver | Functional (fallback) | `AmgclSolver.hpp:77-101`: BiCGSTAB+AMG, header-only |
| SUNDIALS CVODE (small systems) | Functional but capped | `SundialsTimeIntegrator.hpp:210`: `CVodeCreate(CV_BDF)` |
| Physics-discretization separation | Design complete | `DiffusionPhysics.hpp` with MOOSE PhysicsBase gatekeepers |
| Composition gatekeepers | Implemented | `DiffusionPhysics.hpp:523-530`: `shouldCreateTimeDerivative()` |

### 2.2 Critical Bottlenecks

**Bottleneck 1: Dense O(n^3) SUNDIALS cap (showstopper for production)**

`SundialsTimeIntegrator.hpp:166-178`: Caps at 256 dofs, subsamples field state. A 32x32 2D mesh has ~1024 dofs - 4x over the cap. Unpacked result is approximate (`SundialsTimeIntegrator.hpp:283-288`).

**Bottleneck 2: No sparse linear solver in CVODE**

`SundialsTimeIntegrator.hpp:241-251`: Dense LS only for `neq <= 64`; for larger systems, `CVodeSetLinearSolver(cvode_mem, nullptr, nullptr)` - no linear solver at all. Code comment: "if no LS, switch to explicit after first failure."

**Bottleneck 3: No DAE support**

SUNDIALS IDA is linked (`CMakeLists.txt:314`) but never used. Phase 3's clustering reaction networks are stiff DAEs (singular mass matrix for auxiliary species). Only CVODE (pure ODE) is used.

**Bottleneck 4: No nonlinear solve (Picard hand-rolled)**

Phase 1 plan (`phase1.md:670-676`) defaults to Picard. Phase 2 (`phase2.md:56-62`) requires `dD/dC` Jacobian paths but engine has no SNES. `DiffusionEngine::solve()` uses `BiCGSTAB+DSmoother` (`DiffusionEngine.hpp:100-105`), not even the planned HypreBoomerAMG.

**Bottleneck 5: No adjoint/sensitivity infrastructure**

Phase 10 Task 15 requires adjoint-based D(x) inversion. No adjoint solver, no sensitivity analysis, no optimization framework exists in the current stack.
---

## 3. PETSc Capability Analysis

### 3.1 Solver Stack Overlap and Replacement Value

| Capability | PETSc | Current Stack | Gap |
|------------|-------|---------------|-----|
| Sparse iterative solve in time integration | KSP + any PC via TS | Dense only (`SundialsTimeIntegrator.hpp:242`) | **Critical** |
| Newton-Krylov nonlinear solve | SNES (line search, trust region, matrix-free) | Hand-rolled Picard | **Critical** |
| DAE support (singular mass matrix) | TS `TSSetIFunction` F(t,u,du/dt)=G(t,u) | CVODE only; IDA linked but unused | **Critical** |
| IMEX time integration | TSARKIMEX (orders 1-5, 15+ schemes) | None | High |
| Matrix-free Jacobian | `MatCreateShell` + `MatFDColoring` | Not available | High |
| Adjoint/sensitivity | TSAdjoint (forward + adjoint) | None | Critical (Phase 10) |
| Optimization (TAO) | 20+ solvers (Newton-TR, BQNLS, POUNDERS, LCL) | None | Critical (Phase 10) |

**Can PETSc TS replace SUNDIALS CVODE?** Yes. TS supports BDF (orders 1-6, adaptive), native mass matrix handling (`F(t,u,du/dt)=M*du/dt-f(t,u)`), DAEs via the same interface, and IMEX (ARKIMEX). SUNDIALS CVODE is retained for serial small-system efficiency. PETSc TS manual documents "Using Sundials from PETSc" - PETSc can wrap SUNDIALS as a TS backend, providing migration path.

*Source: PETSc TS manual, https://petsc.org/main/docs/manual/ts/, section "Using Sundials from PETSc"*

**MFEM-PETSc Integration (strongest argument):**

MFEM ships native PETSc wrappers at `linalg/petsc.hpp` (1018 lines, author: Stefano Zampini):

| MFEM Class | PETSc Wrap | Purpose |
|------------|-----------|---------|
| `PetscParVector` | `Vec` | Parallel vector, extends `mfem::Vector` |
| `PetscParMatrix` | `Mat` | Converts from `HypreParMatrix` (zero-copy) |
| `PetscLinearSolver` | `KSP` | Linear solver, accepts MFEM operators |
| `PetscNonlinearSolver` | `SNES` | Nonlinear solver, accepts MFEM `Operator` |
| `PetscODESolver` | `TS` | Time integrator, accepts `TimeDependentOperator` |
| `PetscBDDCSolver` | PC BDDC | Balancing domain decomposition |

MFEM PETSc examples at `F:\dev\mfem\examples\petsc\` include ex9p (time-dependent PDE), ex10p (JFNK + matrix-free), ex5p (fieldsplit for block systems). Current MFEM build has `MFEM_USE_MPI` but NOT `MFEM_USE_PETSC` (`_config.hpp`). Requires `PETSC_HAVE_HYPRE` (`petsc.hpp:45-47`).

*Source: MFEM source at `F:\dev\mfem\linalg\petsc.hpp:19-47,155-299`; examples at `F:\dev\mfem\examples\petsc\`*

### 3.2 Scalable Linear Algebra

PETSc `Mat` supports sparse AIJ/BAIJ, nested MATNEST (multi-species), matrix-free MATSHELL. KSP provides 25+ Krylov methods with any preconditioner. `PetscParMatrix` from `HypreParMatrix` uses `MatCreateMPIAIJWithArrays` - zero-copy. `PetscParVector` extends `mfem::Vector` and shares `Memory<real_t>` via `PlaceMemory()`. `PetscODESolver` accepts any `TimeDependentOperator` - `DiffusionEngine` implements this without changing assembly code.

*Source: PETSc Mat manual https://petsc.org/main/docs/manual/mat/; MFEM `linalg/petsc.hpp` class declarations*

### 3.3 Nonlinear Solve Quality

SNES provides: `SNESNEWTONLS` (default, cubic backtracking), `SNESNEWTONTR` (trust region), matrix-free Newton (`-snes_mf`), colored FD Jacobian (`MatFDColoring`), `SNESSetPicard` (Picard mode), `SNESVINEWTONRSLS` (variational inequalities for solid solubility caps), `SNESFAS` (nonlinear multigrid), `SNESNASM` (nonlinear ASM for block preconditioning).

For Phase 2 Fermi diffusion: `FermiDCoef` (`phase2.md:74-84`) serves as the nonlinear residual. SNES computes Jacobian via colored FD - no `dD/dC` derivation needed. `SNESVINEWTONRSLS` handles Phase 2 solid solubility as proper bound constraints instead of ad-hoc clamping (`PhysicsField.hpp:409-413`).

*Source: PETSc SNES manual https://petsc.org/main/docs/manual/snes/, Tables 9-10*

### 3.4 Optimization and Calibration (Phase 10 Tasks 12-15)

**TAO Solvers relevant to Phase 10:**

| TAO Solver | Type | Phase 10 Task |
|------------|------|---------------|
| `TAOPOUNDERS` | Nonlinear least-squares (derivative-free) | Task 13: SIMS profile fitting |
| `TAOBNLS` / `TAOBNTR` | Bounded Newton LS/TR | Task 15: D(x) inversion with bounds |
| `TAOLCL` | PDE-constrained Lagrangian | Task 15: adjoint-based inversion |
| `TAOBQPIP` | Bounded quadratic IP | Task 12: Sobol batch optimization |
| `TAOTRON` | Trust-region Newton (bounds) | Task 15: regularized D(x) inversion |

**Adjoint-based D(x) inversion (Task 15):** PETSc TSAdjoint provides forward + adjoint timestepping. MOOSE's `OptimizeSolve` wraps TAO with 16 solver types. The `ElementOptimizationDiffusionCoefFunctionInnerProduct` pattern maps directly to TAO's gradient callback. PETSc's `TaoSetObjectiveAndGradient()` + `TaoSetHessian()` provide the optimization framework. The adjoint solve uses `TSAdjointSolve()` - no manual adjoint code needed.

**MCMC/Bayesian (Task 13):** PETSc has NO built-in MCMC. MOOSE's `AffineInvariantStretchSampler` (emcee) is a MOOSE-layer feature, not PETSc. For Task 13, PETSc provides the forward model and TAO for MAP estimation, but MCMC must be implemented separately (e.g., via Python wrapper calling PETSc forward solves).

**GP surrogates (Task 14):** PETSc has NO GP. `PetscRegressor` provides linear/logistic regression only. GP must be implemented separately (e.g., via scikit-learn wrapping PETSc forward solves).

*Source: PETSc TAO manual https://petsc.org/main/docs/manual/tao/; TAO solver table https://petsc.org/main/docs/overview/tao_solve_table/; Phase 10 plan `phase10.md:171-186`*

### 3.5 Phase-Field Support (Phase 9)

PETSc TSARKIMEX is well-suited for Allen-Cahn equations. The IMEX formulation splits stiff diffusion (implicit) from non-stiff reaction (explicit). ARKIMEX schemes up to order 5 with stiff accuracy are available. For the coupled thermal + phase-field + diffusion system (Phase 9 Task 7), TS with IMEX and field split (`PCFIELDSPLIT`) provides the multi-physics coupling framework. SUNDIALS CVODE BDF treats everything implicitly - less efficient for Allen-Cahn where the reaction term is non-stiff compared to diffusion.

*Source: PETSc TS manual https://petsc.org/main/docs/manual/ts/, Table 14 (IMEX Runge-Kutta schemes); Phase 9 plan `phase9.md:38-52`*

### 3.6 AMR Integration (Phase 11)

PETSc has `DMForest` (p4est-based AMR) but MFEM's native AMR is superior: MFEM supports h- and p-refinement; PETSc DMForest is h-only. MFEM's `ThresholdRefiner` + `L2ZienkiewiczZhuFluxEstimator` are already planned (`phase11.md:69-75`). MFEM's prolongation/restriction operators enable O(N) field transfer after AMR.

**Recommendation:** Keep MFEM native AMR. PETSc DMForest is not needed.

*Source: Phase 11 plan `phase11.md:65-100`; PETSc DMForest manual pages https://petsc.org/main/docs/manualpages/DMForest/DMFOREST*

### 3.7 MOOSE Alignment

Since the plans extensively reference MOOSE patterns, and MOOSE uses PETSc internally, adopting PETSc makes porting MOOSE patterns significantly easier:

| MOOSE Pattern | PETSc Equivalent | Currently in ViennaPS |
|---------------|-----------------|----------------------|
| `DerivativeMaterialInterface<Kernel>` (auto Jacobian) | `MatFDColoring` (colored FD Jacobian) | Hand-rolled `dD/dC` (`phase2.md:88-98`) |
| `OptimizeSolve` (16 TAO solvers) | TAO direct (same solvers) | None |
| `TransientMultiApp` (sub-cycling) | `TSRHSSplitSetIS` (fast-slow split) | Not implemented |
| `Restartable.h` (checkpoint) | `TSSetSaveTrajectory` + `VecLoad` | Not implemented |

`DerivativeMaterialInterface` is purely a MOOSE-layer feature (AD in C++). PETSc's `MatFDColoring` achieves the same goal (auto Jacobian) via colored finite differences. Both eliminate hand-derived Jacobians.

*Source: Phase 1 plan `phase1.md:670-676`; Phase 10 plan `phase10.md:171-186`; PETSc SNES manual "Finite Difference Jacobian Approximations"*

### 3.8 Build Complexity

| Item | Cost | Detail |
|------|------|--------|
| PETSc library | Zero | Already in vcpkg (`F:\dev\vcpkg\installed\x64-windows\lib\libpetsc.lib`) and prebuilt (`F:\dev\petsc\arch-win32-c-opt\`) |
| MFEM rebuild | Medium | Rebuild with `MFEM_USE_PETSC=ON` - requires `PETSC_HAVE_HYPRE` in PETSc build |
| CMake integration | Low | `find_package(PETSc)` + `target_link_libraries` |
| Header gate | Low | `#ifdef VIENNAPS_HAS_PETSC` following existing pattern |
| MPI dependency | Medium | PETSc requires MPI; can use `PETSC_COMM_SELF` for single-process |

**PETSc-SUNDIALS conflicts:** None. PETSc can use SUNDIALS as a TS backend (documented in TS manual). Both coexist.

**MPI requirement:** PETSc requires MPI, but `PETSC_COMM_SELF` allows single-process execution. Serial code continues to work. MPI parallelism becomes available as a future option.

*Source: PETSc install docs; vcpkg installation verified at `F:\dev\vcpkg\installed\x64-windows\`*

### 3.9 DAE Support

PETSc TS handles DAEs natively via `F(t,u,du/dt) = G(t,u)`:

- **Hessenberg Index-1 DAE** (semi-explicit): `du/dt = f(t,u,z)`, `0 = h(t,u,z)` - directly supported. Covers Phase 3's equilibrium species (charge-state fractions from primary species).
- **Hessenberg Index-2 DAE**: `du/dt = f(t,u,z)`, `0 = h(t,u)` - supported via constraint differentiation.
- **Fully implicit DAE**: `F(t,u,du/dt) = 0` - supported with TSBDF, TSTHETA, TSDIRK.

Phase 3's clustering (dC_311/dt = k_f*C_I^n - k_r*C_311) with equilibrium species is Index-1 DAE. TS handles it directly; SUNDIALS IDA would require a separate code path.

*Source: PETSc TS manual https://petsc.org/main/docs/manual/ts/, section "DAE Formulations"*

### 3.10 Matrix-Free and GPU

PETSc provides `MATSHELL` for matrix-free Jacobian via `MatCreateShell()` + `MatShellSetOperation()`. `MatFDColoring` enables efficient FD Jacobian with graph coloring (O(colors) function evaluations). SNES `-snes_mf` flag activates matrix-free Newton with no matrix storage.

**GPU support:** PETSc provides `VECCUDA`, `VECKOKKOS`, `MATCUDA`, `MATAIJKOKKOS`. MFEM supports GPU via `mfem::MemoryType`. The vcpkg PETSc build includes `petsc_kokkos.hpp`. PETSc's GPU path is more mature than SUNDIALS' (which requires custom `SUNMemory` hooks). For the current `VIENNAPS_USE_GPU` option, PETSc provides a clearer GPU migration path for the solver layer while MFEM handles FEM assembly GPU offload.

*Source: PETSc GPU Support Roadmap https://petsc.org/main/docs/overview/gpu_roadmap/; PETSc SNES manual "Matrix-Free Methods"*
---

## 4. Per-Phase Impact Analysis

### Phase 1: Foundation FEM Diffusion Engine
**PETSc impact: Medium.** The current implicit Euler fallback (`DiffusionEngine.hpp:54-107`) uses `BiCGSTAB+DSmoother`. PETSc KSP would replace this with runtime-configurable solvers (`-ksp_type gmres -pc_type hypre`). The SUNDIALS CVODE path would be replaced by `PetscODESolver` (wrapping TS BDF) for production, eliminating the 256-dof cap. The MFEM assembly code (BilinearForm, LinearForm) remains unchanged - only the solver backend changes.

**Action:** Add `PetscODESolver` as an alternative to `SundialsTimeIntegrator`, gated by `VIENNAPS_HAS_PETSC`. No change to assembly code.

### Phase 2: Fermi/ChargedFermi + Segregation + Solid Solubility
**PETSc impact: High.** Replaces hand-rolled Picard with SNES Newton-Krylov. The `FermiDCoef` coefficient (`phase2.md:74-84`) serves as the nonlinear residual; SNES computes the Jacobian via `MatFDColoring` - eliminating the need for `FermiDdCCoef` (`phase2.md:88-98`). `SNESVINEWTONRSLS` handles solid solubility caps as proper variational inequalities instead of ad-hoc clamping. Segregation BCs (4 Jacobian blocks) map naturally to `MATNEST` + `PCFIELDSPLIT`.

**Action:** Implement `PetscNonlinearSolver` path for Fermi diffusion. Use `MATNEST` for multi-species coupled system.

### Phase 3: CDD + Clustering + Reaction Networks
**PETSc impact: Critical.** This is where PETSc provides the most value:
1. **DAE support:** Clustering with equilibrium species is Index-1 DAE. TS handles it natively; current CVODE cannot.
2. **Stiff nonlinear solve:** SNES with line search handles the stiff reaction networks (I+V recombination, cluster formation/dissociation).
3. **Block preconditioning:** `PCFIELDSPLIT` for the multi-species system (dopant + I + V + clusters) - each species is a field.
4. **SUPG stabilization:** TSARKIMEX treats the convection-dominated pair diffusion with implicit diffusion + explicit advection.

**Action:** Switch to TS for the full CDD system. Use `TSSetIFunction` for the DAE residual. Use `PCFIELDSPLIT` with `MATNEST` for block preconditioning.

### Phase 4: OED + TED + Moving Interface
**PETSc impact: Medium.** Moving boundary mesh treatment uses MFEM's `Mesh::Transform` and subdomain relabeling - PETSc is not directly involved. However, mass-conservative field transfer across mesh topology changes can use PETSc's `DM` intergrid transfer. The dose loss Robin BC maps to a PETSc boundary integral.

**Action:** No PETSc-specific changes. MFEM mesh operations remain.

### Phase 5: Polysilicon Diffusion
**PETSc impact: Low.** Voronoi tessellation and dual mesh are pre-processing steps. The diffusion solve on the dual mesh uses the same KSP/SNES/TS backend.

**Action:** No PETSc-specific changes.

### Phase 6: SiGe/SiGeC + III-V
**PETSc impact: Low.** Bandgap model and strain effects change coefficients (D, n_i), not the solver structure. Ge-B pairing is another reaction network (Phase 3 pattern).

**Action:** No PETSc-specific changes beyond what Phase 3 already requires.

### Phase 7: KMC + Continuum Coupling
**PETSc impact: Medium.** Sub-cycling (KMC dt << FEM dt) maps to TS fast-slow split (`TSRHSSplitSetIS` + `TSARKIMEXSetFastSlowSplit`). Field transfer between KMC bins and FEM mesh uses L2 projection - the mass matrix solve can use KSP. Time-interpolated transfers use `TSGetInterpolationTime`.

**Action:** Use TS fast-slow split for KMC-FEM coupling. KSP for L2 projection solves.

### Phase 8: KMC Lattice Epitaxy
**PETSc impact: Low.** Surface chemistry and coordination-based growth are KMC operations. FEM mesh update during epitaxy is a mesh operation, not a solver operation.

**Action:** No PETSc-specific changes.

### Phase 9: Flash/Laser Anneal
**PETSc impact: High.** Allen-Cahn equations (non-conserved order parameter) are ideal for TSARKIMEX: implicit diffusion + explicit reaction. The coupled thermal + phase-field + diffusion system uses `PCFIELDSPLIT` with separate fields for T, eta_m, eta_c, and dopant. Latent heat coupling is a source term. The multi-physics coupled solve benefits from SNES for the nonlinear coupling. Restart/checkpoint uses `TSSetSaveTrajectory`.

**Action:** Use TSARKIMEX for Allen-Cahn. Use `PCFIELDSPLIT` for multi-physics. Use `TSSetSaveTrajectory` for restart (Phase 9 Task 9).

### Phase 10: PDE API + Results + Calibration
**PETSc impact: Critical for Tasks 12-15.**
- **Task 12 (Sobol/Morris):** PETSc provides batch solve infrastructure. TAO `TAOBQPIP` for optimization. Quasi-random sampling via `PetscDT` (Sobol sequences).
- **Task 13 (Bayesian MCMC):** PETSc provides forward model; MCMC must be external. TAO provides MAP estimate. `PetscDA` (data assimilation) provides ensemble-based methods.
- **Task 14 (GP surrogate):** PETSc has no GP. External library needed. `PetscRegressor` provides linear regression only.
- **Task 15 (Adjoint D(x) inversion):** TSAdjoint provides the adjoint solve. TAO provides optimization (`TAOBNLS`, `TAOLCL`, `TAOTRON`). The gradient computation (inner product) is a TAO gradient callback. This is the single most PETSc-dependent task.

**Action:** Use TSAdjoint + TAO for Task 15. Use `PetscDT` for Task 12 sampling. External MCMC for Task 13.

### Phase 11: Adaptive Mesh Refinement
**PETSc impact: None.** MFEM native AMR is used. PETSc DMForest is not needed (h-only, no p-refinement). Field transfer uses MFEM's prolongation/restriction operators.

**Action:** No PETSc-specific changes.

---

## 5. Cost-Benefit Table

| Dimension | Cost | Benefit |
|-----------|------|---------|
| Build complexity | Medium: rebuild MFEM with `MFEM_USE_PETSC=ON`, add `find_package(PETSc)` to CMake | Eliminates 256-dof cap; enables sparse solvers |
| Maintenance | Low: `VIENNAPS_HAS_PETSC` gate; MFEM wrappers are maintained upstream | PETSc is actively maintained by Argonne (v3.25.3, July 2026) |
| Runtime dependency | Medium: MPI required (can use `PETSC_COMM_SELF` for serial) | MPI parallelism available for free |
| Code changes | Low: new `PetscTimeIntegrator.hpp` + `PetscNonlinearSolver.hpp` wrappers; existing assembly unchanged | SNES replaces Picard; TS replaces CVODE for production; TAO for optimization |
| Performance gain | High: sparse KSP vs dense O(n^3); Newton vs Picard (2-4x fewer iterations) | Production-scale diffusion (>256 dofs) becomes possible |
| Capability gain | Critical: DAE support, adjoint/sensitivity, optimization, IMEX, matrix-free | Phases 3, 9, 10 become feasible |
| Disk size | ~50MB (libpetsc.lib + headers) | Already installed |
| Learning curve | Medium: PETSc options database, callback model | Runtime-configurable solvers (`-ksp_type`, `-snes_type`) |
| MFEM rebuild risk | Medium: must verify `PETSC_HAVE_HYPRE` in vcpkg PETSc | MFEM wrappers are well-tested (9 examples) |
| SUNDIALS conflict | None | PETSc can wrap SUNDIALS as TS backend |

---

## 6. Recommendation

**Adopt PETSc as an optional solver backend, gated by `VIENNAPS_HAS_PETSC`.**

### Specific boundaries:

**DO adopt PETSc for:**
1. Time integration (TS) - replace CVODE for production (>256 dofs); retain CVODE for serial small-system fast path
2. Nonlinear solve (SNES) - replace hand-rolled Picard for Fermi/CDD (Phases 2-3)
3. DAE support (TS) - for Phase 3 clustering with equilibrium species
4. Optimization (TAO) + adjoint (TSAdjoint) - for Phase 10 Task 15
5. IMEX time integration (TSARKIMEX) - for Phase 9 Allen-Cahn
6. Block preconditioning (PCFIELDSPLIT + MATNEST) - for multi-species coupled systems

**DO NOT adopt PETSc for:**
1. AMR - keep MFEM native (h+p refinement > DMForest h-only)
2. MCMC - use external library (emcee via Python)
3. GP surrogates - use external library (scikit-learn/GPy)
4. FEM assembly - keep MFEM BilinearForm/LinearForm
5. Mesh generation - keep MFEM Mesh + LevelSetToMesh
6. Header-only fallback - keep AMGCL for builds without PETSc/MPI

### Coexistence with existing stack:
- `VIENNAPS_HAS_MFEM` + `VIENNAPS_HAS_PETSC`: full production path (PetscODESolver, PetscNonlinearSolver)
- `VIENNAPS_HAS_MFEM` + `VIENNAPS_HAS_SUNDIALS` (no PETSc): current path (CVODE, dense, 256-dof cap)
- `VIENNAPS_HAS_MFEM` only: implicit Euler fallback (BiCGSTAB + DSmoother)
- No MFEM: profile-only mode

---

## 7. Integration Plan

### Step 1: Rebuild MFEM with PETSc

```powershell
# Verify vcpkg PETSc has HYPRE support
# Check: F:\dev\vcpkg\installed\x64-windows\include\petscconf.h for PETSC_HAVE_HYPRE

# Reconfigure MFEM
cmake -B build -DCMAKE_BUILD_TYPE=Release `
  -DMFEM_USE_PETSC=ON `
  -DMFEM_USE_MPI=ON `
  -DPETSC_DIR="F:/dev/vcpkg/installed/x64-windows" `
  -DHYPRE_DIR="F:/dev/vcpkg/installed/x64-windows" `
  ...
cmake --build build --config Release --parallel
```

### Step 2: CMake Integration in ViennaPS

Add to `CMakeLists.txt` after the SUNDIALS block (line ~322):

```cmake
# PETSc for scalable solvers (KSP/SNES/TS/TAO) - optional, requires MPI
find_package(PETSc CONFIG QUIET)
if(PETSc_FOUND)
  message(STATUS "[ViennaPS] Found PETSc ${PETSc_VERSION}")
  target_link_libraries(${PROJECT_NAME} INTERFACE PETSc::petsc)
  target_compile_definitions(${PROJECT_NAME} INTERFACE VIENNAPS_HAS_PETSC)
  if(TARGET mfem AND MFEM_USE_PETSC)
    message(STATUS "[ViennaPS] MFEM has PETSc support - PetscODESolver/PetscNonlinearSolver available")
  else()
    message(WARNING "[ViennaPS] PETSc found but MFEM not built with MFEM_USE_PETSC - rebuild MFEM with -DMFEM_USE_PETSC=ON")
  endif()
else()
  message(STATUS "[ViennaPS] PETSc not found - using SUNDIALS/AMGCL fallback solvers")
endif()
```

### Step 3: New Headers

Create `include/viennaps/fields/PetscTimeIntegrator.hpp`:
- Wraps `mfem::PetscODESolver` (which wraps PETSc TS)
- Implements `TimeDependentOperator` interface
- Gated by `#ifdef VIENNAPS_HAS_PETSC`
- Uses `TSSetType(ts, TSBDF)` for BDF, `TSSetIFunction` for DAE support
- Runtime-configurable via `-ts_*` options

Create `include/viennaps/fields/PetscNonlinearSolver.hpp`:
- Wraps `mfem::PetscNonlinearSolver` (which wraps PETSc SNES)
- Accepts MFEM `Operator` as the nonlinear residual
- Uses `MatFDColoring` for automatic Jacobian
- Gated by `#ifdef VIENNAPS_HAS_PETSC`

### Step 4: Migration Path

1. **Phase 1-2:** Add `PetscTimeIntegrator` as alternative to `SundialsTimeIntegrator`. Engine selects based on `#ifdef VIENNAPS_HAS_PETSC`.
2. **Phase 3:** Switch to TS DAE mode (`TSSetIFunction`). Use `MATNEST` + `PCFIELDSPLIT` for multi-species.
3. **Phase 9:** Use `TSARKIMEX` for Allen-Cahn. Use `PCFIELDSPLIT` for multi-physics.
4. **Phase 10:** Use `TSAdjoint` + TAO for adjoint-based D(x) inversion.

---

## 8. Sources

### Primary Sources (PETSc)
- PETSc User Guide: https://petsc.org/main/docs/manual/
- PETSc TS (ODE/DAE solvers): https://petsc.org/main/docs/manual/ts/
- PETSc SNES (Nonlinear solvers): https://petsc.org/main/docs/manual/snes/
- PETSc TAO (Optimization): https://petsc.org/main/docs/manual/tao/
- PETSc Time Integrator Table: https://petsc.org/main/docs/overview/integrator_table/
- PETSc TAO Solver Table: https://petsc.org/main/docs/overview/tao_solve_table/
- PETSc GPU Support Roadmap: https://petsc.org/main/docs/overview/gpu_roadmap/
- PETSc DMForest: https://petsc.org/main/docs/manualpages/DMForest/DMFOREST

### Primary Sources (MFEM)
- MFEM PETSc wrappers: `F:\dev\mfem\linalg\petsc.hpp` (1018 lines, author: Stefano Zampini)
- MFEM PETSc implementation: `F:\dev\mfem\linalg\petsc.cpp`
- MFEM PETSc examples: `F:\dev\mfem\examples\petsc\` (ex1p-ex11p, including JFNK and matrix-free variants)
- MFEM config: `F:\dev\mfem\build\config\_config.hpp` (confirms `MFEM_USE_MPI` enabled, `MFEM_USE_PETSC` disabled)

### Local Installation
- PETSc in vcpkg: `F:\dev\vcpkg\installed\x64-windows\lib\libpetsc.lib`, headers at `include\petsc*.h`
- PETSc prebuilt: `F:\dev\petsc\arch-win32-c-opt\lib\libpetsc.lib`
- SUNDIALS in vcpkg: linked via `SUNDIALS::cvode`, `SUNDIALS::ida` (`CMakeLists.txt:312-318`)

### ViennaPS Source Code
- `include/viennaps/fields/SundialsTimeIntegrator.hpp` - current CVODE integration (314 lines)
- `include/viennaps/fields/DiffusionEngine.hpp` - FEM assembly + solver (134 lines)
- `include/viennaps/fields/AmgclSolver.hpp` - AMGCL fallback solver (195 lines)
- `include/viennaps/fields/PhysicsField.hpp` - field container with MFEM GridFunction (500 lines)
- `CMakeLists.txt` - build system (454 lines)

### Plan Documents
- Phase 1: `docs/superpowers/plans/2026-07-20-diffusion-phase1.md`
- Phase 2: `docs/superpowers/plans/2026-07-20-diffusion-phase2.md`
- Phase 3: `docs/superpowers/plans/2026-07-20-diffusion-phase3.md`
- Phase 9: `docs/superpowers/plans/2026-07-20-diffusion-phase9.md`
- Phase 10: `docs/superpowers/plans/2026-07-20-diffusion-phase10.md`
- Phase 11: `docs/superpowers/plans/2026-07-20-diffusion-phase11.md`
