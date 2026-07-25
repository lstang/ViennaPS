# PETSc vs Hypre Comparative Analysis for ViennaPS Diffusion Engine

**Date:** 2026-07-22
**Author:** Research Agent
**Scope:** Comparative evaluation of PETSc, Hypre, both, or neither for the 11-phase diffusion engine plan.
**Prior work:** Builds on `2026-07-22-petsc-feasibility.md` (399 lines), which recommended "partial adopt PETSc" but did NOT evaluate Hypre as a standalone option.

---

## 1. Executive Summary

**Recommendation: Adopt Hypre immediately (zero build cost), defer PETSc until Phase 10.**

Hypre is already compiled into the existing MFEM build (`MFEM_HYPRE_VERSION 23200`, `MFEM_USE_MPI`). MFEM's `hypre.hpp` (2356 lines) exposes `HypreBoomerAMG`, `HyprePCG`, `HypreGMRES`, `HypreFGMRES`, `HypreSmoother`, `HypreParaSails`, `HypreEuclid`, `HypreILU`, and more -- all available right now with zero rebuild. Replacing `mfem::DSmoother` with `mfem::HypreBoomerAMG` in `DiffusionEngine.hpp:100` is a one-line change.

Furthermore, SUNDIALS 7.8.0 provides iterative solvers (`SUNLinSol_SPGMR`, `SUNLinSol_SPBCGS`, `SUNLinSol_SPTFQMR`, `SUNLinSol_PCG`) that accept user-defined preconditioners via `SUNLinSolSetPreconditioner()`. A Hypre preconditioner can be wrapped as SUNDIALS `PSetup`/`Psolve` callbacks, eliminating the 256-dof dense cap in `SundialsTimeIntegrator.hpp:167` WITHOUT needing PETSc.

PETSc, by contrast, has a **blocking issue**: the vcpkg PETSc build lacks `PETSC_HAVE_HYPRE` (verified: `petscconf.h` has only 74 `PETSC_HAVE` defines, packages `:blaslapack:mathlib:mpi:`, uses `MPIUNI` stub). MFEM's `petsc.hpp:45-47` requires `PETSC_HAVE_HYPRE` with `#error`. Using MFEM's PETSc wrappers requires rebuilding PETSc with Hypre support AND rebuilding MFEM with `MFEM_USE_PETSC=ON` -- a significant effort.

MFEM also provides native solvers that cover most needs without PETSc: `NewtonSolver` (with Eisenstat-Walker adaptive tolerance), `LBFGSSolver`, IMEX ODE solvers (`IMEXRK2`, `IMEX_DIRK_RK3`), `ARKStepSolver` (SUNDIALS ARKode with IMEX mode), and `KINSolver` (SUNDIALS KINSOL with JFNK support).

**Bottom line:** Hypre fixes 3 of 5 bottlenecks at zero cost. MFEM's native NewtonSolver + IMEX + KINSOL cover Phases 1-9. PETSc's unique value is Phase 10 (TAO optimization, TSAdjoint sensitivity). Defer PETSc until Phase 10 is actually needed.

---

## 2. Hypre Standalone Assessment

### 2.1 What Hypre Provides

Hypre is a parallel linear solver library from LLNL. It provides **linear algebra only** -- no nonlinear solvers, no time integrators, no DAE support, no optimization.

**Solver/preconditioner capabilities (from MFEM `hypre.hpp`):**

| MFEM Class | Hypre Backend | Type | Standalone solver? | Preconditioner? |
|-----------|--------------|------|-------------------|----------------|
| `HypreBoomerAMG` | BoomerAMG | AMG | YES (Mult() solves) | YES |
| `HyprePCG` | PCG | Krylov (CG) | YES | N/A |
| `HypreGMRES` | GMRES | Krylov (GMRES) | YES | N/A |
| `HypreFGMRES` | FGMRES | Krylov (flexible) | YES | N/A |
| `HypreSmoother` | Various | Smoother | NO (single sweep) | YES |
| `HypreParaSails` | ParaSails | Sparse approx inverse | NO | YES |
| `HypreEuclid` | Euclid | Parallel ILU | NO | YES |
| `HypreILU` | ILU | Parallel ILU(k) | NO | YES |
| `HypreAMS` | AMS | Aux-space Maxwell | YES | YES |
| `HypreADS` | ADS | Aux-space divergence | YES | YES |
| `HypreDiagScale` | Jacobi | Diagonal scaling | NO | YES |
| `HypreIdentity` | Identity | No-op | NO | YES |
| `HypreTriSolve` | TriSolve | Triangular solve | YES | YES |
| `HypreLOBPCG` | LOBPCG | Eigenvalue solver | YES (eigen) | N/A |
| `HypreAME` | AME | Maxwell eigenvalue | YES (eigen) | N/A |

*Source: `F:\dev\mfem\linalg\hypre.hpp` (2356 lines, all classes verified)*

**Key observations:**
- `HypreBoomerAMG` can be used as a **standalone solver** (inherits `HypreSolver` -> `Solver`, has `Mult()`) or as a **preconditioner** for `HyprePCG`/`HypreGMRES`/`HypreFGMRES`.
- `HypreBoomerAMG::SetSystemsOptions(int dim)` configures for systems (e.g., vector diffusion).
- `HypreBoomerAMG::SetAdvectiveOptions()` provides AIR-AMG for advection-dominated problems (Hypre >= 2.19.0; MFEM has 2.23.0).
- `HypreParMatrix` provides parallel CSR matrix operations, BC elimination, thresholding, block extraction.

### 2.2 What MFEM's Hypre Wrapper Exposes

```
Solver
  +-- HypreSolver (abstract base)
       +-- HypreBoomerAMG      (AMG - the workhorse)
       +-- HyprePCG            (parallel CG)
       +-- HypreGMRES          (parallel GMRES)
       +-- HypreFGMRES         (parallel flexible GMRES)
       +-- HypreSmoother       (Jacobi, GS, Chebyshev, FIR)
       +-- HypreParaSails      (sparse approximate inverse)
       +-- HypreEuclid         (parallel ILU)
       +-- HypreILU            (native ILU, Hypre >= 2.19.0)
       +-- HypreAMS            (Maxwell auxiliary-space)
       +-- HypreADS            (divergence auxiliary-space)
       +-- HypreDiagScale      (diagonal/Jacobi)
       +-- HypreIdentity       (no-op)
       +-- HypreTriSolve       (triangular solve)
```

Additionally, `HypreParVector` (extends `mfem::Vector`) and `HypreParMatrix` (extends `mfem::Operator`) provide parallel data structures. All available because MFEM was built with `MFEM_USE_MPI` and `MFEM_HYPRE_VERSION 23200`.

*Source: `F:\dev\mfem\linalg\hypre.hpp:67-2356`; MFEM config at `F:\dev\mfem\build\config\_config.hpp`*

### 2.3 Can Hypre Alone Fix the Dense O(n^3) SUNDIALS Bottleneck?

**Not directly.** Hypre provides linear solvers, not time integrators. The SUNDIALS bottleneck (`SundialsTimeIntegrator.hpp:166-178`, 256-dof cap) is in CVODE's linear solver setup. The current code uses `SUNLinSol_Dense` for `neq <= 64` and no linear solver for larger systems.

**But yes, indirectly.** SUNDIALS 7.8.0 provides iterative linear solvers that accept user-defined preconditioners. See Section 4 for the full integration path.

### 2.4 What Hypre Does NOT Provide

| Capability | Hypre | PETSc | MFEM native (no PETSc) |
|-----------|-------|-------|----------------------|
| Nonlinear solve | NO | YES (SNES) | YES (NewtonSolver, KINSolver) |
| Time integration | NO | YES (TS) | YES (ODESolver, CVODE, ARKStep) |
| DAE support | NO | YES (TS) | Partial (IDA linked, unused) |
| IMEX time integration | NO | YES (TSARKIMEX) | YES (IMEXRK2, IMEX_DIRK_RK3, ARKStepSolver::IMEX) |
| Adjoint/sensitivity | NO | YES (TSAdjoint) | NO |
| Optimization | NO | YES (TAO) | Partial (SLBQPOptimizer) |
| Matrix-free Jacobian | NO | YES (MATSHELL) | Partial (Operator with GetGradient()) |

*Sources: `F:\dev\mfem\linalg\solvers.hpp` (NewtonSolver :780, LBFGSSolver :865, SLBQPOptimizer :1036); `F:\dev\mfem\linalg\ode.hpp` (IMEX solvers :1030-1090); `F:\dev\mfem\linalg\sundials.hpp` (KINSolver :896, ARKStepSolver :720)*

---

## 3. PETSc Assessment (Updated from Prior Report)

### 3.1 Critical Finding: PETSc Build Lacks HYPRE

The prior report stated "PETSc is already in vcpkg" and "the marginal build cost is low." This is **incorrect** for the purpose of MFEM integration.

**Verified PETSc configuration** (`F:\dev\vcpkg\installed\x64-windows\include\petscconf.h`):
- `PETSC_HAVE_MPIUNI 1` -- MPI **stub**, not real MPI
- **NO** `PETSC_HAVE_HYPRE` -- Hypre support is **absent**
- `PETSC_HAVE_PACKAGES ":blaslapack:mathlib:mpi:"` -- minimal packages only
- 74 total `PETSC_HAVE` defines (very minimal build)

**MFEM's PETSc wrappers REQUIRE Hypre** (`F:\dev\mfem\linalg\petsc.hpp:45-47`):

```cpp
#if !defined(PETSC_HAVE_HYPRE)
#error "MFEM requires PETSc built with HYPRE support"
#endif
```

**Implication:** The current PETSc build **cannot** be used with MFEM's PETSc wrappers. To use MFEM's `PetscParMatrix`, `PetscLinearSolver`, `PetscNonlinearSolver`, `PetscODESolver`, PETSc must be rebuilt with Hypre support, and then MFEM must be rebuilt with `MFEM_USE_PETSC=ON`.

The same configuration was verified at `F:\dev\petsc\arch-win32-c-opt\include\petscconf.h` -- identical minimal build.

### 3.2 PETSc's Unique Value (Still Valid)

Even with the build issue, PETSc provides capabilities that neither Hypre nor MFEM native solvers offer:

| PETSc capability | Phase needed | Hypre/MFEM alternative | Gap |
|-----------------|-------------|----------------------|-----|
| **TSAdjoint** (adjoint timestepping) | Phase 10 Task 15 | None | **Critical** |
| **TAO** (20+ optimization solvers) | Phase 10 Tasks 12-15 | SLBQPOptimizer (QP only) | **Critical** |
| **TS DAE** (F(t,u,du/dt)=G(t,u)) | Phase 3 | SUNDIALS IDA (linked, unused) | Medium |
| **PCFIELDSPLIT** + **MATNEST** | Phase 3 multi-species | Manual block assembly | Medium |
| **MatFDColoring** (colored FD Jacobian) | Phase 2 | Manual dD/dC derivation | Low |
| **SNESVINEWTONRSLS** (variational inequalities) | Phase 2 solid solubility | Ad-hoc clamping | Low |

### 3.3 PETSc Build Cost (Corrected)

| Item | Prior report estimate | Actual cost |
|------|----------------------|-------------|
| PETSc library | "Zero - already in vcpkg" | **High**: must rebuild PETSc with --with-hypre |
| MFEM rebuild | "Medium" | **High**: must rebuild with MFEM_USE_PETSC=ON |
| MPI dependency | "Medium: can use PETSC_COMM_SELF" | **High**: current PETSc uses MPIUNI stub |
| CMake integration | "Low" | Low (unchanged) |

---

## 4. SUNDIALS + Hypre Integration (Critical Path)

### 4.1 The API Path

SUNDIALS 7.8.0 provides a clean iterative solver + preconditioner API that can wrap Hypre.

**Step 1: Create an iterative SUNLinearSolver**

| Solver | Header | Krylov method | Use case |
|--------|--------|--------------|----------|
| `SUNLinSol_SPGMR` | `sunlinsol_spgmr.h:88` | Scaled Preconditioned GMRES | General nonsymmetric |
| `SUNLinSol_SPBCGS` | `sunlinsol_spbcgs.h:81` | Scaled Preconditioned Bi-CGStab | Moderate size, nonsymmetric |
| `SUNLinSol_SPFGMR` | `sunlinsol_spfgmr.h` | Scaled Preconditioned FGMRES | Variable preconditioner |
| `SUNLinSol_SPTFQMR` | `sunlinsol_sptfqmr.h` | Scaled Preconditioned TFQMR | No restart needed |
| `SUNLinSol_PCG` | `sunlinsol_pcg.h` | Preconditioned CG | SPD systems only |

*Source: `F:\dev\vcpkg\installed\x64-windows\include\sunlinsol\` (all headers verified)*

**Step 2: Set matrix-vector product callback**

```c
// sundials_linearsolver.h:175-176
SUNErrCode SUNLinSolSetATimes(SUNLinearSolver S, void* A_data, SUNATimesFn ATimes);
```

`ATimes` computes `z = A * v`. For the diffusion engine, this is a sparse matvec using `SparseMatrix::Mult()` or `HypreParMatrix::Mult()`.

*Source: `F:\dev\vcpkg\installed\x64-windows\include\sundials\sundials_iterative.h:97-108`*

**Step 3: Set preconditioner callbacks**

```c
// sundials_linearsolver.h:179-180
SUNErrCode SUNLinSolSetPreconditioner(SUNLinearSolver S, void* P_data,
                                       SUNPSetupFn Pset, SUNPSolveFn Psol);
```

- `SUNPSetupFn` (`sundials_iterative.h:120`): `int (*)(void* P_data)` -- builds preconditioner
- `SUNPSolveFn` (`sundials_iterative.h:148`): `int (*)(void* P_data, N_Vector r, N_Vector z, sunrealtype tol, int lr)` -- solves Pz = r

For a Hypre BoomerAMG preconditioner:
- `PSetup`: calls `HypreBoomerAMG::Setup()` (builds AMG hierarchy)
- `Psolve`: calls `HypreBoomerAMG::Mult()` (applies one AMG V-cycle)

**Step 4: Attach to CVODE**

```c
// cvode_ls.h:82-83
int CVodeSetLinearSolver(void* cvode_mem, SUNLinearSolver LS, SUNMatrix A);
// cvode_ls.h:97-98
int CVodeSetPreconditioner(void* cvode_mem, CVLsPrecSetupFn pset, CVLsPrecSolveFn psolve);
```

`CVodeSetLinearSolver` accepts ANY `SUNLinearSolver` type, including iterative solvers. `SUNMatrix A` can be `NULL` for matrix-free (Jacobian applied via `ATimes`).

*Source: `F:\dev\vcpkg\installed\x64-windows\include\cvode\cvode_ls.h:82,97-98`*

### 4.2 Code Changes Required in SundialsTimeIntegrator.hpp

Current code (`SundialsTimeIntegrator.hpp:238-251`):

```cpp
if (neq <= 64) {
    A = SUNDenseMatrix(neq, neq, sunctx);
    LS = SUNLinSol_Dense(y, A, sunctx);
    if (A && LS) CVodeSetLinearSolver(cvode_mem, LS, A);
} else {
    CVodeSetLinearSolver(cvode_mem, nullptr, nullptr);  // NO LINEAR SOLVER
}
```

Proposed replacement (conceptual):

```cpp
// For any neq: use SPBCGS with HypreBoomerAMG preconditioner
LS = SUNLinSol_SPBCGS(y, PREC_LEFT, 5, sunctx);
SUNLinSolSetATimes(LS, &ode_data, cvodeATimes);
SUNLinSolSetPreconditioner(LS, &hypre_amg, cvodePSetup, cvodePSolve);
CVodeSetLinearSolver(cvode_mem, LS, NULL);  // NULL matrix = matrix-free
CVodeSetPreconditioner(cvode_mem, cvodePrecSetup, cvodePrecSolve);
```

The N_Vector to HypreParVector conversion is straightforward: `N_VGetArrayPointer()` returns a `double*` that can be wrapped into `mfem::Vector` via `Vector::SetData()`. For serial execution (MPI_COMM_SELF), `HypreParVector` wraps the same data.

**This eliminates:**
1. The 256-dof cap (line 167) -- no dense matrix needed
2. The subsampling approximation (lines 168-178) -- full state evolves
3. The "no linear solver for neq > 64" fallback (lines 246-251) -- iterative solver works for any size

### 4.3 What Percentage of Bottlenecks Does Hypre-Alone Fix?

| Bottleneck | Hypre alone fixes it? | How |
|-----------|----------------------|-----|
| 1. Dense O(n^3) SUNDIALS cap (256 dofs) | **YES** | SUNDIALS iterative solver + Hypre PC |
| 2. No sparse linear solver in CVODE | **YES** | SUNLinSol_SPBCGS + HypreBoomerAMG PC |
| 3. No DAE support | **NO** | Hypre is linear-only; need IDA or TS |
| 4. No nonlinear solve (Picard hand-rolled) | **PARTIAL** | MFEM NewtonSolver + Hypre LS = Newton-Krylov; but no MatFDColoring |
| 5. No adjoint/sensitivity infrastructure | **NO** | No adjoint in Hypre or MFEM native |

**Score: 2.5 of 5 bottlenecks fixed by Hypre alone (zero build cost).**

With MFEM native NewtonSolver + KINSolver, bottleneck 4 is also addressed (total: 3.5 of 5).

---

## 5. Capability Matrix

| Capability | Hypre alone (via MFEM) | PETSc alone (current build) | PETSc + Hypre (requires rebuild) | Current stack (DSmoother + CVODE dense) |
|-----------|----------------------|---------------------------|----------------------------------|----------------------------------------|
| Sparse AMG preconditioner | YES (HypreBoomerAMG) | NO (no HYPRE in build) | YES (PCHYPRE) | NO (DSmoother) |
| Sparse iterative solve inside time integration | YES (wrap in SUNDIALS) | YES (KSP) | YES | NO (dense, 256 cap) |
| Nonlinear solve (Newton-Krylov) | YES (MFEM NewtonSolver + Hypre LS) | YES (SNES) | YES | NO (Picard hand-rolled) |
| DAE support | NO (linear only) | YES (TS) | YES | NO (CVODE only) |
| IMEX time integration | YES (MFEM IMEX RK or ARKStep) | YES (TSARKIMEX) | YES | NO |
| Adjoint/sensitivity | NO | YES (TSAdjoint) | YES | NO |
| Optimization | NO (SLBQPOptimizer only) | YES (TAO 20+ solvers) | YES | NO |
| Matrix-free Jacobian | Partial (Operator + GetGradient) | YES (MATSHELL + MatFDColoring) | YES | NO |
| Block preconditioning | Partial (HypreParMatrixFromBlocks) | YES (PCFIELDSPLIT + MATNEST) | YES | NO |
| Variational inequalities | NO | YES (SNESVINEWTONRSLS) | YES | NO (ad-hoc clamping) |
| AMR | YES (MFEM native, h+p) | YES (DMForest, h-only) | YES | YES (MFEM native) |
| GPU support | YES (MFEM CUDA + Hypre GPU) | YES (VECCUDA, MATCUDA) | YES | Partial (MFEM CUDA) |

---

## 6. Cost Comparison

| Dimension | Hypre alone | PETSc alone | Hypre + PETSc | Neither (status quo) |
|-----------|-------------|-------------|---------------|----------------------|
| Build cost | **ZERO** (already in MFEM) | HIGH (rebuild PETSc + MFEM) | HIGH (rebuild PETSc + MFEM) | 0 |
| MFEM rebuild needed? | NO | YES (MFEM_USE_PETSC=ON) | YES | NO |
| New runtime deps | None | MPI (real, not MPIUNI) | MPI | None |
| Code changes | Low (swap DSmoother, add SUNDIALS PC wrappers) | Medium (new PetscTimeIntegrator, PetscNonlinearSolver) | Medium | None |
| Bottlenecks fixed | 3.5/5 | 4.5/5 (if built with HYPRE) | 5/5 | 0/5 |
| Maintenance burden | Low (MFEM maintains wrappers) | Medium (PETSc options database) | Medium | None |
| Disk size | 0 (already installed) | ~50MB (already installed but useless without rebuild) | ~50MB | 0 |
| Learning curve | Low (MFEM Solver interface) | Medium (PETSc callback model) | Medium | None |

---

## 7. Per-Phase Recommendation

| Phase | What it does | Hypre sufficient? | Need PETSc? | Recommended solver |
|-------|-------------|-----------------|-------------|-------------------|
| 1 | Foundation: FEM mesh, constant D, CVODE, HypreBoomerAMG | **YES** | NO | Swap DSmoother to HypreBoomerAMG; add SUNDIALS iterative LS + Hypre PC |
| 2 | Fermi/ChargedFermi: concentration-dependent D, segregation, solid solubility | **MOSTLY** | NO for basic | MFEM NewtonSolver + HypreBoomerAMG LS; KINSolver for JFNK |
| 3 | CDD: multi-species reaction networks, clustering, SUPG | **NO for DAE** | YES for DAE | SUNDIALS IDA + Hypre PC for DAE; or PETSc TS if available |
| 4 | OED, TED init, dose loss, moving interface | **YES** | NO | Same as Phase 1 solver; mesh ops are MFEM |
| 5 | Polysilicon: Voronoi, dual mesh, grain boundary segregation | **YES** | NO | Pre-processing; same solver |
| 6 | SiGe/III-V: bandgap, strain, Ge-B pairing | **YES** | NO | Coefficient changes; same solver |
| 7 | KMC atomistic + continuum coupling | **YES** | NO | Sub-cycling; L2 projection with Hypre LS |
| 8 | KMC lattice epitaxy | **YES** | NO | KMC + mesh update |
| 9 | Flash/laser anneal: heat transfer, Allen-Cahn phase field, melt diffusion | **YES** | NO for basic | ARKStepSolver::IMEX or MFEM IMEX_DIRK_RK3 + Hypre LS |
| 10 | PDE API, results extraction, calibration (Sobol, MCMC, GP, adjoint inversion) | **NO** | **YES** (TAO + TSAdjoint) | PETSc TAO + TSAdjoint (rebuild required) |
| 11 | AMR: h/p refinement, ZZ estimator, moving boundary AMR | **YES** | NO | MFEM native AMR (neither Hypre nor PETSc needed) |

**Summary:** Phases 1-2, 4-9, 11: Hypre sufficient (zero build cost). Phase 3: Hypre + SUNDIALS IDA can work for DAE; PETSc TS is nicer but not required. Phase 10: PETSc TAO + TSAdjoint is the only path for adjoint-based optimization.

---

## 8. Phased Adoption Strategy

### Phase A: Immediate (Zero Build Cost) -- Adopt Hypre

**Actions:**
1. Replace `mfem::DSmoother` with `mfem::HypreBoomerAMG` in `DiffusionEngine.hpp:100-103`
2. Add SUNDIALS iterative solver + Hypre preconditioner wrapper to `SundialsTimeIntegrator.hpp`
3. Remove the 256-dof cap and subsampling logic (lines 166-178)
4. Use `HyprePCG` + `HypreBoomerAMG` as the default solver combination for SPD diffusion systems
5. For nonsymmetric systems (advection, SUPG), use `HypreGMRES` + `HypreBoomerAMG`

**Build cost:** Zero. MFEM already has Hypre linked.

**Code changes:** ~100-200 lines (new SUNDIALS PC wrapper + DiffusionEngine solver swap).

### Phase B: Near-Term (Phases 2-9) -- MFEM Native Nonlinear + IMEX

**Actions:**
1. Use `mfem::NewtonSolver` with `HypreBoomerAMG` as the linear solver for Phase 2 Fermi diffusion
2. Use `mfem::KINSolver` (SUNDIALS KINSOL wrapper) for JFNK if NewtonSolver is insufficient
3. Use `mfem::ARKStepSolver(ARKStepSolver::IMEX)` for Phase 9 Allen-Cahn
4. Use `mfem::IMEX_DIRK_RK3` as alternative IMEX if ARKStep is unavailable
5. For Phase 3 DAE: integrate SUNDIALS IDA (already linked in CMakeLists.txt:314) with Hypre PC

**Build cost:** Zero. All solvers are in the existing MFEM build.

**Risk of rewriting solver code twice:** LOW. MFEM's `Solver` interface is abstract. Whether using `HypreBoomerAMG` or `PetscLinearSolver`, the `Operator` and `Solver` interfaces are the same. When PETSc is adopted later, only the solver instantiation changes, not the assembly code.

### Phase C: Deferred (Phase 10) -- Rebuild PETSc with Hypre

**Actions (when Phase 10 is reached):**
1. Rebuild PETSc with Hypre support: `--with-hypre=1` (or vcpkg overlay)
2. Rebuild MFEM with `MFEM_USE_PETSC=ON`
3. Add `PetscTimeIntegrator.hpp` wrapping `mfem::PetscODESolver` (TS)
4. Use `TSAdjoint` for adjoint-based D(x) inversion (Task 15)
5. Use TAO for optimization (Tasks 12-15)
6. Use `PetscParMatrix` from `HypreParMatrix` (zero-copy conversion available)

**Migration path:** `HypreParMatrix` to `PetscParMatrix` is zero-copy (`petsc.hpp` confirms this). Existing assembly code using `BilinearForm` and `SparseMatrix` does not change. Only the solver backend changes.

**Why not do Phase C now?**
1. PETSc rebuild is significant effort (rebuild PETSc + MFEM)
2. Phases 1-9 do not need PETSc (Hypre + MFEM native suffices)
3. Phase 10 is the last phase -- may be months away
4. Premature adoption adds maintenance burden without immediate benefit

---

## 9. Final Recommendation

### Adopt Hypre now. Defer PETSc until Phase 10.

**Rationale:**

1. **Hypre is free.** Zero build cost. Already in MFEM. One-line change to swap `DSmoother` for `HypreBoomerAMG`.

2. **SUNDIALS + Hypre eliminates the 256-dof cap.** SUNDIALS 7.8.0 iterative solvers (`SUNLinSol_SPBCGS`, `SUNLinSol_SPGMR`) accept user-defined preconditioners (`SUNLinSolSetPreconditioner`). Wrap `HypreBoomerAMG` as `PSetup`/`Psolve` callbacks. This fixes the most critical bottleneck without PETSc.

3. **MFEM has native nonlinear solvers.** `NewtonSolver` (with Eisenstat-Walker adaptive tolerance, `solvers.hpp:780`), `KINSolver` (SUNDIALS KINSOL with JFNK, `sundials.hpp:896`), `LBFGSSolver` (`solvers.hpp:865`). No need for PETSc SNES for Phases 2-3.

4. **MFEM has native IMEX.** `IMEXRK2`, `IMEX_DIRK_RK3` (`ode.hpp:1030-1090`), `ARKStepSolver::IMEX` (`sundials.hpp:728`). No need for PETSc TSARKIMEX for Phase 9.

5. **PETSc is blocked by missing HYPRE.** The vcpkg PETSc build lacks `PETSC_HAVE_HYPRE`. MFEM's `petsc.hpp:45-47` errors out without it. Rebuilding PETSc + MFEM is a significant effort that is not justified until Phase 10.

6. **PETSc's unique value is Phase 10 only.** `TSAdjoint` (adjoint timestepping) and `TAO` (20+ optimization solvers) have no Hypre or MFEM equivalent. These are needed for Task 15 (adjoint-based D(x) inversion) and Tasks 12-14 (calibration/optimization).

7. **Migration to PETSc later is clean.** `HypreParMatrix` to `PetscParMatrix` is zero-copy. MFEM's `Operator`/`Solver` interfaces are backend-agnostic. Assembly code does not change. Only solver instantiation changes.

### Specific Boundaries

**DO adopt Hypre now for:**
1. Linear solve in DiffusionEngine -- replace `DSmoother` with `HypreBoomerAMG`
2. SUNDIALS CVODE linear solver -- wrap Hypre as SUNDIALS preconditioner
3. Any sparse linear solve in the diffusion engine

**DO use MFEM native solvers for:**
1. Nonlinear solve -- `NewtonSolver` + `HypreBoomerAMG` (Phases 2-3)
2. IMEX time integration -- `ARKStepSolver::IMEX` or `IMEX_DIRK_RK3` (Phase 9)
3. JFNK -- `KINSolver` with JFNK mode (Phase 2-3)
4. AMR -- MFEM native `ThresholdRefiner` (Phase 11)

**DO adopt PETSc later (Phase 10) for:**
1. Adjoint-based D(x) inversion -- `TSAdjoint` (Task 15)
2. Optimization -- `TAO` solvers (Tasks 12-15)
3. DAE support (if SUNDIALS IDA integration fails) -- `TS` with `TSSetIFunction`

**DO NOT adopt PETSc for:**
1. AMR -- MFEM native (h+p refinement > DMForest h-only)
2. MCMC -- external library (emcee via Python)
3. GP surrogates -- external library (scikit-learn/GPy)
4. FEM assembly -- keep MFEM BilinearForm/LinearForm

### Coexistence Model

```
VIENNAPS_HAS_MFEM + VIENNAPS_HAS_SUNDIALS + Hypre (via MFEM):
  -> Production path for Phases 1-9
  -> HypreBoomerAMG for linear solve
  -> SUNDIALS CVODE + iterative LS + Hypre PC for time integration
  -> MFEM NewtonSolver / KINSolver for nonlinear
  -> ARKStepSolver::IMEX for Allen-Cahn

VIENNAPS_HAS_MFEM + VIENNAPS_HAS_SUNDIALS + VIENNAPS_HAS_PETSC (future):
  -> Production path for Phase 10
  -> PetscODESolver (TS) for adjoint
  -> TAO for optimization
  -> Zero-copy HypreParMatrix -> PetscParMatrix conversion

VIENNAPS_HAS_MFEM only (no SUNDIALS):
  -> Implicit Euler fallback (BiCGSTAB + HypreBoomerAMG)

No MFEM:
  -> Profile-only mode
```

---

## 10. Sources

### Primary Sources (MFEM Source Code)
- MFEM Hypre wrappers: `F:\dev\mfem\linalg\hypre.hpp` (2356 lines)
  - `HypreBoomerAMG` at :1828, `HyprePCG` at :1352, `HypreGMRES` at :1444
  - `HypreFGMRES` at :1523, `HypreSmoother` at :1076, `HypreSolver` (base) at :1238
  - `HypreParMatrix` at :418, `HypreParVector` at :229
  - `HypreParaSails` at :1638, `HypreEuclid` at :1725, `HypreILU` at :1781
  - `HypreAMS` at :1988, `HypreADS` at :2065
- MFEM PETSc wrappers: `F:\dev\mfem\linalg\petsc.hpp` (1018 lines, author: Stefano Zampini)
  - `PETSC_HAVE_HYPRE` requirement at :45-47
- MFEM solvers: `F:\dev\mfem\linalg\solvers.hpp`
  - `NewtonSolver` at :780, `LBFGSSolver` at :865, `SLBQPOptimizer` at :1036
  - `CGSolver` at :626, `GMRESSolver` at :660, `BiCGSTABSolver` at :709, `BlockILU` at :1102
- MFEM ODE solvers: `F:\dev\mfem\linalg\ode.hpp`
  - `IMEXExpImplEuler` at :1030, `IMEXRK2` at :1046, `IMEX_DIRK_RK3` at :1079
- MFEM SUNDIALS wrappers: `F:\dev\mfem\linalg\sundials.hpp`
  - `CVODESolver` at :429, `ARKStepSolver` (IMEX) at :720, `KINSolver` (JFNK) at :896
- MFEM config: `F:\dev\mfem\build\config\_config.hpp`
  - `MFEM_USE_MPI` enabled, `MFEM_HYPRE_VERSION 23200`, `MFEM_USE_SUNDIALS` enabled
  - `MFEM_USE_PETSC` disabled, `MFEM_USE_CUDA` enabled

### Primary Sources (SUNDIALS Headers)
- SUNDIALS version: `F:\dev\vcpkg\installed\x64-windows\include\sundials\sundials_config.h` (7.8.0)
- SUNDIALS linear solver API: `F:\dev\vcpkg\installed\x64-windows\include\sundials\sundials_linearsolver.h`
  - `SUNLinSolSetPreconditioner` at :179-180
  - `SUNLinSolSetATimes` at :175-176
- SUNDIALS iterative solver types: `F:\dev\vcpkg\installed\x64-windows\include\sunlinsol\`
  - `sunlinsol_spgmr.h` (SPGMR, :88), `sunlinsol_spbcgs.h` (BiCGStab, :81)
  - `sunlinsol_spfgmr.h` (FGMRES), `sunlinsol_sptfqmr.h` (TFQMR), `sunlinsol_pcg.h` (PCG)
- SUNDIALS callback types: `F:\dev\vcpkg\installed\x64-windows\include\sundials\sundials_iterative.h`
  - `SUNATimesFn` at :108, `SUNPSetupFn` at :120, `SUNPSolveFn` at :148
- CVODE linear solver interface: `F:\dev\vcpkg\installed\x64-windows\include\cvode\cvode_ls.h`
  - `CVodeSetLinearSolver` at :82, `CVodeSetPreconditioner` at :97-98

### Primary Sources (PETSc Configuration)
- vcpkg PETSc config: `F:\dev\vcpkg\installed\x64-windows\include\petscconf.h`
  - `PETSC_HAVE_MPIUNI 1` (MPI stub, not real MPI)
  - **NO** `PETSC_HAVE_HYPRE` (Hypre support absent)
  - `PETSC_HAVE_PACKAGES ":blaslapack:mathlib:mpi:"` (minimal)
  - 74 total `PETSC_HAVE` defines
- Prebuilt PETSc config: `F:\dev\petsc\arch-win32-c-opt\include\petscconf.h` (identical minimal build)
- PETSc Hypre integration: `F:\dev\vcpkg\installed\x64-windows\include\petscpc.h`
  - `PCHYPRESetType` at :292, `PCHYPREGetType` at :293
- PETSc Hypre type: `F:\dev\vcpkg\installed\x64-windows\include\petscpctypes.h`
  - `PCHYPRE "hypre"` at :52

### ViennaPS Source Code
- `include/viennaps/fields/DiffusionEngine.hpp` (134 lines) -- uses `DSmoother` + `BiCGSTABSolver` at :100-105
- `include/viennaps/fields/SundialsTimeIntegrator.hpp` (314 lines) -- 256-dof cap at :167, dense LS at :241-251
- `include/viennaps/fields/AmgclSolver.hpp` (195 lines) -- standalone AMGCL fallback, not integrated with MFEM

### Prior Research
- PETSc feasibility report: `docs/superpowers/research/2026-07-22-petsc-feasibility.md` (399 lines)

### Plan Documents
- Phase 1-11 plans: `docs/superpowers/plans/2026-07-20-diffusion-phase{1-11}.md`
