# ADR-0003: Adopt Hypre as the Linear Solver Backend

**Date:** 2026-07-22

## Decision

Adopt **Hypre** (via MFEM's existing wrappers) as the linear solver backend for all diffusion engine phases. Defer PETSc until Phase 10 (optimization/adjoint only).

## Context

The diffusion engine plans (Phases 1-11) require scalable sparse linear solvers for:
- Implicit time integration (CVODE/IDA linear solves)
- Nonlinear Newton iterations (concentration-dependent diffusivity)
- L2 projection field transfers (KMC coupling, AMR, moving boundary)
- Block systems (multi-species CDD)

The existing code uses `mfem::BiCGSTABSolver` + `mfem::DSmoother` (no AMG) and SUNDIALS CVODE with dense O(n^3) linear solver capped at 256 dofs. This is a showstopper for production-scale meshes (a 32x32 2D mesh already exceeds the cap).

Three options were evaluated (full analysis in `docs/superpowers/research/2026-07-22-hypre-deep-dive.md`):

1. **Hypre alone** (via MFEM wrappers) -- zero build cost, fixes 3.5 of 5 bottlenecks
2. **PETSc** -- requires rebuilding both PETSc (vcpkg build lacks `PETSC_HAVE_HYPRE`) and MFEM
3. **Hypre now, PETSc later** -- phased adoption

## Rationale

**Hypre is already available.** MFEM is built with `MFEM_USE_MPI` and `MFEM_HYPRE_VERSION 23200` (`_config.hpp:47,201`). MFEM's `linalg/hypre.hpp` (2356 lines) exposes 15 solver classes including `HypreBoomerAMG`, `HyprePCG`, `HypreGMRES`, `HypreFGMRES`, `HypreILU`. No rebuild needed.

**Serial-to-parallel bridge exists.** The engine assembles serial `SparseMatrix` but `HypreBoomerAMG` requires `HypreParMatrix`. The constructor `HypreParMatrix(MPI_COMM_SELF, glob_size, row_starts, row_starts, &sparseMat)` (`hypre.hpp:567-569`) bridges this with zero-copy (when `HYPRE_BIGINT` is undefined, verified `HYPRE_config.h:24`). No `ParFiniteElementSpace` or `ParMesh` needed.

**SUNDIALS + Hypre eliminates the 256-dof cap.** SUNDIALS 7.8.0 provides iterative solvers (`SUNLinSol_SPBCGS`, `SUNLinSol_SPGMR`) that accept user-defined preconditioners via `SUNLinSolSetPreconditioner()` (`sundials_linearsolver.h:179`). `HypreBoomerAMG::Mult()` is wrapped as the `Psolve` callback. MFEM's `CVODESolver` wrapper uses `SUN_PREC_NONE` (`sundials.cpp:905`) so the raw SUNDIALS C API is used (as `SundialsTimeIntegrator.hpp` already does).

**MFEM provides complementary native solvers** -- no PETSc needed for Phases 1-9:
- `NewtonSolver` (`solvers.hpp:780`) with `HypreBoomerAMG` as inner solver -- nonlinear (Phases 2-3)
- `KINSolver` (`sundials.hpp:896`) with JFNK -- matrix-free Newton (Phase 2-3)
- `ARKStepSolver::IMEX` (`sundials.hpp:720,728`) -- Allen-Cahn (Phase 9)
- `IMEX_DIRK_RK3` (`ode.hpp:1079`) -- pure MFEM IMEX alternative

**PETSc is deferred to Phase 10.** The vcpkg PETSc build lacks `PETSC_HAVE_HYPRE` (`petscconf.h`: only 74 defines, `MPIUNI` stub). MFEM's `petsc.hpp:45-47` `#error`s without it. Rebuilding both PETSc and MFEM is not justified until Phase 10 Tasks 12-15 (adjoint D(x) inversion via `TSAdjoint`, optimization via `TAO`) -- the only capabilities with no Hypre or MFEM equivalent.

## What Hypre Fixes

| Bottleneck | Status | How |
|---|---|---|
| Dense O(n^3) SUNDIALS cap at 256 dofs | **Fixed** | `SUNLinSol_SPBCGS` + `HypreBoomerAMG` preconditioner |
| No sparse solver for neq > 64 | **Fixed** | `HyprePCG`/`HypreGMRES` + `HypreBoomerAMG` |
| No nonlinear solve (hand-rolled Picard) | **Partial** | MFEM `NewtonSolver` + `HypreBoomerAMG` as inner solver |
| No DAE support (IDA linked, unused) | **Not fixed** | Requires raw SUNDIALS IDA C API (MFEM has no IDA wrapper) |
| No adjoint/optimization | **Not fixed** | Requires PETSc TAO/TSAdjoint (Phase 10 only) |

## What Hypre Does NOT Provide

- Nonlinear solvers (use MFEM `NewtonSolver` / `KINSolver` instead)
- Time integrators (use SUNDIALS CVODE/IDA/ARKode instead)
- DAE support (use raw SUNDIALS IDA C API instead)
- Optimization (defer to PETSc TAO, Phase 10)
- Adjoint/sensitivity (defer to PETSc TSAdjoint, Phase 10)
- GPU acceleration (vcpkg Hypre is CPU-only: `HYPRE_USING_CUDA` undef, `HYPRE_config.h:90`)

## Sources

- MFEM Hypre wrappers: `F:\dev\mfem\linalg\hypre.hpp` (2356 lines)
- MFEM config: `F:\dev\mfem\build\config\_config.hpp` (`MFEM_USE_MPI:47`, `MFEM_HYPRE_VERSION 23200:201`)
- Hypre config: `F:\dev\vcpkg\installed\x64-windows\include\HYPRE_config.h` (v2.32.0, CPU-only)
- SUNDIALS headers: `F:\dev\vcpkg\installed\x64-windows\include\sundials\` (v7.8.0)
- Full analysis: `docs/superpowers/research/2026-07-22-hypre-deep-dive.md` (837 lines)
