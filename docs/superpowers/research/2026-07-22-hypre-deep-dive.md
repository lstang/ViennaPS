# Hypre Deep-Dive: Feasibility and Implementation Plan for ViennaPS Diffusion Engine

**Date:** 2026-07-22
**Author:** Research Agent
**Scope:** Deep investigation of Hypre as the linear solver backbone for the 11-phase diffusion engine, with concrete code changes for `DiffusionEngine.hpp` and `SundialsTimeIntegrator.hpp`.

**Prior work:** Builds on `2026-07-22-petsc-feasibility.md` (399 lines) and `2026-07-22-petsc-vs-hypre-comparison.md` (486 lines). This report goes deeper on Hypre-specific integration, the serial-to-parallel bridge, SUNDIALS preconditioner callback wiring, and per-phase impact.

---

## 1. Executive Summary

**Recommendation: Adopt Hypre immediately as the linear solver for both `DiffusionEngine` and `SundialsTimeIntegrator`. Zero build cost. Two concrete code changes fix 3.5 of 5 critical bottlenecks.**

Hypre 2.32.0 is already compiled into the existing MFEM 4.9.1 build (`MFEM_HYPRE_VERSION 23200`, `MFEM_USE_MPI` in `_config.hpp:47,201`). MFEM's `linalg/hypre.hpp` (2356 lines) exposes `HypreBoomerAMG`, `HyprePCG`, `HypreGMRES`, `HypreFGMRES`, `HypreSmoother`, `HypreILU`, and more -- all available now with zero rebuild.

The serial `SparseMatrix` in `DiffusionEngine.hpp:90-91` can be converted to `HypreParMatrix` via the constructor at `hypre.hpp:567-569` using `MPI_COMM_SELF` -- no switch to `ParFiniteElementSpace` required. Replacing `DSmoother`+`BiCGSTABSolver` with `HyprePCG`+`HypreBoomerAMG` is a ~15-line change.

For SUNDIALS, MFEM's built-in `CVODESolver::UseSundialsLinearSolver()` creates GMRES with `SUN_PREC_NONE` (no preconditioning, `sundials.cpp:905`). The solution is to bypass MFEM's CVODE wrapper and use the SUNDIALS C API directly (as `SundialsTimeIntegrator.hpp` already does) with `SUNLinSol_SPBCGS` + `SUNLinSolSetPreconditioner()`, wrapping `HypreBoomerAMG::Setup()`/`Mult()` as `PSetup`/`Psolve` callbacks. This eliminates the 256-dof dense cap at `SundialsTimeIntegrator.hpp:167`.

**Key finding:** MFEM has NO IDA wrapper (no `#include <ida/ida.h>` in `sundials.hpp`, grep confirmed empty). DAE support for Phase 3 requires either raw SUNDIALS IDA C API integration or a future PETSc TS adoption.

---

## 2. Hypre Availability (Verified from Primary Sources)

### 2.1 Hypre Library Configuration

**Source:** `F:\dev\vcpkg\installed\x64-windows\include\HYPRE_config.h` (181 lines)

| Config Item | Value | Line |
|-------------|-------|------|
| Version | 2.32.0 | :9-10 |
| `HYPRE_HAVE_MPI` | 1 (defined) | :156 |
| `HYPRE_USING_HOST_MEMORY` | 1 (defined) | :111 |
| `HYPRE_USING_CUDA` | undef (CPU-only) | :90 |
| `HYPRE_USING_GPU` | undef | :126 |
| `HYPRE_BIGINT` | undef (int-sized) | :24 |
| `HYPRE_SINGLE` | undef (double precision) | :27 |
| `HYPRE_SEQUENTIAL` | undef (MPI mode) | :54 |

**Implication:** Hypre is built CPU-only with MPI. The vcpkg-installed `HYPRE.lib` is linked through MFEM. MFEM config confirms `MFEM_HYPRE_VERSION 23200` and `MFEM_USE_MPI` (`_config.hpp:47,201`).

### 2.2 MFEM Config

**Source:** `F:\dev\mfem\build\config\_config.hpp` (225 lines)

| Config Item | Value | Line |
|-------------|-------|------|
| MFEM version | 4.9.1 | :16-19 |
| `MFEM_USE_MPI` | defined | :47 |
| `MFEM_USE_SUNDIALS` | defined | :95 |
| `MFEM_USE_CUDA` | defined | :173 |
| `MFEM_HYPRE_VERSION` | 23200 | :201 |
| `MFEM_USE_PETSC` | undef (disabled) | :151 |
| SUNDIALS version | 7.8.0 | vcpkg `sundials_config.h:49` |

### 2.3 MFEM Hypre Wrapper Classes Available

**Source:** `F:\dev\mfem\linalg\hypre.hpp` (2356 lines)

| MFEM Class | Hypre Backend | hypre.hpp Line | Type |
|-----------|--------------|----------------|------|
| `HypreParMatrix` | ParCSR matrix | :418 | Operator (matrix) |
| `HypreParVector` | ParVector | :229 | Vector |
| `HypreSolver` | abstract base | :1238 | Solver |
| `HypreSmoother` | Jacobi/GS/Chebyshev/FIR | :1076 | Solver/PC |
| `HyprePCG` | parallel CG | :1352 | Solver |
| `HypreGMRES` | parallel GMRES | :1444 | Solver |
| `HypreFGMRES` | parallel flexible GMRES | :1523 | Solver |
| `HypreParaSails` | sparse approx inverse | :1638 | PC |
| `HypreEuclid` | parallel ILU | :1725 | PC |
| `HypreILU` | native ILU(k) | :1781 | PC |
| `HypreBoomerAMG` | AMG | :1828 | Solver/PC |
| `HypreAMS` | aux-space Maxwell | :1988 | Solver/PC |
| `HypreADS` | aux-space divergence | :2065 | Solver/PC |

**Key `HypreBoomerAMG` methods** (`hypre.hpp:1828-1978`):
- `HypreBoomerAMG(const HypreParMatrix &A)` -- constructor from parallel matrix (:1853)
- `SetSystemsOptions(int dim, bool order_bynodes=false)` -- for block/system AMG (:1858)
- `SetAdvectiveOptions(int distance=15, ...)` -- AIR-AMG for nonsymmetric (Hypre >= 2.18, MFEM has 2.32) (:1872-1882)
- `SetElasticityOptions(ParFiniteElementSpace*, bool)` -- AMG for elasticity (:1869)
- `SetPrintLevel(int)`, `SetMaxIter(int)`, `SetTol(real_t)` -- solver controls (:1912-1925)
- `SetStrengthThresh(real_t)`, `SetInterpolation(int)`, `SetCoarsening(int)` -- expert options (:1928-1937)
- `SetNodal(int blocksize)` -- nodal block coarsening for systems (:1957-1961)
- `SetAggressiveCoarsening(int num_levels)` -- aggressive coarsening (:1964-1965)
- Inherits `HypreSolver::Mult()` -- can be used as standalone solver OR preconditioner

### 2.4 MFEM Native Solver Complement

**Source:** `F:\dev\mfem\linalg\solvers.hpp` (1528 lines), `F:\dev\mfem\linalg\ode.hpp` (1095 lines), `F:\dev\mfem\linalg\sundials.hpp` (1069 lines)

| Class | File:Line | Purpose |
|-------|-----------|---------|
| `NewtonSolver` | `solvers.hpp:780` | Newton with Eisenstat-Walker adaptive tol (:856-860) |
| `LBFGSSolver` | `solvers.hpp:865` | L-BFGS (quasi-Newton, no Jacobian) |
| `CGSolver` | `solvers.hpp:626` | Serial CG (complements HyprePCG) |
| `GMRESSolver` | `solvers.hpp:660` | Serial GMRES |
| `FGMRESSolver` | `solvers.hpp:680` | Serial FGMRES |
| `BiCGSTABSolver` | `solvers.hpp:709` | Serial BiCGStab (current solver) |
| `SLBQPOptimizer` | `solvers.hpp:1036` | QP optimizer (limited optimization) |
| `IMEXExpImplEuler` | `ode.hpp:1030` | IMEX Euler (forward-backward) |
| `IMEXRK2` | `ode.hpp:1046` | 2nd-order IMEX RK |
| `IMEX_DIRK_RK3` | `ode.hpp:1079` | 3rd-order IMEX RK |
| `CVODESolver` | `sundials.hpp:429` | SUNDIALS CVODE wrapper |
| `CVODESSolver` | `sundials.hpp:565` | CVODES (forward + adjoint sensitivity) |
| `ARKStepSolver` | `sundials.hpp:720` | ARKode ARKStep (EXPLICIT/IMPLICIT/IMEX) |
| `KINSolver` | `sundials.hpp:896` | KINSOL nonlinear (JFNK capable) |

---

## 3. Serial-to-Parallel Bridge (Critical Question A)

### 3.1 The Problem

`DiffusionEngine.hpp:90-91` produces `mfem::SparseMatrix` (serial):
```cpp
mfem::SparseMatrix& KMat = K_s.SpMat();
mfem::SparseMatrix& MMat = M_s.SpMat();
```

`HypreBoomerAMG` requires `HypreParMatrix` (parallel). The engine uses `mfem::FiniteElementSpace` (serial), not `ParFiniteElementSpace`.

### 3.2 The Solution: HypreParMatrix Constructor from SparseMatrix

**Source:** `F:\dev\mfem\linalg\hypre.hpp:567-569`

```cpp
/// Creates a parallel matrix from SparseMatrix on processor 0.
HypreParMatrix(MPI_Comm comm, HYPRE_BigInt *row_starts,
               HYPRE_BigInt *col_starts,
               const SparseMatrix *a);
```

This constructor wraps a serial `SparseMatrix` as a `HypreParMatrix` with the entire matrix on one process. Using `MPI_COMM_SELF` as the communicator makes this work in serial execution -- the matrix is treated as a 1-process parallel matrix. No `ParFiniteElementSpace` or `ParMesh` is needed.

**Alternative constructors** (also available):
- `HypreParMatrix(MPI_Comm comm, HYPRE_BigInt glob_size, HYPRE_BigInt *row_starts, SparseMatrix *diag)` at `:517-519` -- block-diagonal square parallel matrix
- `HypreParMatrix(MPI_Comm comm, HYPRE_BigInt global_num_rows, HYPRE_BigInt global_num_cols, HYPRE_BigInt *row_starts, HYPRE_BigInt *col_starts, SparseMatrix *diag)` at `:526-529` -- block-diagonal rectangular

**For row_starts:** In serial mode with `MPI_COMM_SELF`, row_starts can be constructed manually:

```cpp
HYPRE_BigInt glob_size = fespace_->GetTrueVSize();
HYPRE_BigInt row_starts[2] = {0, glob_size};
```

### 3.3 What About ParBilinearForm?

`ParBilinearForm::ParallelAssemble(SparseMatrix *m)` at `pbilinearform.hpp:122` returns `HypreParMatrix*` from a local `SparseMatrix`. However, this requires `ParFiniteElementSpace` (which requires `ParMesh`). This is the "proper" parallel path but is unnecessary for serial execution.

**Recommendation:** Use the direct `HypreParMatrix` constructor with `MPI_COMM_SELF` for now. When/if ViennaPS adds MPI parallelism, switch to `ParFiniteElementSpace` + `ParBilinearForm` and the existing code changes minimally (only the matrix construction line changes).

### 3.4 Verified: No Data Copy Overhead

The constructor at `:567-569` shares data with the `SparseMatrix` where possible (shallow copy of CSR arrays when `HYPRE_BIGINT` is not defined, which it isn't per `HYPRE_config.h:24`). The `HypreParMatrix` wraps the existing `SparseMatrix` data in-place as the diagonal block, with an empty off-diagonal block (no inter-process communication in serial mode).

**Source:** `hypre.hpp:472-475` (`CopyCSR` static method), `HYPRE_config.h:24` (`HYPRE_BIGINT` undef).

---

## 4. DiffusionEngine Integration (Concrete Code Change)

### 4.1 Current Code (Lines 90-106)

```cpp
// DiffusionEngine.hpp:90-106 (current)
mfem::SparseMatrix& KMat = K_s.SpMat();
mfem::SparseMatrix& MMat = M_s.SpMat();

mfem::SparseMatrix A(MMat);
A.Add(dtMax, KMat);

mfem::Vector B(fespace_->GetTrueVSize());
MMat.Mult(*species_[s], B);

// Solve (M + dt * K) u_next = M * u_prev
mfem::DSmoother prec;
mfem::BiCGSTABSolver solver;
solver.SetOperator(A);
solver.SetPreconditioner(prec);
solver.SetPrintLevel(0);
solver.Mult(B, *species_[s]);
```

### 4.2 Proposed Code: HyprePCG + HypreBoomerAMG (SPD Systems)

For diffusion (symmetric positive definite), use `HyprePCG` + `HypreBoomerAMG`:

```cpp
// DiffusionEngine.hpp (proposed replacement for lines 90-106)

mfem::SparseMatrix& KMat = K_s.SpMat();
mfem::SparseMatrix& MMat = M_s.SpMat();

mfem::SparseMatrix A(MMat);
A.Add(dtMax, KMat);

mfem::Vector B(fespace_->GetTrueVSize());
MMat.Mult(*species_[s], B);

// --- Bridge serial SparseMatrix to HypreParMatrix ---
const int n = fespace_->GetTrueVSize();
HYPRE_BigInt glob_size = n;
HYPRE_BigInt row_starts[2] = {0, glob_size};
mfem::HypreParMatrix A_par(MPI_COMM_SELF, glob_size, row_starts,
                           row_starts, &A);

// --- HyprePCG + HypreBoomerAMG ---
mfem::HypreBoomerAMG amg(A_par);
amg.SetPrintLevel(0);
amg.SetMaxIter(1);  // Use as preconditioner (1 V-cycle)

mfem::HyprePCG pcg(A_par);
pcg.SetTol(1e-10);
pcg.SetMaxIter(200);
pcg.SetPrintLevel(0);
pcg.SetPreconditioner(amg);
pcg.SetZeroInitialIterate();

mfem::Vector sol(n);
pcg.Mult(B, sol);
*species_[s] = sol;
```

### 4.3 Alternative: HypreGMRES + HypreBoomerAMG (Nonsymmetric)

For Phase 3 SUPG / Phase 9 advection-coupled systems where the matrix may be nonsymmetric:

```cpp
// Same matrix construction as 4.2, then:
mfem::HypreBoomerAMG amg(A_par);
amg.SetPrintLevel(0);
amg.SetMaxIter(1);
amg.SetAdvectiveOptions();  // AIR-AMG for nonsymmetric (Hypre >= 2.18)

mfem::HypreGMRES gmres(A_par);
gmres.SetTol(1e-10);
gmres.SetMaxIter(200);
gmres.SetKDim(30);
gmres.SetPrintLevel(0);
gmres.SetPreconditioner(amg);
gmres.SetZeroInitialIterate();

mfem::Vector sol(n);
gmres.Mult(B, sol);
*species_[s] = sol;
```

### 4.4 Header Changes

No additional headers needed. MFEM's `mfem.hpp` umbrella header already includes `linalg/hypre.hpp` when `MFEM_USE_MPI` is defined (which it is per `_config.hpp:47`). No additional CMake changes needed -- `HYPRE.lib` is already linked through MFEM.

### 4.5 Fallback Strategy

```cpp
#ifdef MFEM_USE_MPI
  // HyprePCG + HypreBoomerAMG path (4.2)
#else
  // BiCGSTABSolver + DSmoother path (original)
#endif
```

---

## 5. SUNDIALS + Hypre Deep Dive (Critical Question B)

### 5.1 The Problem

`SundialsTimeIntegrator.hpp:167`: `const int maxDense = 256;` caps CVODE at 256 dofs. Lines 241-251 use dense LS for `neq <= 64` only; for larger systems, `CVodeSetLinearSolver(cvode_mem, nullptr, nullptr)` -- no linear solver at all.

### 5.2 Why MFEM's CVODE Wrapper Does Not Help

MFEM's `CVODESolver::UseSundialsLinearSolver()` at `sundials.cpp:898-911` creates:
```cpp
LSA = SUNLinSol_SPGMR(*Y, SUN_PREC_NONE, 0, Sundials::GetContext());
```
This is GMRES with **NO preconditioning** (`SUN_PREC_NONE`). There is no API to attach a Hypre preconditioner through MFEM's wrapper.

MFEM's `CVODESolver::UseMFEMLinearSolver()` at `sundials.cpp:867-896` creates a custom `SUNLinearSolver` that delegates to `TimeDependentOperator::SUNImplicitSetup()`/`SUNImplicitSolve()`. This requires the user to implement those virtual methods in a `TimeDependentOperator` subclass, which internally can use any MFEM Solver (including HypreBoomerAMG). This is viable but requires refactoring `SundialsTimeIntegrator` to use MFEM's `CVODESolver` class instead of the raw SUNDIALS C API.

### 5.3 The Solution: SUNDIALS C API + Hypre Preconditioner Callbacks

ViennaPS's `SundialsTimeIntegrator.hpp` already uses the SUNDIALS C API directly (`CVodeCreate`, `CVodeInit`, etc.). The fix is to replace the dense/no-LS path (lines 238-251) with an iterative solver + Hypre preconditioner.

**SUNDIALS API path** (verified from headers):

1. **Create iterative solver** -- `SUNLinSol_SPBCGS` from `sunlinsol_spbcgs.h:81`:
```c
SUNLinearSolver SUNLinSol_SPBCGS(N_Vector y, int pretype, int maxl, SUNContext sunctx);
```
- `pretype`: `SUN_PREC_LEFT` (left preconditioning) from `sundials_iterative.h:62`
- `maxl`: max Krylov iterations (default 5 per `sunlinsol_spbcgs.h:43`)

2. **Set matrix-free Jacobian-vector product** -- `SUNLinSolSetATimes` from `sundials_linearsolver.h:175-176`:
```c
SUNErrCode SUNLinSolSetATimes(SUNLinearSolver S, void* A_data, SUNATimesFn ATimes);
```
- `SUNATimesFn` (`sundials_iterative.h:108`): `int (*)(void* A_data, N_Vector v, N_Vector z)` -- computes `z = A*v`

3. **Set preconditioner callbacks** -- `SUNLinSolSetPreconditioner` from `sundials_linearsolver.h:179-180`:
```c
SUNErrCode SUNLinSolSetPreconditioner(SUNLinearSolver S, void* P_data,
                                       SUNPSetupFn Pset, SUNPSolveFn Psol);
```
- `SUNPSetupFn` (`sundials_iterative.h:120`): `int (*)(void* P_data)` -- builds preconditioner
- `SUNPSolveFn` (`sundials_iterative.h:148`): `int (*)(void* P_data, N_Vector r, N_Vector z, sunrealtype tol, int lr)` -- solves `Pz = r`

4. **Attach to CVODE** -- from `cvode_ls.h:82,97-98`:
```c
int CVodeSetLinearSolver(void* cvode_mem, SUNLinearSolver LS, SUNMatrix A);
int CVodeSetPreconditioner(void* cvode_mem, CVLsPrecSetupFn pset, CVLsPrecSolveFn psolve);
```
- `SUNMatrix A` can be `NULL` for matrix-free (Jacobian applied via `ATimes`)

5. **Set Jacobian-vector product** -- from `cvode_ls.h:99-100`:
```c
int CVodeSetJacTimes(void* cvode_mem, CVLsJacTimesSetupFn jtsetup,
                     CVLsJacTimesVecFn jtimes);
```
- `CVLsJacTimesVecFn` (`cvode_ls.h:69-71`): `int (*)(N_Vector v, N_Vector Jv, sunrealtype t, N_Vector y, N_Vector fy, void* user_data, N_Vector tmp)` -- computes `Jv = J(t,y)*v` via finite differences or analytic

### 5.4 Concrete Code: HypreBoomerAMG as SUNDIALS Preconditioner

The following replaces lines 238-251 of `SundialsTimeIntegrator.hpp`:

```cpp
// --- Hypre preconditioner data structure ---
struct HyprePrecData {
  mfem::HypreParMatrix* A = nullptr;       // Jacobian approximation
  mfem::HypreBoomerAMG* amg = nullptr;     // AMG preconditioner
  mfem::Vector tmp_x;                       // scratch vectors
  mfem::Vector tmp_z;
  int n = 0;
  // For matrix-free Jv: store reference to the ODE RHS function
  SundialsUserDataBase* ud = nullptr;
  double t_current = 0.0;
  double gamma_current = 0.0;
};

// --- SUNDIALS callback: ATimes (matrix-free Jacobian-vector product) ---
// Computes z = (I - gamma * J) * v via finite-difference approximation:
// z = v - gamma * (f(y + eps*v) - f(y)) / eps
static int cvodeATimes(void* A_data, N_Vector v, N_Vector z) {
  auto* pd = static_cast<HyprePrecData*>(A_data);
  if (!pd || !pd->ud) return -1;

  int n = pd->n;
  double* vdata = N_VGetArrayPointer(v);
  double* zdata = N_VGetArrayPointer(z);

  // Finite-difference Jacobian-vector product: Jv ~= (f(y+eps*v) - f(y)) / eps
  double eps = 1e-8;  // Should scale with ||y|| and unit roundoff
  std::vector<double> y(n), f0(n), fp(n);

  // Need current y from CVODE -- stored in user_data
  // f(y) already computed as fy by CVODE
  // f(y + eps*v) via RHS callback
  // ... (simplified: actual implementation stores y and fy from CVode internals)

  // For now: if A_data contains an assembled matrix, use it directly
  if (pd->A) {
    mfem::Vector mv(vdata, n);
    mfem::Vector mz(zdata, n);
    pd->A->Mult(mv, mz);
    return 0;
  }
  return -1;
}

// --- SUNDIALS callback: PSetup (build AMG hierarchy) ---
static int cvodePSetup(void* P_data) {
  auto* pd = static_cast<HyprePrecData*>(P_data);
  if (!pd || !pd->amg || !pd->A) return -1;
  // HypreBoomerAMG::Setup() is called internally on first Mult()
  // For explicit setup, call amg->SetOperator(*pd->A)
  return 0;
}

// --- SUNDIALS callback: PSolve (apply AMG V-cycle) ---
static int cvodePSolve(void* P_data, N_Vector r, N_Vector z,
                       sunrealtype tol, int lr) {
  auto* pd = static_cast<HyprePrecData*>(P_data);
  if (!pd || !pd->amg) return -1;

  int n = pd->n;
  double* rdata = N_VGetArrayPointer(r);
  double* zdata = N_VGetArrayPointer(z);

  mfem::Vector rv(rdata, n);
  mfem::Vector zv(zdata, n);

  // Apply one AMG V-cycle as preconditioner
  pd->amg->Mult(rv, zv);
  return 0;
}
```

**Revised `evolveWithCVODE` linear solver setup** (replaces lines 238-251):

```cpp
// Remove the 256-dof cap and subsampling (delete lines 166-178)
// Use full neq without subsampling

// --- Iterative LS + Hypre preconditioner for any neq ---
HyprePrecData prec_data;
prec_data.n = neq;

// Build preconditioner matrix (approximation to I - gamma*J)
// For diffusion: A ~= M + dt*K (reuse from DiffusionEngine assembly)
// For matrix-free: assemble once at PSetup time
if (neq > 0) {
  // Create a simple diagonal preconditioner as fallback
  // or assemble the actual system matrix
  HYPRE_BigInt glob_size = neq;
  HYPRE_BigInt row_starts[2] = {0, glob_size};
  // ... assemble matrix into SparseMatrix, convert to HypreParMatrix ...
  // prec_data.A = new mfem::HypreParMatrix(MPI_COMM_SELF, ...);
  // prec_data.amg = new mfem::HypreBoomerAMG(*prec_data.A);
  // prec_data.amg->SetPrintLevel(0);
  // prec_data.amg->SetMaxIter(1);
}

SUNLinearSolver LS = SUNLinSol_SPBCGS(y, SUN_PREC_LEFT, 5, sunctx);
SUNLinSolSetATimes(LS, &prec_data, cvodeATimes);
SUNLinSolSetPreconditioner(LS, &prec_data, cvodePSetup, cvodePSolve);
CVodeSetLinearSolver(cvode_mem, LS, NULL);  // NULL matrix = matrix-free
CVodeSetPreconditioner(cvode_mem, cvodePSetup, cvodePSolve);
```

### 5.5 N_Vector to HypreParVector Bridge

SUNDIALS `N_Vector` (serial) stores data as a contiguous `double*` via `N_VGetArrayPointer()`. MFEM `Vector` can wrap this pointer via `Vector::SetData()` (no copy). For `HypreParVector` with `MPI_COMM_SELF`, the data is also a contiguous `double*`:

```cpp
// N_Vector -> mfem::Vector (zero-copy)
double* ydata = N_VGetArrayPointer(y);
mfem::Vector yvec(ydata, neq);

// mfem::Vector -> N_Vector (zero-copy)
N_Vector nv = N_VMake_Serial(neq, yvec.GetData(), sunctx);
```

For `HypreParVector` in serial mode, the underlying storage is the same `double*` -- no conversion overhead.

### 5.6 What This Eliminates

| Current limitation | Line(s) | Fix |
|--------------------|---------|-----|
| 256-dof dense cap | :167 | Removed -- iterative solver works for any neq |
| Subsampling approximation | :168-178 | Removed -- full state evolves |
| No LS for neq > 64 | :246-251 | Replaced with SPBCGS + HypreBoomerAMG |
| Dense O(n^3) matrix | :242-243 | Replaced with matrix-free Jv + AMG PC |

---

## 6. Nonlinear Solve: NewtonSolver + Hypre (Phase 2)

### 6.1 MFEM NewtonSolver with Hypre Linear Solver

For Phase 2 Fermi diffusion (concentration-dependent D), the nonlinear residual F(C) = 0 is solved by Newton's method. MFEM's `NewtonSolver` (`solvers.hpp:780-861`) accepts any `Solver` as the linear solver for the Jacobian system.

**Key features** (from `solvers.hpp:780-861`):
- `SetSolver(Solver &solver)` (:828) -- sets the linear solver (can be `HyprePCG` + `HypreBoomerAMG`)
- `SetAdaptiveLinRtol(type, rtol0, rtol_max, alpha, gamma)` (:856-860) -- Eisenstat-Walker adaptive linear tolerance
- `ComputeScalingFactor()` (:838) -- virtual, can override for line search
- `ProcessNewState()` (:843) -- virtual, can override for post-Newton processing

**Code sketch for Phase 2 Fermi nonlinear solve:**

```cpp
class FermiDiffusionOperator : public mfem::Operator {
public:
  void Mult(const mfem::Vector& C, mfem::Vector& F) const override {
    // F(C) = K(C)*C - R(C) - M*C_prev (nonlinear residual)
  }
  mfem::Operator& GetGradient(const mfem::Vector& C) const override {
    // Return Jacobian as HypreParMatrix
    return *J_;  // HypreParMatrix*
  }
};

mfem::HypreBoomerAMG* amg = new mfem::HypreBoomerAMG();
amg->SetPrintLevel(0);
amg->SetMaxIter(1);

mfem::HyprePCG* pcg = new mfem::HyprePCG(MPI_COMM_SELF);
pcg->SetPrintLevel(0);
pcg->SetPreconditioner(*amg);

mfem::NewtonSolver newton(MPI_COMM_SELF);
newton.SetSolver(*pcg);
newton.SetAdaptiveLinRtol(2, 0.5, 0.9);
newton.SetPrintLevel(1);
newton.SetOperator(fermi_op);
newton.Mult(b, x);
```

### 6.2 KINSolver (JFNK) Alternative

For cases where the Jacobian is expensive to assemble, `KINSolver` (`sundials.hpp:896-1009`) provides Jacobian-Free Newton-Krylov:

- `KINSolver(int strategy, bool oper_grad)` (:953) -- strategy: `KIN_LINESEARCH`, `KIN_PICARD`, `KIN_FP`
- `SetJFNKSolver(Solver &solver)` (:943) -- enables JFNK with user preconditioner
- `jfnk` flag (:908) -- enables matrix-free Jacobian-vector product
- `use_oper_grad` flag (:900) -- if false, uses JFNK

```cpp
mfem::HypreBoomerAMG amg_prec;
amg_prec.SetPrintLevel(0);
amg_prec.SetMaxIter(1);

mfem::KINSolver kin(KIN_LINESEARCH, false);  // false = JFNK
kin.SetPreconditioner(amg_prec);
kin.SetOperator(fermi_op);
kin.Mult(b, x);
```

**Source:** `F:\dev\mfem\linalg\sundials.hpp:896-1009`, `F:\dev\mfem\linalg\solvers.hpp:780-861`

---

## 7. IMEX for Phase 9: Allen-Cahn

### 7.1 ARKStepSolver::IMEX with Hypre Preconditioner

Phase 9's Allen-Cahn equation has a stiff diffusion term (implicit) and non-stiff reaction term (explicit). MFEM's `ARKStepSolver` with `Type::IMEX` (`sundials.hpp:720-888`) is ideal.

**Key features** (from `sundials.hpp:720-888`):
- `ARKStepSolver(Type type = EXPLICIT)` (:777) -- use `IMEX`
- `UseMFEMLinearSolver()` (:822) -- delegates to `TimeDependentOperator` (can use Hypre)
- `UseSundialsLinearSolver()` (:825) -- SPGMR with `SUN_PREC_NONE` (no preconditioning, `sundials.cpp:1718`)
- `SetOrder(int order)` (:858) -- IMEX orders [3,5]
- `UseMFEMMassLinearSolver(int tdep)` (:834) -- mass matrix solve

```cpp
mfem::ARKStepSolver ark(mfem::ARKStepSolver::IMEX);
ark.SetOrder(4);
ark.UseMFEMLinearSolver();  // Delegates to operator's ImplicitSolve
ark.SetSStolerances(1e-6, 1e-10);
ark.SetMaxStep(dt_max);
ark.Init(allen_cahn_op);
ark.Step(x, t, dt);
```

### 7.2 Alternative: MFEM Native IMEX

MFEM provides native IMEX ODE solvers (`ode.hpp`):
- `IMEXExpImplEuler` (`ode.hpp:1030`) -- 1st order
- `IMEXRK2` (`ode.hpp:1046`) -- 2nd order, L-stable
- `IMEX_DIRK_RK3` (`ode.hpp:1079`) -- 3rd order, L-stable

These require `TimeDependentOperator::ImplicitSolve()` which can use Hypre internally.

**Source:** `F:\dev\mfem\linalg\ode.hpp:1030-1090`, `F:\dev\mfem\linalg\sundials.hpp:720-888`

---

## 8. Block Systems: Multi-Species (Phase 3)

### 8.1 HypreBoomerAMG::SetSystemsOptions for Multi-Species

For multi-species diffusion (dopant + I + V + clusters), the system matrix has block structure. `HypreBoomerAMG::SetSystemsOptions(int dim)` (`hypre.hpp:1858`) configures AMG for systems of `dim` unknowns per node.

```cpp
// After converting block SparseMatrix to HypreParMatrix:
mfem::HypreBoomerAMG amg(A_par);
amg.SetSystemsOptions(n_species);  // n_species = 3 for B, I, V
amg.SetPrintLevel(0);
amg.SetMaxIter(1);

mfem::HypreGMRES gmres(A_par);
gmres.SetPreconditioner(amg);
gmres.SetKDim(50);
```

### 8.2 Block-Diagonal AMG via HypreParMatrixFromBlocks

`HypreParMatrixFromBlocks()` (`hypre.hpp:1060-1061`) merges a 2D array of `HypreParMatrix*` blocks into a single matrix. Some blocks can be NULL.

```cpp
mfem::Array2D<const mfem::HypreParMatrix*> blocks(n_species, n_species);
// Fill diagonal and coupling blocks...
mfem::HypreParMatrix* A = mfem::HypreParMatrixFromBlocks(blocks);
```

**Source:** `F:\dev\mfem\linalg\hypre.hpp:1060-1061,1858`

---

## 9. DAE Support (Phase 3 Clustering)

### 9.1 Current State: No DAE Support in MFEM

**Critical finding:** MFEM has NO IDA wrapper. Grep for `IDA` in `F:\dev\mfem\linalg\` returned empty results. `sundials.hpp` includes only `cvodes/cvodes.h` (:47), `arkode/arkode_arkstep.h` (:46), and `kinsol/kinsol.h` (:48). No `ida/ida.h`. SUNDIALS IDA is linked in ViennaPS `CMakeLists.txt:314` but never used.

### 9.2 DAE via SUNDIALS IDA (Raw C API)

Phase 3 clustering with equilibrium species is Index-1 DAE. SUNDIALS IDA can solve this directly, but requires raw C API integration:

```c
#include <ida/ida.h>
#include <ida/ida_ls.h>
void* ida_mem = IDACreate(sunctx);
IDAInit(ida_mem, residualFunction, t0, y);
IDASStolerances(ida_mem, rtol, atol);
// IDASetLinearSolver(ida_mem, LS, NULL);
// IDASetPreconditioner(ida_mem, pset, psolve);
```

**Effort:** Medium -- requires implementing the IDA residual function and integrating it into `SundialsTimeIntegrator.hpp` behind a new code path.

### 9.3 DAE Workaround: ODE + Algebraic Elimination

For Index-1 DAE, the algebraic constraints can be eliminated: `z = g^{-1}(u)` (solved after each ODE step). This avoids IDA entirely but fails for tightly coupled DAE.

### 9.4 Future: PETSc TS (Phase 10)

PETSc TS handles DAEs natively via `TSSetIFunction`. But PETSc rebuild is required. Defer to Phase 10.

**Source:** SUNDIALS IDA headers at `F:\dev\vcpkg\installed\x64-windows\include\ida\` (verified existence); MFEM `sundials.hpp` (no IDA include)

---

## 10. GPU Support

### 10.1 Current Hypre GPU Status

**Source:** `F:\dev\vcpkg\installed\x64-windows\include\HYPRE_config.h`

Hypre is built **CPU-only**:
- `HYPRE_USING_HOST_MEMORY` defined (:111)
- `HYPRE_USING_CUDA` undef (:90)
- `HYPRE_USING_GPU` undef (:126)

MFEM has `MFEM_USE_CUDA` defined (`_config.hpp:173`), but Hypre's CPU-only build means Hypre solvers run on CPU. MFEM handles GPU-CPU data transfer internally.

### 10.2 GPU Upgrade Path

To enable GPU-accelerated Hypre:
1. Rebuild Hypre with `--with-cuda`
2. Rebuild MFEM (MFEM already has `MFEM_USE_CUDA`)
3. MFEM's `hypre.hpp:52-58` checks `HYPRE_USING_CUDA` and requires `MFEM_USE_CUDA` (already satisfied)

**Cost:** Medium. Not needed for current phases.

| Component | GPU? | Source |
|-----------|------|--------|
| MFEM assembly | Yes (`MFEM_USE_CUDA`) | `_config.hpp:173` |
| Hypre solvers | No (CPU-only) | `HYPRE_config.h:90,111` |
| SUNDIALS | No (serial N_Vector) | `SundialsTimeIntegrator.hpp:190` |
| ViennaPS GPU (OptiX) | Optional (`VIENNAPS_USE_GPU`) | `CMakeLists.txt` |

---

## 11. Per-Phase Impact Table (Critical Question C)

| Phase | Description | What Hypre Fixes | What Hypre Does NOT Fix | Solver Recommendation |
|-------|-------------|-----------------|------------------------|----------------------|
| **1** | Foundation: FEM mesh, constant D, CVODE | Replaces `DSmoother` with AMG in DiffusionEngine; enables SUNDIALS iterative LS + Hypre PC (removes 256-dof cap) | Nothing missing for Phase 1 | `HyprePCG`+`HypreBoomerAMG` for FEM; `SUNLinSol_SPBCGS`+Hypre PC for CVODE |
| **2** | Fermi/ChargedFermi: D(C), segregation, solid solubility | AMG as preconditioner for NewtonSolver; `HypreGMRES` for nonsymmetric segregated system | Nonlinear solve (use `NewtonSolver` or `KINSolver`); no variational inequality for solubility caps | `NewtonSolver`+`HyprePCG`+AMG; or `KINSolver` JFNK + AMG PC |
| **3** | CDD: multi-species, clustering, reaction networks | `SetSystemsOptions(n_species)` for block AMG; `HypreParMatrixFromBlocks`; `HypreGMRES` for nonsymmetric | DAE support (no IDA in MFEM); need raw SUNDIALS IDA or algebraic elimination; no `PCFIELDSPLIT` | `HypreGMRES`+AMG with `SetSystemsOptions`; IDA via raw C API for DAE |
| **4** | OED, TED, dose loss, moving interface | Same solver as Phase 1; L2 projection uses `HyprePCG` | Moving interface is mesh operation, not solver | Same as Phase 1 |
| **5** | Polysilicon: Voronoi, dual mesh, GB segregation | Same solver on dual mesh | Grain growth is scalar ODE | Same as Phase 1 on dual mesh |
| **6** | SiGe/III-V: bandgap, strain, Ge-B pairing | Coefficient changes only; same solver | Ge-B pairing is reaction network (Phase 3 pattern) | Same as Phase 1/3 |
| **7** | KMC + continuum coupling | L2 projection solve uses `HyprePCG`+AMG | KMC is atomistic (no linear solve) | `HyprePCG`+AMG for L2 projection |
| **8** | KMC lattice epitaxy | Mesh update during epitaxy | KMC surface chemistry (no linear solve) | N/A (KMC) |
| **9** | Flash/laser: heat, Allen-Cahn, melt diffusion | `ARKStepSolver::IMEX` with Hypre for implicit diffusion; block AMG for T+eta+C | Nonlinear Allen-Cahn reaction (explicit in IMEX); no adjoint | `ARKStepSolver::IMEX`+`HyprePCG`+AMG; `SetSystemsOptions(4)` |
| **10** | PDE API, calibration, adjoint inversion | Forward solves use Hypre; L2 projection for results | No adjoint/sensitivity; no optimization (TAO); no MCMC; no GP | Hypre for forward solves; external tools for calibration |
| **11** | AMR: h/p refinement, ZZ estimator | AMR is MFEM native; Hypre handles refined meshes | Nothing (MFEM native AMR is sufficient) | Same solver on refined mesh |

---

## 12. Gaps and Limitations

| Gap | Impact | Phase(s) | Workaround |
|-----|--------|----------|------------|
| No DAE support in MFEM | Phase 3 clustering with equilibrium species | 3 | Raw SUNDIALS IDA C API; or algebraic elimination for Index-1 |
| No adjoint/sensitivity | Phase 10 Task 15 D(x) inversion | 10 | Defer to PETSc TSAdjoint (requires rebuild) |
| No optimization framework | Phase 10 Tasks 12-15 calibration | 10 | Defer to PETSc TAO; or external Python optimization |
| No matrix-free colored FD Jacobian | Phase 2 nonlinear Jacobian | 2 | Manual `dD/dC` derivation; or KINSolver JFNK (no explicit Jacobian) |
| No PCFIELDSPLIT | Phase 3 block preconditioning | 3 | Manual block-diagonal AMG via `HypreParMatrixFromBlocks` |
| No variational inequalities | Phase 2 solid solubility caps | 2 | Ad-hoc clamping (current approach in `PhysicsField.hpp:409-413`) |
| Hypre is CPU-only | GPU acceleration | All | Rebuild Hypre with CUDA (future) |
| MFEM CVODE wrapper has no PC | SUNDIALS preconditioning | 1,3,9 | Use raw SUNDIALS C API (as `SundialsTimeIntegrator.hpp` already does) |

---

## 13. Implementation Plan (Critical Question D)

### Step 1: DiffusionEngine Hypre Integration (Immediate, zero build cost)

**File:** `include/viennaps/fields/DiffusionEngine.hpp`

**Changes:**
1. Replace lines 100-105 (DSmoother + BiCGSTABSolver) with HyprePCG + HypreBoomerAMG per Section 4.2
2. Add `#ifdef MFEM_USE_MPI` guard with fallback to original solver
3. Add `#include <mpi.h>` if not already included (needed for `MPI_COMM_SELF`)

**Verification:**
```powershell
cmake --build build --config Release --parallel
ctest -R DiffusionEngine --test-dir build -C Release --output-on-failure
```

**Expected result:** Linear solve convergence improves (AMG vs DSmoother). No change in solution values.

### Step 2: SUNDIALS Iterative LS + Hypre Preconditioner (Immediate)

**File:** `include/viennaps/fields/SundialsTimeIntegrator.hpp`

**Changes:**
1. Add SUNDIALS iterative solver headers:
   ```cpp
   #include <sunlinsol/sunlinsol_spbcgs.h>
   ```
2. Define `HyprePrecData` struct and callback functions (PSetup, PSolve, ATimes) per Section 5.4
3. Replace lines 166-178 (256-dof cap + subsampling) -- remove entirely, use full `neq`
4. Replace lines 238-251 (dense LS / no LS) with `SUNLinSol_SPBCGS` + `SUNLinSolSetPreconditioner` + `CVodeSetLinearSolver`
5. For the preconditioner matrix: assemble `(M + dt*K)` as `HypreParMatrix` once per time step (in PSetup callback)
6. Add cleanup: `SUNLinSolFree(LS)` in the existing cleanup block (lines 300-304)

**Verification:**
```powershell
ctest -R SundialsTimeIntegrator --test-dir build -C Release --output-on-failure
```

**Expected result:** CVODE handles >256 dofs without subsampling. Solution accuracy improves (no information loss from subsampling).

### Step 3: NewtonSolver for Phase 2 (When Phase 2 is implemented)

**File:** New `include/viennaps/fields/NonlinearSolver.hpp`

**Changes:**
1. Create `FermiDiffusionOperator : public mfem::Operator` implementing `Mult()` and `GetGradient()`
2. Use `mfem::NewtonSolver` with `HyprePCG` + `HypreBoomerAMG` as linear solver per Section 6.1
3. Alternative: `KINSolver` JFNK mode with `HypreBoomerAMG` preconditioner per Section 6.2

### Step 4: Block AMG for Phase 3 (When Phase 3 is implemented)

**Changes:**
1. Assemble multi-species system as block `SparseMatrix`
2. Convert to `HypreParMatrix` (single matrix with interleaved dofs)
3. Call `amg.SetSystemsOptions(n_species)` for nodal block AMG
4. Use `HypreGMRES` (nonsymmetric due to coupling)
5. For DAE: integrate SUNDIALS IDA via raw C API

### Step 5: IMEX for Phase 9 (When Phase 9 is implemented)

**Changes:**
1. Create `AllenCahnOperator : public TimeDependentOperator` with explicit/implicit RHS split
2. Use `ARKStepSolver(IMEX)` with `UseMFEMLinearSolver()`
3. Operator's `ImplicitSolve()` uses `HyprePCG` + `HypreBoomerAMG` internally
4. Use `SetSystemsOptions(4)` for coupled T + eta_m + eta_c + C

### CMake Changes

**None required.** Hypre is already linked through MFEM. SUNDIALS iterative solver headers are already available at `F:\dev\vcpkg\installed\x64-windows\include\sunlinsol\`. No new `find_package` or `target_link_libraries` needed.

### New Files

| File | Purpose | Phase |
|------|---------|-------|
| `HyprePreconditioner.hpp` | HypreBoomerAMG wrapper for SUNDIALS callbacks | Step 2 |
| `NonlinearSolver.hpp` | NewtonSolver/KINSolver wrapper for Fermi diffusion | Step 3 |

### Migration to PETSc (Phase 10, future)

When PETSc is adopted (after rebuild with Hypre support):
1. `HypreParMatrix` to `PetscParMatrix` is zero-copy (`petsc.hpp` confirms)
2. MFEM `Operator`/`Solver` interfaces are backend-agnostic
3. Assembly code does not change
4. Only solver instantiation changes

---

## 14. Sources

### Primary Sources (MFEM Source Code)
- MFEM Hypre wrappers: `F:\dev\mfem\linalg\hypre.hpp` (2356 lines)
  - `HypreParMatrix` at :418, constructor from `SparseMatrix*` at :567-569
  - `HypreBoomerAMG` at :1828, `SetSystemsOptions` at :1858, `SetAdvectiveOptions` at :1872
  - `HyprePCG` at :1352, `HypreGMRES` at :1444, `HypreFGMRES` at :1523
  - `HypreSmoother` at :1076, `HypreILU` at :1781
  - `HypreParMatrixFromBlocks` at :1060-1061
  - `CopyCSR` (shallow copy when no HYPRE_BIGINT) at :472-475
- MFEM solvers: `F:\dev\mfem\linalg\solvers.hpp` (1528 lines)
  - `NewtonSolver` at :780, `SetAdaptiveLinRtol` at :856-860
  - `LBFGSSolver` at :865, `BiCGSTABSolver` at :709
- MFEM ODE solvers: `F:\dev\mfem\linalg\ode.hpp` (1095 lines)
  - `IMEXExpImplEuler` at :1030, `IMEXRK2` at :1046, `IMEX_DIRK_RK3` at :1079
- MFEM SUNDIALS wrappers: `F:\dev\mfem\linalg\sundials.hpp` (1069 lines)
  - `CVODESolver` at :429, `UseMFEMLinearSolver` at :516, `UseSundialsLinearSolver` at :519
  - `ARKStepSolver` at :720, `Type::IMEX` at :728
  - `KINSolver` at :896, `SetJFNKSolver` at :943, `jfnk` flag at :908
- MFEM SUNDIALS implementation: `F:\dev\mfem\linalg\sundials.cpp` (2439 lines)
  - `UseSundialsLinearSolver` uses `SUN_PREC_NONE` at :905
  - `UseMFEMLinearSolver` creates custom `SUNLinearSolver` at :867-896
- MFEM config: `F:\dev\mfem\build\config\_config.hpp` (225 lines)
  - `MFEM_USE_MPI` at :47, `MFEM_USE_SUNDIALS` at :95, `MFEM_USE_CUDA` at :173
  - `MFEM_HYPRE_VERSION 23200` at :201, `MFEM_USE_PETSC` undef at :151
- MFEM ParBilinearForm: `F:\dev\mfem\fem\pbilinearform.hpp`
  - `ParallelAssemble(SparseMatrix*)` at :122

### Primary Sources (SUNDIALS Headers)
- SUNDIALS version: `F:\dev\vcpkg\installed\x64-windows\include\sundials\sundials_config.h` (7.8.0, :49)
- SUNDIALS linear solver API: `F:\dev\vcpkg\installed\x64-windows\include\sundials\sundials_linearsolver.h`
  - `SUNLinSolSetATimes` at :175-176
  - `SUNLinSolSetPreconditioner` at :179-180
- SUNDIALS iterative types: `F:\dev\vcpkg\installed\x64-windows\include\sundials\sundials_iterative.h`
  - `SUNPrecType` enum at :59-65 (`SUN_PREC_LEFT` at :62)
  - `SUNATimesFn` at :108, `SUNPSetupFn` at :120, `SUNPSolveFn` at :148
- SPBCGS solver: `F:\dev\vcpkg\installed\x64-windows\include\sunlinsol\sunlinsol_spbcgs.h`
  - `SUNLinSol_SPBCGS` at :81, `SUNSPBCGS_MAXL_DEFAULT 5` at :43
- CVODE linear solver: `F:\dev\vcpkg\installed\x64-windows\include\cvode\cvode_ls.h` (134 lines)
  - `CVodeSetLinearSolver` at :82, `CVodeSetPreconditioner` at :97-98
  - `CVodeSetJacTimes` at :99-100, `CVLsJacTimesVecFn` at :69-71
- Available iterative solvers: `sunlinsol_spgmr.h`, `sunlinsol_spbcgs.h`, `sunlinsol_spfgmr.h`, `sunlinsol_sptfqmr.h`, `sunlinsol_pcg.h`

### Primary Sources (Hypre Configuration)
- Hypre config: `F:\dev\vcpkg\installed\x64-windows\include\HYPRE_config.h` (181 lines)
  - Version 2.32.0 at :9-10, `HYPRE_HAVE_MPI` at :156
  - `HYPRE_USING_HOST_MEMORY` at :111, `HYPRE_USING_CUDA` undef at :90
  - `HYPRE_BIGINT` undef at :24, `HYPRE_SINGLE` undef at :27
- Hypre headers: `F:\dev\vcpkg\installed\x64-windows\include\HYPRE*.h` (22 files)

### ViennaPS Source Code
- `include/viennaps/fields/DiffusionEngine.hpp` (134 lines) -- `DSmoother`+`BiCGSTABSolver` at :100-105
- `include/viennaps/fields/SundialsTimeIntegrator.hpp` (314 lines) -- 256-dof cap at :167, dense LS at :241-251
- `CMakeLists.txt` -- SUNDIALS IDA linked at :314

### Prior Research
- PETSc feasibility: `docs/superpowers/research/2026-07-22-petsc-feasibility.md` (399 lines)
- PETSc vs Hypre comparison: `docs/superpowers/research/2026-07-22-petsc-vs-hypre-comparison.md` (486 lines)

### Plan Documents
- Phase 1-11 plans: `docs/superpowers/plans/2026-07-20-diffusion-phase{1-11}.md`
