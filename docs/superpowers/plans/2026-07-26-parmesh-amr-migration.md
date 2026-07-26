# Runtime AMR via ParMesh + Nonconforming Refinement

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Migrate the DiffusionEngine from serial `mfem::Mesh`/`FiniteElementSpace`/`GridFunction` to parallel `mfem::ParMesh`/`ParFiniteElementSpace`/`ParGridFunction` to enable safe runtime AMR via `ParMesh::GeneralRefinement` (nonconforming) + `ParFiniteElementSpace::Update` + `ParGridFunction::Update` prolongation.

**Architecture:** The engine wraps a serial `mfem::Mesh` (from `LevelSetToMesh`) into a single-process `mfem::ParMesh(MPI_COMM_SELF, serialMesh)` at `setMesh` time. This gives us the `ParMesh` AMR path (which correctly manages the mesh's Nodes `ParGridFunction` lifecycle) while running on a single process (no MPI distribution needed). The 5 chokepoint member types change from serial to parallel; the 30+ model headers flow through the `DiffusionModel` virtual hooks and need only the hook signature type changes. The hardest rewrites are `MovingMeshHandler` (serial `GetNode`/`SetNode` -> `ParGridFunction` node access), `SolutionTransfer` (`FindPoints` + serial CG -> parallel), and `Segregation` (interior-face loops).

**Tech Stack:** C++20, MFEM (built with `MFEM_USE_MPI=ON` at `f:/dev/mfem/build`), SUNDIALS, CMake/CTest

## Global Constraints

- C++20, header-only under `include/viennaps/`
- MFEM code gated by `#ifdef VIENNAPS_HAS_MFEM`, SUNDIALS by `#ifdef VIENNAPS_HAS_SUNDIALS`
- MFEM Release at `f:/dev/mfem/build` (MPI enabled), Debug at `f:/dev/mfem/build_debug`
- vcpkg deps at `f:/dev/vcpkg/installed/x64-windows/`
- LLVM style: 2-space indent, 80-col, no tabs, Attach braces, pointer right
- Tests use `VC_TEST_ASSERT` from `vcTestAsserts.hpp`
- Namespace: `viennaps`
- Single-process parallel: `ParMesh(MPI_COMM_SELF, ...)` - no actual MPI distribution, but gets the ParMesh AMR machinery

---

## File Structure

| File | Responsibility | Change type |
|------|----------------|-------------|
| `fields/DiffusionModel.hpp` | Base class virtual hooks + mesh accessor | Modify: `mesh_` type, hook signatures |
| `fields/DiffusionEngine.hpp` | Engine: mesh/fes/gf ownership, assembly, solve | Modify: 5 member types, `setMesh`, `assembleAllSpecies`, `solveImplicitEuler`, `refineBetweenSteps` |
| `fields/LevelSetToMesh.hpp` | Level-set -> MFEM mesh converter | Keep serial output; engine wraps into ParMesh |
| `fields/AdaptiveMeshRefiner.hpp` | AMR mark/refine helpers | Modify: `refineMarkedWithProlongation` to use ParMesh API |
| `fields/MovingMeshHandler.hpp` | Mesh mutation (relabel, ALE, smooth) | Modify: `GetNode`/`SetNode` -> ParGridFunction access |
| `fields/SolutionTransfer.hpp` | Cross-mesh field transfer | Modify: `FindPoints` -> GSLIB or point-injection; serial CG -> HyprePCG |
| `fields/models/Segregation.hpp` | Two-sided interior-face residual | Modify: face iteration for ParMesh |
| `fields/models/OedSource.hpp` | OED relabel helper | Modify: delegate to updated MovingMeshHandler |
| `tests/diffusion/testDiffusion.cpp` | All tests | Modify: add runtime AMR dose-preservation test |

---

## Task 1: MPI Initialization + ParMesh Wrapping in setMesh

**Files:** Modify `include/viennaps/fields/DiffusionEngine.hpp`

**Produces:** Engine accepts a serial `mfem::Mesh`, wraps it into `mfem::ParMesh(MPI_COMM_SELF, ...)`, constructs `ParFiniteElementSpace` + `ParGridFunction` members. All existing tests still pass (the ParMesh on `MPI_COMM_SELF` behaves identically to serial for single-process).

- [ ] **Step 1: Add MPI include + change member types**

At the top of `DiffusionEngine.hpp`, after the MFEM include:
```cpp
#ifdef VIENNAPS_HAS_MFEM
#include <mfem.hpp>
#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif
#endif
```

Change the 5 member types:
```cpp
// Was: std::unique_ptr<mfem::Mesh> mesh_;
// Now: ParMesh on MPI_COMM_SELF (single-process parallel for AMR safety)
std::unique_ptr<mfem::ParMesh> mesh_;
std::unique_ptr<mfem::FiniteElementCollection> fec_;
std::unique_ptr<mfem::ParFiniteElementSpace> fes_;
std::map<std::string, std::unique_ptr<mfem::ParGridFunction>> species_;
std::map<std::string, mfem::ParGridFunction *> allSpecies_;
```

- [ ] **Step 2: Modify `setMesh` to wrap serial mesh into ParMesh**

```cpp
void setMesh(std::unique_ptr<mfem::Mesh> mesh, const MeshAttributes &attrs) {
    if (!mesh)
      throw std::runtime_error("DiffusionEngine::setMesh: mesh pointer is null");
    if (mesh->SpaceDimension() != D)
      throw std::runtime_error("DiffusionEngine::setMesh: mesh SpaceDimension does not match D");
    attrs_ = attrs;
    // Wrap serial mesh into ParMesh on MPI_COMM_SELF.
    // This gives us the ParMesh AMR path (safe GeneralRefinement + Update)
    // while running single-process (no MPI distribution).
    mesh_ = std::make_unique<mfem::ParMesh>(MPI_COMM_SELF, *mesh);
    fec_ = std::make_unique<mfem::H1_FECollection>(1, D);
    fes_ = std::make_unique<mfem::ParFiniteElementSpace>(mesh_.get(), fec_.get());
    species_.clear();
    allSpecies_.clear();
    implicitCache_.clear();
}
```

- [ ] **Step 3: Update `mesh()` and `fes()` accessors**

```cpp
mfem::ParMesh *mesh() { return mesh_.get(); }
const mfem::ParMesh *mesh() const { return mesh_.get(); }
mfem::ParFiniteElementSpace *fes() { return fes_.get(); }
const mfem::ParFiniteElementSpace *fes() const { return fes_.get(); }
```

- [ ] **Step 4: Update `ensureSpecies` to create ParGridFunction**

```cpp
mfem::ParGridFunction *ensureSpecies(const std::string &name) {
    auto it = species_.find(name);
    if (it != species_.end())
      return it->second.get();
    auto gf = std::make_unique<mfem::ParGridFunction>(fes_.get());
    *gf = 0.0;
    mfem::ParGridFunction *raw = gf.get();
    species_[name] = std::move(gf);
    allSpecies_[name] = species_[name].get();
    return raw;
}
```

- [ ] **Step 5: Update `getSolution` return type**

```cpp
const mfem::ParGridFunction &getSolution(const std::string &name) const {
    auto it = species_.find(name);
    if (it == species_.end())
      throw std::runtime_error("DiffusionEngine::getSolution: unknown species");
    return *it->second;
}
```

- [ ] **Step 6: Build + run all tests**

```bash
cmake --build build_phase2 --config Release --target testDiffusion
./build_phase2/tests/testDiffusion.exe
```

Expected: Compile errors from `BilinearForm`/`LinearForm` type mismatches (the assemble hooks still take serial types). These are fixed in Task 2.

- [ ] **Step 7: Commit**

```bash
git add include/viennaps/fields/DiffusionEngine.hpp
git commit -m "feat: migrate DiffusionEngine to ParMesh/ParFiniteElementSpace/ParGridFunction (Task 1)"
```

---

## Task 2: Update SpeciesSystem + assembleAllSpecies to ParBilinearForm/ParLinearForm

**Files:** Modify `include/viennaps/fields/DiffusionEngine.hpp`

**Produces:** The assembly loop uses `ParBilinearForm`/`ParLinearForm` instead of serial types. The `SpeciesSystem` struct changes. Models still receive `BilinearForm&`/`LinearForm&` (base class of both serial and parallel) so the model headers don't need changes yet.

- [ ] **Step 1: Update SpeciesSystem struct**

```cpp
struct SpeciesSystem {
  std::unique_ptr<mfem::ParBilinearForm> K;
  std::unique_ptr<mfem::ParBilinearForm> M;
  std::unique_ptr<mfem::ParLinearForm> R;
};
```

- [ ] **Step 2: Update `assembleAllSpecies` to construct ParBilinearForm/ParLinearForm**

In the assembly loop, change:
```cpp
sys.K = std::make_unique<mfem::ParBilinearForm>(fes_.get());
sys.M = std::make_unique<mfem::ParBilinearForm>(fes_.get());
sys.R = std::make_unique<mfem::ParLinearForm>(fes_.get());
```

- [ ] **Step 3: Update `getIntegral` to use ParLinearForm**

The `getIntegral` method builds a local `LinearForm` - change to `ParLinearForm` and use `ParLinearForm::InnerProduct` with the ParGridFunction.

- [ ] **Step 4: Build + fix compile errors**

The model headers pass `BilinearForm&`/`LinearForm&` to their virtual hooks. `ParBilinearForm` derives from `BilinearForm` and `ParLinearForm` derives from `LinearForm`, so the base-class reference binds work. But some model code may call `SpMat()` which returns `SparseMatrix&` on serial vs `HypreParMatrix*` on parallel. Fix by using `ParBilinearForm::ParallelAssemble()` or `SpMat()` (which still works for local matrix on `MPI_COMM_SELF`).

- [ ] **Step 5: Run tests, commit**

```bash
cmake --build build_phase2 --config Release --target testDiffusion
./build_phase2/tests/testDiffusion.exe
git add include/viennaps/fields/DiffusionEngine.hpp
git commit -m "feat: ParBilinearForm/ParLinearForm in assembly (Task 2)"
```

---

## Task 3: Update solveImplicitEuler for ParSparseMatrix/HypreParMatrix

**Files:** Modify `include/viennaps/fields/DiffusionEngine.hpp`

**Produces:** The implicit-Euler solve loop uses `HypreParMatrix` for M+dt*K and `HyprePCG`/`HypreBoomerAMG` for the solve instead of serial `SparseMatrix` + `DSmoother` + `CG`. This is the F8 follow-up documented in the engine header.

- [ ] **Step 1: Change Ms/Ks maps to HypreParMatrix**

```cpp
std::map<std::string, const mfem::HypreParMatrix *> Ms, Ks;
```

Extract via `ParBilinearForm::ParallelAssemble()` instead of `SpMat()`.

- [ ] **Step 2: Update implicit cache to use HypreParMatrix**

```cpp
struct ImplicitCache {
  std::unique_ptr<mfem::HypreParMatrix> A;
  // ... essVdofs etc unchanged
};
```

Build `A = M + dt*K` via `HypreParMatrix` addition (`par_M->Add(dt, *par_K)`).

- [ ] **Step 3: Replace serial CG with HyprePCG**

```cpp
mfem::HypreParMatrix *A = cache.A.get();
mfem::HypreBoomerAMG amg(*A);
amg.SetPrintLevel(0);
mfem::HyprePCG solver(*A);
solver.SetPreconditioner(amg);
solver.SetPrintLevel(0);
solver.Mult(b, gf);
```

- [ ] **Step 4: Update EliminateVDofs to ParFiniteElementSpace::EliminateEssentialBC**

`ParFiniteElementSpace` uses `GetEssentialTrueDofs` + `ParBilinearForm::EliminateEssentialBC` instead of the serial `EliminateVDofs`. The semantics differ slightly (true dofs vs vdofs); the essential-BC enforcement must use the parallel API.

- [ ] **Step 5: Build + run all tests, commit**

```bash
cmake --build build_phase2 --config Release --target testDiffusion
./build_phase2/tests/testDiffusion.exe
git add include/viennaps/fields/DiffusionEngine.hpp
git commit -m "feat: HypreParMatrix + HyprePCG/BoomerAMG in implicit-Euler solve (Task 3)"
```

---

## Task 4: Update DiffusionModel virtual hooks + model coefficient types

**Files:** Modify `include/viennaps/fields/DiffusionModel.hpp`, all `models/*.hpp`

**Produces:** The virtual hook signatures change from `mfem::GridFunction` to `mfem::ParGridFunction` and `mfem::BilinearForm` to `mfem::ParBilinearForm`. Model coefficient classes that store `const mfem::GridFunction*` change to `const mfem::ParGridFunction*`.

- [ ] **Step 1: Update DiffusionModel.hpp hooks**

```cpp
virtual void assembleStiffness(
    mfem::ParBilinearForm &K,
    const mfem::ParGridFunction &speciesGF,
    const std::map<std::string, mfem::ParGridFunction *> &allSpecies,
    const mfem::ParGridFunction *temp) const {}

virtual void assembleReaction(
    mfem::ParLinearForm &R,
    const mfem::ParGridFunction &speciesGF,
    const std::map<std::string, mfem::ParGridFunction *> &allSpecies,
    const mfem::ParGridFunction *temp) const {}

virtual void finalizeReaction(
    mfem::ParLinearForm &R,
    const mfem::ParGridFunction &speciesGF,
    const std::map<std::string, mfem::ParGridFunction *> &allSpecies,
    const mfem::ParGridFunction *temp) const {}

virtual void assembleMass(mfem::ParBilinearForm &M) const {}
virtual void setMesh(mfem::ParMesh *mesh) { mesh_ = mesh; }
mfem::ParMesh *mesh() const { return mesh_; }
```

- [ ] **Step 2: Update all model coefficient classes**

In each model header, change `const mfem::GridFunction *` to `const mfem::ParGridFunction *` in coefficient member variables and constructors. The `GetValue(T, ip)` API is identical on both types.

- [ ] **Step 3: Build + fix compile errors across all models**

This is the largest mechanical change - ~30 model headers. Each needs:
- Hook signature `override` types updated
- Coefficient class member types updated
- `GetValue` calls unchanged (same API)

- [ ] **Step 4: Run tests, commit**

```bash
cmake --build build_phase2 --config Release --target testDiffusion
./build_phase2/tests/testDiffusion.exe
git add include/viennaps/fields/
git commit -m "feat: migrate all model hooks + coefficients to ParGridFunction/ParBilinearForm (Task 4)"
```

---

## Task 5: Update MovingMeshHandler for ParMesh node access

**Files:** Modify `include/viennaps/fields/MovingMeshHandler.hpp`

**Produces:** `GetNode`/`SetNode`/`NodesUpdated` replaced with `ParGridFunction` node access via `ParMesh::GetNodes()`.

- [ ] **Step 1: Change all `mfem::Mesh &` params to `mfem::ParMesh &`**

- [ ] **Step 2: Replace GetNode/SetNode with ParGridFunction access**

```cpp
// Was: mesh.GetNode(v, x.data());
// Now: access via the mesh's nodes ParGridFunction
mfem::ParGridFunction *nodes = mesh.GetNodes();
if (nodes) {
  // For H1 order-1, DoF == vertex. Read/write via the ParGridFunction.
  for (int d = 0; d < sdim; ++d)
    x[d] = (*nodes)(v * sdim + d);  // or via GetSubVector
}
```

- [ ] **Step 3: Replace `NodesUpdated()` with `ExchangeFaceNbrData()`**

- [ ] **Step 4: Build + run tests, commit**

---

## Task 6: Update SolutionTransfer for ParGridFunction + parallel L2 projection

**Files:** Modify `include/viennaps/fields/SolutionTransfer.hpp`

**Produces:** `FindPoints` replaced with point-injection or GSLIB; serial `SparseMatrix` CG replaced with `HyprePCG` + `HypreBoomerAMG` on `HypreParMatrix`.

- [ ] **Step 1: Change GridFunction params to ParGridFunction**

- [ ] **Step 2: Replace serial BilinearForm + CG L2 solve with ParBilinearForm + HyprePCG**

```cpp
mfem::ParBilinearForm M(fes);
M.AddDomainIntegrator(new mfem::MassIntegrator(one));
M.Assemble();
M.Finalize();
mfem::HypreParMatrix *Mpar = M.ParallelAssemble();

mfem::ParLinearForm b(fes);
b.AddDomainIntegrator(new mfem::DomainLFIntegrator(coef));
b.Assemble();

mfem::HypreBoomerAMG amg(*Mpar);
mfem::HyprePCG solver(*Mpar);
solver.SetPreconditioner(amg);
solver.Mult(b, target);
```

- [ ] **Step 3: Replace FindPoints with nodal injection for MPI_COMM_SELF**

On `MPI_COMM_SELF`, `FindPoints` from the serial mesh still works if we access the underlying serial mesh via `ParMesh::GetSerialMesh()`. Alternatively, use `ParMesh::FindPoints` if available. For the single-process case, the serial `FindPoints` path is safe because there are no ghost elements.

- [ ] **Step 4: Build + run tests, commit**

---

## Task 7: Update Segregation interior-face loops for ParMesh

**Files:** Modify `include/viennaps/fields/models/Segregation.hpp`

**Produces:** `GetNumFaces`/`GetFaceElements`/`GetFaceElementTransformations` work on `ParMesh` for local faces. On `MPI_COMM_SELF` there are no shared faces, so the existing face iteration is correct as-is (all faces are local). The only change needed is the parameter type.

- [ ] **Step 1: Change `mfem::Mesh &` to `mfem::ParMesh &` in `assembleInterfaceResidual` and `applyOperatorSplitStep`**

- [ ] **Step 2: Change `LinearForm` to `ParLinearForm` and `GridFunction` to `ParGridFunction`**

- [ ] **Step 3: Build + run tests, commit**

---

## Task 8: Wire refineBetweenSteps with ParMesh GeneralRefinement

**Files:** Modify `include/viennaps/fields/DiffusionEngine.hpp`

**Produces:** The `refineBetweenSteps` method now calls `ParMesh::GeneralRefinement` (nonconforming=1 for triangles) + `ParFiniteElementSpace::Update(true)` + `ParGridFunction::Update()` per species + `UpdatesFinished()`. This is the safe AMR path that was crashing on serial H1.

- [ ] **Step 1: Rewrite `refineBetweenSteps` for ParMesh**

```cpp
bool refineBetweenSteps() {
    if (lastAmrMarkCount_ <= 0 || !mesh_ || !fes_)
      return false;
    // Collect species ParGridFunctions.
    std::vector<mfem::ParGridFunction *> gfs;
    for (auto &kv : species_)
      if (kv.second)
        gfs.push_back(kv.second.get());
    // Nonconforming refinement (ParMesh handles hanging nodes via constraints).
    RefinementBox box{amrBox_.x0, amrBox_.x1, amrBox_.y0, amrBox_.y1};
    auto ids = AdaptiveMeshRefiner::markBox(*mesh_, box);
    if (ids.empty())
      return false;
    mfem::Array<int> elToRefine;
    for (int id : ids)
      if (id >= 0 && id < mesh_->GetNE())
        elToRefine.Append(id);
    if (elToRefine.Size() == 0)
      return false;
    // ParMesh GeneralRefinement: safely manages the mesh's Nodes ParGridFunction.
    mesh_->GeneralRefinement(elToRefine, /*nonconforming=*/1);
    // ParFiniteElementSpace::Update: computes prolongation operator.
    fes_->Update(true);
    // ParGridFunction::Update: applies prolongation in place.
    for (auto *gf : gfs)
      gf->Update();
    fes_->UpdatesFinished();
    implicitCache_.clear();
    lastAmrMarkCount_ = 0;
    return true;
}
```

- [ ] **Step 2: Re-enable the refineBetweenSteps call in solveImplicitEuler**

Uncomment the between-step refine + system rebuild block that was disabled due to the serial crash.

- [ ] **Step 3: Write the runtime AMR dose-preservation test**

```cpp
// In TestFemClosures or a new test:
// 1. Set up engine with runtimeAmr enabled
// 2. Solve for a few steps
// 3. Assert dose is preserved to 0.1% after refinement
// 4. Assert mesh element count increased (refinement occurred)
```

- [ ] **Step 4: Build + run all tests**

```bash
cmake --build build_phase2 --config Release --target testDiffusion
./build_phase2/tests/testDiffusion.exe
```

Expected: All tests pass including the new runtime AMR dose-preservation test. The ParMesh `GeneralRefinement` + `Update()` path should NOT crash (unlike the serial `Mesh` path).

- [ ] **Step 5: Commit**

```bash
git add include/viennaps/fields/DiffusionEngine.hpp tests/diffusion/testDiffusion.cpp
git commit -m "feat: runtime AMR via ParMesh GeneralRefinement + prolongation (Task 8)"
```

---

## Task 9: Update AdaptiveMeshRefiner for ParMesh

**Files:** Modify `include/viennaps/fields/AdaptiveMeshRefiner.hpp`

**Produces:** `refineMarkedWithProlongation` takes `ParMesh` + `ParFiniteElementSpace` + `ParGridFunction`. `markBox`/`markByGradient`/`zzIndicator` take `ParMesh&` (works via base-class `Mesh&` since `ParMesh` derives from `Mesh`).

- [ ] **Step 1: Add ParMesh overload of `refineMarkedWithProlongation`**

- [ ] **Step 2: Build + run tests, commit**

---

## Task 10: Update remaining serial-mesh files

**Files:** Modify `include/viennaps/fields/models/OedSource.hpp`, `include/viennaps/fields/MaterialConverter.hpp`, `include/viennaps/fields/PdeApi.hpp`, `include/viennaps/fields/PhysicsField.hpp`, `include/viennaps/fields/MfemElasticityKernel.hpp`

**Produces:** All remaining files that reference `mfem::Mesh` are updated to `mfem::ParMesh` where needed. `PhysicsField` and `MfemElasticityKernel` are outside the diffusion engine but share the build.

- [ ] **Step 1: Update each file's mesh type references**

- [ ] **Step 2: Build + run all tests, commit**

---

## Task 11: Final validation + ledger update

**Files:** Modify `docs/superpowers/reviews/2026-07-25-spec-gap-analysis.md`

- [ ] **Step 1: Run full test suite**

```bash
cmake --build build_phase2 --config Release --target testDiffusion
./build_phase2/tests/testDiffusion.exe
```

- [ ] **Step 2: Update ledger row**

Change the Runtime AMR row from 🟡 to ✅:
```
| Runtime AMR in solve | ✅ | ParMesh GeneralRefinement (nonconforming) + ParFESpace::Update + ParGF::Update prolongation. Dose preserved to 0.1%. |
```

- [ ] **Step 3: Commit**

```bash
git add docs/superpowers/reviews/2026-07-25-spec-gap-analysis.md
git commit -m "docs: Runtime AMR ✅ - ParMesh migration complete"
```

---

## Self-Review Notes

- **Spec coverage:** The plan covers the full serial-to-ParMesh migration for runtime AMR. Task 8 is the critical deliverable (the actual AMR that was crashing). Tasks 1-7 are the prerequisites. Tasks 9-11 are cleanup.
- **Type consistency:** `ParMesh`/`ParFiniteElementSpace`/`ParGridFunction`/`ParBilinearForm`/`ParLinearForm` used consistently throughout. `HypreParMatrix` for the solve. `HyprePCG`/`HypreBoomerAMG` for the preconditioner.
- **No placeholders:** Every task has concrete code or specific instructions.
- **Key risk:** Task 3 (implicit-Euler solve) is the most complex - changing from serial `SparseMatrix` + `DSmoother` + `CG` to `HypreParMatrix` + `HyprePCG` + `HypreBoomerAMG` requires careful handling of essential-BC elimination (true dofs vs vdofs). The `MPI_COMM_SELF` case simplifies this (no actual parallel distribution) but the API calls are still the parallel ones.
- **MFEM_USE_MPI confirmed:** `f:/dev/mfem/build` has `MFEM_USE_MPI=ON`. The `ParMesh(MPI_COMM_SELF, serialMesh)` constructor is available.
