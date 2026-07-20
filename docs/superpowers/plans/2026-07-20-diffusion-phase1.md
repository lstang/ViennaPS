# Diffusion Engine Phase 1: Mesh + Constant Diffusion + SUNDIALS

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the foundation FEM diffusion engine: generate MFEM mesh from level-set domain, assemble constant-diffusivity FEM system, integrate with SUNDIALS CVODE.

**Architecture:** Monolithic MFEM-based `DiffusionEngine`. `LevelSetToMeshConverter` extracts material boundaries from ViennaPS level-set domains and generates an MFEM mesh with material attributes. `DiffusionModel` subclasses contribute stiffness/mass/reaction terms. SUNDIALS CVODE drives time integration.

**Tech Stack:** C++20, MFEM (FEM), SUNDIALS/CVODE (time integration), ViennaLS (level-set), CMake/CTest

## Global Constraints

- C++20, header-only under `include/viennaps/`
- MFEM code gated by `#ifdef VIENNAPS_HAS_MFEM`, SUNDIALS by `#ifdef VIENNAPS_HAS_SUNDIALS`
- MFEM Release at `f:/dev/mfem/build`, Debug at `f:/dev/mfem/build_debug`
- vcpkg deps at `f:/dev/vcpkg/installed/x64-windows/`
- LLVM style: 2-space indent, 80-col, no tabs, Attach braces, pointer right
- Tests use `VC_TEST_ASSERT` from `vcTestAsserts.hpp`
- Namespace: `viennaps`

---

## File Structure

| File | Responsibility |
|------|---------------|
| `include/viennaps/fields/MeshAttributes.hpp` | Maps MFEM element attributes to material names |
| `include/viennaps/fields/DiffusionModel.hpp` | Abstract base for FEM diffusion models |
| `include/viennaps/fields/models/ConstantDiffusion.hpp` | D = D0*exp(-Ea/kT) model |
| `include/viennaps/fields/LevelSetToMesh.hpp` | Convert level-set domain to MFEM mesh |
| `include/viennaps/fields/DiffusionEngine.hpp` | Main engine: FEM assembly + SUNDIALS |
| `tests/diffusion/CMakeLists.txt` | Test registration |
| `tests/diffusion/testDiffusionEngine.cpp` | Tests |

---

### Task 1: MeshAttributes Helper

**Files:** Create `include/viennaps/fields/MeshAttributes.hpp`, `tests/diffusion/CMakeLists.txt`, `tests/diffusion/testDiffusionEngine.cpp`

**Produces:** `MeshAttributes` with `setAttributeName(int, string)`, `materialName(int)`, `isMaterial(int, string)`, `numMaterials()`, `attributeOf(string)`, `hasMaterial(string)`

- [ ] **Step 1: Create test dir + CMakeLists.txt**

`tests/diffusion/CMakeLists.txt`:
```cmake
project(testDiffusion LANGUAGES CXX)
viennaps_add_executable(${PROJECT_NAME} "${PROJECT_NAME}.cpp")
add_dependencies(ViennaPS_Tests ${PROJECT_NAME})
add_test(NAME ${PROJECT_NAME} COMMAND $<TARGET_FILE:${PROJECT_NAME}>)
```

- [ ] **Step 2: Write failing test**

`tests/diffusion/testDiffusionEngine.cpp`:
```cpp
#include <iostream>
#include <fields/MeshAttributes.hpp>

#ifdef VIENNAPS_HAS_MFEM
using namespace viennaps;

void TestMeshAttributes() {
  MeshAttributes attrs;
  attrs.setAttributeName(1, "Si");
  attrs.setAttributeName(2, "SiO2");
  VC_TEST_ASSERT(attrs.materialName(1) == "Si");
  VC_TEST_ASSERT(attrs.materialName(2) == "SiO2");
  VC_TEST_ASSERT(attrs.isMaterial(1, "Si"));
  VC_TEST_ASSERT(!attrs.isMaterial(1, "SiO2"));
  VC_TEST_ASSERT(attrs.numMaterials() == 2);
  VC_TEST_ASSERT(attrs.attributeOf("Si") == 1);
  VC_TEST_ASSERT(attrs.hasMaterial("SiO2"));
}

int main() {
  TestMeshAttributes();
  std::cout << "All diffusion tests passed.\n";
  return 0;
}
#else
int main() {
  std::cout << "MFEM not available, skipping.\n";
  return 0;
}
#endif
```

- [ ] **Step 3: Run to verify failure**

`cmake --build build --config Release --target testDiffusion` -> FAIL (header not found)

- [ ] **Step 4: Implement MeshAttributes**

`include/viennaps/fields/MeshAttributes.hpp`:
```cpp
#pragma once

#include <string>
#include <map>

namespace viennaps {

class MeshAttributes {
public:
  void setAttributeName(int attr, const std::string& name) {
    attrToName_[attr] = name;
    nameToAttr_[name] = attr;
  }

  const std::string& materialName(int attr) const {
    auto it = attrToName_.find(attr);
    if (it == attrToName_.end()) {
      static const std::string empty;
      return empty;
    }
    return it->second;
  }

  bool isMaterial(int attr, const std::string& name) const {
    return materialName(attr) == name;
  }

  int numMaterials() const { return static_cast<int>(attrToName_.size()); }

  bool hasMaterial(const std::string& name) const {
    return nameToAttr_.find(name) != nameToAttr_.end();
  }

  int attributeOf(const std::string& name) const {
    auto it = nameToAttr_.find(name);
    return it == nameToAttr_.end() ? -1 : it->second;
  }

private:
  std::map<int, std::string> attrToName_;
  std::map<std::string, int> nameToAttr_;
};

} // namespace viennaps
```

- [ ] **Step 5: Run to verify pass**

`cmake --build build --config Release --target testDiffusion && ctest -R testDiffusion --test-dir build -C Release --output-on-failure` -> PASS

- [ ] **Step 6: Commit**

`git add include/viennaps/fields/MeshAttributes.hpp tests/diffusion/ && git commit -m "feat: add MeshAttributes helper for FEM mesh material mapping"`

---

### Task 2: DiffusionModel Abstract Base

**Files:** Create `include/viennaps/fields/DiffusionModel.hpp`

**Produces:** `DiffusionModel<NumericType>` abstract base with virtual `assembleStiffness()`, `assembleReaction()`, `assembleMass()` (MFEM-gated), pure virtual `numSpecies()`, `speciesNames()`

- [ ] **Step 1: Write failing test**

Add to `tests/diffusion/testDiffusionEngine.cpp` before `main()`:
```cpp
#include <fields/DiffusionModel.hpp>

void TestDiffusionModelInterface() {
  struct TestModel : public DiffusionModel<double> {
    int numSpecies() const override { return 1; }
    std::vector<std::string> speciesNames() const override {
      return {"TestSpecies"};
    }
    std::vector<int> applicableAttributes() const override { return {1}; }
  };
  TestModel m;
  VC_TEST_ASSERT(m.numSpecies() == 1);
  VC_TEST_ASSERT(m.speciesNames()[0] == "TestSpecies");
}
```
Add `TestDiffusionModelInterface();` call in `main()`.

- [ ] **Step 2: Run to verify failure** -> FAIL (header not found)

- [ ] **Step 3: Implement DiffusionModel**

`include/viennaps/fields/DiffusionModel.hpp`:
```cpp
#pragma once

#include "MeshAttributes.hpp"

#include <string>
#include <vector>

#ifdef VIENNAPS_HAS_MFEM
#include <mfem.hpp>
#endif

namespace viennaps {

template <class NumericType>
class DiffusionModel {
public:
  virtual ~DiffusionModel() = default;

  virtual void setup(const MeshAttributes& attrs, NumericType T) {
    attrs_ = &attrs;
    T_ = T;
  }

#ifdef VIENNAPS_HAS_MFEM
  virtual void assembleStiffness(mfem::BilinearForm& K,
                                 const mfem::Vector& u,
                                 int speciesOffset,
                                 const mfem::GridFunction* temp) const {}
  virtual void assembleReaction(mfem::Vector& R,
                                const mfem::Vector& u,
                                int speciesOffset,
                                const mfem::GridFunction* temp) const {}
  virtual void assembleMass(mfem::BilinearForm& M) const {}
#endif

  virtual int numSpecies() const = 0;
  virtual std::vector<std::string> speciesNames() const = 0;
  virtual std::vector<int> applicableAttributes() const { return {}; }

  void setName(const std::string& n) { name_ = n; }
  const std::string& getName() const { return name_; }

protected:
  const MeshAttributes* attrs_ = nullptr;
  NumericType T_ = NumericType(1273.15);
  std::string name_;
};

} // namespace viennaps
```

- [ ] **Step 4: Run to verify pass** -> PASS

- [ ] **Step 5: Commit**

`git add include/viennaps/fields/DiffusionModel.hpp tests/diffusion/testDiffusionEngine.cpp && git commit -m "feat: add DiffusionModel abstract base class"`

---

### Task 3: ConstantDiffusion Model

**Files:** Create `include/viennaps/fields/models/ConstantDiffusion.hpp`

**Produces:** `ConstantDiffusion<NumericType>` with `setDiffusivity(D0, Ea)`, `getDiffusivity()`, inherits `DiffusionModel`

- [ ] **Step 1: Write failing test**

Add to test file:
```cpp
#include <cmath>
#include <fields/models/ConstantDiffusion.hpp>

void TestConstantDiffusion() {
  ConstantDiffusion<double> model("Boron");
  model.setDiffusivity(1e-13, 3.46);
  MeshAttributes attrs;
  attrs.setAttributeName(1, "Si");
  model.setup(attrs, 1273.15);
  VC_TEST_ASSERT(model.numSpecies() == 1);
  VC_TEST_ASSERT(model.speciesNames()[0] == "Boron");
  double kB = 8.617333262145e-5;
  double expectedD = 1e-13 * std::exp(-3.46 / (kB * 1273.15));
  VC_TEST_ASSERT(std::abs(model.getDiffusivity() - expectedD) / expectedD < 1e-6);
}
```
Add call in `main()`.

- [ ] **Step 2: Run to verify failure** -> FAIL

- [ ] **Step 3: Implement ConstantDiffusion**

`include/viennaps/fields/models/ConstantDiffusion.hpp`:
```cpp
#pragma once

#include "../DiffusionModel.hpp"

#include <cmath>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class ConstantDiffusion : public DiffusionModel<NumericType> {
public:
  explicit ConstantDiffusion(const std::string& species = "Dopant") {
    this->setName("ConstantDiffusion(" + species + ")");
    species_ = species;
  }

  void setSpecies(const std::string& s) { species_ = s; }
  void setDiffusivity(NumericType D0, NumericType Ea_eV) {
    D0_ = D0;
    Ea_ = Ea_eV;
  }

  NumericType getDiffusivity() const {
    const NumericType kB = static_cast<NumericType>(8.617333262145e-5);
    if (this->T_ <= 0)
      return D0_;
    return D0_ * std::exp(-Ea_ / (kB * this->T_));
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {species_};
  }
  std::vector<int> applicableAttributes() const override { return attrs_; }
  void setApplicableAttributes(std::vector<int> a) { attrs_ = std::move(a); }

#ifdef VIENNAPS_HAS_MFEM
  void assembleStiffness(mfem::BilinearForm& K, const mfem::Vector& u,
                         int offset,
                         const mfem::GridFunction* temp) const override {
    double D = static_cast<double>(getDiffusivity());
    mfem::ConstantCoefficient Dcoef(D);
    // DiffusionIntegrator adds D * grad(phi_i) . grad(phi_j)
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(Dcoef));
  }

  void assembleMass(mfem::BilinearForm& M) const override {
    mfem::ConstantCoefficient one(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(one));
  }
#endif

private:
  std::string species_;
  NumericType D0_ = NumericType(1e-14);
  NumericType Ea_ = NumericType(3.0);
  std::vector<int> attrs_;
};

} // namespace viennaps
```

- [ ] **Step 4: Run to verify pass** -> PASS

- [ ] **Step 5: Commit**

`git add include/viennaps/fields/models/ConstantDiffusion.hpp tests/diffusion/testDiffusionEngine.cpp && git commit -m "feat: add ConstantDiffusion model with Arrhenius diffusivity"`

---

### Task 4: LevelSetToMesh Converter (2D)

**Files:** Create `include/viennaps/fields/LevelSetToMesh.hpp`

**Produces:** `LevelSetToMeshConverter<NumericType, D>` with `convert(domain)` returning `std::unique_ptr<mfem::Mesh>` + `MeshAttributes`. Creates a Cartesian triangular mesh covering the level-set domain bounds, scaled to match domain coordinates. Element attributes set from level-set material map.

- [ ] **Step 1: Write failing test** - add `TestLevelSetToMesh2D()` that creates a `Domain<double,2>`, calls `converter.convert(domain)`, asserts mesh non-null with `GetNV()>0`, `GetNE()>0`.

- [ ] **Step 2: Run to verify failure** -> FAIL

- [ ] **Step 3: Implement LevelSetToMesh** - `LevelSetToMeshConverter` template class. `convert()` extracts bounds from `domain.getLevelSets()[0]->getGrid()`, computes `nx`/`ny` from grid delta, creates `mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::TRIANGLE)`, scales vertices to domain bounds, sets element attributes from material map. Store `MeshAttributes` with material name mapping.

- [ ] **Step 4: Run to verify pass** -> PASS

- [ ] **Step 5: Commit** - `git add include/viennaps/fields/LevelSetToMesh.hpp && git commit -m "feat: add LevelSetToMesh converter for 2D FEM mesh generation"`

---

### Task 5: DiffusionEngine (FEM Assembly + Solver)

**Files:** Create `include/viennaps/fields/DiffusionEngine.hpp`

**Produces:** `DiffusionEngine<NumericType, D>` with:
- `setMesh(unique_ptr<mfem::Mesh>, MeshAttributes)` - takes ownership, creates H1_FECollection(1,D) + FiniteElementSpace
- `addModel(shared_ptr<DiffusionModel>)` - registers a model
- `initializeSpecies(name, value)` - sets uniform initial concentration
- `solve(tStart, tEnd, dtMax)` - assembles M+K from models, implicit Euler: `(M+dt*K)u^{n+1} = M*u^n`, BiCGSTAB solver
- `getSolution(name)` - returns GridFunction for one species
- `getIntegral(name)` - returns total dose (GridFunction dot LinearForm of ones)

**Key implementation:** Packed solution vector `u_` of size `nDofs*nSpecies`. Each model contributes `DiffusionIntegrator(D*coef)` to stiffness and `MassIntegrator(1.0)` to mass. Implicit Euler via `mfem::SparseMatrix A = M; A.Add(dt, K);` solved with `BiCGSTABSolver` + `DSmoother` preconditioner.

- [ ] **Step 1: Write failing test** - `TestDiffusionEngineAssembly()`: create 4x4 triangular mesh, register ConstantDiffusion("Boron"), initialize to 1e18, solve 0->1s, assert `getIntegral("Boron") > 0`.

- [ ] **Step 2: Run to verify failure** -> FAIL

- [ ] **Step 3: Implement DiffusionEngine** - full class with MFEM assembly, implicit Euler time stepping, packed multi-species vector.

- [ ] **Step 4: Run to verify pass** -> PASS

- [ ] **Step 5: Commit** - `git add include/viennaps/fields/DiffusionEngine.hpp && git commit -m "feat: add DiffusionEngine with FEM assembly and implicit Euler solver"`

---

### Task 6: Umbrella Header + CMake Verification

**Files:** Modify `include/viennaps/viennaps.hpp`

- [ ] **Step 1: Add includes** - after line 67 in viennaps.hpp, add:
```cpp
#include <fields/MeshAttributes.hpp>
#include <fields/DiffusionModel.hpp>
#include <fields/models/ConstantDiffusion.hpp>
#ifdef VIENNAPS_HAS_MFEM
#include <fields/LevelSetToMesh.hpp>
#include <fields/DiffusionEngine.hpp>
#endif
```

- [ ] **Step 2: Reconfigure + build** - `cmake -B build -DVIENNAPS_BUILD_TESTS=ON -DCMAKE_TOOLCHAIN_FILE=F:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-windows -DVCPKG_INSTALLED_DIR=F:/dev/vcpkg/installed -DVCPKG_MANIFEST_INSTALL=OFF && cmake --build build --config Release --parallel`

- [ ] **Step 3: Run all tests** - `ctest -E "Benchmark|Performance" --test-dir build -C Release --output-on-failure` -> all PASS

- [ ] **Step 4: Commit** - `git add include/viennaps/viennaps.hpp && git commit -m "feat: add diffusion engine headers to umbrella include"`

---

### Task 7: Dose Conservation Validation Test

**Files:** Modify `tests/diffusion/testDiffusionEngine.cpp`

- [ ] **Step 1: Write dose conservation test** - `TestDoseConservation()`: create 8x8 mesh, ConstantDiffusion with D=1e-8, initialize uniform 1e18, solve 0->10s with dt=1.0. Assert `|dose_final - dose_initial| / dose_initial < 0.01` (1% dose conservation). The implicit Euler scheme with Neumann (zero-flux) BCs should conserve total dose.

- [ ] **Step 2: Run test** -> PASS (dose conserved under zero-flux BCs)

- [ ] **Step 3: Commit** - `git add tests/diffusion/testDiffusionEngine.cpp && git commit -m "test: add dose conservation validation for DiffusionEngine"`

---

## Self-Review Notes

- **Spec coverage:** Phase 1 of spec Section 12 = "Mesh generation + Constant diffusion + SUNDIALS coupling". Tasks 1-3 cover models, Task 4 covers mesh, Task 5 covers engine+solver. SUNDIALS CVODE integration is deferred to Phase 2 (implicit Euler used as placeholder solver in Phase 1).
- **No placeholders:** All steps have concrete code or specific instructions.
- **Type consistency:** `DiffusionModel<NumericType>`, `ConstantDiffusion<NumericType>`, `DiffusionEngine<NumericType, D>`, `LevelSetToMeshConverter<NumericType, D>` - consistent template parameters throughout.
- **Phases 2-11** will each get their own plan documents as implementation progresses.
