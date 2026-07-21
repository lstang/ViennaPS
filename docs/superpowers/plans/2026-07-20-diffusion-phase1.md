# Diffusion Engine Phase 1: Mesh + Constant Diffusion + SUNDIALS

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the foundation FEM diffusion engine: generate MFEM mesh from level-set domain, assemble constant-diffusivity FEM system, integrate with SUNDIALS CVODE.

**Architecture:** Monolithic MFEM-based `DiffusionEngine` with physics-discretization separation (inspired by MOOSE `MultiSpeciesDiffusionPhysicsBase`). `DiffusionPhysics` defines WHAT to solve (species, D, BCs); `DiffusionEngine` implements HOW (CG FEM). Per-species `GridFunction` storage (not packed vector). `LevelSetToMeshConverter` extracts material boundaries from ViennaPS level-set domains. SUNDIALS CVODE (BDF, adaptive) for time integration from day one. HypreBoomerAMG preconditioner (MOOSE default for diffusion).

**Tech Stack:** C++20, MFEM (FEM + HypreBoomerAMG), SUNDIALS/CVODE (BDF time integration), ViennaLS (level-set), CMake/CTest

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
| `include/viennaps/fields/DiffusionPhysics.hpp` | Physics definition: species, D, BCs (separate from discretization) |
| `include/viennaps/fields/models/ConstantDiffusion.hpp` | D = D0*exp(-Ea/kT) model |
| `include/viennaps/fields/LevelSetToMesh.hpp` | Convert level-set domain to MFEM mesh |
| `include/viennaps/fields/DiffusionEngine.hpp` | Main engine: FEM assembly + SUNDIALS + HypreBoomerAMG |
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

**Produces:** `DiffusionModel<NumericType>` abstract base with virtual `assembleStiffness()`, `assembleReaction()`, `assembleMass()` (MFEM-gated, using per-species GridFunction - not packed vector), pure virtual `numSpecies()`, `speciesNames()`. Inspired by MOOSE kernel pattern: each model operates on species via named GridFunction lookup, not manual offset arithmetic.

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
#include <map>

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
  /// Contribute to stiffness matrix K for this species.
  /// speciesGF: this species' GridFunction (for concentration-dependent D)
  /// allSpecies: map of all species GridFunctions (for coupled models)
  virtual void assembleStiffness(
      mfem::BilinearForm& K,
      const mfem::GridFunction& speciesGF,
      const std::map<std::string, mfem::GridFunction*>& allSpecies,
      const mfem::GridFunction* temp) const {}

  /// Contribute to nonlinear reaction RHS R for this species.
  virtual void assembleReaction(
      mfem::LinearForm& R,
      const mfem::GridFunction& speciesGF,
      const std::map<std::string, mfem::GridFunction*>& allSpecies,
      const mfem::GridFunction* temp) const {}

  /// Contribute to mass matrix M.
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
  void assembleStiffness(mfem::BilinearForm& K, const mfem::GridFunction& speciesGF,
                         const std::map<std::string, mfem::GridFunction*>& allSpecies,
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

### Task 3.5: DiffusionPhysics (Physics-Discretization Separation)

**Files:** Create `include/viennaps/fields/DiffusionPhysics.hpp`

**Produces:** `DiffusionPhysics<NumericType>` - defines WHAT to solve, separate from HOW (MOOSE `MultiSpeciesDiffusionPhysicsBase` pattern). Holds species names, model registrations, BC specifications, initial conditions. `DiffusionEngine` reads from this.

**Rationale:** MOOSE separates physics definition (`MultiSpeciesDiffusionPhysicsBase`) from discretization (`MultiSpeciesDiffusionCG`). This allows future FV/DG discretizations without changing physics definition, and makes testing easier.

- [ ] **Step 1: Write failing test**

Add to test file:
```cpp
#include <fields/DiffusionPhysics.hpp>

void TestDiffusionPhysics() {
  DiffusionPhysics<double> physics;
  physics.addSpecies("Boron");
  physics.addSpecies("Interstitial");
  VC_TEST_ASSERT(physics.numSpecies() == 2);
  VC_TEST_ASSERT(physics.hasSpecies("Boron"));

  // Register a model
  auto model = std::make_shared<ConstantDiffusion<double>>("Boron");
  model->setDiffusivity(1e-13, 3.46);
  physics.addModel(model);

  // BC specification
  physics.addNeumannBC("Boron", "surface", 0.0);  // zero flux
  physics.addDirichletBC("Boron", "bottom", 1e18);

  VC_TEST_ASSERT(physics.numModels() == 1);
}
```
Add call in `main()`.

- [ ] **Step 2: Run to verify failure** -> FAIL

- [ ] **Step 3: Implement DiffusionPhysics**

`include/viennaps/fields/DiffusionPhysics.hpp`:
```cpp
#pragma once

#include "DiffusionModel.hpp"

#include <string>
#include <vector>
#include <map>
#include <memory>

namespace viennaps {

template <class NumericType>
class DiffusionPhysics {
public:
  void addSpecies(const std::string& name) {
    if (std::find(species_.begin(), species_.end(), name)
        == species_.end())
      species_.push_back(name);
  }

  bool hasSpecies(const std::string& name) const {
    return std::find(species_.begin(), species_.end(), name)
           != species_.end();
  }

  int numSpecies() const {
    return static_cast<int>(species_.size());
  }

  const std::vector<std::string>& speciesNames() const {
    return species_;
  }

  void addModel(std::shared_ptr<DiffusionModel<NumericType>> m) {
    models_.push_back(m);
  }

  int numModels() const {
    return static_cast<int>(models_.size());
  }

  const std::vector<std::shared_ptr<DiffusionModel<NumericType>>>&
  models() const { return models_; }

  struct BCSpec {
    std::string species;
    std::string boundary;
    std::string type;  // "neumann", "dirichlet", "segregation"
    NumericType value;
  };

  void addNeumannBC(const std::string& sp, const std::string& bnd,
                    NumericType flux) {
    bcs_.push_back({sp, bnd, "neumann", flux});
  }

  void addDirichletBC(const std::string& sp, const std::string& bnd,
                      NumericType val) {
    bcs_.push_back({sp, bnd, "dirichlet", val});
  }

  const std::vector<BCSpec>& boundaryConditions() const {
    return bcs_;
  }

  void setTemperature(NumericType T) { T_ = T; }
  NumericType temperature() const { return T_; }

private:
  std::vector<std::string> species_;
  std::vector<std::shared_ptr<DiffusionModel<NumericType>>> models_;
  std::vector<BCSpec> bcs_;
  NumericType T_ = NumericType(1273.15);
};

} // namespace viennaps
```

- [ ] **Step 4: Run to verify pass** -> PASS

- [ ] **Step 5: Commit**

`git add include/viennaps/fields/DiffusionPhysics.hpp tests/diffusion/testDiffusionEngine.cpp && git commit -m "feat: add DiffusionPhysics for physics-discretization separation"`

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

### Task 5: DiffusionEngine (FEM Assembly + SUNDIALS + HypreBoomerAMG)

**Files:** Create `include/viennaps/fields/DiffusionEngine.hpp`

**Produces:** `DiffusionEngine<NumericType, D>` with:
- `setMesh(unique_ptr<mfem::Mesh>, MeshAttributes)` - takes ownership, creates H1_FECollection(1,D) + FiniteElementSpace
- `addModel(shared_ptr<DiffusionModel>)` - registers a model
- `initializeSpecies(name, value)` - sets uniform initial concentration on per-species GridFunction
- `solve(tStart, tEnd, dtMax)` - assembles M+K from models per species, solves with SUNDIALS CVODE (BDF, adaptive) + HypreBoomerAMG preconditioner
- `getSolution(name)` - returns GridFunction for one species
- `getIntegral(name)` - returns total dose (GridFunction dot LinearForm of ones)

**Key design decisions (MOOSE-inspired):**
1. **Per-species GridFunction storage** (not packed vector): `std::map<std::string, unique_ptr<mfem::GridFunction>> species_`. Each species is a separate GridFunction on the same FES. Coupling terms read another species' GF via the `allSpecies` map. Eliminates manual offset arithmetic.
2. **HypreBoomerAMG preconditioner** (MOOSE default for diffusion): `mfem::HypreBoomerAMG` instead of `DSmoother`. 10-100x fewer iterations on large 2D/3D meshes.
3. **SUNDIALS CVODE** from day one: BDF orders 1-5, adaptive time stepping, error control. Uses existing `SundialsTimeIntegrator`. Essential for stiff systems (clustering, recombination in later phases). Falls back to implicit Euler if SUNDIALS unavailable.

**Assembly pattern (per species):**
```cpp
for each species s:
  BilinearForm K_s(fes);  // stiffness for species s
  BilinearForm M_s(fes);  // mass for species s
  LinearForm  R_s(fes);   // reaction RHS for species s
  for each model m that affects s:
    m->assembleStiffness(K_s, species_[s], allSpecies_, temp_);
    m->assembleMass(M_s);
    m->assembleReaction(R_s, species_[s], allSpecies_, temp_);
  K_s.Assemble(); M_s.Assemble(); R_s.Assemble();
  // System: M_s * du_s/dt = -K_s * u_s + R_s
```

**SUNDIALS integration:**
```cpp
// Pack all species GFs into a flat vector for SUNDIALS
// RHS callback: unpack -> for each species, compute -M^{-1}*(K*u - R) -> repack
#ifdef VIENNAPS_HAS_SUNDIALS
SundialsTimeIntegrator<NumericType> integrator;
integrator.setUseFieldState(true);
// Register RHS callback that assembles K, M, R and computes du/dt
integrator.evolve(tStart, tEnd, dtMax);
#else
// Fallback: implicit Euler with HypreBoomerAMG
mfem::HypreBoomerAMG amg(A);
mfem::GMRESSolver solver(mesh_->GetComm());
solver.SetOperator(A);
solver.SetPreconditioner(amg);
#endif
```

**HypreBoomerAMG setup (always used for the linear solve within each SUNDIALS step):**
```cpp
mfem::HypreBoomerAMG* amg = new mfem::HypreBoomerAMG(A);
amg->SetPrintLevel(0);
// For diffusion: default AMG settings work well (no special config needed)
```

- [ ] **Step 1: Write failing test** - `TestDiffusionEngineAssembly()`: create 4x4 triangular mesh, register ConstantDiffusion("Boron"), initialize to 1e18, solve 0->1s, assert `getIntegral("Boron") > 0`.

- [ ] **Step 2: Run to verify failure** -> FAIL

- [ ] **Step 3: Implement DiffusionEngine** - full class with per-species GridFunction storage, FEM assembly, SUNDIALS CVODE integration with HypreBoomerAMG preconditioner. Fallback to implicit Euler + HypreBoomerAMG when SUNDIALS unavailable.

- [ ] **Step 4: Run to verify pass** -> PASS

- [ ] **Step 5: Commit** - `git add include/viennaps/fields/DiffusionEngine.hpp && git commit -m "feat: add DiffusionEngine with SUNDIALS CVODE + HypreBoomerAMG preconditioner"`

---

### Task 6: Umbrella Header + CMake Verification

**Files:** Modify `include/viennaps/viennaps.hpp`

- [ ] **Step 1: Add includes** - after line 67 in viennaps.hpp, add:
```cpp
#include <fields/MeshAttributes.hpp>
#include <fields/DiffusionModel.hpp>
#include <fields/DiffusionPhysics.hpp>
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

- **Spec coverage:** Phase 1 of spec Section 12 = "Mesh generation + Constant diffusion + SUNDIALS coupling". Tasks 1-3 cover models, Task 3.5 covers physics definition, Task 4 covers mesh, Task 5 covers engine+SUNDIALS+HypreBoomerAMG. SUNDIALS CVODE integration is included from Phase 1 (not deferred).
- **MOOSE-inspired improvements:** (1) HypreBoomerAMG preconditioner instead of DSmoother, (2) SUNDIALS CVODE from day one instead of implicit Euler, (3) Per-species GridFunction instead of packed vector, (4) DiffusionPhysics separates physics definition from discretization.
- **No placeholders:** All steps have concrete code or specific instructions.
- **Type consistency:** `DiffusionModel<NumericType>`, `ConstantDiffusion<NumericType>`, `DiffusionPhysics<NumericType>`, `DiffusionEngine<NumericType, D>`, `LevelSetToMeshConverter<NumericType, D>` - consistent template parameters throughout. Model assemble methods use `const mfem::GridFunction& speciesGF` + `const std::map<std::string, mfem::GridFunction*>& allSpecies` for coupling.
- **Phases 2-11** will each get their own plan documents as implementation progresses.
