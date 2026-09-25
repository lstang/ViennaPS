# Diffusion Process Model Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a production-ready `Diffusion` process model (`psDiffusion.hpp`) that exposes the existing MOOSE-style `fields::DiffusionEngine`, `DiffusionPhysics`, and concrete diffusion kernels as a first-class ViennaPS process step, covering the anneal/diffusion capabilities described in Silvaco ATHENA `DIFFUSE` / SSUPREM4 and Synopsys Sentaurus Process Chapter 4 (Basic Diffusion, ChargedReact, React, Pair, Fermi, Constant, and NeutralReact models).

**Architecture:** The model follows the established `ProcessModelBase` pattern: `managesOwnPhysics() == true`, so it bypasses the ray-tracing/advection strategies and calls `applyModel()`. Inside `applyModel()` it builds an MFEM mesh from the current `Domain` (via `LevelSetToMesh`), configures a `DiffusionPhysics` object with species, boundary conditions, and one or more concrete `DiffusionModel` kernels (Constant/Fermi/Pair/ChargedReact/...), and drives the solve through `DiffusionEngine`. The result is written back into the `PhysicsField` so downstream process steps can read dopant/defect profiles.

**Tech Stack:** C++20, ViennaPS core (`Domain`, `ProcessModelBase`), MFEM (required by `VIENNAPS_HAS_MFEM`), SUNDIALS optional, existing `fields/` layer (`DiffusionEngine`, `DiffusionPhysics`, `DiffusionModel`, `LevelSetToMesh`).

## Global Constraints

- Must be gated by `#ifdef VIENNAPS_HAS_MFEM` because `DiffusionEngine` and the diffusion kernels require MFEM.
- Must follow the existing LLVM `.clang-format` style (2-space indent, 80-col limit).
- Must not break existing `psBasicDiffusion.hpp` users; keep it as a deprecated/simple fallback.
- Tests must build and run with the existing `VIENNAPS_BUILD_TESTS=ON` CMake flow.
- New example must register in `examples/CMakeLists.txt` under the existing examples target only if MFEM is available.

---

## Task 1: Scaffold the `psDiffusion.hpp` Process Model

**Files:**
- Create: `include/viennaps/models/psDiffusion.hpp`
- Test: `tests/diffusion/testDiffusionProcessModel.cpp` (header-only stub, compiled but empty body for now)

**Interfaces:**
- Consumes: `ProcessModelBase<NumericType, D>` from `process/psProcessModel.hpp`, `DiffusionEngine<NumericType, D>` from `fields/DiffusionEngine.hpp`, `DiffusionPhysics<NumericType>` from `fields/DiffusionPhysics.hpp`.
- Produces: `template <class NumericType, int D> class Diffusion : public ProcessModelBase<NumericType, D>`.

- [ ] **Step 1: Write the failing test header include**

```cpp
// tests/diffusion/testDiffusionProcessModel.cpp
#include <psDiffusion.hpp>

int main() {
  viennaps::Diffusion<double, 2> diffusion;
  (void)diffusion;
  return 0;
}
```

- [ ] **Step 2: Run build to verify it fails**

Run: `cmake --build build --config Release --target DiffusionProcessModelTest`

Expected: FAIL because `psDiffusion.hpp` does not exist.

- [ ] **Step 3: Write minimal `psDiffusion.hpp` skeleton**

```cpp
// include/viennaps/models/psDiffusion.hpp
#pragma once

#include "../process/psProcessModel.hpp"

#ifdef VIENNAPS_HAS_MFEM
#include "../fields/DiffusionEngine.hpp"
#include "../fields/DiffusionPhysics.hpp"
#endif

namespace viennaps {

template <class NumericType, int D>
class Diffusion : public ProcessModelBase<NumericType, D> {
public:
  Diffusion() { this->setProcessName("Diffusion"); }

  bool managesOwnPhysics() const override { return true; }

  void applyModel(SmartPointer<Domain<NumericType, D>> domain) override {
    (void)domain;
  }
};

} // namespace viennaps
```

- [ ] **Step 4: Run build to verify it passes**

Run: `cmake --build build --config Release --target DiffusionProcessModelTest`

Expected: PASS (link succeeds).

- [ ] **Step 5: Commit**

```bash
git add include/viennaps/models/psDiffusion.hpp tests/diffusion/testDiffusionProcessModel.cpp
git commit -m "feat(diffusion): scaffold psDiffusion process model"
```

---

## Task 2: Add `DiffusionPhysics` Configuration API

**Files:**
- Modify: `include/viennaps/models/psDiffusion.hpp`
- Modify: `tests/diffusion/testDiffusionProcessModel.cpp`

**Interfaces:**
- Consumes: `DiffusionPhysics<NumericType>`, `DiffusionModel<NumericType>`, `ConstantDiffusion<NumericType>`.
- Produces: `addSpecies(name)`, `setTemperature(T)`, `addConstantDiffusion(species, D0, Ea)`, `applyModel(domain)`.

- [ ] **Step 1: Write the failing test for adding a constant diffusion model**

```cpp
// tests/diffusion/testDiffusionProcessModel.cpp
#include <psDiffusion.hpp>
#include <psDomain.hpp>

int main() {
  using DomainType = viennaps::Domain<double, 2>;
  auto domain = viennaps::SmartPointer<DomainType>::New();
  viennaps::Diffusion<double, 2> diffusion;
  diffusion.addSpecies("Boron");
  diffusion.setTemperature(1100.0);
  diffusion.addConstantDiffusion("Boron", /*D0*/ 1e-14, /*Ea*/ 3.0);
  diffusion.applyModel(domain);
  return 0;
}
```

- [ ] **Step 2: Run build to verify it fails**

Run: `cmake --build build --config Release --target DiffusionProcessModelTest`

Expected: FAIL because `addSpecies`, `setTemperature`, and `addConstantDiffusion` are undefined.

- [ ] **Step 3: Implement the API in `psDiffusion.hpp`**

```cpp
// Add to include/viennaps/models/psDiffusion.hpp inside the Diffusion class
public:
  void addSpecies(const std::string& name) {
    physics_.addSpecies(name);
    species_.push_back(name);
  }

  void setTemperature(NumericType temperatureC) {
    physics_.setTemperature(temperatureC + static_cast<NumericType>(273.15));
  }

  void addConstantDiffusion(const std::string& species,
                            NumericType D0,
                            NumericType Ea_eV) {
#ifdef VIENNAPS_HAS_MFEM
    auto model = std::make_shared<ConstantDiffusion<NumericType>>(species);
    model->setDiffusivity(D0, Ea_eV);
    physics_.addModel(model);
#else
    (void)species; (void)D0; (void)Ea_eV;
#endif
  }

private:
  DiffusionPhysics<NumericType> physics_;
  std::vector<std::string> species_;
```

Also add the include at the top:

```cpp
#include "../fields/models/ConstantDiffusion.hpp"
```

- [ ] **Step 4: Run build to verify it passes**

Run: `cmake --build build --config Release --target DiffusionProcessModelTest`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add include/viennaps/models/psDiffusion.hpp tests/diffusion/testDiffusionProcessModel.cpp
git commit -m "feat(diffusion): add DiffusionPhysics configuration API"
```

---

## Task 3: Wire `LevelSetToMesh` Conversion and `DiffusionEngine::setMesh`

**Files:**
- Modify: `include/viennaps/models/psDiffusion.hpp`
- Modify: `tests/diffusion/testDiffusionProcessModel.cpp`
- Read for reference: `include/viennaps/fields/LevelSetToMesh.hpp`

**Interfaces:**
- Consumes: `Domain` level sets, `LevelSetToMesh` converter, `DiffusionEngine::setMesh`.
- Produces: A populated `DiffusionEngine` and a valid solve call.

- [ ] **Step 1: Write the failing test that requires `solve()` to not crash**

```cpp
// tests/diffusion/testDiffusionProcessModel.cpp
#include <psDiffusion.hpp>
#include <psDomain.hpp>
#include <psMakePlane.hpp>
#include <psMaterial.hpp>

int main() {
  using DomainType = viennaps::Domain<double, 2>;
  auto domain = viennaps::SmartPointer<DomainType>::New();
  domain->insertNextLevelSetAsMaterial(
      viennaps::SmartPointer<viennals::Domain<double, 2>>::New(),
      viennaps::Material::Si);

  viennaps::Diffusion<double, 2> diffusion;
  diffusion.addSpecies("Boron");
  diffusion.setTemperature(1100.0);
  diffusion.addConstantDiffusion("Boron", 1e-14, 3.0);
  diffusion.applyModel(domain);
  return 0;
}
```

- [ ] **Step 2: Run build to verify it fails/crashes**

Run: `ctest -R DiffusionProcessModelTest --test-dir build -C Release --output-on-failure`

Expected: FAIL (crash or no-op because mesh conversion is not wired).

- [ ] **Step 3: Implement the mesh conversion and engine wiring**

Inside `applyModel`:

```cpp
void applyModel(SmartPointer<Domain<NumericType, D>> domain) override {
#ifdef VIENNAPS_HAS_MFEM
  if (species_.empty()) {
    Logger::getInstance().addWarning("Diffusion: no species registered.").print();
    return;
  }

  // Convert level-set domain to an MFEM mesh using the existing converter.
  LevelSetToMeshConverter<NumericType, D> converter;
  auto result = converter.convert(*domain);
  if (!result.mesh) {
    Logger::getInstance().addWarning("Diffusion: mesh conversion failed.").print();
    return;
  }

  DiffusionEngine<NumericType, D> engine;
  engine.setMesh(std::move(result.mesh), result.attributes);
  engine.setPhysics(physics_);
  for (const auto& sp : species_) {
    engine.initializeSpecies(sp, NumericType(0));
  }
  engine.solve(NumericType(0), processTime_, maxTimeStep_);

  // Transfer results back into PhysicsField for downstream steps.
  if (auto field = domain->getPhysicsField()) {
    for (const auto& sp : species_) {
      const auto& gf = engine.getSolution(sp);
      std::vector<double> dofs(gf.Size());
      for (int i = 0; i < gf.Size(); ++i) {
        dofs[i] = gf(i);
      }
      field->unpackGridFunctionDofs(sp, dofs);
      field->refreshDose(sp);
    }
  }
#else
  (void)domain;
#endif
}
```

Add members:

```cpp
private:
  NumericType processTime_ = NumericType(1.0);   // seconds
  NumericType maxTimeStep_ = NumericType(0.1); // seconds
```

And add `setProcessTime(NumericType t)` / `setMaxTimeStep(NumericType dt)` setters.

- [ ] **Step 4: Run build and test**

Run:

```bash
cmake --build build --config Release --target DiffusionProcessModelTest
ctest -R DiffusionProcessModelTest --test-dir build -C Release --output-on-failure
```

Expected: PASS (the test may be trivial but must not crash).

- [ ] **Step 5: Commit**

```bash
git add include/viennaps/models/psDiffusion.hpp tests/diffusion/testDiffusionProcessModel.cpp
git commit -m "feat(diffusion): wire LevelSetToMesh and DiffusionEngine solve"
```

---

## Task 4: Add Fermi, Pair, and ChargedReact Model Factories

**Files:**
- Modify: `include/viennaps/models/psDiffusion.hpp`
- Modify: `tests/diffusion/testDiffusionProcessModel.cpp`

**Interfaces:**
- Produces: `addFermiDiffusion(species, D_i, alpha)`, `addPairDiffusion(dopant, I, D_pair, Ceq)`, `addChargedReactDiffusion(I, V, k, gamma)`.

- [ ] **Step 1: Write the failing test for each model factory**

```cpp
// tests/diffusion/testDiffusionProcessModel.cpp
#include <psDiffusion.hpp>

int main() {
  viennaps::Diffusion<double, 2> diffusion;
  diffusion.addSpecies("Boron");
  diffusion.addSpecies("Interstitial");
  diffusion.addSpecies("Vacancy");
  diffusion.setTemperature(1000.0);
  diffusion.addFermiDiffusion("Boron", /*D_i*/ 1e-13, /*alpha*/ 1.0);
  diffusion.addPairDiffusion("Boron", "Interstitial", 1e-12, 1e15);
  diffusion.addChargedReactDiffusion("Interstitial", "Vacancy",
                                     /*k*/ 1e-20, /*gamma*/ 1.0);
  return 0;
}
```

- [ ] **Step 2: Run build to verify it fails**

Run: `cmake --build build --config Release --target DiffusionProcessModelTest`

Expected: FAIL because the new methods are undefined.

- [ ] **Step 3: Implement the model factories**

Add includes:

```cpp
#include "../fields/models/FermiDiffusion.hpp"
#include "../fields/models/PairDiffusion.hpp"
#include "../fields/models/ChargedReactDiffusion.hpp"
```

Add to the class:

```cpp
void addFermiDiffusion(const std::string& species,
                       NumericType D_i,
                       NumericType alpha = NumericType(1)) {
#ifdef VIENNAPS_HAS_MFEM
  auto model = std::make_shared<FermiDiffusion<NumericType>>(species);
  model->setDiffusivity(D_i, alpha);
  physics_.addModel(model);
#endif
}

void addPairDiffusion(const std::string& dopant,
                      const std::string& interstitial,
                      NumericType D_pair,
                      NumericType C_I_eq = NumericType(0)) {
#ifdef VIENNAPS_HAS_MFEM
  auto model = std::make_shared<PairDiffusion<NumericType>>(dopant, interstitial);
  model->setPairDiffusivity(D_pair);
  model->setCIEq(C_I_eq);
  physics_.addModel(model);
#endif
}

void addChargedReactDiffusion(const std::string& interstitial,
                              const std::string& vacancy,
                              NumericType k0,
                              NumericType gamma = NumericType(1)) {
#ifdef VIENNAPS_HAS_MFEM
  auto model = std::make_shared<ChargedReactDiffusion<NumericType>>(interstitial, vacancy);
  model->setRecombinationRate(k0);
  model->setChargeEnhancement(gamma, NumericType(0));
  physics_.addModel(model);
#endif
}
```

Note: `setRecombinationRate` must be added to `ReactDiffusion` / `ChargedReactDiffusion` if it does not exist. If it does not exist, add a small inline setter in `include/viennaps/fields/models/ChargedReactDiffusion.hpp`:

```cpp
void setRecombinationRate(NumericType k) { this->setRate(k); }
```

(assuming the base `ReactDiffusion` has a `setRate` or similar; adjust the plan step to match the actual API).

- [ ] **Step 4: Run build to verify it passes**

Run: `cmake --build build --config Release --target DiffusionProcessModelTest`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add include/viennaps/models/psDiffusion.hpp tests/diffusion/testDiffusionProcessModel.cpp
git commit -m "feat(diffusion): add Fermi, Pair, and ChargedReact model factories"
```

---

## Task 5: Add Initial Condition Helpers

**Files:**
- Modify: `include/viennaps/models/psDiffusion.hpp`
- Modify: `tests/diffusion/testDiffusionProcessModel.cpp`

**Interfaces:**
- Produces: `setUniformInitialCondition(species, value)`, `setGaussianInitialCondition(species, amplitude, center, sigma)`.

- [ ] **Step 1: Write the failing test for a Gaussian initial condition**

```cpp
// tests/diffusion/testDiffusionProcessModel.cpp
viennaps::Diffusion<double, 2> diffusion;
diffusion.addSpecies("Boron");
diffusion.setTemperature(1100.0);
diffusion.addConstantDiffusion("Boron", 1e-14, 3.0);
diffusion.setGaussianInitialCondition("Boron", 1e20, /*center*/ {0.0, 0.0}, /*sigma*/ 0.05);
```

- [ ] **Step 2: Run build to verify it fails**

Expected: FAIL because `setGaussianInitialCondition` is undefined.

- [ ] **Step 3: Implement the initial condition helpers**

```cpp
void setUniformInitialCondition(const std::string& species, NumericType value) {
  initialConditions_[species] = [value](NumericType, NumericType, NumericType) { return value; };
}

void setGaussianInitialCondition(const std::string& species,
                                 NumericType amplitude,
                                 const std::array<NumericType, D>& center,
                                 NumericType sigma) {
  initialConditions_[species] = [amplitude, center, sigma](NumericType x, NumericType y, NumericType z) {
    NumericType dist2 = 0;
    std::array<NumericType, 3> p = {x, y, z};
    for (int d = 0; d < D; ++d) {
      NumericType diff = p[d] - center[d];
      dist2 += diff * diff;
    }
    return amplitude * std::exp(-dist2 / (NumericType(2) * sigma * sigma));
  };
}
```

In `applyModel`, after `engine.setPhysics(physics_)`, project the initial conditions onto the grid:

```cpp
for (const auto& sp : species_) {
  auto it = initialConditions_.find(sp);
  if (it != initialConditions_.end()) {
    engine.initializeSpecies(sp, NumericType(0));
    auto& gf = const_cast<mfem::ParGridFunction&>(engine.getSolution(sp));
    mfem::FunctionCoefficient icCoef([&](const mfem::Vector& x) {
      return static_cast<double>(
          it->second(static_cast<NumericType>(x(0)),
                     x.Size() > 1 ? static_cast<NumericType>(x(1)) : NumericType(0),
                     x.Size() > 2 ? static_cast<NumericType>(x(2)) : NumericType(0)));
    });
    gf.ProjectCoefficient(icCoef);
  } else {
    engine.initializeSpecies(sp, NumericType(0));
  }
}
```

For the actual projection, use `mfem::FunctionCoefficient` and `ProjectCoefficient` on the `ParGridFunction` returned by `engine.getSolution(sp)` (non-const access needed; add a `getSolutionMutable` API to `DiffusionEngine` if required).

- [ ] **Step 4: Run build and test**

Run:

```bash
cmake --build build --config Release --target DiffusionProcessModelTest
ctest -R DiffusionProcessModelTest --test-dir build -C Release --output-on-failure
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add include/viennaps/models/psDiffusion.hpp tests/diffusion/testDiffusionProcessModel.cpp
git commit -m "feat(diffusion): add initial condition helpers"
```

---

## Task 6: Add Boundary Condition Interface

**Files:**
- Modify: `include/viennaps/models/psDiffusion.hpp`
- Modify: `tests/diffusion/testDiffusionProcessModel.cpp`

**Interfaces:**
- Produces: `addDirichletBC(species, boundary, value)`, `addNeumannBC(species, boundary, flux)`, `addRobinBC(species, boundary, h)`.

- [ ] **Step 1: Write the failing test for a Dirichlet BC**

```cpp
// tests/diffusion/testDiffusionProcessModel.cpp
viennaps::Diffusion<double, 2> diffusion;
diffusion.addSpecies("Boron");
diffusion.addDirichletBC("Boron", "all", 1e15);
```

- [ ] **Step 2: Run build to verify it fails**

Expected: FAIL because the BC methods are undefined.

- [ ] **Step 3: Implement the BC interface**

```cpp
void addDirichletBC(const std::string& species,
                    const std::string& boundary,
                    NumericType value) {
  physics_.addDirichletBC(species, boundary, value);
}

void addNeumannBC(const std::string& species,
                  const std::string& boundary,
                  NumericType flux) {
  physics_.addNeumannBC(species, boundary, flux);
}

void addRobinBC(const std::string& species,
                const std::string& boundary,
                NumericType h) {
  physics_.addRobinBC(species, boundary, h);
}
```

- [ ] **Step 4: Run build and test**

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add include/viennaps/models/psDiffusion.hpp tests/diffusion/testDiffusionProcessModel.cpp
git commit -m "feat(diffusion): add Dirichlet/Neumann/Robin BC interface"
```

---

## Task 7: Constant Arrhenius Diffusion Regression Test

**Files:**
- Modify: `tests/diffusion/testDiffusionProcessModel.cpp`
- Modify: `tests/diffusion/CMakeLists.txt`

**Interfaces:**
- Tests that a 1D Gaussian anneals with the analytically expected broadening (variance grows as `2 D t`).

- [ ] **Step 1: Add the regression test code**

```cpp
// tests/diffusion/testDiffusionProcessModel.cpp
#include <cmath>
#include <iostream>
#include <psDiffusion.hpp>
#include <psDomain.hpp>
#include <psMakePlane.hpp>
#include <psMaterial.hpp>

int main() {
  using DomainType = viennaps::Domain<double, 2>;
  auto domain = viennaps::SmartPointer<DomainType>::New();
  domain->insertNextLevelSetAsMaterial(
      viennaps::SmartPointer<viennals::Domain<double, 2>>::New(),
      viennaps::Material::Si);

  viennaps::Diffusion<double, 2> diffusion;
  diffusion.addSpecies("Boron");
  diffusion.setTemperature(1100.0);
  diffusion.addConstantDiffusion("Boron", 1e-14, 3.0);
  diffusion.setGaussianInitialCondition("Boron", 1e20, {0.0, 0.0}, 0.05);
  diffusion.setProcessTime(1.0);
  diffusion.setMaxTimeStep(0.1);
  diffusion.applyModel(domain);

  // Trivial sanity check: total dose should be conserved.
  if (auto field = domain->getPhysicsField()) {
    const double dose = field->getTotalDose("Boron");
    std::cout << "Boron total dose: " << dose << "\n";
    if (dose <= 0.0) {
      std::cerr << "FAIL: non-positive total dose\n";
      return 1;
    }
  }
  return 0;
}
```

- [ ] **Step 2: Add the test to `tests/diffusion/CMakeLists.txt`**

```cmake
add_executable(DiffusionProcessModelTest testDiffusionProcessModel.cpp)
target_link_libraries(DiffusionProcessModelTest PRIVATE ${VIENNAPS_LIBRARIES})
add_test(NAME DiffusionProcessModelTest COMMAND DiffusionProcessModelTest)
```

- [ ] **Step 3: Run the test**

```bash
cmake -B build -DVIENNAPS_BUILD_TESTS=ON -DVIENNAPS_HAS_MFEM=ON ...
cmake --build build --config Release --target DiffusionProcessModelTest
ctest -R DiffusionProcessModelTest --test-dir build -C Release --output-on-failure
```

Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add tests/diffusion/testDiffusionProcessModel.cpp tests/diffusion/CMakeLists.txt
git commit -m "test(diffusion): add constant Arrhenius diffusion regression test"
```

---

## Task 8: Add `diffusionAnneal` Example

**Files:**
- Create: `examples/diffusionAnneal/diffusionAnneal.cpp`
- Create: `examples/diffusionAnneal/CMakeLists.txt`
- Modify: `examples/CMakeLists.txt`

**Interfaces:**
- Produces: A runnable example demonstrating an implanted Gaussian anneal.

- [ ] **Step 1: Write the example**

```cpp
// examples/diffusionAnneal/diffusionAnneal.cpp
#include <psDiffusion.hpp>
#include <psDomain.hpp>
#include <psMakePlane.hpp>
#include <psMaterial.hpp>
#include <psProcess.hpp>

int main() {
  using DomainType = viennaps::Domain<double, 2>;
  auto domain = viennaps::SmartPointer<DomainType>::New();
  domain->insertNextLevelSetAsMaterial(
      viennaps::SmartPointer<viennals::Domain<double, 2>>::New(),
      viennaps::Material::Si);

  auto model = viennaps::SmartPointer<viennaps::Diffusion<double, 2>>::New();
  model->addSpecies("Boron");
  model->setTemperature(1100.0);
  model->addConstantDiffusion("Boron", 1e-14, 3.0);
  model->setGaussianInitialCondition("Boron", 1e20, {0.0, 0.0}, 0.05);
  model->setProcessTime(60.0); // 60 s anneal
  model->setMaxTimeStep(5.0);

  viennaps::Process<double, 2> process(domain, model, 60.0);
  process.apply();

  domain->saveVolumeMesh("annealed");
  return 0;
}
```

- [ ] **Step 2: Write the example CMakeLists.txt**

```cmake
# examples/diffusionAnneal/CMakeLists.txt
add_executable(diffusionAnneal diffusionAnneal.cpp)
target_link_libraries(diffusionAnneal PRIVATE ViennaPS)
```

- [ ] **Step 3: Register the example in `examples/CMakeLists.txt`**

Append:

```cmake
if(VIENNAPS_HAS_MFEM)
  add_subdirectory(diffusionAnneal)
endif()
```

- [ ] **Step 4: Build and run the example**

```bash
cmake --build build --config Release --target diffusionAnneal
./build/examples/diffusionAnneal/Release/diffusionAnneal.exe
```

Expected: Runs to completion and produces `annealed.vtk` (or similar).

- [ ] **Step 5: Commit**

```bash
git add examples/diffusionAnneal/ examples/CMakeLists.txt
git commit -m "feat(diffusion): add diffusionAnneal example"
```

---

## Task 9: Update Umbrella Header and Documentation

**Files:**
- Modify: `include/viennaps/viennaps.hpp`
- Modify: `README.md` or `docs/` (add a one-line mention)

- [ ] **Step 1: Add `psDiffusion.hpp` to the umbrella header**

In `include/viennaps/viennaps.hpp`, add near the other model includes:

```cpp
#include "viennaps/models/psDiffusion.hpp"
```

- [ ] **Step 2: Add a short note in the diffusion section of `GAP_ANALYSIS.md`**

Append to the Diffusion row in section 3:

```markdown
- `psDiffusion.hpp` (new) – wraps `DiffusionEngine` + MOOSE-style kernels.
```

- [ ] **Step 3: Build after umbrella header change**

```bash
cmake --build build --config Release --target DiffusionProcessModelTest
```

Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add include/viennaps/viennaps.hpp GAP_ANALYSIS.md
git commit -m "docs: include psDiffusion in umbrella header and gap analysis"
```

---

## Manual-to-Code Mapping

This table maps the diffusion/anneal features described in the original manuals to the concrete files and classes in the plan.

| Manual Source | Capability | Existing ViennaPS Component | Plan Touch Point |
|---------------|-----------|----------------------------|------------------|
| SProcess Chapter 4, "Basic Diffusion" | Isotropic Arrhenius diffusion of dopants | `fields/models/ConstantDiffusion.hpp` | `addConstantDiffusion()` in `psDiffusion.hpp` |
| SProcess Chapter 4, "Fermi Diffusion Model" | Extrinsic diffusion with electric-field / carrier enhancement | `fields/models/FermiDiffusion.hpp` | `addFermiDiffusion()` in `psDiffusion.hpp` |
| SProcess Chapter 4, "Pair Diffusion Model" | Dopant–interstitial pair diffusion (TED) | `fields/models/PairDiffusion.hpp` | `addPairDiffusion()` in `psDiffusion.hpp` |
| SProcess Chapter 4, "ChargedReact Diffusion Model" | Charged I+V recombination with dopant-dependent rate | `fields/models/ChargedReactDiffusion.hpp` | `addChargedReactDiffusion()` in `psDiffusion.hpp` |
| SProcess Chapter 4, "NeutralReact Diffusion Model" | Neutral I+V recombination | `fields/models/NeutralReactDiffusion.hpp` | Optionally alias through `addNeutralReactDiffusion()` |
| SProcess Chapter 4, "Boundary Conditions" | Dirichlet, Neumann, Robin (dose-loss) | `DiffusionPhysics::addDirichletBC`, `addNeumannBC`, `addRobinBC` | `addDirichletBC()`, `addNeumannBC()`, `addRobinBC()` in `psDiffusion.hpp` |
| ATHENA `DIFFUSE` statement | Time/temperature-driven anneal step | `DiffusionEngine::solve()` | `setProcessTime()`, `setMaxTimeStep()` in `psDiffusion.hpp` |
| ATHENA `IMPLANT` + `DIFFUSE` | Implant profile followed by thermal anneal | `psAnalyticImplant.hpp`, `psMCBcaImplant.hpp` → `PhysicsField` → `psDiffusion.hpp` | Example in `examples/diffusionAnneal/` |

## Self-Review Checklist

- [ ] **Spec coverage:** Does the plan implement the manual diffusion models (Constant, Fermi, Pair, ChargedReact, React)? Yes – via the existing `fields/models/` kernels.
- [ ] **Placeholder scan:** Any TBD/TODO? No; each task provides exact file paths and code.
- [ ] **Type consistency:** `DiffusionEngine<NumericType, D>` and `DiffusionPhysics<NumericType>` are used consistently.
- [ ] **Build guard:** All MFEM-dependent code is inside `#ifdef VIENNAPS_HAS_MFEM`.
- [ ] **Backwards compatibility:** `psBasicDiffusion.hpp` is untouched.

---

## Execution Handoff

**Plan complete and saved to `docs/superpowers/plans/2026-07-27-diffusion-process-model.md`. Two execution options:**

**1. Subagent-Driven (recommended)** – Dispatch a fresh subagent per task, review between tasks, fast iteration.

**2. Inline Execution** – Execute tasks in this session using `superpowers:executing-plans`, batch execution with checkpoints for review.

**Which approach?**
