# Diffusion Engine Phase 2: Fermi + ChargedFermi + Segregation + Solid Solubility

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans.

**Goal:** Add Fermi and ChargedFermi transport models (charge-state-dependent diffusivity), segregation boundary conditions at material interfaces, and solid solubility deactivation.

**Architecture:** Extends Phase 1 `DiffusionModel` base. Fermi model computes D_eff from local electron/hole concentration (n,p) and intrinsic carrier concentration n_i. Segregation implemented as MFEM boundary integrator on interface elements. Solid solubility caps active concentration, pushing excess into cluster species.

**Tech Stack:** C++20, MFEM, SUNDIALS, existing Phase 1 DiffusionEngine

**Depends on:** Phase 1 (DiffusionModel, DiffusionEngine, MeshAttributes)

## Global Constraints
Same as Phase 1. All MFEM-gated. Namespace `viennaps`. LLVM style.

## File Structure

| File | Responsibility |
|------|---------------|
| `include/viennaps/fields/models/FermiDiffusion.hpp` | Fermi model: D = D_i*(1+alpha*n/ni) + D_v*(1+beta*p/ni) |
| `include/viennaps/fields/models/ChargedFermiDiffusion.hpp` | Full charge-state decomposition |
| `include/viennaps/fields/models/SolidSolubility.hpp` | Caps active C, excess -> cluster |
| `include/viennaps/fields/models/Segregation.hpp` | Interface BC: C_2 = m(T)*C_1 |
| `include/viennaps/fields/IntrinsicCarrier.hpp` | n_i(T, material) calculator |
| `tests/diffusion/testDiffusionEngine.cpp` | Extended with Phase 2 tests |

---

### Task 1: IntrinsicCarrier Calculator

**Files:** Create `include/viennaps/fields/IntrinsicCarrier.hpp`

**Produces:** `IntrinsicCarrier<NumericType>` with `ni(T, material)`, `electronConcentration(C_dopant, T, material)`, `holeConcentration(C_dopant, T, material)`

- [ ] **Step 1: Write failing test** - `TestIntrinsicCarrier()`: assert `ni(300, "Si")` ~ 1e10 cm^-3, `ni(1273, "Si")` > 1e10 (increases with T). Verify `electronConcentration(1e17, 300, "Si")` ~ 1e17 (n-type).

- [ ] **Step 2: Run to verify failure** -> FAIL

- [ ] **Step 3: Implement** - `ni(T) = sqrt(Nc*Nv) * exp(-Eg/(2*kB*T))`. Nc, Nv, Eg from `MaterialPropertySystem`. `electronConcentration` = max(C_dopant, ni) for n-type. `holeConcentration` = ni^2 / electronConcentration.

- [ ] **Step 4: Run to verify pass** -> PASS

- [ ] **Step 5: Commit** - `"feat: add IntrinsicCarrier calculator for Fermi-level-dependent diffusion"`

---

### Task 2: FermiDiffusion Model

**Files:** Create `include/viennaps/fields/models/FermiDiffusion.hpp`

**Produces:** `FermiDiffusion<NumericType>` extending `DiffusionModel`. Computes per-element D_eff using `IntrinsicCarrier`. Uses `mfem::GridFunctionCoefficient` to evaluate local concentration -> n,p -> D_eff.

- [ ] **Step 1: Write failing test** - `TestFermiDiffusion()`: create model with D_i=1e-13, alpha=1.0. At high dopant (1e20, extrinsic), assert `getDiffusivity(C=1e20, T=1273)` > `getDiffusivity(C=1e15, T=1273)` (extrinsic enhancement).

- [ ] **Step 2: Run to verify failure** -> FAIL

- [ ] **Step 3: Implement** - `D_eff = D_i * (1 + alpha * n/ni)`. Override `assembleStiffness` to use `mfem::GridFunctionCoefficient` wrapping the concentration-dependent D. The coefficient evaluates D at each quadrature point from local u value.

Key assembly code pattern:
```cpp
void assembleStiffness(mfem::BilinearForm& K, const mfem::Vector& u,
                       int offset, const mfem::GridFunction* temp) const override {
  // Create a GridFunction view of this species' block in u
  mfem::GridFunction conc(fes_, const_cast<double*>(u.GetData() + offset));
  // Custom coefficient: D(conc) = D_i * (1 + alpha * n(conc) / ni(T))
  FermiDCoef coef(D_i_, alpha_, ni_, T_);
  coef.SetConcentrationField(&conc);
  K.AddDomainIntegrator(new mfem::DiffusionIntegrator(coef));
}
```

- [ ] **Step 4: Run to verify pass** -> PASS

- [ ] **Step 5: Commit** - `"feat: add FermiDiffusion model with charge-state-dependent diffusivity"`

---

### Task 3: ChargedFermiDiffusion Model

**Files:** Create `include/viennaps/fields/models/ChargedFermiDiffusion.hpp`

**Produces:** `ChargedFermiDiffusion<NumericType>` - full charge-state decomposition: `D = sum_z D^z * f^z(n, p, T)` where f^z are charge-state fractions computed from Fermi level.

- [ ] **Step 1: Write failing test** - `TestChargedFermi()`: verify D_eff transitions smoothly between intrinsic and extrinsic regimes. At C=ni, D ~ D_intrinsic. At C=1e20 n-type, D significantly different.

- [ ] **Step 2: Run to verify failure** -> FAIL

- [ ] **Step 3: Implement** - Compute Fermi level E_F from charge neutrality. Compute charge-state fractions f^z = exp(-z*(E_F-E_i)/kT) / sum. D = sum_z D^z * f^z. Uses `IntrinsicCarrier` for n_i, E_i.

- [ ] **Step 4: Run to verify pass** -> PASS

- [ ] **Step 5: Commit** - `"feat: add ChargedFermiDiffusion with full charge-state decomposition"`

---

### Task 4: SolidSolubility Model

**Files:** Create `include/viennaps/fields/models/SolidSolubility.hpp`

**Produces:** `SolidSolubility<NumericType>` - reaction model that caps active concentration. `C_active = min(C, C_ss(T))`. Excess pushed to a cluster species. Contributes to `assembleReaction`.

- [ ] **Step 1: Write failing test** - `TestSolidSolubility()`: set C_ss=1e20 at T=1000C. Initialize C=1e21. After one reaction step, assert C_active <= 1e20 and cluster species gained the excess.

- [ ] **Step 2: Run to verify failure** -> FAIL

- [ ] **Step 3: Implement** - `assembleReaction` iterates elements, for each: `C_excess = max(0, C - C_ss(T))`. R[species] -= C_excess/dt. R[cluster] += C_excess/dt. C_ss(T) = C_ss0 * exp(-Ea_ss/(kB*T)).

- [ ] **Step 4: Run to verify pass** -> PASS

- [ ] **Step 5: Commit** - `"feat: add SolidSolubility deactivation model"`

---

### Task 5: Segregation Boundary Condition

**Files:** Create `include/viennaps/fields/models/Segregation.hpp`

**Produces:** `Segregation<NumericType>` - MFEM boundary integrator on interface elements between two material attributes. Enforces `C_2 = m(T) * C_1` via penalty method.

- [ ] **Step 1: Write failing test** - `TestSegregation()`: 2-material mesh (Si + SiO2). Initialize Boron=1e18 in Si. After diffusion with segregation coefficient m=0.1, assert C_SiO2 ~ 0.1 * C_Si at interface.

- [ ] **Step 2: Run to verify failure** -> FAIL

- [ ] **Step 3: Implement** - Uses `mfem::BoundaryIntegrator` with penalty formulation. The segregation BC adds a penalty term: `penalty * (C_2 - m*C_1)` on interface boundary elements. The penalty coefficient is large (1e6) to enforce the constraint. `m(T)` from `ParameterDatabase`.

- [ ] **Step 4: Run to verify pass** -> PASS

- [ ] **Step 5: Commit** - `"feat: add Segregation boundary condition for material interfaces"`

---

### Task 6: Integration Test - Fermi Diffusion with Segregation

**Files:** Modify `tests/diffusion/testDiffusionEngine.cpp`

- [ ] **Step 1: Write integration test** - `TestFermiWithSegregation()`: 2-material mesh, FermiDiffusion for Boron in Si, ConstantDiffusion in SiO2, Segregation BC at interface. Solve 0->30s at 1000C. Verify dose conservation across both materials and concentration ratio at interface matches m.

- [ ] **Step 2: Run** -> PASS

- [ ] **Step 3: Commit** - `"test: add Fermi diffusion with segregation integration test"`
