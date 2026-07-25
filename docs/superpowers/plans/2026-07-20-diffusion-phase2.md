# Diffusion Engine Phase 2: Fermi + ChargedFermi + Segregation + Solid Solubility

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans.

**Goal:** Add Fermi and ChargedFermi transport models (charge-state-dependent diffusivity), segregation boundary conditions at material interfaces, and solid solubility deactivation.

**Architecture:** Extends Phase 1 `DiffusionModel` base. Fermi model computes D_eff from local electron/hole concentration (n,p) and intrinsic carrier concentration n_i. Uses MOOSE `MatDiffusionBase` pattern: `DiffusivityMaterial` computes D(C,T) at quadrature points via `GridFunctionCoefficient`, Jacobian gets dD/dC automatically. Segregation implemented as `SegregationBC` using MOOSE `BinaryRecombinationBC` pattern (dynamic rate K_seg*C_1 - K_deseg*C_2, equilibrium gives m(T) = K_seg/K_deseg). Solid solubility caps active concentration, pushing excess into cluster species.

**Tech Stack:** C++20, MFEM, SUNDIALS, existing Phase 1 DiffusionEngine

**Depends on:** Phase 1 (DiffusionModel, DiffusionEngine, MeshAttributes)

## Global Constraints
Same as Phase 1. All MFEM-gated. Namespace `viennaps`. LLVM style.

## File Structure

| File | Responsibility |
|------|---------------|
| `include/viennaps/fields/IntrinsicCarrier.hpp` | n_i(T, material) calculator |
| `include/viennaps/fields/DiffusivityMaterial.hpp` | MOOSE MatDiffusion pattern: computes D(C,T) at quadrature points |
| `include/viennaps/fields/models/FermiDiffusion.hpp` | Fermi model: D = D_i*(1+alpha*n/ni) + D_v*(1+beta*p/ni) |
| `include/viennaps/fields/models/ChargedFermiDiffusion.hpp` | Full charge-state decomposition |
| `include/viennaps/fields/models/SolidSolubility.hpp` | Caps active C, excess -> cluster |
| `include/viennaps/fields/models/Segregation.hpp` | Interface BC: dynamic rate K_seg*C1 - K_deseg*C2 (BinaryRecombinationBC pattern) |
| `tests/diffusion/testDiffusionEngine.cpp` | Extended with Phase 2 tests |

---

### Task 1: IntrinsicCarrier Calculator

**Files:** Create `include/viennaps/fields/IntrinsicCarrier.hpp`

**Produces:** `IntrinsicCarrier<NumericType>` with `ni(T, material)`, `electronConcentration(C_dopant, T, material)`, `holeConcentration(C_dopant, T, material)`, and a hook for non-ideal activity coefficients.

**MOOSE reference for activity coefficients (verified):** `modules/chemical_reactions/include/kernels/CoupledDiffusionReactionSub.h` carries `_gamma_u`, `_gamma_v[]`, `_gamma_eq` (activity coefficients) — the canonical MOOSE idiom for non-ideal thermodynamics. At degenerate doping (C >> N_c), Boltzmann statistics break down and n ≠ C; the Fermi-Dirac activity coefficient γ_n(C,T) must be applied. Expose this hook now so Phase 2 Task 3 (`ChargedFermiDiffusion`) and Phase 4 Task 4 (`ChargedEquilibriumDiffusion`) can use it without re-architecting.

- [ ] **Step 1: Write failing test** - `TestIntrinsicCarrier()`: assert `ni(300, "Si")` ~ 1e10 cm^-3, `ni(1273, "Si")` > 1e10 (increases with T). Verify `electronConcentration(1e17, 300, "Si")` ~ 1e17 (n-type, Boltzmann regime). **Add:** at C=1e21 (degenerate), `electronConcentration` with `useFermiDirac=true` returns a value *less than* the Boltzmann-statistics value (activity γ < 1).

- [ ] **Step 2: Run to verify failure** -> FAIL

- [ ] **Step 3: Implement** - `ni(T) = sqrt(Nc*Nv) * exp(-Eg/(2*kB*T))`. Nc, Nv, Eg from `MaterialPropertySystem`. `electronConcentration` = max(C_dopant, ni) for n-type (Boltzmann). `holeConcentration` = ni^2 / electronConcentration. **Add** `activity(C, T, material, statistics)` returning γ ∈ (0,1] for Fermi-Dirac (default γ=1 for Boltzmann). Provide `electronConcentration(C, T, material, useFermiDirac)` overload that divides C by γ.

- [ ] **Step 4: Run to verify pass** -> PASS

- [ ] **Step 5: Commit** - `"feat: add IntrinsicCarrier calculator with Fermi-Dirac activity hook"`

---

### Task 2: DiffusivityMaterial + FermiDiffusion Model

**Files:** Create `include/viennaps/fields/DiffusivityMaterial.hpp`, `include/viennaps/fields/models/FermiDiffusion.hpp`

**Produces:** `DiffusivityMaterial<NumericType>` (MOOSE `MatDiffusionBase` pattern) + `FermiDiffusion<NumericType>` extending `DiffusionModel`. DiffusivityMaterial computes D(C,T) at quadrature points via `GridFunctionCoefficient`. FermiDiffusion uses it for concentration-dependent D.

**MOOSE reference:** `MatDiffusionBase::precomputeQpResidual()` = `_diffusivity[_qp] * _grad_v[_qp]`. Material property evaluated at quadrature points. Jacobian: `precomputeQpJacobian()` adds `dD/dC * phi * grad_v`. **Critically**, MOOSE supplies `dD/dC` automatically via `DerivativeMaterialInterface<Kernel>` (verified `framework/include/materials/DerivativeMaterialInterface.h`); MFEM does not.

**Jacobian strategy (DECISION — affects every nonlinear model in Phases 2-9):** Per the Phase 1 Task 5 default, the engine uses Picard iteration for nonlinear models unless the model explicitly opts in to providing `dD/dC`. For FermiDiffusion at extrinsic doping the Picard lag is significant, so this task MUST supply a `dD/dC` path. Pick one and note it in the commit:
- **(b1) QuadratureFunction Jacobian** — model exposes `evalDdC(trans, ip)` returning `dD/dC` evaluated at quadrature points; engine assembles the chain-rule term via `mfem::MixedGradGradIntegrator` with a `QuadratureFunctionCoefficient` updated each Newton step.
- **(b2) Picard with frequent reassembly** — re-evaluate `FermiDCoef` from current C at each Newton iteration; converge slowly (factor 2-4 more iterations than (b1)) but trivially correct.

**Phase 2 default: (b1)** - extrinsic diffusion is the regime where Picard hurts most.

**Hypre-backed Newton (recommended over both (a) and (b)):** MFEM's `NewtonSolver` (`solvers.hpp:780`) with `HypreBoomerAMG` as the inner linear solver eliminates hand-derived Jacobians entirely. The engine assembles the nonlinear residual `F(u) = K(u)*u - R(u)`, and Newton's method computes the Jacobian via `MatFDColoring`-style finite differences (one residual evaluation per Jacobian-vector product). `HypreBoomerAMG` (constructed from the `HypreParMatrix` bridge per Phase 1 Task 5) serves as the AMG preconditioner inside Newton's Krylov iteration. This gives quadratic convergence (vs Picard's linear) without `dD/dC` derivation. Set `NewtonSolver::SetAdaptiveLinRtol()` (`solvers.hpp:856`) for Eisenstat-Walker adaptive tolerance.

- [ ] **Step 1: Write failing test** - `TestFermiDiffusion()`: create model with D_i=1e-13, alpha=1.0. At high dopant (1e20, extrinsic), assert `getDiffusivity(C=1e20, T=1273)` > `getDiffusivity(C=1e15, T=1273)` (extrinsic enhancement). **Add:** `evalDdC(C=1e20, T=1273)` returns positive value consistent with analytic `d/dC [D_i*(1+alpha*n/ni)]`.

- [ ] **Step 2: Run to verify failure** -> FAIL

- [ ] **Step 3: Implement DiffusivityMaterial** - wraps D(C,T) computation as MFEM coefficient AND exposes `dD/dC` for Jacobian path (b1):
```cpp
// DiffusivityMaterial: evaluates D and dD/dC at quadrature points from local C.
// Mirrors MOOSE MatDiffusionBase + DerivativeMaterialInterface<Kernel>:
//   precomputeQpResidual()  = _diffusivity[_qp] * _grad_v[_qp]
//   precomputeQpJacobian() += dD/dC[_qp] * phi * grad_v
class FermiDCoef : public mfem::Coefficient {
  const mfem::GridFunction* conc_;
  NumericType D_i_, alpha_, T_, ni_;
public:
  void SetConcentrationField(const mfem::GridFunction* c) { conc_ = c; }
  double Eval(mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip) override {
    double C = conc_->GetValue(T, ip);
    double n = std::max(C, (double)ni_);  // n-type approximation
    return D_i_ * (1.0 + alpha_ * n / ni_);
  }
};

// Companion coefficient for the Jacobian chain-rule term (dD/dC * phi * grad_v).
// Engine uses this with mfem::MixedGradGradIntegrator when assembling Jacobian.
class FermiDdCCoef : public mfem::Coefficient {
  const mfem::GridFunction* conc_;
  NumericType D_i_, alpha_, ni_;
public:
  double Eval(mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip) override {
    double C = conc_->GetValue(T, ip);
    double n = std::max(C, (double)ni_);
    // d/dC [D_i*(1 + alpha*n/ni)] = D_i*alpha/ni * dn/dC
    // dn/dC = 1 when C > ni (extrinsic), 0 otherwise.
    return (C > ni_) ? D_i_ * alpha_ / ni_ : 0.0;
  }
};
```

- [ ] **Step 4: Implement FermiDiffusion** - uses DiffusivityMaterial in `assembleStiffness`. **Also override `assembleStiffnessJacobian`** (new virtual on `DiffusionModel` from Phase 1) supplying the `FermiDdCCoef` term:
```cpp
void assembleStiffness(mfem::BilinearForm& K,
                       const mfem::GridFunction& speciesGF,
                       const std::map<std::string, mfem::GridFunction*>& allSpecies,
                       const mfem::GridFunction* temp) const override {
  FermiDCoef coef(D_i_, alpha_, ni_, T_);
  coef.SetConcentrationField(&speciesGF);
  K.AddDomainIntegrator(new mfem::DiffusionIntegrator(coef));
}

// New: Jacobian chain-rule term. Engine calls this when assembling Newton J
// instead of (or in addition to) the Picard-lagged K above.
void assembleStiffnessJacobian(mfem::MixedBilinearForm& dKdC,
                               const mfem::GridFunction& speciesGF) const override {
  FermiDdCCoef dcoef(D_i_, alpha_, ni_);
  dcoef.SetConcentrationField(&speciesGF);
  // dK/dC contributes: integral of (dD/dC * phi_j) * grad(phi_i) . grad(test)
  dKdC.AddDomainIntegrator(new mfem::MixedGradGradIntegrator(dcoef));
}
```

- [ ] **Step 5: Run to verify pass** -> PASS

- [ ] **Step 6: Commit** - `"feat: add DiffusivityMaterial + FermiDiffusion with analytic dD/dC Jacobian (MOOSE MatDiffusion + DerivativeMaterialInterface pattern)"`

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

### Task 5: Segregation Interface Condition (InterfaceReaction Pattern)

**Files:** Create `include/viennaps/fields/models/Segregation.hpp`

**Produces:** `Segregation<NumericType>` - **two-sided interface condition** between material 1 (e.g. Si) and material 2 (e.g. SiO2), using the MOOSE `InterfaceReaction` pattern. The rate is `kf*C_1 - kb*C_2`; at equilibrium `C_2/C_1 = kf/kb = m(T)` (the segregation coefficient). This is more physical than a static penalty `penalty*(C_2 - m*C_1)` because the dynamic-rate form handles transients correctly.

**MOOSE reference (PRIMARY):** `framework/include/interfacekernels/InterfaceReaction.{h,C}` — verified. The residual is two-sided with opposite sign on each side of the interface, and crucially the **full 4-block Jacobian** is assembled explicitly:
```cpp
// InterfaceReaction.C lines 30-72 (verified):
// Residual:
//   Element  side: r =  _test         * (kf*u - kb*v)
//   Neighbor side: r = -_test_neighbor * (kf*u - kb*v)   // sign flip
// Jacobian (4 blocks):
//   ElementElement:    _test         *  kf * _phi
//   NeighborNeighbor: -_test_neighbor * -kb * _phi_neighbor
//   NeighborElement:  -_test_neighbor *  kf * _phi
//   ElementNeighbor:   _test         * -kb * _phi_neighbor
```
Mapping: `u` = concentration on material 1 side, `v` = concentration on material 2 side, `kf` = `K_seg`, `kb` = `K_deseg`.

**MOOSE reference (SECONDARY, surface loss only):** `modules/scalar_transport/include/bcs/BinaryRecombinationBC.h` models `A + B -> C` *at a boundary*, not at an interior interface. It is the right pattern for **surface** dose-loss where one species is consumed by another at the gas surface (Phase 4), but it is the **wrong** reference for Si/SiO2 segregation because it is one-sided. Use `InterfaceReaction` for segregation.

**Why two-sided matters in MFEM:** A plain `mfem::BoundaryIntegrator` is one-sided — it only contributes the Element residual and the Element-Element / Element-Neighbor Jacobian blocks. The other two blocks (Neighbor residual, Neighbor-Neighbor / Neighbor-Element Jacobian) must be assembled by a second integrator on the neighboring submesh, OR by using MFEM's `InterfaceIntegrator` machinery (`mfem::InterfaceSubmesh`, `mfem::L2 restricted to interface`). Skipping the neighbor side breaks mass conservation across the interface — the most common segregation bug.

- [ ] **Step 1: Write failing test** - `TestSegregation()`: 2-material mesh (Si + SiO2). Initialize Boron=1e18 in Si, 0 in SiO2. Solve 0->1s with m=0.1, kf=1e-3, kb=1e-2. Assert (a) C_SiO2/C_Si at interface ~ 0.1 within 5%, (b) **total dose across both materials conserved** within 1% (`|dose_final - dose_initial|/dose_initial < 0.01`). The dose-conservation assertion is what catches the one-sided-integrator bug.

- [ ] **Step 2: Run to verify failure** -> FAIL

- [ ] **Step 3: Implement SegregationCondition** - two-sided MFEM interface integrator. Implement BOTH sides in one class; do not ship a one-sided version:

```cpp
// Segregation interface condition — mirrors MOOSE InterfaceReaction exactly.
// rate = kf * C_1 - kb * C_2;  equilibrium C_2/C_1 = kf/kb = m(T)
//   kf = m(T) * k0      (forward: mat1 -> mat2)
//   kb = k0             (backward: mat2 -> mat1)
//
// Two-sided: must be applied to BOTH submeshes meeting at the interface,
// OR use mfem::InterfaceSubmesh + mfem::InterfaceIntegrator to handle
// both sides in one assembly. The 4 Jacobian blocks below must all be
// present or dose is not conserved across the interface.
class SegregationCondition {
public:
  void setSegregationCoefficient(NumericType m, NumericType k0) {
    kf_ = m * k0;   // forward rate (mat1 -> mat2)
    kb_ = k0;       // backward rate (mat2 -> mat1)
  }

  // Element-side (material 1) residual contribution at quadrature point:
  //   R_elem += test * (kf * C1 - kb * C2)
  // Element-side Jacobian:
  //   dR_elem/dC1      += test * kf * phi           (ElementElement)
  //   dR_elem/dC2      += test * (-kb) * phi_nbr    (ElementNeighbor)
  void assembleElementSide(mfem::LinearForm& R_elem,
                           mfem::DenseMatrix& K_elem_elem,
                           mfem::DenseMatrix& K_elem_nbr,
                           const mfem::GridFunction& C1,
                           const mfem::GridFunction& C2,
                           const mfem::FiniteElement& fe_elem,
                           const mfem::FiniteElement& fe_nbr,
                           const mfem::IntegrationRule& ir,
                           mfem::ElementTransformation& trans_elem,
                           mfem::ElementTransformation& trans_nbr);

  // Neighbor-side (material 2) residual contribution — SIGN FLIP:
  //   R_nbr += -test_nbr * (kf * C1 - kb * C2)
  // Neighbor-side Jacobian:
  //   dR_nbr/dC2       += -test_nbr * (-kb) * phi_nbr  (NeighborNeighbor)
  //   dR_nbr/dC1       += -test_nbr * kf * phi          (NeighborElement)
  void assembleNeighborSide(/* mirrors above with sign flip */);

private:
  NumericType kf_, kb_;
};
```

**Implementation options (pick one and note it in commit message):**
1. **Two integrators** — register `SegregationElementIntegrator` on the material-1 submesh boundary faces and `SegregationNeighborIntegrator` on the material-2 submesh boundary faces. Simpler, but requires keeping C1/C2 grid functions accessible across submeshes.
2. **mfem::InterfaceIntegrator** — uses `mfem::InterfaceSubmesh` to expose the interior interface as a first-class object with element + neighbor DoF on both sides. Cleaner long-term; this is what the Phase 4 multi-material coupling should standardize on.

Either way, the dose-conservation test in Step 1 will fail if any of the four Jacobian blocks is missing — that's the safety net.

- [ ] **Step 4: Run to verify pass** -> PASS (interface ratio within 5%, dose conserved within 1%)

- [ ] **Step 5: Commit** - `"feat: add Segregation two-sided interface condition (MOOSE InterfaceReaction pattern)"`

---

### Task 6: Integration Test - Fermi Diffusion with Segregation

**Files:** Modify `tests/diffusion/testDiffusionEngine.cpp`

- [ ] **Step 1: Write integration test** - `TestFermiWithSegregation()`: 2-material mesh, FermiDiffusion for Boron in Si, ConstantDiffusion in SiO2, Segregation BC at interface. Solve 0->30s at 1000C. Verify dose conservation across both materials and concentration ratio at interface matches m.

- [ ] **Step 2: Run** -> PASS

- [ ] **Step 3: Commit** - `"test: add Fermi diffusion with segregation integration test"`
