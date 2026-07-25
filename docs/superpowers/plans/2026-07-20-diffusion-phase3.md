# Diffusion Engine Phase 3: CDD + React/Pair/Charged Variants + Clustering

**Goal:** Add the CDD (Classical Dopant Diffusion) full-coupled model, React/ChargedReact/Pair/ChargedPair transport models, and point-defect clustering ({311}, VC, BIC, dislocation loops). Refactors `DiffusionModel` to composable `KernelTerm` pattern (MOOSE kernel composition: each physics term is a separate tiny kernel, independently testable).

**Depends on:** Phase 2 (Fermi, IntrinsicCarrier)

## File Structure

| File | Responsibility |
|------|---------------|
| `fields/KernelTerm.hpp` | Composable kernel base (MOOSE Kernel pattern: one term per kernel) |
| `fields/KernelTerms.hpp` | Built-in terms: DiffusionTerm, ReactionTerm, CoupledForceTerm, SourceTerm |
| `models/ReactDiffusion.hpp` | I+V recombination: dC_I/dt -= k*C_I*C_V |
| `models/ChargedReactDiffusion.hpp` | React + charge-state-dependent rates |
| `models/PairDiffusion.hpp` | D = D_pair * (C_I/C_I_eq), TED-like |
| `models/ChargedPairDiffusion.hpp` | Pair + Fermi-level coupling |
| `models/CddDiffusion.hpp` | Full coupled: composes DiffusionTerm + CoupledForceTerm + ReactionTerm |
| `models/NeutralReactDiffusion.hpp` | Neutral defect reactions |
| `models/Cluster311.hpp` | {311} interstitial cluster growth/dissociation |
| `models/VacancyCluster.hpp` | Vacancy cluster (VC) model |
| `models/ImpurityCluster.hpp` | Boron-interstitial clustering (BIC) |
| `models/DislocationLoop.hpp` | Loop growth from I supersaturation |
| `PointDefectEquilibrium.hpp` | C_I^eq, C_V^eq calculators |

## Tasks

### Task 0: KernelTerm Composable Base (MOOSE Kernel Pattern)

**Files:** Create `include/viennaps/fields/KernelTerm.hpp`, `include/viennaps/fields/KernelTerms.hpp`

**Produces:** `KernelTerm` abstract base + built-in terms. Each term handles ONE physics contribution (like MOOSE kernels: `MatDiffusion`, `Reaction`, `CoupledForce`, `BodyForce`). `DiffusionModel` becomes a container of `KernelTerm` objects.

**MOOSE reference:** `MatDiffusion` (30 lines), `Reaction` (50 lines), `CoupledForce` (68 lines) - each tiny, focused, independently testable.

**Primary/equilibrium species split (MOOSE `chemical_reactions` pattern — verified):** The canonical MOOSE multi-species coupling example lives at `modules/chemical_reactions/include/kernels/`. Kernels come in pairs:
- `PrimaryDiffusion` / `PrimaryConvection` / `PrimaryTimeDerivative` — operate on **primary** species (the N species being solved as unknowns).
- `CoupledDiffusionReactionSub` / `CoupledConvectionReactionSub` / `CoupledBEEquilibriumSub` — operate on **equilibrium (secondary)** species, derived as AuxVariables from primary species via stoichiometry + log_k. `CoupledDiffusionReactionSub` carries `_weight`, `_log_k`, `_sto_u`, `std::vector<Real> _sto_v` (stoichiometric coefs), plus `_gamma_u`, `_gamma_v[]`, `_gamma_eq` activity coefficients.

**Why this matters for the diffusion plan:** `ChargedFermiDiffusion` (Phase 2 Task 3) decomposes a dopant into charge states B⁰, B⁺, B⁻. Solving all three as primary species is redundant — they are in instantaneous equilibrium via the Fermi level. Instead, solve 1 primary (total B) and compute 3 equilibrium AuxVariables (the charge-state fractions) the way `CoupledBEEquilibriumSub` does. **Add to this task:** an `EquilibriumSpeciesAuxKernel` concept that the engine evaluates *after* each Newton step to update AuxVariables from primary species + equilibrium constants. This avoids redundant solves and is the only correct formulation when charge-state equilibrium is fast compared to diffusion.

- [ ] **Step 1: Write failing test** - `TestKernelTerm()`: create `DiffusionTerm` with constant D, `ReactionTerm` with rate k, `CoupledForceTerm` coupling species A to B. Verify each produces correct residual contribution. **Add:** `TestEquilibriumSpeciesAuxKernel()` — define a primary species "B_total" and an equilibrium species "B_active = K_eq(T) * B_total"; after evaluating the aux kernel, verify B_active = K_eq * B_total within tolerance.

- [ ] **Step 2: Run to verify failure** -> FAIL

- [ ] **Step 3: Implement KernelTerm base + built-in terms:**
```cpp
// Base: one physics term per KernelTerm (MOOSE Kernel pattern)
class KernelTerm {
public:
  virtual ~KernelTerm() = default;
#ifdef VIENNAPS_HAS_MFEM
  virtual void assembleResidual(mfem::LinearForm& R,
      const std::map<std::string, mfem::GridFunction*>& species,
      const mfem::GridFunction* temp) const {}
  virtual void assembleStiffness(mfem::BilinearForm& K,
      const std::map<std::string, mfem::GridFunction*>& species,
      const mfem::GridFunction* temp) const {}
  virtual void assembleMass(mfem::BilinearForm& M) const {}
#endif
  virtual std::string targetSpecies() const = 0;
};

// DiffusionTerm: ∇·(D∇C) — MOOSE MatDiffusion equivalent
class DiffusionTerm : public KernelTerm { ... };

// ReactionTerm: λC — MOOSE Reaction equivalent
class ReactionTerm : public KernelTerm { ... };

// CoupledForceTerm: -σ*v — MOOSE CoupledForce equivalent (coupling)
class CoupledForceTerm : public KernelTerm { ... };

// SourceTerm: f(x,t) — MOOSE BodyForce equivalent
class SourceTerm : public KernelTerm { ... };
```

- [ ] **Step 4: Run to verify pass** -> PASS

- [ ] **Step 5: Commit** - `"feat: add KernelTerm composable base with MOOSE kernel pattern"`

---

### Task 1: PointDefectEquilibrium Calculator
- Create `PointDefectEquilibrium<NumericType>` with `C_I_eq(T, material)`, `C_V_eq(T, material)` using Arrhenius from parameter DB
- Test: assert C_I_eq(1273, "Si") ~ 1e10-1e12 range
- Commit: `"feat: add PointDefectEquilibrium calculator"`

### Task 2: ReactDiffusion Model
- `ReactDiffusion<NumericType>` - 2 species (I, V). Stiffness: D_I*grad(C_I), D_V*grad(C_V). Reaction: R_I -= k*C_I*C_V, R_V -= k*C_I*C_V
- Test: initialize I=1e15, V=1e15. After 1s, both decrease. Verify mass action.
- Commit: `"feat: add ReactDiffusion I+V recombination model"`

### Task 3: ChargedReactDiffusion Model
- Extends React with charge-state-dependent k: `k = k0 * (1 + gamma * n/ni)`
- Test: at high dopant, recombination faster than intrinsic
- Commit: `"feat: add ChargedReactDiffusion model"`

### Task 4: PairDiffusion Model

- `PairDiffusion<NumericType>` - dopant+I pair. D_eff = D_pair * (C_I / C_I_eq). Consumes I, enhances dopant mobility (TED).

**Stabilization (REQUIRED — pair diffusion is convection-dominated):** `D_eff = D_pair * C_I/C_I_eq` is mathematically a *concentration-dependent diffusion* that, when linearized around a sharp C_I front (the typical TED initial condition), behaves like advection along the C_I gradient. Plain Galerkin will oscillate. MOOSE ships ready-made SUPG (Streamline Upwind Petrov-Galerkin) at `modules/level_set/include/kernels/LevelSetAdvectionSUPG.{h,C}` (verified):
- `LevelSetAdvectionSUPG` is `ADKernelGrad`; `precomputeQpResidual()` returns `tau * v * (v · ∇u)` with the classic Hughes-Brooks parameter **`tau = hmin / (2 * ||v||)`**.
- `LevelSetTimeDerivativeSUPG` provides the corresponding SUPG-stabilized time derivative `tau * v * u_dot`.
- `LevelSetForcingFunctionSUPG` covers the source term.

**Action:** port these to MFEM as `SupgAdvectionTerm` / `SupgTimeDerivativeTerm` `KernelTerm` subclasses (Phase 3 Task 0). The velocity field `v` is the effective drift = `D_pair * grad(C_I) / C_I_eq` — a coupled-species-derived coefficient, exactly the case `LevelSetAdvectionSUPG` handles (velocity is a `ADVectorVariableValue`). The same SUPG terms are needed by Phase 4 `OedSource` and Phase 9 `MeltDiffusion`.

- Test: inject excess I (step profile), verify dopant diffusion enhanced vs constant-D baseline AND verify no spurious oscillations in the dopant profile (SUPG on vs off comparison; SUPG-off should show >5% overshoot, SUPG-on <0.1%).
- Commit: `"feat: add PairDiffusion model for TED with SUPG stabilization (MOOSE LevelSetAdvectionSUPG pattern)"`

### Task 5: ChargedPairDiffusion Model
- Pair + Fermi coupling: D_pair depends on charge state
- Test: verify D_eff changes with dopant concentration
- Commit: `"feat: add ChargedPairDiffusion model"`

### Task 6: Cluster311 Model
- `Cluster311<NumericType>` - adds C_311 species. Reaction: dC_311/dt = k_f*C_I^n - k_r*C_311. Consumes I
- Test: inject excess I, verify 311 grows, I decreases. At long time, 311 dissociates
- Commit: `"feat: add {311} cluster model"`

### Task 7: VacancyCluster Model
- Analogous to Cluster311 for vacancies: dC_VC/dt = k_f*C_V^m - k_r*C_VC
- Test: inject excess V, verify VC grows
- Commit: `"feat: add vacancy cluster model"`

### Task 8: ImpurityCluster (BIC) Model
- `ImpurityCluster<NumericType>` - B + I -> BIC. dC_BIC/dt = k_f*C_B*C_I - k_r*C_BIC
- Test: high B + I -> BIC forms, active B decreases
- Commit: `"feat: add boron-interstitial cluster (BIC) model"`

### Task 9: DislocationLoop Model
- `DislocationLoop<NumericType>` - loop growth from I supersaturation. dC_loop/dt = k * (C_I/C_I_eq - 1)^p
- Test: sustained I supersaturation -> loop grows
- Commit: `"feat: add dislocation loop growth model"`

### Task 10: CDD (Classical Dopant Diffusion) Model - Composable

**Files:** Create `include/viennaps/fields/models/CddDiffusion.hpp`

**Produces:** `CddDiffusion<NumericType>` - composes `KernelTerm` objects into full coupled system. Each physics term is a separate tiny kernel (MOOSE pattern). This is the "kitchen sink" model.

**MOOSE alignment reference (verified):** `modules/chemical_reactions/include/physics/ReactionNetworkPhysicsBase.{h,C}` is the actual MOOSE physics base for coupled reaction networks. It holds `std::vector<VariableName> _solver_species` (primary) separate from `std::vector<AuxVariableName> _aux_species` (equilibrium), parses a `reactions` string into `ReactionNetworkUtils::Reaction` objects, and supports per-equation `equation_scaling`. Aligning the CDD API with this pattern (rather than hand-composing KernelTerms) makes the model declaratively reaction-string-driven — closer to how TCAD engineers specify these networks.

**Action:** make `CddDiffusion` accept a small reaction-specification format (list of `(reactants, products, rate_constant, equilibrium_constant)` tuples) and translate each into the appropriate `KernelTerm` composition internally. Keep the explicit `addTerm(...)` API from below for power users.

**Also:** per Phase 1 Task 3.5, `CddDiffusion` MUST call `physics.shouldCreateTimeDerivative(species, *this)` for each species it composes terms on, and skip the `TimeDerivativeTerm` when denied. This is the mechanism that prevents a double-`dC/dt` bug when CDD is composed with FermiDiffusion on the same species.

- [ ] **Step 1: Write failing test** - implant B -> anneal -> verify TED (transient enhancement), 311 formation, dose retention. **Add:** verify that composing CDD alongside FermiDiffusion on Boron does NOT produce a double dC/dt (use the Phase 1 Task 3.5 `shouldCreateTimeDerivative` gate).

- [ ] **Step 2: Run to verify failure** -> FAIL

- [ ] **Step 3: Implement CDD as composition of KernelTerms:**
```cpp
CddDiffusion() {
  // Dopant diffusion with pair enhancement
  addTerm(std::make_shared<PairDiffusionTerm>("Boron", D_pair, "Interstitial"));
  // Interstitial diffusion
  addTerm(std::make_shared<DiffusionTerm>("Interstitial", D_I));
  // Vacancy diffusion
  addTerm(std::make_shared<DiffusionTerm>("Vacancy", D_V));
  // I+V recombination (CoupledForceTerm pattern)
  addTerm(std::make_shared<CoupledForceTerm>("Interstitial", "Vacancy", -k_recomb));
  addTerm(std::make_shared<CoupledForceTerm>("Vacancy", "Interstitial", -k_recomb));
  // {311} cluster formation (ReactionTerm + CoupledForce)
  addTerm(std::make_shared<Cluster311FormationTerm>("311", "Interstitial"));
  // BIC formation
  addTerm(std::make_shared<BicFormationTerm>("BIC", "Boron", "Interstitial"));
  // Time derivative for each species
  addTerm(std::make_shared<TimeDerivativeTerm>("Boron"));
  addTerm(std::make_shared<TimeDerivativeTerm>("Interstitial"));
  addTerm(std::make_shared<TimeDerivativeTerm>("Vacancy"));
  addTerm(std::make_shared<TimeDerivativeTerm>("311"));
  addTerm(std::make_shared<TimeDerivativeTerm>("BIC"));
}
```

- [ ] **Step 4: Run to verify pass** -> PASS

- [ ] **Step 5: Commit** - `"feat: add CDD as composable KernelTerm composition"`

**Block preconditioning with Hypre (CDD multi-species system):** The CDD system has 5+ coupled species (Boron, Interstitial, Vacancy, 311, BIC). Two paths:
  - **Monolithic + `SetSystemsOptions`:** Pack all species into one `HypreParMatrix` and call `amg.SetSystemsOptions(numSpecies)` to configure BoomerAMG for block systems (nodal ordering). Simple but may struggle with strongly coupled off-diagonal blocks.
  - **Block-diagonal:** Use `mfem::BlockOperator` + separate `HypreBoomerAMG` per species block. Each species gets its own AMG cycle. Better for loosely coupled species (dopant diffusion is weakly coupled to cluster kinetics). This is the MFEM equivalent of MOOSE's `PhysicsBase::addPreconditioning()` with field-split.
  Prefer block-diagonal for CDD; switch to monolithic if convergence stalls.

**DAE support (critical for Phase 3):** Clustering with equilibrium species (charge-state fractions from `CoupledBEEquilibriumSub` pattern, Task 0) forms an Index-1 DAE. MFEM has **no IDA wrapper** (grep of `sundials.hpp` for `ida` is empty). Two options:
  - **(A) Raw SUNDIALS IDA C API** (recommended): use `IDACreate()` + `IDASInit()` with the same `HypreBoomerAMG` preconditioner wrapping as Phase 1 CVODE (`SUNLinSol_SPBCGS` + `IDASSetLinearSolver` + `IDASSetPreconditioner`). IDA handles the singular mass matrix natively.
  - **(B) DAE-to-ODE reformulation:** eliminate algebraic constraints (equilibrium species) analytically, solving only primary species as ODE. Simpler but fragile when equilibrium constants change.
  Prefer (A). The raw IDA API follows the same integration pattern as the CVODE + Hypre path from Phase 1.

### Task 11: NeutralReactDiffusion Model
- Neutral defect reactions without charge coupling
- Test: verify basic recombination
- Commit: `"feat: add NeutralReactDiffusion model"`

### Task 12: Integration Test - CDD TED Sequence
- Full sequence: implant damage profile -> CDD anneal -> verify dopant profile matches expected TED behavior (transient enhancement then relaxation)
- Commit: `"test: add CDD TED integration test"`
