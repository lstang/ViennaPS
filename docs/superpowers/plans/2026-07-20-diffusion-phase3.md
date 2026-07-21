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

- [ ] **Step 1: Write failing test** - `TestKernelTerm()`: create `DiffusionTerm` with constant D, `ReactionTerm` with rate k, `CoupledForceTerm` coupling species A to B. Verify each produces correct residual contribution.

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
- `PairDiffusion<NumericType>` - dopant+I pair. D_eff = D_pair * (C_I / C_I_eq). Consumes I, enhances dopant mobility (TED)
- Test: inject excess I, verify dopant diffusion enhanced vs constant-D baseline
- Commit: `"feat: add PairDiffusion model for TED"`

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

- [ ] **Step 1: Write failing test** - implant B -> anneal -> verify TED (transient enhancement), 311 formation, dose retention

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

### Task 11: NeutralReactDiffusion Model
- Neutral defect reactions without charge coupling
- Test: verify basic recombination
- Commit: `"feat: add NeutralReactDiffusion model"`

### Task 12: Integration Test - CDD TED Sequence
- Full sequence: implant damage profile -> CDD anneal -> verify dopant profile matches expected TED behavior (transient enhancement then relaxation)
- Commit: `"test: add CDD TED integration test"`
