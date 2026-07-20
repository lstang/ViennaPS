# Diffusion Engine Phase 3: CDD + React/Pair/Charged Variants + Clustering

**Goal:** Add the CDD (Classical Dopant Diffusion) full-coupled model, React/ChargedReact/Pair/ChargedPair transport models, and point-defect clustering ({311}, VC, BIC, dislocation loops).

**Depends on:** Phase 2 (Fermi, IntrinsicCarrier)

## File Structure

| File | Responsibility |
|------|---------------|
| `models/ReactDiffusion.hpp` | I+V recombination: dC_I/dt -= k*C_I*C_V |
| `models/ChargedReactDiffusion.hpp` | React + charge-state-dependent rates |
| `models/PairDiffusion.hpp` | D = D_pair * (C_I/C_I_eq), TED-like |
| `models/ChargedPairDiffusion.hpp` | Pair + Fermi-level coupling |
| `models/CddDiffusion.hpp` | Full coupled: dopant+I+V with all terms |
| `models/NeutralReactDiffusion.hpp` | Neutral defect reactions |
| `models/Cluster311.hpp` | {311} interstitial cluster growth/dissociation |
| `models/VacancyCluster.hpp` | Vacancy cluster (VC) model |
| `models/ImpurityCluster.hpp` | Boron-interstitial clustering (BIC) |
| `models/DislocationLoop.hpp` | Loop growth from I supersaturation |
| `PointDefectEquilibrium.hpp` | C_I^eq, C_V^eq calculators |

## Tasks

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

### Task 10: CDD (Classical Dopant Diffusion) Model
- `CddDiffusion<NumericType>` - composes Pair + React + clustering into one model. Registers dopant, I, V, 311, VC, BIC, loop species. Full coupled system
- Test: implant B -> anneal -> verify TED (transient enhancement), 311 formation, dose retention
- Commit: `"feat: add CDD full-coupled diffusion model"`

### Task 11: NeutralReactDiffusion Model
- Neutral defect reactions without charge coupling
- Test: verify basic recombination
- Commit: `"feat: add NeutralReactDiffusion model"`

### Task 12: Integration Test - CDD TED Sequence
- Full sequence: implant damage profile -> CDD anneal -> verify dopant profile matches expected TED behavior (transient enhancement then relaxation)
- Commit: `"test: add CDD TED integration test"`
