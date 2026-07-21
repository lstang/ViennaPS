# Diffusion Engine Phase 4: OED + TED + Dose Loss + Interface Physics

**Goal:** Add oxidation-enhanced diffusion (OED), transient-enhanced diffusion (TED) initial conditions, dose loss (evaporation), and advanced interface physics.

**Depends on:** Phase 3 (CDD, clustering)

## File Structure

| File | Responsibility |
|------|---------------|
| `models/OedSource.hpp` | I injection at Si/SiO2 interface during oxidation |
| `models/TedInitializer.hpp` | Initialize C_I from implant damage profile |
| `models/DoseLossBC.hpp` | Surface evaporation: -D*dC/dn = h*C |
| `models/ChargedEquilibriumDiffusion.hpp` | Equilibrium concentrations with charge |
| `models/CarbonDiffusion.hpp` | Carbon I-trapping, B suppression |
| `models/NitrogenDiffusion.hpp` | Nitrogen diffusion model |
| `models/CopperDiffusion.hpp` | Cu diffusion + ion-pairing |
| `models/MobileImpurity.hpp` | Mobile impurities / ion-pairing general |

## Tasks

### Task 0: Moving-Interface Idiom for the Si/SiO2 Oxidation Boundary

**Required before Task 1.** The Si/SiO2 interface moves as oxidation proceeds. Phase 4 must pick one of three MOOSE-verified idioms and name the class it mirrors. The current plan leaves this unspecified.

- **(A) Subdomain relabeling — RECOMMENDED for sharp Si↔SiO2 transitions.** `framework/include/meshmodifiers/ElementSubdomainModifierBase.h` (verified lines 56-75) — elements flip from Si to SiO2 subdomain as oxide grows; the moving boundary is auto-tracked via `_moving_boundaries` (`map<SubdomainPair, BoundaryID>`); fields on newly-oxidized elements are reinitialized via `ReinitStrategy { IC, POLYNOMIAL_NEIGHBOROR, POLYNOMIAL_WHOLE, POLYNOMIAL_NEARBY, NONE }`. **`POLYNOMIAL_NEARBY` (KD-tree extrapolation from active neighbors) is exactly what reinitializing dopant concentration on a newly-oxidized element needs.** Concrete subclass to mirror: `ThresholdElementSubdomainModifier.h`.
- **(B) MFEM-side ALE for small smooth oxide growth.** `framework/include/mfem/mesh/MFEMMesh.h::displace()` + `framework/include/mfem/problem/MFEMProblem.h` methods `displaceMesh()` / `updateFESpaces()` / `updateGridFunctions()`. Use when oxide growth is small enough that no element inversion occurs and no topology change is needed.
- **(C) Boundary-node ALE.** `modules/tensor_mechanics/include/kernels/ALEKernel.h` (mechanics residual on undisplaced mesh) + `modules/navier_stokes/include/bcs/INSADDisplaceBoundaryBC.h` (velocity·dt boundary-node displacement) — the template for an oxidation-rate-driven boundary BC.

**Phase 4 default: (A) subdomain relabeling** because oxidation is a sharp Si→SiO2 transition with topology change. Use `INSADDisplaceBoundaryBC` pattern (C) as a secondary for the oxide free surface if a smooth-surface ALE is wanted on top.

- [ ] **Step 1: Document the choice** in a short ADR (or this plan file): which idiom, which MOOSE class, and why.
- [ ] **Step 2: No standalone code in Task 0** — the choice propagates into Tasks 1 and 9 (Integration Test).

---

### Task 1: OedSource Model

- `OedSource<NumericType>` - boundary flux BC on C_I at Si/SiO2 interface. Flux = k_ox * dx_ox/dt (interstitial injection proportional to oxidation rate). Requires coupling to `psOxidation` for oxide growth rate.

**Moving-interface coupling (per Task 0):** since the interface is tracked by `ElementSubdomainModifier` (idiom A) or MFEM ALE (idiom B), `OedSource` is *not* a fixed-boundary BC. It must:
- Query the current interface location from the modifier's `_moving_boundaries` (idiom A) or from the displaced MFEM mesh (idiom B) each time step.
- Use SUPG stabilization on the resulting advection-like I-injection (per Phase 3 Task 4 `LevelSetAdvectionSUPG` pattern) when the interface velocity is high.

- Test: simulate oxidation step, verify C_I increases near interface AND the injection tracks the moving interface (not a fixed boundary).
- Commit: `"feat: add OED interstitial injection source with moving-interface coupling"`

### Task 2: TedInitializer

- `TedInitializer<NumericType>` - reads MCBcaImplant damage output (I/V profiles) and sets initial conditions for CDD diffusion. `initialize(field, implantResult)`.

**MOOSE references (verified):**
- `framework/include/ics/SolutionIC.h` — the canonical MOOSE pattern for "initialize a field from another solver's output." `TedInitializer` is conceptually a `SolutionIC` whose source is the implant Monte Carlo result.
- `framework/include/ics/IntegralPreservingFunctionIC.h` — **critical for dose conservation.** When projecting a damage profile onto the FEM mesh, the *total interstitial count* must be preserved, not the pointwise values. `IntegralPreservingFunctionIC` rescales the projected field so `∫C_projected = ∫C_source`. Without this, mesh-resolution dependence in the projection silently changes the I dose.
- `framework/include/ics/FunctionIC.h` / `FunctorIC.h` — for analytic damage profiles (e.g., PearsonIV from implant moments).

**Action:** make `TedInitializer` accept either a `MCBcaImplantResult` (discrete damage points → `SolutionIC` path) or an analytic profile function (`FunctionIC` path). In both cases, after projection, apply the `IntegralPreservingFunctionIC` rescaling so total interstitial count is invariant under mesh refinement.

- Test: provide mock implant damage with known total I count. After projection onto coarse and fine FEM meshes, verify total I count matches to within 0.1% on BOTH meshes (the conservation invariant).
- Commit: `"feat: add TED initializer from implant damage with IntegralPreservingFunctionIC dose conservation"`

### Task 3: DoseLossBC Model

- `DoseLossBC<NumericType>` - surface boundary condition: `-D*dC/dn = h*C`. Implemented as MFEM boundary integrator with Robin BC.

**MOOSE references (verified):**
- `framework/include/bcs/ADRobinBC.h` — Robin BC `du/dn = coef * u`. **Note:** MOOSE's Robin BC does NOT include the diffusion coefficient D in `coef`; ViennaPS's `DoseLossBC` form `−D·dC/dn = h·C` must explicitly fold D into `coef = h/D`. Document this in the header.
- `modules/scalar_transport/include/bcs/DissociationFluxBC.{h,C}` — complement to `BinaryRecombinationBC`: residual `-_test * Kd * v`. Use this when one species is lost at the surface by transforming into another (e.g., B_evap → B_gas). Provide a `DissociationLossBC` variant alongside `DoseLossBC`.
- `framework/include/bcs/VacuumBC.h` — Marshak-style `du/dn = coef*u` with `coef = 1/(4·v_th)`. More physical than generic Robin for evaporation to vacuum.

- Test: high h → dose decreases over time; h=0 → dose conserved. **Add:** the time-integrated dose loss equals `∫h·C·dt` to within 0.1% (mass-balance check on the surface flux).
- Commit: `"feat: add DoseLossBC (Robin) + DissociationLossBC + VacuumBC variants (MOOSE ADRobinBC + DissociationFluxBC + VacuumBC pattern)"`

### Task 4: ChargedEquilibriumDiffusion Model
- `ChargedEquilibriumDiffusion<NumericType>` - computes equilibrium concentrations with charge states. D = D0 * f_eq(T, n, p)
- Test: verify D_eff matches equilibrium at given T
- Commit: `"feat: add ChargedEquilibrium diffusion model"`

### Task 5: CarbonDiffusion Model
- `CarbonDiffusion<NumericType>` - C traps interstitials: C + I -> C-I complex. Suppresses B transient diffusion. Adds C_I_trap species
- Test: inject C + excess I -> I trapped, B TED suppressed
- Commit: `"feat: add Carbon diffusion with I trapping"`

### Task 6: NitrogenDiffusion Model
- `NitrogenDiffusion<NumericType>` - N-specific transport, affects oxidation rate
- Test: verify N diffuses, basic transport
- Commit: `"feat: add Nitrogen diffusion model"`

### Task 7: CopperDiffusion Model
- `CopperDiffusion<NumericType>` - Cu diffusion + ion-pairing with charged dopants. D_Cu affected by local electric field
- Test: verify Cu diffusion faster in doped regions (ion-pairing)
- Commit: `"feat: add Copper diffusion with ion-pairing"`

### Task 8: MobileImpurity Model
- `MobileImpurity<NumericType>` - general mobile impurity framework with ion-pairing. Substrate for Cu, Na, etc.
- Test: verify basic mobile impurity transport
- Commit: `"feat: add MobileImpurity framework with ion-pairing"`

### Task 9: Integration Test - Full OED Sequence
- Implant B -> oxidize (with OED) -> anneal -> verify enhanced diffusion near Si/SiO2 interface vs bulk
- Commit: `"test: add OED integration test"`
