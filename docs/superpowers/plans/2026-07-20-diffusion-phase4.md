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

### Task 1: OedSource Model
- `OedSource<NumericType>` - boundary flux BC on C_I at Si/SiO2 interface. Flux = k_ox * dx_ox/dt (interstitial injection proportional to oxidation rate). Requires coupling to `psOxidation` for oxide growth rate
- Test: simulate oxidation step, verify C_I increases near interface
- Commit: `"feat: add OED interstitial injection source"`

### Task 2: TedInitializer
- `TedInitializer<NumericType>` - reads MCBcaImplant damage output (I/V profiles) and sets initial conditions for CDD diffusion. `initialize(field, implantResult)`
- Test: provide mock implant damage, verify C_I initialized correctly
- Commit: `"feat: add TED initializer from implant damage"`

### Task 3: DoseLossBC Model
- `DoseLossBC<NumericType>` - surface boundary condition: `-D*dC/dn = h*C`. Implemented as MFEM boundary integrator with Robin BC
- Test: high h -> dose decreases over time. h=0 -> dose conserved
- Commit: `"feat: add dose loss (evaporation) boundary condition"`

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
