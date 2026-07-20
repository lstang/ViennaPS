# Diffusion Engine Phase 6: SiGe/SiGeC + III-V Compound Diffusion

**Goal:** Add SiGe interdiffusion, boron diffusion in SiGe (bandgap effects), carbon suppression of B TED, Ge-B pairing, and III-V compound semiconductor diffusion.

**Depends on:** Phase 3 (CDD)

## File Structure

| File | Responsibility |
|------|---------------|
| `models/SiGeDiffusion.hpp` | SiGe interdiffusion + B in SiGe |
| `models/SiGeCDiffusion.hpp` | Carbon I-trapping, B suppression |
| `BandgapModel.hpp` | E_g(x_Ge, strain) -> n_i |
| `models/IIIVDiffusion.hpp` | III-V compound (GaAs, InP) diffusion |
| `MaterialConverter.hpp` | Convert mesh material (Si -> GaAs) |

## Tasks

### Task 1: BandgapModel
- `BandgapModel<NumericType>` - E_g(x_Ge, T, strain). Affects n_i which affects Fermi-dependent D. Bowing parameter for SiGe
- Test: E_g(Si, 300K) ~ 1.12 eV, E_g(Ge, 300K) ~ 0.66 eV, E_g(Si0.5Ge0.5) between
- Commit: `"feat: add BandgapModel for SiGe"`

### Task 2: SiGe Interdiffusion
- `SiGeDiffusion<NumericType>` - D_SiGe(x_Ge, T) with concentration-dependent terms. Low-to-high Ge and low-Ge models. Si and Ge interdiffuse
- Test: Si/Ge interface -> intermixing over time. Verify profile smooths
- Commit: `"feat: add SiGe interdiffusion model"`

### Task 3: Boron Diffusion in SiGe
- B diffusivity modified by Ge content via bandgap and point-defect parameters. D_B(SiGe) != D_B(Si)
- Test: verify D_B in SiGe differs from D_B in Si at same T
- Commit: `"feat: add boron diffusion in SiGe with bandgap effects"`

### Task 4: SiGeC - Carbon Suppression
- `SiGeCDiffusion<NumericType>` - C traps interstitials, suppresses B transient diffusion. C + I -> C-I complex. Extends CarbonDiffusion from Phase 4 with SiGe context
- Test: with C present, B TED suppressed vs without C
- Commit: `"feat: add SiGeC carbon suppression of B TED"`

### Task 5: Ge-B Pairing
- Germanium-boron pairing model. Ge + B -> Ge-B pair (immobile). Reduces active B. Cluster initialization from Ge-B pairs
- Test: high Ge + B -> pairing reduces mobile B
- Commit: `"feat: add Ge-B pairing model"`

### Task 6: Strain Effects
- SiGe strain affects dopant activation and point-defect equilibrium. Strain from lattice mismatch. Modifies D via strain-dependent activation energies
- Test: strained SiGe has different D than relaxed
- Commit: `"feat: add strain effects on diffusion in SiGe"`

### Task 7: MaterialConverter
- `MaterialConverter` - converts mesh elements from Si to GaAs/InP/etc. Changes material properties, attribute names
- Test: convert Si mesh to GaAs, verify attribute change
- Commit: `"feat: add MaterialConverter for III-V"`

### Task 8: III-V Diffusion
- `IIIVDiffusion<NumericType>` - species-specific D in GaAs/InP. Different I/V equilibrium (Ga/As sublattice vacancies)
- Test: verify basic diffusion in GaAs
- Commit: `"feat: add III-V compound semiconductor diffusion"`

### Task 9: Integration Test
- Deposit SiGe layer -> implant B -> anneal -> verify B profile in SiGe differs from Si. With C co-implant, verify TED suppression
- Commit: `"test: add SiGe/SiGeC diffusion integration test"`
