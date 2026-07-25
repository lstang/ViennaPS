# Diffusion Engine Phase 9: Flash/Laser Anneal

**Goal:** Add flash and laser anneal simulation: heat transfer, melting phase field, crystallinity phase field, dopant diffusion in melt, and intensity models. Phase fields use Allen-Cahn equation (MOOSE `ADAllenCahn` + `ACInterface` pattern) for non-conserved order parameters.

**Depends on:** Phase 3 (CDD)

**Manual sources (verified via `iconv` UTF-16 decode of `Manual/sprocess_ug.md`):**
- SProcess §"Flash or Laser Anneal Model", body lines **20802–24080** (TOC p.201–241)
- Key sub-sections: §"Heat Transfer Equation" line 20926; §"Transfer Matrix Method" line 21237; §"FDTD" line 21972; §"Phase Field Method" line 20597; §"Melting Phase Field Equation" line 22575; §"Crystallinity Phase Field Equation" line 23069; §"Intensity Models for Flash Anneal" TOC p.217; §"Intensity Model for Scanning Laser" TOC p.219.

## Manual Equation References (SProcess §4.201–223)

**Heat transfer equation (SProcess eq. 213 at line 20928):**

```
ρ·c_P·∂T/∂t = ∇·(κ(T)·∇T) + G + ρ·L·∂φ/∂t                                 (SProcess 213)
```

- `κ` = thermal conductivity (phase-dependent: amorphous / liquid / crystalline, SProcess lines 21115–21124)
- `ρ`, `c_P` = mass density, specific heat (also phase-dependent)
- `G` = volumetric heat source from optical absorption (Beer–Lambert below, or TMM)
- `L` = latent heat (SProcess Eq. 244); `∂φ/∂t` couples to the melting phase field

**Beer–Lambert absorption (SProcess line 21231):**

```
G(z) = α(λ, T, material) · I₀(t) · exp(−∫₀^z α(z') dz')
```

`α` set via `pdbSet <material> Absorptivity {<expression>}` (SProcess line 21232). For multilayer / interference, replace with TMM (next).

**Transfer Matrix Method (SProcess §"TMM", line 21237):** for thin-film stacks where layer thickness ≈ wavelength, Beer–Lambert fails; use Fresnel + transfer matrices. Each layer `j` has complex wave impedance `Z_j = Z_0/n_j` (refractive index). Forward/backward amplitudes related by 2×2 transfer matrices at each interface. Result: per-segment heat generation `G(x, z, t)` accounting for thin-film interference. SProcess splits 2D structures into vertical segments (`Minimum.Angle.Between.Segments`, `Minimum.Segment.Width`, line 21280).

**FDTD (SProcess §"FDTD", line 21972):** for sub-wavelength features, solve Maxwell's equations on a tensor mesh (SProcess uses Sentaurus Mesh + EMW). Shares the complex refractive index with TMM. Local vs global temperature modes (line 22032).

**Melting phase field (SProcess §"Melting Phase Field Equation", line 22575):** Allen-Cahn for liquid/solid order parameter `φ_m ∈ [−1, +1]`:

```
∂φ_m/∂t = −M_φ · δF/δφ_m
F[φ_m, T] = ∫ [ (κ_φ/2)·|∇φ_m|² + (φ_m²−1)²/4 − λ·(T − T_m)·φ_m ] dV
δF/δφ_m = −κ_φ·∇²φ_m + (φ_m³ − φ_m) − λ·(T − T_m)
```

- `M_φ` = phase-field mobility (controls interface thickness vs width)
- `T_m` = melting point (1685 K for Si)
- Latent heat `L·∂φ_m/∂t` feeds back into the heat equation (SProcess 213 above)

**Crystallinity phase field (SProcess §"Crystallinity Phase Field Equation", line 23069):** Allen-Cahn for amorphous/crystalline `φ_c ∈ [−1, +1]`, driven by SPER velocity:

```
δF/δφ_c = −κ_φ·∇²φ_c + (φ_c³ − φ_c) − λ·v_SPER(T)·φ_c
v_SPER(T) = v_0 · exp(−Ea_SPER / kT)                                       (SPER velocity, SProcess §4.197)
```

Coupled to heat equation (line 23159): crystallization releases latent heat too.

## File Structure

| File | Responsibility |
|------|---------------|
| `models/AllenCahnTerm.hpp` | Allen-Cahn kernel: dη/dt = L*(κ∇²η - df/dη) (MOOSE ADAllenCahn pattern) |
| `models/ACInterfaceTerm.hpp` | Gradient energy term: κ∇η·∇test (MOOSE ACInterface pattern) |
| `models/HeatTransfer.hpp` | Thermal FEM: rho*cp*dT/dt = div(k*grad(T)) + Q (MOOSE ADHeatConduction pattern) |
| `models/LaserIntensity.hpp` | Gaussian, table lookup, scanning laser intensity |
| `models/TransferMatrix.hpp` | Optical absorption in multilayer (thin film optics) |
| `models/MeltingPhaseField.hpp` | Liquid/solid Allen-Cahn: f(η,T) = (η²-1)²/4 - λ(T-T_m)*η |
| `models/CrystallinityPhaseField.hpp` | Amorphous/crystalline Allen-Cahn: f(η,T) driven by SPER velocity |
| `models/MeltDiffusion.hpp` | Liquid-phase dopant diffusion (very fast D, activated where η_melt>0.5) |
| `models/FlashLaserAnneal.hpp` | Orchestrator: thermal + phase + diffusion |

## Tasks

### Task 1: HeatTransfer Model
- `HeatTransfer<NumericType>` - FEM thermal solver. rho*c_p*dT/dt = div(k*grad(T)) + Q(x,t). Temperature-dependent k, c_p. Uses same MFEM mesh as diffusion. Adds T as a species in the unknown vector
- Test: uniform Q -> T rises uniformly. Q=0, boundary T=1000 -> steady-state gradient
- Commit: `"feat: add HeatTransfer FEM thermal solver"`

### Task 2: TransferMatrix - Optical Absorption
- `TransferMatrix<NumericType>` - thin film optics for multilayer stack. Computes absorbed energy Q(x) from laser intensity, layer refractive indices, absorption coefficients. Fresnel equations at each interface
- Test: single Si layer, known absorption coefficient -> verify exponential decay of Q with depth
- Commit: `"feat: add TransferMatrix optical absorption model"`

### Task 3: LaserIntensity Models
- `LaserIntensity<NumericType>` - Gaussian (flash): I(r,t) = I0*exp(-r^2/2sigma^2)*pulse(t). Scanning laser: I(x,t) = I0*exp(-(x-v*t)^2/2sigma^2). Table lookup: user-specified I(x,t). User-specified: via PDE API callback
- Test: Gaussian peak at center, scanning laser moves with v*t
- Commit: `"feat: add LaserIntensity models (Gaussian, scanning, table)"`

### Task 3.5: AllenCahnTerm + ACInterfaceTerm (MOOSE Phase Field Pattern)
- Create `models/AllenCahnTerm.hpp` + `models/ACInterfaceTerm.hpp`
- **MOOSE reference:** `ADAllenCahn::computeDFDOP()` returns `_dFdEta[_qp]` (bulk driving force). `ACInterface::computeQpResidual()` = `_grad_u * kappaNablaLPsi()` (gradient energy).
- `AllenCahnTerm`: bulk driving force -L*df/dη. `ACInterfaceTerm`: gradient energy L*κ*∇η·∇test. Combined: dη/dt = L*(κ∇²η - df/dη)
- Test: double-well free energy f=(η²-1)²/4. Starting η=0.1 -> evolves to η=1. Starting η=-0.1 -> evolves to η=-1
- Commit: `"feat: add AllenCahnTerm + ACInterfaceTerm (MOOSE ADAllenCahn + ACInterface pattern)"`

**IMEX time integration with Hypre (Phase 9 solver strategy):** Allen-Cahn equations are ideal for IMEX (implicit diffusion, explicit reaction). The diffusion term is stiff and implicit; the reaction term is non-stiff and explicit.
  - **Recommended: `mfem::ARKStepSolver` (`sundials.hpp:720`) with `Type::IMEX` (`:728`).** Wraps SUNDIALS ARKode IMEX mode. The implicit solve uses the same `SUNLinSol_SPBCGS` + `HypreBoomerAMG` preconditioner path from Phase 1.
  - **Alternative: `mfem::IMEX_DIRK_RK3` (`ode.hpp:1079`).** Pure MFEM, no SUNDIALS dependency for the time integrator. Third-order accurate.
  - For the multi-physics coupled solve (Task 7 FlashLaserAnneal): pack T, eta_m, eta_c, and dopant into a block system. Use `mfem::BlockOperator` + block-diagonal `HypreBoomerAMG` (one AMG per field). This mirrors the Phase 3 CDD block preconditioning approach.

### Task 4: MeltingPhaseField (Allen-Cahn Equation)
- `MeltingPhaseField<NumericType>` - Allen-Cahn equation for liquid/solid order parameter η_m in [-1,1] (-1=solid, 1=liquid)
- Free energy: f(η,T) = (η²-1)²/4 - λ(T-T_melt)*η (double-well + thermal driving force)
- df/dη = η³ - η - λ(T-T_melt). When T>T_melt, driving force pushes η->1 (liquid)
- Latent heat release couples back to heat equation via L*∂η/∂t term in HeatTransfer
- Composed using AllenCahnTerm + ACInterfaceTerm + TimeDerivativeTerm
- Test: T > T_melt -> η_m -> 1 (liquid). T < T_melt -> η_m -> -1 (solid). Sharp interface at T_melt
- Commit: `"feat: add MeltingPhaseField using Allen-Cahn equation"`

### Task 5: CrystallinityPhaseField (Allen-Cahn Equation)
- `CrystallinityPhaseField<NumericType>` - Allen-Cahn for amorphous/crystalline η_c in [-1,1] (-1=amorphous, 1=crystalline)
- Free energy: f(η,T) = (η²-1)²/4 - λ*v_SPER(T)*η. Driving force from SPER velocity
- df/dη = η³ - η - λ*v_SPER(T). SPER velocity v = v0*exp(-Ea/kT)
- During melt (η_m>0): η_c resets to 0 (liquid has no crystallinity) via coupling term
- Coupled with existing SPERKernel for interface velocity
- Test: amorphous layer -> anneal -> η_c advances from crystalline seed. Melt -> η_c resets
- Commit: `"feat: add CrystallinityPhaseField using Allen-Cahn equation"`

### Task 6: MeltDiffusion
- `MeltDiffusion<NumericType>` - liquid-phase dopant diffusion. D_liquid >> D_solid (orders of magnitude). Activated where η_m > 0 (liquid phase field from Allen-Cahn). D_eff = D_solid + (D_liquid - D_solid) * (1+η_m)/2. Solute transport in liquid Si. Resolidification traps dopant at solidification front
- Test: melt region -> dopant diffuses much faster than solid. Resolidify -> trap
- Commit: `"feat: add MeltDiffusion for liquid-phase dopant transport"`

### Task 7: FlashLaserAnneal Orchestrator
- `FlashLaserAnneal<NumericType>` - ties together: laser intensity -> Q(x,t) -> heat transfer -> T(x,t) -> melting phase field -> crystallinity phase field -> dopant diffusion (solid + melt). Multi-physics coupled solve
- Test: flash pulse -> T rises -> melts -> dopant redistributes -> cools -> resolidifies
- Commit: `"feat: add FlashLaserAnneal multi-physics orchestrator"`

### Task 8: FDTD for Sub-Wavelength Features
- Optional: FDTD simulation of electromagnetic field for sub-wavelength structures. Computes Q(x,y,z) from Maxwell's equations. Used when feature size ~ wavelength
- Test: simple 1D layered -> verify matches transfer matrix
- Commit: `"feat: add FDTD for sub-wavelength laser absorption"`

### Task 9: Saving Thermal Profile

- Save the full thermal-and-phase state for use in subsequent diffusion steps; load for non-thermal-coupled diffusion.

**Scope (must be specified — the original task was too thin):**
1. **What is saved:** T(x), the phase fields η_m (melt) and η_c (crystallinity) from Tasks 4-5, AND any latent-heat accumulators / internal state the time integrator carries. Saving only T silently corrupts the restart if a subsequent step resumes from a partially-molten state.
2. **Format:** binary `dataStore` stream per MOOSE convention (verified `framework/include/restart/Backup.h` lines 25-33 — two stringstreams for header+data, plus `vector<pair<string,string>> mesh_files` sidecar). Use VTK/ParaView output only for inspection, not restart.
3. **Mesh-topology survival:** Phase 11 AMR changes mesh topology between save and load. The `mesh_files` sidecar exists precisely for this — save the mesh alongside the fields. A restart across an AMR cycle without mesh sidecar silently misinterprets DoF indices.

**MOOSE idiom (verified):** `framework/include/restart/Restartable.h::declareRestartableData<T>(name, args...)` (lines 126-127) — *individual fields* opt into restart by name; the engine's `save()`/`load()` walks the declared set. This avoids the "did I remember to serialize μ_accumulator?" bug class. Adopt this pattern: every model member that should survive a save/load declares itself with a name on construction.

**Action:** make every Phase 9 model (`MeltingPhaseField`, `CrystallinityPhaseField`, `MeltDiffusion`, latent-heat accumulator inside `HeatTransfer`) declare its persistent state via a `declareRestartable<T>("name")` API on the engine. `save()`/`load()` serialize the declared set + the mesh sidecar.

- Test: (a) save → load → verify all fields identical bit-for-bit; (b) save → AMR cycle (Phase 11) → load → verify fields identical after topology change (uses mesh sidecar); (c) save → engine destroyed → new engine → load → verify identical (no in-memory leakage).
- Commit: `"feat: add full thermal+phase state save/load via declareRestartableData pattern + mesh sidecar (MOOSE Restartable + Backup pattern)"`

### Task 10: Integration Test - Flash Anneal
- Implant B -> flash anneal (melt) -> verify: (1) T reaches T_melt, (2) dopant redistributes in melt, (3) resolidification traps dopant, (4) final profile is shallower/broader than solid-state anneal
- Commit: `"test: add flash anneal integration test"`

### Task 11: Integration Test - Laser Scan
- Scanning laser across surface -> verify T profile follows laser, melting zone moves with laser, dopant redistribution follows melt zone
- Commit: `"test: add scanning laser anneal integration test"`
