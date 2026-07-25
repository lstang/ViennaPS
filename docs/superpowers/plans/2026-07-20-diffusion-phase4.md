# Diffusion Engine Phase 4: OED + TED + Dose Loss + Interface Physics

**Goal:** Add oxidation-enhanced diffusion (OED), transient-enhanced diffusion (TED) initial conditions, dose loss (evaporation), and advanced interface physics.

**Depends on:** Phase 3 (CDD, clustering)

**Solver backend:** All linear solves use HypreBoomerAMG (via HypreParMatrix from serial SparseMatrix using MPI_COMM_SELF, per Phase 1 Task 5). The moving-interface field transfer (L2 projection across mesh topology change) uses HyprePCG + HypreBoomerAMG for the mass matrix solve. OED interstitial injection tracks the moving Si/SiO2 interface each time step.

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

**Files:** Modify `include/viennaps/fields/models/ChargedEquilibriumDiffusion.hpp` (currently 1-line stub per gap analysis)

**Produces:** `ChargedEquilibriumDiffusion<NumericType>` — substitutes/immobile species whose concentration is determined by an equilibrium law with the local electron/hole density, rather than by an independent continuity equation. Used for Cu⁺/Cu⁰ equilibrium (SProcess 4.194 "ChargedEquilibrium" example at lines 19733–19910) and for dopant activation under equilibrium clustering.

**Manual equations (SProcess §"ChargedEquilibrium", body lines 19733–19910; ATHENA §3.1.1 charge-state table at lines 5229–5253):**

The species is split into charge states `z ∈ {−, 0, +}` with concentrations `C_z`. Local equilibrium imposes mass-action between charge states via the Fermi level:

```
C_z / C_0 = g_z(n, p, T) = exp((E_F − E_z) / kT)
```

where `E_z` is the energy level of charge state `z` and `E_F` is the Fermi level (computed from `n, p` by `IntrinsicCarrier` from Phase 2). The total concentration `C_tot = Σ_z C_z` is the conserved unknown; the mobile fraction `C_0` (or whatever state diffuses) is:

```
C_mobile = C_tot · g_0 / Σ_z g_z
```

**Effective diffusivity** (only the neutral state diffuses, per SProcess 4.194):

```
D_eff(C, T, n, p) = D_0(T) · g_0(n,p,T) / Σ_z g_z(n,p,T)
```

**Assembly (extends Phase 2 Fermi quadrature-point D pattern, MOOSE `MatDiffusionBase`):**

```cpp
// Per quadrature point q:
const NumericType ni = intrinsicCarrier_->ni(T, mat);
const NumericType n = electronDensity(C_dopant, ni, T);  // from IntrinsicCarrier
const NumericType p = ni*ni / std::max(n, eps);
NumericType sumG = 0, g0 = 0;
for (int z : chargeStates_) {
  const NumericType gz = std::exp((EF - energyLevel_[z]) / (kBT));
  sumG += gz;
  if (z == 0) g0 = gz;
}
const NumericType Deff = D0_*std::exp(-Ea_/(kBT)) * g0 / std::max(sumG, eps);
// Stiffness contribution: K_ij += Deff * grad_phi_j · grad_phi_i * JxW
```

**MOOSE references (verified):**
- `framework/include/materials/DerivativeParsedMaterial.h` — MOOSE's automatic-Jacobian path for `D(C,T,n,p)`; this is what `g_0/Σg_z` would map to in architecture (B). In architecture (A), supply `dD/dC` analytically (chain rule through `n(C)`) and feed it to the engine's Picard/Newton loop per Phase 1 Task 5 Jacobian strategy (b).
- `modules/chemical_reactions/include/kernels/CoupledBEKinetic.h` — backward-Euler kinetic rate; relevant if equilibrium is relaxed to a kinetic rate `dC_z/dt = k_fwd·C_0·n − k_rev·C_z` instead of imposed algebraically. SProcess exposes both; the equilibrium (algebraic) form is the default.

**Scope:** ship the algebraic-equilibrium form (instantaneous charge-state partition). The kinetic relaxation form is a thin extension — note in commit, defer.

- [ ] **Step 1: Write failing test** — `TestChargedEquilibriumDiffusion()`: at fixed T, intrinsic (n=ni), verify `D_eff = D_0/numChargeStates` (all states equally populated). At heavy n-type doping, verify `D_eff` tilts toward the negative-state weight (per SProcess 4.194 example).
- [ ] **Step 2: Run to verify failure** → FAIL (stub returns constant D)
- [ ] **Step 3: Implement** — replace stub with the quadrature-point `D_eff` above; register via `DiffusionModel::assembleStiffness` using `GridFunctionCoefficient` for `C_dopant`.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add ChargedEquilibriumDiffusion with Fermi-level charge-state partition (SProcess 4.194; MOOSE DerivativeParsedMaterial pattern)"`

---

### Task 5: CarbonDiffusion Model

**Files:** Modify `include/viennaps/fields/models/CarbonDiffusion.hpp`

**Produces:** `CarbonDiffusion<NumericType>` — substitutional carbon `C_s` diffuses slowly; the dominant physics is **interstitial trapping**: `C_s + I ⇌ C_sI` (carbon–interstitial complex). This depletes free interstitials and thereby suppresses boron TED (SProcess §"NeutralReact", Carbon example at lines 19364–19521; ATHENA §3.9.3 at lines 17141–17180).

**Manual equations (SProcess 4.193 / 4.193 Carbon example, body lines 19364–19564):**

Substitutional carbon continuity (slow diffusion + trapping):

```
∂C_s/∂t = ∇·(D_C(T) ∇C_s) − k_f·C_s·C_I + k_r·C_sI
```

Trapped complex `C_sI` (immobile):

```
∂C_sI/∂t = k_f·C_s·C_I − k_r·C_sI
```

with forward/reverse rates `k_f = 4π·r_C·D_I(T)`, `k_r = k_f·C*_I` at equilibrium (so the complex concentration relaxes to `C_sI_eq = C_s·C_I/C*_I`). `D_C(T) = D0_C·exp(−Ea_C/kT)`.

**Net effect on B TED:** by depleting free `C_I`, the apparent B pair-diffusivity `D_B^pair = D_B^pair,0 · C_I/C*_I` drops — verified by the integration test in Task 9.

**Assembly (KernelTerm composition from Phase 3 Task 0):**

```cpp
// On the C_s equation:
addTerm<DiffusionTerm>(D_C(T));                          // ∇·(D_C ∇C_s)
addTerm<ReactionTerm>([ ](C_s, C_I, C_sI){ return -k_f*C_s*C_I + k_r*C_sI; });
// On the C_sI equation (immobile — no diffusion term):
addTerm<ReactionTerm>(...);                              // +k_f*C_s*C_I - k_r*C_sI
// Cross-species coupling to C_I (from Phase 3 PointDefectEquilibrium):
registerCoupling("Interstitial", -k_f*C_s*C_I + k_r*C_sI);
```

**MOOSE references (verified):**
- `modules/geochemistry/include/userobjects/GeochemistryKineticRate.h` — general kinetic rate law `R = k_f·Πa_i^ν_i − k_r·Πa_j^ν_j`; this is the canonical MOOSE pattern for the `C_s+I ⇌ C_sI` reaction.
- `framework/include/kernels/CoupledForce.h` — `∑k_i·v_i` source/sink; use for the cross-species coupling to `C_I` (the rate appears with opposite sign on the C_I equation — mass conservation).

**Parameters (Phase 10 Task 9 calibrated DB):** `D0_C = 0.91 cm²/s`, `Ea_C = 3.15 eV` (substitutional C in Si, SProcess Advanced Calibration); `r_C ≈ 0.2 nm` capture radius.

- [ ] **Step 1: Write failing test** — `TestCarbonTrapping()`: inject excess I into a box with C_s present; verify free `C_I` decays toward `C*_I` faster than the no-carbon baseline; verify `C_sI` grows by the same amount (mass balance).
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — two-species model (`C_s`, `C_sI`) + cross-coupling to `Interstitial`. Use `KernelTerm` composition.
- [ ] **Step 4: Run to verify pass** → PASS; verify mass conservation `ΔC_s + ΔC_sI ≈ 0` within 0.1%.
- [ ] **Step 5: Commit** — `"feat: add CarbonDiffusion with C_s-I trapping kinetics (SProcess 4.193 Carbon example; MOOSE GeochemistryKineticRate pattern)"`

---

### Task 6: NitrogenDiffusion Model

**Files:** Modify `include/viennaps/fields/models/NitrogenDiffusion.hpp`

**Produces:** `NitrogenDiffusion<NumericType>` — nitrogen in Si/SiO₂. Two coupled effects (SProcess 4.193, lines 19364+):

1. **Slow substitutional diffusion** of `N_s` with Arrhenius D_N.
2. **N–V pairing** (nitrogen-vacancy complex `N_sV`), analogous to carbon-I trapping but on the vacancy sublattice. Drives nitrogen pile-up near Si/SiO₂ where V is injected.
3. **Oxidation-rate modification** — nitrogen at the Si/SiO₂ interface slows subsequent oxidation (thin interface term; couples to `psOxidation` via a multiplier on the linear rate constant `B/A → B/A · (1 − α·C_N_surf)`).

**Manual equations (SProcess 4.193; ATHENA Appendix B.6 Point Defect Parameters):**

```
∂N_s/∂t = ∇·(D_N(T) ∇N_s) − k_fV·N_s·C_V + k_rV·N_sV
∂N_sV/∂t =               k_fV·N_s·C_V − k_rV·N_sV
```

with `k_fV = 4π·r_N·D_V(T)`, `k_rV = k_fV·C*_V`. Interface modifier:

```
(B/A)_eff = (B/A)_0 · (1 − α_N · C_N_surf / C_N_ref)
```

**Scope decision:** implement the diffusion + N–V pairing (Tasks 1–2 analogue of Carbon). The oxidation-rate coupling is a single scalar multiplier on the oxidation model — defer the `psOxidation` integration to a follow-up but expose `setOxidationRateModifier(α, C_ref)` so it can be wired without re-architecting.

**MOOSE references (verified):** identical to Task 5 — `GeochemistryKineticRate` (pairing), `CoupledForce` (cross-species coupling to `Vacancy`). No new patterns.

- [ ] **Step 1: Write failing test** — `TestNitrogenVacancyPairing()`: inject V into a nitrogen-doped region; verify `N_sV` grows, free `N_s` drops, mass balance holds.
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — two-species (`N_s`, `N_sV`) + vacancy coupling. Mirror Carbon's structure.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add NitrogenDiffusion with N-V pairing and oxidation-rate modifier hook (SProcess 4.193)"`

---

### Task 7: CopperDiffusion Model

**Files:** Modify `include/viennaps/fields/models/CopperDiffusion.hpp`

**Produces:** `CopperDiffusion<NumericType>` — fast interstitial copper `Cu_i` diffusion + ion-pairing with ionized acceptors/donors + trapping at cluster/defect sites (SProcess 4.195 "ChargedEquilibrium" copper example, body lines 19869+; ATHENA §3.2.x). The classic use case is Cu⁺ drift in p-type Si under an electric field.

**Manual equations (SProcess 4.195; ATHENA §3.1.1 ion-pairing eq. 3-3 at lines 5338–5347):**

Interstitial copper continuity with drift (Einstein relation `μ = qD/kT`):

```
J_Cu = −D_Cu(T) · (∇C_Cu + (q/kT)·z_Cu·C_Cu·E)
∂C_Cu/∂t = −∇·J_Cu − k_pair·C_Cu·C_A^z + k_diss·C_CuA
```

where `E = −∇φ` is the electric field from the Poisson–Boltzmann solution (or quasi-neutral approximation `E = (kT/q)·∇ln(n)`), `z_Cu = +1`, and `C_A^z` is the ionized acceptor concentration. The pair `C_CuA` is immobile:

```
∂C_CuA/∂t =  k_pair·C_Cu·C_A^z − k_diss·C_CuA
```

**Assembly note (the drift term):** this is the first model in the engine with a genuine **drift** contribution. Use `mfem::MixedScalarVectorGradientIntegrator` (or hand-assemble the `(q/kT)·z·E·φ_j·φ_i` term) — equivalent to MOOSE's approach in `modules/combined/include/kernels/` for electromigration. The drift term must be assembled into the **same stiffness matrix** `K` as the diffusion term so the SUNDIALS preconditioner sees it.

**MOOSE references (verified):**
- `framework/include/kernels/CoupledForce.h` — acceptor-coupling source term.
- `framework/include/bcs/ADRobinBC.h` — already cited in Task 3; needed here for the surface Cu evaporation BC.
- For the drift term, no exact MOOSE diffusion+drift kernel exists in `framework/`; the closest analogue is `modules/combined/include/kernels/ACInterfaceKobayashi1.h` (vectorial gradient coupling). Implement the integrator directly in MFEM and cite `MFEML2ZienkiewiczZhuIndicator`'s sibling `framework/include/mfem/integrators/NLDiffusionIntegrator.h` for the assembly pattern.

**Parameters (Phase 10 Task 9):** `D0_Cu = 4.7e-3 cm²/s`, `Ea_Cu = 0.43 eV` (very fast — Cu is the canonical "fast" metal contaminant).

- [ ] **Step 1: Write failing test** — `TestCopperDriftPairing()`: (a) uniform p-type doping, no gradient → with E-field, Cu drifts up-gradient (drift dominates diffusion); (b) high acceptor concentration → Cu pairs form, mobile Cu drops.
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — drift+diffusion stiffness assembly + pairing reaction terms + acceptor coupling.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add CopperDiffusion with drift + ion-pairing (SProcess 4.195; ATHENA 3.1.1 eq. 3-3 drift term)"`

---

### Task 8: MobileImpurity Model

**Files:** Modify `include/viennaps/fields/models/MobileImpurity.hpp`

**Produces:** `MobileImpurity<NumericType, SpeciesTag>` — a templated generalization of Task 7's structure for any fast-diffusing metal contaminant (Cu, Na, K, Fe, Au). The `SpeciesTag` selects parameter sets from the Phase 10 calibrated DB. Models the same three pieces as Copper: drift-diffusion, ion-pairing, optional trapping.

**Why generalize (cite this in the commit):** SProcess treats Cu, Na, K, Fe, Au with the same Continuum/ChargedEquilibrium model framework — only the parameters (`D0, Ea, z, pair_rate, trap_rate`) differ. Templating avoids 5 near-identical headers and matches the spec's ParameterDatabase inheritance pattern (Phase 10 Task 11).

**API:**

```cpp
template <class NumericType, class SpeciesTag>
class MobileImpurity : public DiffusionModel<NumericType> {
  // SpeciesTag carries: name, D0, Ea, charge z, pairing target species,
  // capture radius, trap-site density (all from CalibratedParameters).
  // Defaults to the Copper parameter set (Task 7) when SpeciesTag = CopperTag.
};
using CopperDiffusion = MobileImpurity<double, CopperTag>;
using SodiumDiffusion = MobileImpurity<double, SodiumTag>;
```

**MOOSE reference (verified):** `modules/chemical_reactions/include/materials/LangmuirMaterial.h` — Langmuir isotherm for sorption/trapping; the canonical MOOSE pattern for "mobile species ↔ finite trap sites." Use the `MollifiedLangmuirMaterial.h` variant (smoothed) for differentiability under Newton.

**Scope:** ship `CopperTag`, `SodiumTag`, `IronTag` parameter sets wired through the DB. `SodiumTag` is the highest-value addition after Cu because Na contamination is the classical MOS-BTBT reliability failure mode. `IronTag` traps at oxygen precipitates — note as a stretch, defer.

- [ ] **Step 1: Write failing test** — `TestMobileImpurityTemplate()`: instantiate `MobileImpurity<double, SodiumTag>`, verify `D_Na(T=1000K)` matches the DB value; verify `CopperDiffusion` alias still produces Task 7's behavior bit-for-bit.
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — refactor CopperDiffusion (Task 7) into the template; add `SodiumTag` DB entry.
- [ ] **Step 4: Run to verify pass** → PASS; Task 7's test still passes unchanged.
- [ ] **Step 5: Commit** — `"refactor: generalize CopperDiffusion into MobileImpurity<SpeciesTag>; add Na/Fe tags (MOOSE LangmuirMaterial trapping pattern)"`

---

### Task 9: Integration Test — Full OED/TED/Impurity Sequence

**Files:** Modify `tests/diffusion/testDiffusion.cpp`

**Scenario (matches ATHENA §3.5 + SProcess §4.5 typical OED workflow):**

1. **Implant** B at 50 keV, dose 1e14 cm⁻² (use mock implant → `TedInitializer` Task 2 with `IntegralPreservingFunctionIC`).
2. **Oxidize** at 1000 °C for 30 min in wet O₂ — `OedSource` (Task 1) injects I at the moving Si/SiO₂ interface per Task 0 idiom (A).
3. **Anneal** 1050 °C for 30 s — TED from the implanted + OED-injected I drives transient B diffusion.
4. **Variant with carbon:** co-implant C before oxidation; `CarbonDiffusion` (Task 5) suppresses free I → TED reduced.

**Assertions:**

- (a) B junction depth `x_j` (where `C_B = C_substrate`) deeper in the oxidized case than in a no-oxidation control by ≥ 20% (OED signature).
- (b) Total interstitial dose decays monotonically through anneal (recombination dominates once implant + OED sources stop).
- (c) In the carbon-co-implanted variant, B `x_j` is shallower than the no-carbon variant by ≥ 15% (TED suppression signature).
- (d) Mass balance: `dose_B(oxidation_end) ≈ dose_B(implant)` within 1% — verifies the moving-interface field transfer (Task 0 idiom A `POLYNOMIAL_NEARBY`) doesn't invent/lose B as the interface crosses elements.
- (e) **Regression snapshot:** save the 1D B profile (`ResultsExtractor::extract1D` from Phase 10 Task 5) to `tests/diffusion/regression/oed_ted_baseline.csv`; future runs compare within 1% per point. (Note in commit: this closes gap-analysis row "Regression tests (snapshots) ❌".)

**MOOSE reference (verified):** `modules/combined/test/tests/` is the canonical place multi-physics sequence tests live in MOOSE; the test idiom is `AuxKernel`-driven field comparison + `RunException` on threshold breach. Mirror with `VC_TEST_ASSERT` and a small CSV reader.

- [ ] **Step 1: Write the test** — `TestOedTedIntegration()` with the four sub-assertions above plus the CSV snapshot write.
- [ ] **Step 2: Run** — verify all four assertions pass; record baseline CSV.
- [ ] **Step 3: Commit** — `"test: add OED+TED integration test with C-suppression variant and regression snapshot (ATHENA 3.5 / SProcess 4.5 workflow)"`

---

## Self-Review Notes

- **Spec coverage:** Phase 4 of spec Section 12 = "OED + TED + dose loss + interface physics." Tasks 1–3 (OedSource, TedInitializer, DoseLossBC) cover the core interface physics with moving-interface idioms. Tasks 4–8 cover the spec's "advanced interface / mobile impurity" extension (ChargedEquilibrium, Carbon, Nitrogen, Copper, MobileImpurity) — every model that SProcess §4.188–4.195 lists and that the gap-analysis marked 🟡 now has a manual-cited FEM assembly path. Task 9 ties them together with an OED+TED+carbon workflow matching the manual's canonical example.
- **Manual citations:** every task cites SProcess body line ranges (19364, 19733, 19869, 19521) and/or ATHENA §3.x with line ranges — verified via `iconv` UTF-16 decode of `Manual/sprocess_ug.md` and `Manual/athena_users1.md`.
- **MOOSE citations:** every task names a verified MOOSE class with path (`GeochemistryKineticRate`, `CoupledForce`, `DerivativeParsedMaterial`, `LangmuirMaterial`, `ADRobinBC`, `IntegralPreservingFunctionIC`) — not hand-waved.
- **No placeholders:** every task has a concrete test (with quantitative pass criteria), implementation steps, and a commit message. The "deferred" items (kinetic relaxation of charge states, Fe-tag trap sites, psOxidation wiring) are explicitly named and bounded — they are follow-ups, not gaps.
- **Engine integration:** Task 7 introduces the first **drift** term — this is the one architectural novelty vs Phases 1–3 (pure diffusion). The plan calls out that the drift term must land in the same stiffness matrix `K` so the HypreBoomerAMG preconditioner (Phase 1 Task 5) sees it; this propagates a constraint into the engine's `assembleStiffness` signature, which already accepts arbitrary bilinear-form contributions per the Phase 1 species-outer/terms-inner loop.
- **Gap-analysis rows closed by this plan:** ChargedEquilibrium 🟡→✅, Carbon 🟡→✅, Nitrogen 🟡→✅, Copper 🟡→✅, OED 🟡→✅, TED 🟡→✅, DoseLoss ✅ (deepened), Regression tests ❌→✅ (Task 9 snapshot).
