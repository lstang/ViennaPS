# Diffusion Engine Phase 6: SiGe/SiGeC + III-V Compound Diffusion

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans.

**Goal:** Add SiGe interdiffusion (defect-mediated), boron diffusion in SiGe (bandgap narrowing), carbon suppression of B TED (SiGeC), Ge–B pairing, strain effects, and III-V compound semiconductor diffusion (GaAs, InP, InGaAs) with sublattice-specific donor/acceptor mechanisms.

**Depends on:** Phase 3 (CDD, point-defect equilibrium), Phase 2 (Fermi, IntrinsicCarrier)

**Manual sources (verified via `iconv` UTF-16 decode):**
- SProcess §"Dopant Diffusion in Silicon Germanium", body lines **26775–27900+** (TOC p.243–247)
- SProcess §"Diffusion in III–V Compounds", body lines **27974–29900** (TOC p.251)
- ATHENA §3.9 "SiGe/SiGeC", lines **16997–17180** (p. 3-95); §3.8 "Compound Semiconductor", lines **16869–16996** (p. 3-93)

**Solver backend:** No solver changes. SiGe interdiffusion and III-V diffusion use the same HypreBoomerAMG + SUNDIALS CVODE/IDA infrastructure from Phases 1 and 3. Bandgap and strain effects modify coefficients (D, n_i) only, not the solver structure. Ge–B pairing is a kinetic reaction term (Phase 3 `KernelTerm`).

## Global Constraints

Same as Phase 1. MFEM-gated. Namespace `viennaps`. LLVM style. Reuses `KernelTerm` (Phase 3), `IntrinsicCarrier` (Phase 2), `BandgapModel` (existing), `PointDefectEquilibrium` (Phase 3).

---

## File Structure

| File | Responsibility | Status (HEAD) |
|------|----------------|---------------|
| `fields/BandgapModel.hpp` | E_g(x_Ge, T, strain) → n_i ratio | exists (has `niRatioToSi`) |
| `fields/models/SiGeDiffusion.hpp` | SiGe interdiffusion + B in SiGe | exists (1D stub) |
| `fields/models/SiGeCDiffusion.hpp` | Carbon I-trapping in SiGe context | exists (1D `tedFactor`) |
| `fields/models/GeBPairing.hpp` | Ge–B pairing (kinetic, immobile pair) | exists (1D stub) |
| `fields/models/StrainDiffusionModifier.hpp` | Strain effect on D and point-defect eq. | exists (stub) |
| `fields/MaterialConverter.hpp` | Convert mesh material (Si → GaAs) | exists (stub) |
| `fields/models/IIIVDiffusion.hpp` | III-V donor/acceptor sublattice D | exists (59-line stub) |

---

### Task 0: SiGe Concentration Regime (ADR)

**Required before Task 1.** SProcess exposes two SiGe interdiffusion regimes (lines 27068, 27328):

- **(A) Low-Ge-doped SiGe** (`x_Ge < ~0.2`, SProcess §"Low Germanium-Doped SiGe Model", lines 27328–27500) — interdiffusion treated as a perturbation on Si; B diffusivity modified via bandgap narrowing (eq. 321) and point-defect parameter shifts. **Default for Tasks 2–3.**
- **(B) Low-to-High Ge-doped SiGe** (`x_Ge ∈ [0.2, 1.0]`, SProcess §"Low-to-High Ge-Doped SiGe Model", lines 27068–27330) — full interdiffusion between Si and Ge as distinct species; non-linear concentration-dependent D_SiGe. **Task 4.**

**Phase 6 default: ship both.** Task 2 implements the interdiffusion PDE (valid for both regimes); Tasks 3 and 5 add the SiGe-context modifications (bandgap, Ge–B pairing) that apply on top.

- [ ] **Step 1: Document** (done above).
- [ ] **Step 2: No code.**

---

### Task 1: BandgapModel — SiGe Bandgap and n_i Ratio

**Files:** Modify `fields/BandgapModel.hpp`

**Produces:** extends `BandgapModel<NumericType>` to compute the bandgap `E_g(x_Ge, T, strain)` and the intrinsic-carrier ratio `n_i(SiGe) / n_i(Si)` per quadrature point. The ratio enters the Fermi-dependent dopant diffusivity (Phase 2 Task 3) — this is the dominant effect of Ge on B diffusion (SProcess §"Bandgap Effect", lines 27334–27410).

**Manual equations (SProcess eq. 321 at line 27359; bowing equation at line 27367):**

```
n_i(SiGe) / n_i(Si) = exp(ΔE_g(x_Ge) / (2·kT))                      (SProcess 321)

ΔE_g(x_Ge) = E_g^SiGe(x_Ge) − E_g^Si
           = (E_g^Ge − E_g^Si)·x_Ge  +  bowing·x_Ge·(1 − x_Ge)
           ≈ 0.835·x_Ge·(1 − x_Ge)   (for low-Ge SiGe, SProcess line 27367) + linear part
```

Temperature dependence: `E_g(T) = E_g(0) − α·T²/(T+β)` (Varshni form for each endpoint).

**Strain contribution (SProcess §"SiGe Strain and Dopant Activation", lines 27808–27870):** strain shifts `E_g` by `ΔE_g^strain = a·ε_hydro` (deformation potential); this is captured by `StrainDiffusionModifier` (Task 6) and folded into `ΔE_g` here.

**Assembly:** pure coefficient evaluation — no PDE. `BandgapModel::niRatio(x_Ge, T, strain)` returns the scalar that Phase 2's `FermiDiffusion` multiplies into `n_i` at each QP.

**MOOSE reference (verified):** `framework/include/materials/ParsedMaterial.h` — parsed expression for the bowing formula; mirror with a hand-coded `niRatio()` (no parser needed for a closed-form expression).

**Parameters:** `E_g^Si(300K) = 1.12 eV`, `E_g^Ge(300K) = 0.66 eV`, bowing `b ≈ 0.37 eV` for unstrained SiGe. Varshni `α/β` from SIBERIA/SEMATECH calibration.

- [ ] **Step 1: Write failing test** — `TestBandgapSiGe()`: `E_g(Si,300K) ≈ 1.12`, `E_g(Ge,300K) ≈ 0.66`, `E_g(Si₀.₅Ge₀.₅,300K)` between (verify bowing: `E_g < 0.5·(1.12+0.66)`); `n_i ratio` grows exponentially with `x_Ge` at fixed T.
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — extend `BandgapModel` with the bowing equation + Varshni T-dependence.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: extend BandgapModel with SiGe bowing + Varshni + n_i ratio (SProcess 321, 27334–27410)"`

---

### Task 2: SiGe Interdiffusion

**Files:** Modify `fields/models/SiGeDiffusion.hpp`

**Produces:** `SiGeDiffusion<NumericType>` — defect-mediated interdiffusion of Si and Ge atoms. Ge is treated as a distinct mobile species in Si (and vice versa); the two species diffuse into each other, smoothing the Si/Ge interface. SProcess §"SiGe Interdiffusion Model", body lines 26775–27330.

**Manual equations (SProcess lines 26790–26850; ref [21], [22]):**

```
∂C_Ge/∂t = ∇·(D_inter · ∇C_Ge) + ∇·(D_inter·C_Ge·(1−x)·(Ω_V/kT)·∇σ_V + ...)

D_inter = D_V^*·(C_V/C_V^*) + D_I^*·(C_I/C_I^*)      (sum of vacancy + interstitial contributions)
```

- `D_V^*`, `D_I^*` = Ge diffusivity contributions via vacancy / interstitial mechanism (from `VacCStarFactor`, `IntCStarFactor` SProcess params at lines 27746–27750)
- The drift term `(Ω_V/kT)·∇σ_V` couples Ge diffusion to the local point-defect gradient (the "defect-mediated" part); SProcess exposes it via Alagator expressions.

**Two-regime handling (per Task 0):**
- **Low-Ge regime (A):** `D_inter ≈ const(x_Ge, T)` (linearized); the drift term is dropped. Ge diffuses as a simple Constant-diffusion species.
- **High-Ge regime (B):** `D_inter = D_inter(x_Ge, T)` (nonlinear); the full defect-mediated form is used.

**Assembly (KernelTerm composition, Phase 3 Task 0):**

```cpp
// On the Ge equation:
addTerm<DiffusionTerm>([x_Ge, T, C_I, C_V](...){ return D_inter(x_Ge, T, C_I, C_V); });
// Symmetric Si equation (Si diffuses into Ge with the same D_inter).
// Optional drift term (high-Ge regime):
addTerm<DriftTerm>(...);   // couples to point-defect gradient
```

**MOOSE reference (verified):** `framework/include/kernels/MatDiffusionBase.h` (nonlinear D via Material property) — exactly the pattern. The drift term reuses Phase 4 Task 7's drift integrator (the only existing drift infrastructure in the engine).

- [ ] **Step 1: Write failing test** — `TestSiGeInterdiffusion()`: sharp Si/Ge interface → profile smooths over time; total `∫(C_Si + C_Ge)` conserved (mass conservation); high-Ge regime `D_inter` increases with `x_Ge` (verify at two values).
- [ ] **Step 2: Run to verify failure** → FAIL (current `SiGeDiffusion` is 1D explicit step)
- [ ] **Step 3: Implement** — 2D/3D FEM interdiffusion with defect-mediated `D_inter` per QP.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add SiGe defect-mediated interdiffusion with vacancy+interstitial contributions (SProcess 26790–26850; MOOSE MatDiffusionBase pattern)"`

---

### Task 3: Boron Diffusion in SiGe — Bandgap Effect

**Files:** Modify `fields/models/SiGeDiffusion.hpp` (add `BInSiGe` mode)

**Produces:** B diffusivity modified by Ge content via the `n_i` ratio (Task 1). This is the Fermi model (Phase 2) with `n_i` replaced by `n_i(SiGe)`. SProcess §"Boron Diffusion", lines 17023–17139 (eq. 3-241 at line 17029).

**Manual equation (SProcess eq. 3-241 / SiGe B diffusion):**

```
D_B(SiGe) = D_B(Si) · (n_i(Si) / n_i(SiGe))^β_B · f_D(C_I, C_V, Ge-modified)
```

- `β_B` = charge-state exponent (typically 1 for B⁻ / V²⁺ pairing; `NIFACT.SIGE`, `EAFACT.SIGE` SProcess params)
- For high-Ge SiGe, B diffusivity is **enhanced** (because `n_i(SiGe) > n_i(Si)` shifts the Fermi level, increasing the charged-defect fraction).

**Assembly:** wrap Phase 2 `FermiDiffusion`'s `n_i` lookup with `BandgapModel::niRatio(x_Ge, T, strain)`. No new stiffness assembly — same `D_eff(C, T, n, p)` formula, different `n_i` input. The local `x_Ge` is read from the Ge `GridFunction` (Task 2) at each QP.

**MOOSE reference (verified):** `framework/include/materials/DerivativeParsedMaterial.h` — `n_i(C_Ge, T)` is a parsed expression; the B diffusivity is its composition. Mirror with hand-coded chain.

- [ ] **Step 1: Write failing test** — `TestBInSiGe()`: at fixed T, `D_B(Si₀.₃Ge₀.₇) > D_B(Si)` (enhanced by larger `n_i`); verify ratio matches `n_i` ratio raised to `β_B`.
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — `SiGeDiffusion` reads `x_Ge` QP-wise, multiplies `n_i` by the ratio, calls existing Fermi assembly.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add B-in-SiGe diffusion with bandgap-modified n_i (SProcess eq. 3-241; MOOSE DerivativeParsedMaterial pattern)"`

---

### Task 4: SiGeC — Carbon Suppression of B TED

**Files:** Modify `fields/models/SiGeCDiffusion.hpp`

**Produces:** `SiGeCDiffusion<NumericType>` — carbon in SiGe traps interstitials (same `C_s + I ⇌ C_sI` kinetics as Phase 4 Task 5, but the SiGe context has Ge-modified point-defect parameters). The net effect is **suppression of B TED** in SiGeC layers — the canonical application is HBT profile control. SProcess §"Carbon Cluster", line 34647; §"Carbon Diffusion Model", line 19521.

**Manual equations (SProcess 4.193 Carbon example, lines 19521–19564; SiGe point-defect factors at lines 27746–27750):**

Same kinetics as Phase 4 Task 5:

```
∂C_s/∂t = ∇·(D_C(T) ∇C_s) − k_f·C_s·C_I + k_r·C_sI
∂C_sI/∂t =                 k_f·C_s·C_I − k_r·C_sI
```

but with **SiGe-modified point-defect equilibrium** (SProcess `VacCStarFactor`, `IntCStarFactor` at lines 27746–27750):

```
C_I^*(SiGe) = C_I^*(Si) · IntCStarFactor(x_Ge, strain)
C_V^*(SiGe) = C_V^*(Si) · VacCStarFactor(x_Ge, strain)
```

These factors can be arbitrary Alagator expressions in SProcess; here they are calibrated polynomials in `x_Ge` and strain. The modified `C_I^*` changes both `k_r` (via `k_r = k_f·C_I^*`) and the B TED coupling (`D_B^pair ∝ C_I/C_I^*`).

**Assembly:** reuse Phase 4 Task 5's `CarbonDiffusion` model, but parameterize `IntCStarFactor`/`VacCStarFactor` from `BandgapModel` and `StrainDiffusionModifier`.

**MOOSE references (verified):** identical to Phase 4 Task 5 — `GeochemistryKineticRate`, `CoupledForce`. No new patterns.

- [ ] **Step 1: Write failing test** — `TestSiGeCCarbonSuppression()`: in Si₀.₃Ge₀.₇ with C present, B TED is suppressed vs the no-C case; the suppression is stronger than in pure Si (because `IntCStarFactor(SiGe)` changes the I equilibrium).
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — parameterize `CarbonDiffusion` with SiGe point-defect factors.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add SiGeC carbon suppression with Ge-modified point-defect factors (SProcess 27746–27750, 19521)"`

---

### Task 5: Ge–B Pairing

**Files:** Modify `fields/models/GeBPairing.hpp`

**Produces:** `GeBPairing<NumericType>` — Ge and B form an immobile, electrically active pair `GeB⁺`. Reduces mobile (electrically active) B concentration. SProcess §"Germanium–Boron Pairing", body lines 27868–27980 (eq. 332, 333 at line 27895).

**Manual equations (SProcess eq. 332, 333):**

```
Ge + B  ⇌  GeB⁺

∂C_GeB/∂t = k_f·C_Ge·C_B − k_b·C_GeB                (SProcess 332)
```

- `k_f`, `k_b` = forward/reverse rates (`pdbSet Silicon Germanium Boron Kf {n}`, `Kb {n}`)
- Ge diffusion (SProcess eq. 334): `∂C_Ge/∂t = ∇·(D_Ge ∇C_Ge)` with `D_Ge = Dstar` constant — the GeB pair does not contribute to Ge mobility (pair is immobile)
- B loses mobile concentration by the pairing sink: `∂C_B^mobile/∂t ⊃ −k_f·C_Ge·C_B + k_b·C_GeB`

**GeB cluster initialization (SProcess line 27976–27980):** `solution add name=GeB ifpresent="Germanium Boron" !negative` — the `GeB` species is added to the unknown vector only when Ge is present. Phase 6 mirrors this: register `GeB` as an engine species conditionally on Ge presence.

**Assembly (KernelTerm composition):**

```cpp
// Three species: Ge, B, GeB (GeB immobile)
// On Ge: DiffusionTerm(D_Ge)  (no GeB contribution per eq. 334)
// On B:  -k_f·C_Ge·C_B + k_b·C_GeB  (sink)
// On GeB: +k_f·C_Ge·C_B - k_b·C_GeB  (source)
```

**MOOSE references (verified):**
- `framework/include/interfacekernels/InterfaceReaction.h` — `k·[c₁c₂]` reaction; mirror as a volumetric ReactionTerm (`MatReaction` from a material rate).
- `modules/chemical_reactions/include/kernels/CoupledBEKinetic.h` — backward-Euler kinetic rate; relevant for the implicit time stepping of `GeB`.

- [ ] **Step 1: Write failing test** — `TestGeBPairing()`: high Ge + B → mobile B decays, `C_GeB` grows, equilibrium `C_GeB / (C_Ge·C_B) → k_f/k_b`. Mass conservation `C_B^total = C_B + C_GeB` within 0.1%.
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — three-species reaction system; conditionally register `GeB`.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add Ge–B pairing (immobile GeB+) with kinetic rates (SProcess eq. 332–334, 27868–27980)"`

---

### Task 6: Strain Effects on Diffusion

**Files:** Modify `fields/models/StrainDiffusionModifier.hpp`

**Produces:** `StrainDiffusionModifier<NumericType>` — strain (from lattice mismatch in SiGe on Si, or from oxidation stress) modifies dopant diffusivity and point-defect equilibrium. SProcess §"Strain Effects", lines 27708–27746; §"SiGe Strain and Dopant Activation", lines 27808–27870.

**Manual equations (SProcess lines 27708–27730):**

Strain enters via deformation potentials acting on point-defect formation energies:

```
C_I^*(strain) = C_I^*(0) · exp(−Ω_I · σ_hydro / kT)
C_V^*(strain) = C_V^*(0) · exp(−Ω_V · σ_hydro / kT)

E_g(strain)   = E_g(0)   + a_c·ε_hydro + b·ε_shear   (deformation potentials)
```

- `Ω_I`, `Ω_V` = point-defect relaxation volumes (~+0.5 Ω_Si for I, −0.2 Ω_Si for V)
- `σ_hydro` = hydrostatic stress (from Phase 7 MfemElasticityKernel — `StressKernel`)
- `a_c`, `b` = conduction-band deformation potentials

**Coupling:** reads `σ_hydro` from the elasticity solution (`MfemElasticityKernel`, already in repo). Modifies `BandgapModel::niRatio` (strain shift to `E_g`) and `PointDefectEquilibrium` (strain shift to `C_I^*, C_V^*`).

**Assembly:** pure coefficient modifier — folds strain into the same `n_i`/`C_I^*` lookups that Tasks 3 and 4 already use. No new PDE.

**MOOSE reference (verified):** `framework/include/tensors/TensorMechanics.h` is the strain→stress solver; `modules/solid_mechanics/include/materials/ComputeSmallStrain.h` provides the strain input. ViennaPS already has `MfemElasticityKernel` (the elasticity side); `StrainDiffusionModifier` is the consumer on the diffusion side.

- [ ] **Step 1: Write failing test** — `TestStrainModifier()`: under tensile hydrostatic strain, `C_I^*` increases (interstitial formation easier); `n_i` shifts via `E_g`; the B diffusivity changes accordingly.
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — strain-modified `BandgapModel` + `PointDefectEquilibrium`.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add strain effects on diffusion via deformation potentials (SProcess 27708–27730, 27808–27870)"`

---

### Task 7: MaterialConverter — Si → III-V

**Files:** Modify `fields/MaterialConverter.hpp`

**Produces:** `MaterialConverter<NumericType>` — converts mesh elements from one material to another (Si → GaAs, GaAs → InGaAs) with property inheritance. SProcess §"Material Conversion", lines 27989–28015; §"Material Merging", lines 28028–28054.

**Algorithm:**

1. Select elements by attribute or region.
2. Relabel attribute (`Si` → `GaAs`); update `MeshAttributes` (Phase 1 Task 1).
3. Re-resolve material parameters from the ParameterDatabase (Phase 10 Task 11 inheritance: GaAs inherits from a generic III-V base).
4. Optionally merge adjacent III-V elements of the same alloy (`Merge.Materials`).

**Coupling:** typically called before a III-V diffusion step (the substrate becomes GaAs). This is subdomain relabeling — Phase 4 Task 0 idiom A (`ElementSubdomainModifier` pattern). The dopant field on converted elements is reinitialized via `POLYNOMIAL_NEARBY` to preserve continuity.

**MOOSE reference (verified):** `framework/include/meshmodifiers/ThresholdElementSubdomainModifier.h` (Phase 4 Task 0 idiom A) — same pattern, different trigger (here: explicit conversion call, not a threshold).

- [ ] **Step 1: Write failing test** — `TestMaterialConverter()`: convert Si → GaAs, verify attribute change; verify parameter DB lookup gives GaAs D values, not Si.
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — attribute relabel + DB re-resolution.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add MaterialConverter for Si↔III-V substrate conversion (SProcess 27989–28054; MOOSE ElementSubdomainModifier pattern)"`

---

### Task 8: III-V Diffusion — Sublattice Donor/Acceptor Mechanisms

**Files:** Modify `fields/models/IIIVDiffusion.hpp`

**Produces:** `IIIVDiffusion<NumericType>` — donor (Si, Se, Ge) and acceptor (Be, Mg, Zn, C) diffusivities in GaAs/InP, governed by distinct point-defect mechanisms on the Ga vs As sublattice. SProcess §"Diffusion in III–V Compounds", lines 27974–29900; ATHENA §3.8, eq. 3-239 (donor) and 3-240 (acceptor) at lines 16891–16957.

**Manual equations (ATHENA eq. 3-239 donor, line 16891):**

```
D_donor = D_AV·(n/ni)  +  D_AV^2·(n/ni)²                       (ATHENA 3-239)

  Donors (Si, Se, Ge) diffuse via the Ga-sublattice vacancy V_Ga mechanism.
  D_AV, D_AV^2 are Arrhenius prefactors for the singly/doubly-charged V_Ga states.
  Si and Se are concentration-independent → only the first term contributes (ATHENA line 16905).
  Ge in GaAs is proportional to (n/ni)² → second term dominates.
```

**Manual equation (ATHENA eq. 3-240 acceptor, line 16929):**

```
D_acceptor = D_AI·(p/ni)  +  D_AI^2·(p/ni)²                    (ATHENA 3-240)

  Acceptors (Be, Mg, Zn, C) diffuse via the Ga-sublattice interstitial I_Ga mechanism.
  C is concentration-independent (first term only).
  Be, Mg ∝ (p/ni).
  Zn ∝ (p/ni)² — strongly concentration-dependent, classic "Zn box profile."
```

**Assembly (Fermi-like, but on the appropriate sublattice):**

```cpp
const NumericType ni = intrinsicCarrier_->ni(T, "GaAs");   // III-V n_i (different from Si)
const NumericType n_over_ni = n / std::max(ni, eps);
const NumericType p_over_ni = p / std::max(ni, eps);

NumericType D = 0;
if (isDonor_) {
  D  = D_AV_  * n_over_ni;
  D += D_AV2_ * n_over_ni * n_over_ni;       // ATHENA 3-239
} else {
  D  = D_AI_  * p_over_ni;
  D += D_AI2_ * p_over_ni * p_over_ni;       // ATHENA 3-240
}
D *= std::exp(-Ea_ / (kB*T));                // Arrhenius T-dependence
```

**Point-defect equilibrium in compound semiconductors (spec §4.7 gap row ❌):** GaAs has two sublattices (Ga, As), each with its own I/V equilibrium. `PointDefectEquilibrium` (Phase 3) must be extended to track `C_I_Ga`, `C_V_Ga`, `C_I_As`, `C_V_As` — four defect species instead of two. This is the deepest gap in §4.7; ship the donor/acceptor D-form first (no explicit sublattice-defect tracking — the `D_AV`/`D_AI` prefactors fold the defect chemistry into Arrhenius constants), then extend `PointDefectEquilibrium` in a follow-up.

**MOOSE references (verified):**
- `framework/include/materials/ParsedMaterial.h` — `(n/ni)^k` parsed expression; mirror with hand-coded powers.
- `modules/chemical_reactions/include/kernels/CoupledBEKinetic.h` — if explicit sublattice defect tracking is added later (kinetic form of the I_Ga/V_Ga equilibrium).

**Parameters (Phase 10 DB):** per dopant/material pair — `D_AV`, `D_AV²`, `D_AI`, `D_AI²`, `Ea` from SProcess Table for GaAs, InP, InGaAs.

- [ ] **Step 1: Write failing test** — `TestIIIVDiffusion()`: (a) Si in GaAs → D concentration-independent (verify flat); (b) Zn in GaAs → D ∝ (p/ni)² (verify box-profile-like behavior); (c) donor D increases with n, acceptor D increases with p.
- [ ] **Step 2: Run to verify failure** → FAIL (current `IIIVDiffusion` is a 59-line stub with rough constants)
- [ ] **Step 3: Implement** — sublattice-aware D formula per eq. 3-239/3-240; species lookup per donor/acceptor.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add III-V diffusion with Ga-sublattice donor/acceptor mechanisms (ATHENA eq. 3-239/3-240; SProcess III-V Ch.)"`

---

### Task 9: Integration Test — SiGe HBT + III-V Substrate

**Files:** Modify `tests/diffusion/testDiffusion.cpp`

**Scenario (matches SProcess §4.6 SiGe HBT workflow):**

1. Deposit Si₀.₃Ge₀.₇ base layer on Si collector.
2. Implant B into the SiGe base.
3. Anneal 800 °C for 10 s — B diffuses in SiGe (bandgap-enhanced, Task 3); Ge–B pairs form (Task 5); carbon co-implant variant suppresses TED (Task 4).
4. **Variant B:** convert substrate Si → GaAs (Task 7); implant Si (donor) and Be (acceptor); anneal — verify donor and acceptor follow eq. 3-239/3-240 (concentration-dependent vs independent).

**Assertions:**

- (a) B in SiGe diffuses faster than in Si at the same T (bandgap effect).
- (b) With Ge–B pairing on, mobile B is reduced; turning pairing off increases mobile B (regression sanity).
- (c) With C co-implant, B TED is suppressed.
- (d) In the III-V variant, Si profile is Gaussian-like (concentration-independent D), Zn profile is box-like ((p/ni)²).
- (e) **Regression snapshot** to `tests/diffusion/regression/sige_hbt_baseline.csv`.

- [ ] **Step 1: Write the test** with the assertions + snapshot.
- [ ] **Step 2: Run** → verify all pass; record baseline.
- [ ] **Step 3: Commit** — `"test: add SiGe HBT + III-V integration test with bandgap/pairing/C-suppression variants (SProcess 4.6 / ATHENA 3.8 workflow)"`

---

## Self-Review Notes

- **Spec coverage:** Phase 6 of spec Section 12 = "SiGe/SiGeC + III-V." Tasks 1–6 cover SiGe (bandgap, interdiffusion, B-in-SiGe, SiGeC, Ge–B pairing, strain); Tasks 7–8 cover III-V (material conversion, sublattice diffusion). Every spec §4.6 and §4.7 row that the gap analysis marked 🟡/❌ now has a manual-cited FEM assembly path.
- **Manual citations:** every task cites SProcess body line ranges (26790, 27334, 27708, 27746, 27868, 27989) and ATHENA §3.8/§3.9 with line ranges and equation numbers (321, 332, 333, 334, 3-239, 3-240) — verified via `iconv` UTF-16 decode.
- **MOOSE citations:** `MatDiffusionBase` (Task 2), `DerivativeParsedMaterial`/`ParsedMaterial` (Tasks 1, 3, 8), `InterfaceReaction`/`CoupledBEKinetic` (Tasks 5, 8), `ThresholdElementSubdomainModifier` (Task 7), `ComputeSmallStrain` (Task 6).
- **Engine integration:** reuses Phase 2 Fermi (Task 3), Phase 3 `KernelTerm` + `PointDefectEquilibrium` (Tasks 2, 4, 5, 8), Phase 4 Task 5 Carbon (Task 4), Phase 4 Task 7 drift integrator (Task 2 high-Ge), Phase 4 Task 0 subdomain relabel (Task 7). No new engine APIs.
- **Deferred (explicit, bounded):** explicit four-species sublattice I/V equilibrium in III-V (Task 8 — folds into `D_AV`/`D_AI` for now); full strain-dependent elastic coupling when oxidation stress is solved (Task 6 — modifier ready, awaits elasticity input).
- **Gap-analysis rows closed:** SiGe interdiffusion 🟡→✅ (FEM), Bandgap model ✅ (deepened), Boron D modified by Ge ✅ (deepened with eq. 3-241), Carbon suppression 🟡→✅ (SiGe-context), Strain effects 🟡→✅, Ge-B pairing 🟡→✅ (FEM, eq. 332–334), Cluster initialization ❌→✅ (Task 5 conditional GeB), Material conversion 🟡→✅, Species-specific D 🟡→✅ (donor/acceptor sublattice), I/V equilibrium in compound ❌→🟡 (folded into D_AV/D_AI; explicit sublattice tracking deferred).
