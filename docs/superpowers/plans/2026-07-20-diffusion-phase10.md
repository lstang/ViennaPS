# Diffusion Engine Phase 10: PDE API + Results Extraction + Parameter DB Calibration

**Goal:** Add C++ function-object PDE specification API (Alagator equivalent), results extraction tools (1D cuts, dose, sheet resistance, fitting), and calibrated parameter database from manuals.

**Depends on:** Phase 2 (Fermi, IntrinsicCarrier)

**Manual sources (verified via `iconv` UTF-16 decode of `Manual/sprocess_ug.md`):**
- SProcess **Chapter 7 "Alagator"**, body lines **69515–71957** (TOC p.819–844) — PDE specification language
- SProcess **Chapter 14 "Extracting Results"**, body lines **90238–91219** (TOC p.871–887) — 1D cuts, dose, fitting
- SProcess §"Parameter Database" body lines **2810–2940** — `pdbSet`/`pdbGet`/Arrhenius builtin, inheritance, blend

**Solver stack note:** Phases 1-9 use Hypre (via MFEM wrappers, zero build cost) for all linear solves: `HypreBoomerAMG` as preconditioner inside SUNDIALS CVODE/IDA SPBCGS, and inside MFEM `NewtonSolver` for nonlinear problems. PETSc/TAO is deferred to Phase 10 only, where `TSAdjoint` (adjoint-based D(x) inversion, Task 15) and `TAO` (optimization, Tasks 12-15) have no Hypre or MFEM equivalent. See `docs/superpowers/research/2026-07-22-hypre-deep-dive.md` for the full solver strategy.

## Manual Equation References (SProcess Ch.7 + Ch.14)

**Alagator PDE form (SProcess §"Basics of Specifying PDEs", line 69779; example at line 69866):**

```
pdbSetString Silicon CX Equation "ddt(CX) - [Arrhenius 0.138 1.37]*grad(CX)"
```

- `ddt(u)` = ∂u/∂t (SProcess line 69700)
- `grad(u)` = ∇u; applied in a divergence context, `grad` auto-becomes `∇·(D∇u)` (SProcess line 69714 — "auto-divergence")
- `[Arrhenius D0 Ea]` = `D0·exp(−Ea/kT)` builtin (SProcess line 2836)
- Interface equations: `Equation_<mat>` form (line 69929) — distinct PDE on each material side of an interface
- Subexpression reuse (line 70001): name a sub-term, reference it in multiple equations

This is the symbolic-string PDE spec that Phase 10 Task 1 explicitly does NOT replicate (decision: C++ `std::function` lambdas instead). Cite the Alagator form as the conceptual reference; the lambda API is strictly richer (arbitrary C++ computations, not just parsed expressions).

**TR-BDF2 time discretization (SProcess §"Time Integration", line 91613; default since N-2017.09):**

```
TR step (t_n → t_n+γ):    u_γ = u_n + γ·Δt/2 · (f(u_n) + f(u_γ))           (trapezoidal)
BDF2 step (t_n → t_{n+1}): u_{n+1} = (1/(γ(2−γ)))·u_γ − ((1−γ)²/(γ(2−γ)))·u_n + (Δt·(1−γ)/(2−γ))·f(u_{n+1})
```

with `γ = 2 − √2`. Local truncation error estimated via Milne's device (SProcess line 91634) drives adaptive time stepping. SProcess uses this for all PDE time integration; ViennaPS uses SUNDIALS CVODE/BDF (Phase 1) which is mathematically equivalent (BDF orders 1–5).

**Results extraction (SProcess Ch.14):**
- `select` / `select list` (line 90308) — choose fields for output
- `print.1d` / `plot.1d` / `slice` (line 90390, 90394, 90442) — 1D data cuts
- `layers` (line 90529) — returns top/bottom/integral/material per layer (dose computation)
- `interpolate` (TOC p.876) — level crossings
- `extract` (TOC p.876) — values during diffuse step
- `SheetResistance` (cmd ref p.1193) — R_s calculation
- `FitArrhenius` (p.878), `FitLine` (p.878), `FitPearson` (p.879), `FitPearsonFloor` (p.879) — least-squares fitting utilities

**Parameter database (SProcess §"Parameter Database", line 2810):**
- `pdbSet <mat> <species> <param> {value}` — write (line 2815)
- `[Arrhenius D0 Ea]` builtin (line 2836) — uses global Tcl temp var; `SetTemp` (line 2844)
- `DiffLimit` (line 2857) — floor on diffusivity for numerical stability
- Inheritance: SiGe inherits from Si (line 2875); "like materials" blend via interpolation (line 2940)



## File Structure

| File | Responsibility |
|------|---------------|
| `PdeTerm.hpp` | Base class for PDE terms |
| `PdeEquation.hpp` | Collection of terms for one species |
| `PdeBC.hpp` | Boundary condition types (Dirichlet, Neumann, segregation, flux) |
| `ResultsExtractor.hpp` | 1D cuts, dose, level crossings, sheet resistance |
| `Fitting.hpp` | FitArrhenius, FitLine, FitPearson, FitPearsonFloor |
| `CalibratedParameters.hpp` | Calibrated D0, Ea, cluster params from manuals |

## Tasks

### Task 1: PdeTerm Base + Built-in Terms

- `PdeTerm<NumericType>` abstract base. Built-in: `DiffusionTerm` (D*grad(C)), `ReactionTerm` (R(C)), `GrowthTerm` (velocity*phi). Each term contributes to stiffness/reaction/mass.

**Coefficient type decision (must be stated — Alagator clarification):** the original Alagator (from the Stanford process-simulation lineage) is a *symbolic* PDE-specification language that generates weak forms from strings. **MOOSE does not have Alagator** — verified: `framework/include/functions/` (33 headers) contains only the Fparser-based family (`MooseParsedFunction.h`, `MooseParsedGradFunction.h`, `MooseParsedVectorFunction.h`) which gives `value(t,pt)` + auto-derived `gradient`/`timeDerivative` for *expressions*, but NOT symbolic weak-form assembly.

**Phase 10 decision:** `PdeTerm` coefficients are **C++ `std::function` lambdas** (e.g., `DiffusionTerm(std::function<NumericType(NumericType C, NumericType T)> D)`), not parsed strings. This is richer than MOOSE's Fparser approach and matches ViennaPS's template-heavy style. If a parsed-string surface is wanted later, wrap Fparser (`MooseParsedFunctionWrapper` reference) — note in commit that this gives free `gradient()` and `timeDerivative()` via AD, useful for `GrowthTerm` and Jacobian assembly. **True Alagator parity (symbolic weak-form generation from a string) is explicitly out of scope.**

- Test: create DiffusionTerm with constant D, Verify it produces same stiffness as ConstantDiffusion.
- Commit: `"feat: add PdeTerm base class with lambda coefficients (Alagator parity out of scope; Fparser wrap optional)"`

### Task 1.5: PdeIC — Initial Conditions for the PDE API

**Files:** Create `include/viennaps/fields/PdeIC.hpp`

**Produces:** initial-condition types attachable to a `PdeEquation`. The PDE API Tasks 2-3 have `.diffusion/.reaction/.growth/.boundary` but no way to set initial conditions — this task fills that gap.

**MOOSE references (verified — three canonical IC patterns):**
- `framework/include/ics/ConstantIC.h` — uniform value. `ConstantIC(value)`.
- `framework/include/ics/FunctionIC.h` / `FunctorIC.h` — analytic function. `FunctionIC(lambda<x,T>)`.
- `framework/include/ics/SolutionIC.h` — initialize from another solver's output (e.g., implant Monte Carlo result, prior diffusion step's GridFunction). This is the pattern Phase 4 Task 2 `TedInitializer` already uses.
- `framework/include/ics/IntegralPreservingFunctionIC.h` — **critical** — rescales so `∫C_projected = ∫C_source`. Required when mesh projection would otherwise silently change the dose.
- `framework/include/ics/RandomIC.h` — for stochastic initial perturbations (e.g., phase-field nucleation).

**Action:** define `PdeIC` abstract base + `ConstantIC`, `FunctionIC`, `FromFileIC` (wraps `SolutionIC` pattern), `PreservingIC` (wraps `IntegralPreservingFunctionIC` pattern), `RandomIC`. Attach via `PdeEquation::initialCondition(species, PdeIC)`.

- Test: each IC type sets a known field; verify GridFunction matches expected values; for `PreservingIC`, verify dose invariance under mesh refinement.
- Commit: `"feat: add PdeIC types (MOOSE ConstantIC/FunctionIC/SolutionIC/IntegralPreservingFunctionIC/RandomIC pattern)"`

### Task 2: PdeBC - Boundary Conditions
- `PdeBC<NumericType>` types: `DirichletBC(C0)`, `NeumannBC(flux)`, `SegregationBC(m, targetMaterial)`, `FluxBC(h)` (evaporation/Robin), `PeriodicBC()`. Each adds appropriate boundary integrator
- Test: DirichletBC holds value. SegregationBC enforces ratio
- Commit: `"feat: add PdeBC boundary condition types"`

### Task 3: PdeEquation - Composable PDE
- `PdeEquation<NumericType>` - builder pattern: `.diffusion(lambda)`, `.reaction(lambda)`, `.growth(lambda)`, `.boundary(name, BC)`. Composes terms into a DiffusionModel-equivalent that the engine can register
- Test: build custom equation, register in engine, solve, verify dose conservation
- Commit: `"feat: add PdeEquation composable PDE specification"`

### Task 4: Modifying Built-in Equations
- `engine.addEquationTerm("Dopant", extraTerm)`, `engine.subEquationTerm("Dopant", term)`, `engine.multiplyTerm("Dopant", factor)`. Allows extending/modifying built-in models without rewriting
- Test: add extra reaction term to Fermi model, verify modified behavior
- Commit: `"feat: add equation term modification API"`

### Task 5: ResultsExtractor - 1D Cuts

- `ResultsExtractor::extract1D(field, startPoint, endPoint, nPoints)` - linear interpolation through FEM solution. Returns depth-concentration profile.

**MOOSE reference (verified):** `framework/include/vectorpostprocessors/LineValueSampler.h` — MOOSE's exact analog. Static helper `generatePointsAndIDs(start, end, num_points, points, ids)` (lines 31-35) is directly portable. Also provides `spatialValue(p)` for transfers. Related samplers in the same dir: `LineFunctionSampler.h`, `PointValueSampler.h`, `SideValueSampler.h`, `NodalValueSampler.h`, `ElementValueSampler.h`.

**Design note:** the free-function API (`ResultsExtractor::extract1D`) is fine for ViennaPS's C++/Python binding style. If a future refactor wants *named, chainable* extraction objects (e.g., `extract1D` → `fitPearson` → `calcSheetResistance` as a DAG), the MOOSE reference is `framework/include/reporters/Reporter.h` lines 81-138 (`declareValue<T>` API with parallel modes). Out of scope for this task; just note it.

- Test: uniform field → flat 1D profile. Linear gradient → linear 1D profile. **Add:** non-trivial mesh with curved elements — verify the sampler walks element-by-element, not by global bounding-box interpolation.
- Commit: `"feat: add 1D data cut extraction (MOOSE LineValueSampler::generatePointsAndIDs pattern)"`

### Task 6: ResultsExtractor - Dose and Level Crossings

- `integrateField(field, region)` - area integral of concentration (dose). `findLevel(field, threshold, direction)` - depth where C = threshold. `findJunction(field)` - depth where C_dopant = C_substrate.

**MOOSE reference (verified):** `framework/include/vectorpostprocessors/LineMaterialRealSampler.h` + `LineMaterialSamplerBase.h` — sample material properties (mobility, D) along the same line as `extract1D`. **The non-obvious design win:** pull mobility(x) AND n(x) on the *same line points* as C(x), so sheet resistance (Task 7) integrates them consistently — no interpolation mismatch between fields sampled at different points.

Also: `modules/level_set/include/postprocessors/LevelSetVolume.h` — `ElementVariablePostprocessor` with threshold for partial-volume integrals. Useful for "integrate C only where C > C_ss" or "dose in the Si subdomain only."

- Test: known uniform C in known area → verify dose. Threshold crossing at known depth. **Add:** dose integrated only over Si subdomain (use MeshAttributes from Phase 1 Task 1) equals expected value.
- Commit: `"feat: add dose integration and level crossing extraction (MOOSE LineMaterialRealSampler + LevelSetVolume pattern)"`

### Task 7: ResultsExtractor - Sheet Resistance

- `calcSheetResistance(field, mobilityModel)` - R_s = 1 / integral(q*mu*n) dx. Mobility model: Caughey-Thomas with doping and temperature dependence.

**MOOSE reference (verified):** compose `LineValueSampler` (for n(x)) + `LineMaterialRealSampler` (for μ(x)) on the *same line points*, then numerically integrate `1/(q·μ·n)` along the line. Pulling μ(x) and n(x) at identical sample points avoids interpolation inconsistency — this is the MOOSE idiom and the right way to do sheet resistance.

**Uncertainty propagation (note for later):** when inputs come from a calibrated posterior (see new Tasks 13-15 below), the sheet resistance should report a confidence interval. MOOSE reference: `modules/stochastic_tools/include/reporters/StatisticsReporter.h` with `ReporterStatisticsContext::finalize` + `BootstrapCalculator`. Optional/stretch — flag in commit.

- Test: uniform n=1e16, μ=1000 → R_s = 1/(q·μ·n·thickness). **Add:** non-uniform n(x) from a Gaussian profile, verify R_s matches direct numerical integration within 0.1%.
- Commit: `"feat: add sheet resistance calculation (compose LineValueSampler + LineMaterialRealSampler on same line points)"`

### Task 8: Fitting Utilities

- `FitArrhenius(T[], D[])` -> (D0, Ea). `FitLine(x[], y[])` -> (slope, intercept). `FitPearson(depth[], conc[])` -> (Rp, dRp, skewness, kurtosis). `FitPearsonFloor` with floor parameter.

**MOOSE reference (verified):** `framework/include/vectorpostprocessors/LeastSquaresFit.h` + `LeastSquaresFitHistory.h` — polynomial LSQ fit with order parameter, truncate-order flag, x/y scale+shift, sample-output vectors (lines 41-87). Exact MOOSE analog of `FitLine`/`FitPearson`. If a future refactor wants shared code, port this class. For tabular parameter lookup, `framework/include/functions/PiecewiseLinearFromVectorPostprocessor.h` + `PiecewiseConstantFromCSV.h` are the references (relevant for Tasks 9-11 calibrated DB).

- Test: generate data from known Arrhenius, fit, verify recovery of D0, Ea.
- Commit: `"feat: add fitting utilities (MOOSE LeastSquaresFit pattern)"`

### Task 9: CalibratedParameters - Dopant Diffusivities
- `CalibratedParameters<NumericType>` - literature values from manuals: D0, Ea for B, P, As, Sb, In, Ga in Si, SiO2, Si3N4, PolySi. I/V equilibrium and diffusivity. Cluster binding energies. Segregation coefficients
- Test: verify D_B(Si, 1000C) ~ 1e-13 cm^2/s (literature value)
- Commit: `"feat: add calibrated dopant diffusivity parameters"`

### Task 10: CalibratedParameters - Cluster and Interface
- {311} binding energy, loop growth rate, BIC binding. Segregation m(T) for B/P/As at Si/SiO2, Si/Si3N4. Oxidation rate constants
- Test: verify {311} binding energy ~ 2.0 eV (literature)
- Commit: `"feat: add calibrated cluster and interface parameters"`

### Task 11: Parameter DB Inheritance and Blend
- Extend `ParameterDatabase` with like-materials interpolation: Si_{1-x}Ge_x blends Si and Ge parameters with x-dependent weighting. Inheritance chain: PolySi -> Si, DopedSiO2 -> SiO2
- Test: Si0.5Ge0.5 D_B should be between Si and Ge values
- Commit: `"feat: add parameter inheritance and like-materials blending"`

### Task 12: Global Sensitivity Analysis (Sobol / Morris)

**Goal:** rank which Arrhenius parameters (D0, Ea for B/P/As; cluster binding energies; segregation coefficients) dominate junction depth / sheet resistance variance. Answers "which parameter do I need to calibrate most carefully?"

**MOOSE references (verified — `stochastic_tools` module):**
- `stochastic_tools/include/samplers/SobolSampler.h` — Sobol quasi-random sampling for total-order sensitivity indices.
- `stochastic_tools/include/samplers/MorrisSampler.h` — elementary-effects screening (cheaper than Sobol; good for first pass with many parameters).
- `stochastic_tools/include/reporters/SobolReporter.h` + `MorrisReporter.h` — publish first-order, total-order, and elementary-effects indices.
- `stochastic_tools/include/samplers/LatinHypercubeSampler.h` — alternative sampling when Sobol's structured grid is too rigid.
- `stochastic_tools/include/multiapps/SamplerFullSolveMultiApp.h` — wraps the diffusion engine as a sub-app; pushes one parameter sample row to each sub-app via command-line args. Supports `batch` mode reusing one sub-app for memory efficiency.

**Action:** wrap `DiffusionEngine` as a `SamplerFullSolveMultiApp`-style driver. Sample the 8-12 dominant parameters (D0/Ea for B, I, V; {311} binding; BIC binding; m_seg). For each sample, run the engine to a fixed anneal endpoint, extract junction depth + sheet resistance (Tasks 5-7). Compute Sobol indices. Output a ranked table.

- Test: on a synthetic engine where only D_B contributes to output variance, Sobol total-order index for D_B ≈ 1.0 and for other parameters ≈ 0.
- Commit: `"feat: add Sobol/Morris global sensitivity analysis (MOOSE stochastic_tools SamplerFullSolveMultiApp + SobolReporter pattern)"`

### Task 13: Bayesian Calibration of D0/Ea from SIMS Depth Profiles

**Goal:** turn the manual `CalibratedParameters.hpp` from Task 9 into *actually calibrated* parameters with posterior credible intervals, fit to measured SIMS depth profiles.

**MOOSE references (verified — `stochastic_tools` module):**
- `stochastic_tools/include/samplers/PMCMCBase.h` — parallel MCMC.
- `stochastic_tools/include/samplers/AffineInvariantStretchSampler.h` — emcee (Foreman-Mackey); robust to correlated parameters, the standard choice for Arrhenius calibration.
- `stochastic_tools/include/samplers/IndependentGaussianMH.h` — Metropolis-Hastings baseline.
- `stochastic_tools/include/reporters/PMCMCDecision.h` / `AffineInvariantStretchDecision.h` / `IndependentMHDecision.h` — accept/reject decisions with `computeEvidence`, `computeTransitionVector`, `nextSamples`.
- `stochastic_tools/include/likelihoods/Gaussian.h` / `TruncatedGaussian.h` — likelihood of observed SIMS profile given simulated profile, assuming per-point Gaussian error.
- `stochastic_tools/include/reporters/StatisticsReporter.h` — posterior mean / STD / percentiles with bootstrap CIs.

**Action:** define the inverse problem: parameters θ = {D0, Ea, ...}, observed data = SIMS depth profile (with measurement uncertainty), forward model = diffusion engine. Run emcee for ~10⁴-10⁵ steps. Posterior gives joint credible intervals on D0/Ea accounting for parameter correlations. Compare Task 9 literature values against posterior medians — disagreements flag stale literature.

- Test: generate synthetic SIMS data from known D0/Ea with known noise; recover D0/Ea within posterior 95% CI; verify CI coverage on a sweep of noise levels.
- Commit: `"feat: add Bayesian calibration of D0/Ea from SIMS via emcee (MOOSE AffineInvariantStretchSampler + Gaussian likelihood pattern)"`

### Task 14: Gaussian-Process Surrogate of the Diffusion Engine

**Goal:** after calibration, train a GP surrogate so parameter sweeps and design-of-experiments become O(ms) instead of O(min). The GP also returns predicted std — the basis for active-learning calibration.

**MOOSE references (verified — `stochastic_tools` module):**
- `stochastic_tools/include/surrogates/GaussianProcessSurrogate.h` + `trainers/GaussianProcessTrainer.h` + `utils/GaussianProcess.h`.
- Covariance kernels: `utils/covariance/SquaredExponentialCovariance.h`, `MaternHalfIntCovariance.h`, `ExponentialCovariance.h`, `LMC.h` (multi-output).
- Hyperparameter tuning: `tuneHyperParamsAdam` with `GPOptimizerOptions` (`num_iter`, `batch_size`, `learning_rate`, `b1/b2/eps`).
- `utils/Standardizer.h` — input/output standardization.
- Also available: `PolynomialChaos.h` (with `computeSobolIndex`, `computeDerivative` — closes the loop with Task 12), `PolynomialRegressionSurrogate.h`, `NearestPointSurrogate.h`, `PODReducedBasisSurrogate.h`, libtorch neural nets (`LibtorchANNSurrogate.h`).
- The base class `SurrogateModel.h::evaluate(x, Real& std)` (verified lines 33-63) returns predicted value AND predicted std — exactly what active-learning calibration needs.

**Action:** train a GP on the (parameters, junction depth/sheet resistance) pairs from Task 12's sampler runs. Cross-validate with leave-one-out. The surrogate then replaces the engine for fast parameter sweeps; the predicted std drives `BayesianActiveLearner` (next-sample selection) for Task 13 follow-ups.

- Test: on the synthetic engine from Task 12's test, GP leave-one-out RMSE < 5% of output range; predicted std correlates with actual error.
- Commit: `"feat: add GP surrogate of diffusion engine (MOOSE GaussianProcessTrainer + Adam hyperparameter tuning pattern)"`

### Task 15 (stretch): Adjoint-Based Spatially-Varying Diffusivity Inversion

**Goal:** fit a spatially-varying D(x,T) field directly to measured profiles with analytic gradients via the adjoint. This is the most powerful form of calibration — instead of fitting scalar D0/Ea, fit the diffusivity field itself.

**MOOSE references (verified — `optimization` module):**
- `optimization/include/executioners/OptimizeSolve.h` — PETSc/TAO wrapper; 16 solvers including `NEWTON_TRUST_REGION`, `BOUNDED_NEWTON_LINE_SEARCH`, `QUASI_NEWTON`, `NELDER_MEAD` (lines 154-171).
- `optimization/include/optimizationreporters/OptimizationReporterBase.h` — `computeObjective()`, `computeGradient()`, `updateParameters()`, Tikhonov regularization (`_tikhonov_coeff`), bounds, equality/inequality constraints.
- **`optimization/include/vectorpostprocessors/ElementOptimizationDiffusionCoefFunctionInnerProduct.h`** — verified — the inner-product VPP specifically for inverting a spatially-varying diffusion coefficient field. This is literally the gradient of misfit w.r.t. D(x) at every element.
- `optimization/include/functions/ParameterMeshFunction.h` — parameter as an FE field (not a scalar vector), enabling spatial regularization.
- `optimization/include/userobjects/AdjointSolutionUserObject.h` + `executioners/TransientAndAdjoint.h` — adjoint transient solve paired with forward.
- Related inner-product VPPs (same module): `ElementOptimizationReactionFunctionInnerProduct.h`, `ElementOptimizationSourceFunctionInnerProduct.h`, `SideOptimizationFunctionInnerProduct.h`.

**Action:** this is a Phase 10+ stretch goal. Flag it for after Tasks 12-14 land. The reference pattern is: forward diffusion solve → misfit to SIMS → adjoint solve → `ElementOptimizationDiffusionCoefFunctionInnerProduct` gives ∂misfit/∂D(x) → TAO updates D(x) → repeat. Useful when the literature diffusivity is known to be wrong in specific device regions (e.g., strain, damage).

- Test: synthetic SIMS from known D(x) = piecewise-constant; verify inversion recovers D(x) within regularization tolerance.
- Commit: `"feat: add adjoint-based D(x,T) inversion (MOOSE OptimizeSolve + ElementOptimizationDiffusionCoefFunctionInnerProduct pattern)"`

### Task 16: Integration Test - Full PDE API Workflow
- Define custom diffusion equation via PDE API with diffusion + reaction + segregation BC. Solve. Extract 1D profile. Fit Pearson. Calculate sheet resistance. Verify end-to-end workflow
- Commit: `"test: add full PDE API workflow integration test"`
