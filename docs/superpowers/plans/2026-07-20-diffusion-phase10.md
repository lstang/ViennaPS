# Diffusion Engine Phase 10: PDE API + Results Extraction + Parameter DB Calibration

**Goal:** Add C++ function-object PDE specification API (Alagator equivalent), results extraction tools (1D cuts, dose, sheet resistance, fitting), and calibrated parameter database from manuals.

**Depends on:** Phase 2 (Fermi, IntrinsicCarrier)

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
- `PdeTerm<NumericType>` abstract base. Built-in: `DiffusionTerm` (D*grad(C)), `ReactionTerm` (R(C)), `GrowthTerm` (velocity*phi). Each term contributes to stiffness/reaction/mass
- Test: create DiffusionTerm with constant D, verify it produces same stiffness as ConstantDiffusion
- Commit: `"feat: add PdeTerm base class with built-in diffusion/reaction/growth terms"`

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
- `ResultsExtractor::extract1D(field, startPoint, endPoint, nPoints)` - linear interpolation through FEM solution. Returns depth-concentration profile
- Test: uniform field -> flat 1D profile. Linear gradient -> linear 1D profile
- Commit: `"feat: add 1D data cut extraction"`

### Task 6: ResultsExtractor - Dose and Level Crossings
- `integrateField(field, region)` - area integral of concentration (dose). `findLevel(field, threshold, direction)` - depth where C = threshold. `findJunction(field)` - depth where C_dopant = C_substrate
- Test: known uniform C in known area -> verify dose. Threshold crossing at known depth
- Commit: `"feat: add dose integration and level crossing extraction"`

### Task 7: ResultsExtractor - Sheet Resistance
- `calcSheetResistance(field, mobilityModel)` - R_s = 1 / integral(q*mu*n) dx. Mobility model: Caughey-Thomas with doping and temperature dependence
- Test: uniform n=1e16, mu=1000 -> R_s = 1/(q*mu*n*thickness)
- Commit: `"feat: add sheet resistance calculation"`

### Task 8: Fitting Utilities
- `FitArrhenius(T[], D[])` -> (D0, Ea). `FitLine(x[], y[])` -> (slope, intercept). `FitPearson(depth[], conc[])` -> (Rp, dRp, skewness, kurtosis). `FitPearsonFloor` with floor parameter
- Test: generate data from known Arrhenius, fit, verify recovery of D0, Ea
- Commit: `"feat: add fitting utilities (Arrhenius, Line, Pearson)"`

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

### Task 12: Integration Test - Full PDE API Workflow
- Define custom diffusion equation via PDE API with diffusion + reaction + segregation BC. Solve. Extract 1D profile. Fit Pearson. Calculate sheet resistance. Verify end-to-end workflow
- Commit: `"test: add full PDE API workflow integration test"`
