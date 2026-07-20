# Diffusion Engine Phase 11: Adaptive Mesh Refinement

**Goal:** Add adaptive mesh refinement (AMR) during diffusion: refine/coarsen based on solution gradients, dose error, and user-specified criteria. Uses MFEM's built-in AMR.

**Depends on:** Phase 1 (DiffusionEngine, LevelSetToMesh)

## File Structure

| File | Responsibility |
|------|---------------|
| `AdaptiveMeshRefiner.hpp` | AMR controller: criteria evaluation, refine/coarsen |
| `RefinementBox.hpp` | User-specified static refinement region |
| `MeshQualityEstimator.hpp` | Skewness, Jacobian, aspect ratio checks |

## Tasks

### Task 1: RefinementBox - Static Refinement
- `RefinementBox<NumericType, D>` - user-specified rectangular region with max element size. Before diffusion, refine mesh inside box to specified resolution. Uses MFEM `Mesh::GeneralRefinement` with attribute-based region selection
- Test: 8x8 mesh, refine box in corner -> verify corner elements smaller than rest
- Commit: `"feat: add RefinementBox for static mesh refinement"`

### Task 2: MeshQualityEstimator
- `MeshQualityEstimator` - computes per-element quality metrics: skewness, Jacobian determinant, aspect ratio. Identifies elements needing refinement or remeshing
- Test: regular mesh -> all quality ~ 1.0. Distorted mesh -> low quality elements flagged
- Commit: `"feat: add MeshQualityEstimator for element quality assessment"`

### Task 3: Relative Difference Refinement Criterion
- `RelativeDifferenceCriterion<NumericType>` - refine elements where |C_new - C_old| / max(|C|, eps) > threshold. Coarsen where below. Tracks solution change between time steps
- Test: sharp profile -> refine near front. Uniform profile -> no refinement
- Commit: `"feat: add relative difference AMR criterion"`

### Task 4: Gradient Refinement Criterion
- `GradientCriterion<NumericType>` - refine elements where |grad(C)| > threshold. Captures steep concentration gradients (implant profiles, junctions)
- Test: error function profile -> refine near steepest gradient
- Commit: `"feat: add gradient-based AMR criterion"`

### Task 5: Local Dose Error Criterion
- `LocalDoseErrorCriterion<NumericType>` - refine elements where integral of |C| change exceeds threshold. Ensures dose accuracy in critical regions
- Test: high-dose region -> refined
- Commit: `"feat: add local dose error AMR criterion"`

### Task 6: Logarithmic and Asinh Criteria
- `LogarithmicCriterion` - |log(C_new/C_old)| > threshold. Better for concentrations spanning many orders of magnitude
- `AsinhCriterion` - inverse hyperbolic sine difference. Handles sign changes and wide dynamic range
- Test: C from 1e10 to 1e20 -> log criterion refines where orders of magnitude change
- Commit: `"feat: add logarithmic and asinh AMR criteria"`

### Task 7: AdaptiveMeshRefiner Controller
- `AdaptiveMeshRefiner<NumericType, D>` - orchestrates AMR. Collects active criteria, evaluates per-element, marks for refine/coarsen, calls MFEM refinement, transfers solution via L2 projection. Configurable: max refinement level, min element size, coarsen threshold
- Test: multiple criteria active -> combined marking. Solution transferred correctly after refinement
- Commit: `"feat: add AdaptiveMeshRefiner controller"`

### Task 8: AMR During Diffusion
- Integrate AMR into `DiffusionEngine::solve()`. After every N time steps (configurable), evaluate criteria, refine/coarsen, transfer solution, continue. Avoids refining every step (expensive)
- Test: implant profile (sharp) -> AMR refines near junction, coarse away. Solution accuracy maintained
- Commit: `"feat: integrate AMR into DiffusionEngine time stepping"`

### Task 9: AMR During Moving Boundary
- For oxidation: after mesh deformation, check quality. If degraded, trigger remesh + AMR. Refine near moving interface. Solution transfer via L2 projection
- Test: oxidation step -> mesh deforms -> quality drops -> remesh with AMR near interface
- Commit: `"feat: add AMR for moving boundary problems"`

### Task 10: Interface-Aligned Refinement
- Automatic refinement at material interfaces (where element attribute changes between neighbors). Densifies mesh near Si/SiO2, Si/Si3N4 boundaries for accurate segregation
- Test: 2-material mesh -> verify finer elements near interface
- Commit: `"feat: add interface-aligned mesh refinement"`

### Task 11: Uniform Mesh Scaling
- Global mesh refinement/coarsening by factor. Useful for mesh convergence studies. `engine.refineGlobally(factor)`, `engine.coarsenGlobally(factor)`
- Test: 8x8 -> refine 2x -> 16x16. Verify solution converges
- Commit: `"feat: add uniform mesh scaling for convergence studies"`

### Task 12: Mesh Convergence Test
- Solve same diffusion problem on progressively finer meshes. Verify solution (1D profile, dose) converges. Richardson extrapolation for error estimate
- Test: 4x4, 8x8, 16x16, 32x32 -> profiles converge. Dose error decreases
- Commit: `"test: add mesh convergence validation test"`
