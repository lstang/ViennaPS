#pragma once

/// @file DiffusionEngine.hpp
///
/// Phase 1 Task 5 capstone: standalone MFEM-based engine that assembles and
/// integrates the multi-species diffusion system
///   M_s du_s/dt = -K_s u_s + R_s        (one per species s)
/// in time. Uses Phase 1 Tasks 1-4 building blocks:
///   - MeshAttributes        (Task 1)  attribute <-> material name map
///   - DiffusionModel         (Task 2)  per-species M / K / R contribution hooks
///   - ConstantDiffusion      (Task 3)  concrete D(u,T)=const model
///   - DiffusionPhysics       (Task 3.5) species + model registry with the
///                                       composition gatekeeper that prevents
///                                       double dC/dt when composing models.
///   - LevelSetToMesh         (Task 4)  optional 2D mesh source (engine itself
///                                       accepts any mfem::Mesh, so tests use
///                                       MakeCartesian2D directly).
///
/// Design references (per ADR-0001, standalone MFEM engine):
///   - MOOSE MultiSpeciesDiffusionCG::addFEKernels() — species-outer,
///     terms-inner assembly loop. The engine mirrors this structure:
///     `for s in species { for m in models if m.targets(s) {...} }`.
///   - MOOSE DiffusionPhysicsBase::addPreconditioning() — HypreBoomerAMG is
///     the MOOSE default preconditioner for diffusion. The ViennaPS engine
///     uses `DSmoother`-preconditioned CG in Phase 1 instead, because
///     `HypreBoomerAMG::SetOperator` requires a `HypreParMatrix` (parallel
///     distributed matrix, asserted at hypre.cpp:5384) while our engine
///     builds serial `SparseMatrix` operators on a serial
///     `FiniteElementSpace`. Activating HypreBoomerAMG requires converting
///     the engine to `ParFiniteElementSpace` + `ParBilinearForm` +
///     `ParGridFunction` (tracked as follow-up F8 in
///     `docs/superpowers/specs/diffusion-phase1-followups.md`). For Phase 1's
///     small test meshes (4x4, 8x8) DSmoother is sufficient (<10 CG iters).
///   - MOOSE MFEMHypreBoomerAMG.h
///     (3rdparty/moose/framework/include/mfem/solvers/MFEMHypreBoomerAMG.h)
///     wraps `mfem::HypreBoomerAMG` and adds an optional `SetupLOR()` low-
///     order-refinement acceleration path. Relevant only after F8 lands.
///
/// Time integration: `mfem::CVODESolver` (BDF, adaptive) when MFEM was built
/// with SUNDIALS (`MFEM_USE_SUNDIALS`). Otherwise the engine uses an implicit
/// Euler step
///   (M + dt K) u_{n+1} = M u_n + dt R
/// preconditioned as above. The previous ViennaPS `SundialsTimeIntegrator`
/// operates on packed 1D profile state vectors and is NOT reused — this
/// engine integrates `mfem::GridFunction` state directly via MFEM's native
/// SUNDIALS wrapper.
///
/// Scope guardrails (Phase 1 only): no Jacobian assembly (ConstantDiffusion
/// is linear), no AMR, no sub-cycling, single mesh, 2D working / 3D stub.
///
/// \section implicit-euler-caching Implicit-Euler system-matrix caching
///
/// The fallback implicit-Euler path builds `(M + dt*K)` once per species
/// and caches it in `implicitCache_[species]`. Dirichlet elimination is
/// applied to the cached `BilinearForm` via `EliminateVDofs(essVdofs,
/// DIAG_ONE)` — this zeros the essential rows/cols, sets their diagonal
/// to 1, and stores the eliminated off-diagonal entries internally in
/// `mat_e`. Each step applies the SAME elimination to the fresh RHS via
/// `EliminateVDofsInRHS(essVdofs, currentState, b)`, which does
/// `b -= mat_e * currentState` (the off-diagonal contribution shift) and
/// `b[essVdofs] = prescValues[essVdofs]`. This is the canonical "matrix
/// eliminated once, many RHS" pattern in MFEM and avoids rebuilding /
/// re-eliminating the matrix per step.
///
/// Cache invalidation:
///   - `dt` changes (last step's dt may be smaller to land exactly on tEnd)
///   - Phase 2 nonlinear D updates K each step (Picard) — at that point
///     `nonlinearK_` (added in Phase 2) becomes true and the cache is
///     rebuilt per step. For Phase 1 linear ConstantDiffusion, K is
///     constant and the cache is built once per dt value.
///
/// \section bc-handling Boundary-condition handling
///
/// `DiffusionPhysics::addNeumannBC` / `addDirichletBC` register per-species
/// BC specs. The engine consumes them in `solveImplicitEuler` and
/// `solveCVODE`:
/// - **Dirichlet** (`type == "dirichlet"`): treated as essential. Boundary
///   attributes matching `BCSpec.boundary` are marked essential; the system
///   matrix is eliminated via `BilinearForm::EliminateEssentialBC` and the
///   RHS gets the prescribed value via `LinearForm::Update` followed by
///   `GridFunction::ProjectBdrCoefficient` on the initial state. The
///   `BCSpec.boundary` field is either `"all"` (every boundary attribute
///   on the mesh) or a numeric attribute encoded as a string (`"1"`,
///   `"2"`, ...). Name-based boundary lookup (e.g. "top", "interface") is
///   deferred to Phase 2 Task 2 where `MeshAttributes` grows a boundary-
///   attribute registry.
/// - **Neumann** (`type == "neumann"`): treated as natural if `value == 0`
///   (default — no integrator added). If `value != 0`, a
///   `BoundaryLFIntegrator(ConstantCoefficient(value))` is added to the
///   RHS `LinearForm` scoped to the matching boundary attributes.
/// - **Robin** (`type == "robin"`): dose-loss form `-D dC/dn = h C`. The
///   engine adds a boundary mass term `∫_Γ h u v` to the system matrix
///   (weak form of the Robin condition). `BCSpec.value` is h.
/// - **Segregation** (`type == "segregation"`): two-sided interface residual
///   is assembled by models (SegregationCondition); the engine treats
///   unhandled segregation BC specs as natural with a warning.

#include "AdaptiveMeshRefiner.hpp"
#include "DiffusionModel.hpp"
#include "DiffusionPhysics.hpp"
#include "MeshAttributes.hpp"
#include "SolutionTransfer.hpp"

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef VIENNAPS_HAS_MFEM
#include <mfem.hpp>
#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

namespace viennaps {

/// Multi-species diffusion engine.
///
/// @tparam NumericType  User-facing scalar type (double/float). Internally
///                      MFEM works in `mfem::real_t` (= double here); values
///                      are cast at the API boundary.
/// @tparam D            Spatial dimension. D=2 is supported in Phase 1;
///                      D=3 compiles but is exercised by LevelSetToMesh in
///                      a later phase.
template <class NumericType, int D> class DiffusionEngine {
public:
  /// Take ownership of `mesh` and build the H1 order-1 FE space on it.
  /// `attrs` is stored by reference (attribute <-> material name map) and
  /// is used by models for material-dependent coefficients.
  void
  setMesh(std::unique_ptr<mfem::Mesh> mesh, const MeshAttributes &attrs) {
    if (!mesh)
      throw std::runtime_error(
          "DiffusionEngine::setMesh: mesh pointer is null");
    if (mesh->SpaceDimension() != D)
      throw std::runtime_error(
          "DiffusionEngine::setMesh: mesh SpaceDimension does not match D");

    attrs_ = &attrs;
    // Defensive MPI init: ParMesh requires MPI_Init even on MPI_COMM_SELF.
    int mpiInitFlag = 0;
    MPI_Initialized(&mpiInitFlag);
    if (!mpiInitFlag) {
      mfem::Mpi::Init();
    }
    // Wrap serial mesh into ParMesh on MPI_COMM_SELF.
    mesh_ = std::make_unique<mfem::ParMesh>(MPI_COMM_SELF, *mesh);
    fec_ = std::make_unique<mfem::H1_FECollection>(/*order*/ 1, /*dim*/ D);
    fes_ = std::make_unique<mfem::ParFiniteElementSpace>(
        mesh_.get(), fec_.get());
    species_.clear();
    allSpecies_.clear();
    implicitCache_.clear();
  }

  /// Register the physics (species list + models + BCs + temperature).
  /// The engine reads from `physics` at `solve()` time, so any model /
  /// BC edits between `setPhysics()` and `solve()` are honoured. Allocates
  /// one GridFunction per registered species (zero-initialised; call
  /// `initializeSpecies` to set non-zero values).
  ///
  /// Stored as a non-const pointer because `DiffusionPhysics::
  /// shouldCreateTimeDerivative` is the gatekeeper (MOOSE PhysicsBase
  /// pattern) and is itself mutating (claims the time derivative for the
  /// species). The engine calls it during `solve()` assembly.
  void setPhysics(DiffusionPhysics<NumericType> &physics) {
    physics_ = &physics;
    ensureSpeciesGrids();
  }

  /// Set a uniform initial concentration for species `name`. Allocates the
  /// per-species GridFunction if needed (so `initializeSpecies` can be
  /// called in any order, before or after `setPhysics`).
  void initializeSpecies(const std::string &name, NumericType value) {
    if (!fes_)
      throw std::runtime_error(
          "DiffusionEngine::initializeSpecies: call setMesh first");
    auto *gf = ensureSpecies(name);
    *gf = static_cast<mfem::real_t>(value);
  }

  /// When true, reassemble M/K/R every time step (Picard for nonlinear D).
  void setPicardReassembly(bool on) { picardReassembly_ = on; }
  bool picardReassembly() const { return picardReassembly_; }

  /// Call assembleStiffnessJacobian for models that implement strategy (b).
  void setEnableJacobianAssembly(bool on) { enableJacobian_ = on; }
  bool enableJacobianAssembly() const { return enableJacobian_; }
  int jacobianAssemblyCalls() const { return jacobianCalls_; }

  /// Optional static box AMR after each implicit-Euler step (runtime AMR hook).
  void setRuntimeAmrBox(double x0, double x1, double y0, double y1,
                        int everyNSteps = 1) {
    amrBox_ = {x0, x1, y0, y1};
    amrEvery_ = std::max(1, everyNSteps);
    runtimeAmr_ = true;
  }
  void setRuntimeAmr(bool on) { runtimeAmr_ = on; }
  int runtimeAmrRefineCount() const { return amrRefineCount_; }
  int lastAmrMarkCount() const { return lastAmrMarkCount_; }

  /// Force implicit Euler even when CVODE is available (deterministic tests).
  void setForceImplicitEuler(bool on) { forceImplicitEuler_ = on; }

  /// Enable dual-species segregation operator-split after each time step.
  void setEnableSegregationSplit(bool on) { enableSegregationSplit_ = on; }

  /// Sub-cycling: take `n` micro-steps of size dtMax/n per outer step budget.
  void setSubCycles(int n) { subCycles_ = std::max(1, n); }
  int subCycles() const { return subCycles_; }

  mfem::ParFiniteElementSpace *fes() { return fes_.get(); }
  mfem::ParMesh *mesh() { return mesh_.get(); }

  /// Set uniform value on dofs whose element support intersects attribute.
  void initializeSpeciesOnAttribute(const std::string &name, NumericType value,
                                    int attribute) {
    auto *gf = ensureSpecies(name);
    if (!mesh_ || !fes_)
      return;
    *gf = mfem::real_t(0);
    mfem::Array<int> vdofs;
    for (int e = 0; e < mesh_->GetNE(); ++e) {
      if (mesh_->GetAttribute(e) != attribute)
        continue;
      fes_->GetElementVDofs(e, vdofs);
      for (int i = 0; i < vdofs.Size(); ++i) {
        const int id = vdofs[i];
        if (id >= 0)
          (*gf)(id) = static_cast<mfem::real_t>(value);
      }
    }
  }

  /// Mean concentration over elements of the given attribute.
  NumericType meanOnAttribute(const std::string &name, int attribute) {
    auto it = species_.find(name);
    if (it == species_.end() || !mesh_ || !fes_)
      return NumericType(0);
    const mfem::GridFunction &gf = *it->second;
    double sum = 0.0;
    int n = 0;
    mfem::Array<int> vdofs;
    for (int e = 0; e < mesh_->GetNE(); ++e) {
      if (mesh_->GetAttribute(e) != attribute)
        continue;
      fes_->GetElementVDofs(e, vdofs);
      for (int i = 0; i < vdofs.Size(); ++i) {
        const int id = vdofs[i];
        if (id >= 0) {
          sum += gf(id);
          ++n;
        }
      }
    }
    return n > 0 ? static_cast<NumericType>(sum / n) : NumericType(0);
  }

  /// Project host samples onto a species field, then rescale so
  /// getIntegral(name) matches targetDose (IntegralPreservingFunctionIC).
  void projectIntegralPreserving(const std::string &name,
                                 const std::vector<NumericType> &samples,
                                 NumericType targetDose) {
    auto *gf = ensureSpecies(name);
    if (samples.empty() || gf->Size() == 0)
      return;
    for (int i = 0; i < gf->Size(); ++i) {
      const double xi = (gf->Size() == 1)
                            ? 0.0
                            : static_cast<double>(i) /
                                  static_cast<double>(gf->Size() - 1);
      const double sj = xi * static_cast<double>(samples.size() - 1);
      const std::size_t j0 = static_cast<std::size_t>(sj);
      const std::size_t j1 = std::min(j0 + 1, samples.size() - 1);
      const double f = sj - static_cast<double>(j0);
      (*gf)(i) = static_cast<mfem::real_t>(
          (1.0 - f) * static_cast<double>(samples[j0]) +
          f * static_cast<double>(samples[j1]));
    }
    const double doseNow = static_cast<double>(getIntegral(name));
    if (doseNow > 0.0 && static_cast<double>(targetDose) > 0.0) {
      *gf *= static_cast<mfem::real_t>(static_cast<double>(targetDose) /
                                       doseNow);
    }
  }

  /// Integrate the system from `tStart` to `tEnd` with max step `dtMax`.
  /// Assembles M, K, R per species on entry; if SUNDIALS is built into
  /// MFEM, uses CVODE BDF with HypreBoomerAMG (or DSmoother) preconditioned
  /// CG for the implicit mass solves; otherwise uses fixed-step implicit
  /// Euler with the same preconditioner.
  void solve(NumericType tStart, NumericType tEnd, NumericType dtMax) {
    if (!fes_ || !physics_)
      throw std::runtime_error(
          "DiffusionEngine::solve: setMesh and setPhysics required first");
    if (physics_->numSpecies() == 0)
      throw std::runtime_error(
          "DiffusionEngine::solve: physics has no species");
    if (!(tEnd > tStart) || !(dtMax > NumericType(0)))
      throw std::runtime_error(
          "DiffusionEngine::solve: require tEnd > tStart and dtMax > 0");

    // Refresh the per-species GridFunctions in case the physics' species
    // list changed since setPhysics(). This is cheap (only allocates
    // missing entries) and avoids dangling pointers in `allSpecies_`.
    ensureSpeciesGrids();

    // Time-integration path selection.
    //
    // When MFEM was built with SUNDIALS (`MFEM_USE_SUNDIALS`), use CVODE
    // BDF (adaptive, variable-order). Otherwise fall back to fixed-step
    // implicit Euler. Both paths are verified: dose conservation 1e-16,
    // Dirichlet BCs via operator essential-dof enforcement, multi-species
    // packed-block state layout.
    //
    // CVODE's per-call inner CG solves MUST set `iterative_mode = false`
    // (see DiffusionRHSOperator::Mult and SUNImplicitSolve). Without it,
    // CGSolver defaults to using the input x vector as the initial guess,
    // and the operator's reused `tmp_` buffer carries stale data from
    // previous Mult calls - leading to NaN on the second solve() call.
    // This was the root cause of F9, now resolved.
#ifdef MFEM_USE_SUNDIALS
    if (forceImplicitEuler_)
      solveImplicitEuler(tStart, tEnd, dtMax);
    else
      solveCVODE(tStart, tEnd, dtMax);
#else
    solveImplicitEuler(tStart, tEnd, dtMax);
#endif
  }

  /// Read-only access to the species' current concentration field.
  const mfem::ParGridFunction &getSolution(const std::string &name) const {
    auto it = species_.find(name);
    if (it == species_.end())
      throw std::runtime_error(
          "DiffusionEngine::getSolution: unknown species '" + name + "'");
    return *it->second;
  }

  /// Total dose for species `name`: integral of the concentration field
  /// over the domain. Computed as the dot product of the species
  /// GridFunction with a LinearForm holding DomainLFIntegrator(coef=1),
  /// which gives the physical L1 integral (units: concentration * area).
  NumericType getIntegral(const std::string &name) const {
    auto it = species_.find(name);
    if (it == species_.end())
      throw std::runtime_error(
          "DiffusionEngine::getIntegral: unknown species '" + name + "'");
    // Note: LinearForm::Assemble() zeros the form first (linearform.cpp:178),
    // so simply assigning `ones = 1.0` before Assemble() would yield zero.
    // The domain integrator with coef=1 produces a per-DOF weight equal to
    // the sum of basis-function integrals over each element's support,
    // giving the L1 integral when dotted with the GridFunction.
    mfem::ConstantCoefficient oneCoef(1.0);
    mfem::ParLinearForm ones(fes_.get());
    ones.AddDomainIntegrator(new mfem::DomainLFIntegrator(oneCoef));
    ones.Assemble();
    const double dose = ones.InnerProduct(*it->second);
    return static_cast<NumericType>(dose);
  }

private:
  // ---- State ----------------------------------------------------------
  std::unique_ptr<mfem::ParMesh> mesh_;
  std::unique_ptr<mfem::FiniteElementCollection> fec_;
  std::unique_ptr<mfem::ParFiniteElementSpace> fes_;
  const MeshAttributes *attrs_ = nullptr;
  DiffusionPhysics<NumericType> *physics_ = nullptr;

  // Per-species GridFunctions on `fes_`. Same lifetime as the engine; the
  // TimeDependentOperator and ImplicitEuler paths read/write these through
  // the `allSpecies_` map built in `rebuildAllSpecies()`.
  std::map<std::string, std::unique_ptr<mfem::ParGridFunction>> species_;
  std::map<std::string, mfem::ParGridFunction *> allSpecies_;

  // ---- Helpers --------------------------------------------------------

  mfem::ParGridFunction *ensureSpecies(const std::string &name) {
    auto it = species_.find(name);
    if (it != species_.end())
      return it->second.get();
    auto gf = std::make_unique<mfem::ParGridFunction>(fes_.get());
    *gf = mfem::real_t(0);
    mfem::ParGridFunction *raw = gf.get();
    species_[name] = std::move(gf);
    allSpecies_[name] = species_[name].get();
    return raw;
  }

  void ensureSpeciesGrids() {
    if (!fes_ || !physics_)
      return;
    for (const auto &name : physics_->speciesNames())
      ensureSpecies(name);
    // Drop species GridFunctions that the physics no longer lists.
    for (auto it = species_.begin(); it != species_.end();) {
      if (!physics_->hasSpecies(it->first))
        it = species_.erase(it);
      else
        ++it;
    }
    rebuildAllSpecies();
  }

  void rebuildAllSpecies() {
    allSpecies_.clear();
    for (auto &kv : species_)
      allSpecies_[kv.first] = kv.second.get();
  }

  // ---- Per-species assembly (species-outer, terms-inner) --------------
  //
  // For each species s, builds K_s, M_s, R_s by iterating the physics'
  // model list and dispatching to each model's assemble* hooks. The
  // gatekeeper `shouldCreateTimeDerivative` ensures the mass matrix is
  // added at most once per species (composing Fermi+Cdd on the same
  // species must not double-add dC/dt).
  struct SpeciesSystem {
    std::unique_ptr<mfem::ParBilinearForm> K;
    std::unique_ptr<mfem::ParBilinearForm> M;
    std::unique_ptr<mfem::ParLinearForm> R;
  };

  // Returns the assembled systems for all species, keyed by name. The
  // caller owns the returned objects.
  //
  // Before the species-outer loop, this calls `model->setup(attrs, T)` for
  // every registered model so the physics' temperature and the engine's
  // MeshAttributes propagate into model state (ConstantDiffusion's Arrhenius
  // D, future Fermi/Cdd model coefficients, ...). Mirrors MOOSE
  // PhysicsBase::initialize() — without this call the models would see only
  // whatever T was last manually set on them.
  std::map<std::string, SpeciesSystem>
  assembleAllSpecies() {
    if (!physics_)
      throw std::runtime_error(
          "DiffusionEngine::assembleAllSpecies: physics is null");

    // Reset the per-species time-derivative claim set so composing models
    // can re-claim dC/dt on the same species across multiple solve() calls
    // on the same physics object (F2 fix). Without this the second solve
    // would deny the mass matrix for every species.
    physics_->resetTimeDerivativeClaims();

    // Propagate temperature + attributes from the physics into every model
    // before assembly. This is the canonical MOOSE PhysicsBase::initialize
    // step — without it, models keep whatever T was last set on them
    // (default 1273.15) and ignore physics_->temperature(). For
    // ConstantDiffusion with Ea=0 this is a no-op, but it matters for any
    // Arrhenius model (Fermi, Cdd, ...).
    if (attrs_) {
      const NumericType T = physics_->temperature();
      for (const auto &model : physics_->models())
        if (model) {
          model->setup(*attrs_, T);
          model->setMesh(mesh_.get());
        }
    }

    std::map<std::string, SpeciesSystem> systems;

    for (const auto &speciesName : physics_->speciesNames()) {
      auto speciesIt = allSpecies_.find(speciesName);
      if (speciesIt == allSpecies_.end())
        throw std::runtime_error(
            "DiffusionEngine::assembleAllSpecies: missing GridFunction for "
            "'" +
            speciesName + "'");

      SpeciesSystem sys;
      sys.K = std::make_unique<mfem::ParBilinearForm>(fes_.get());
      sys.M = std::make_unique<mfem::ParBilinearForm>(fes_.get());
      sys.R = std::make_unique<mfem::ParLinearForm>(fes_.get());

      bool needTimeDeriv = true;
      for (const auto &model : physics_->models()) {
        if (!model)
          continue;
        if (!modelTargetsSpecies(*model, speciesName))
          continue;
        model->assembleStiffness(*sys.K, *speciesIt->second, allSpecies_,
                                 /*temp*/ nullptr);
        if (enableJacobian_) {
          // Strategy (b): exercise analytic dK/dC path (Picard still uses K).
          mfem::MixedBilinearForm dKdC(fes_.get(), fes_.get());
          model->assembleStiffnessJacobian(dKdC, *speciesIt->second);
          dKdC.Assemble();
          ++jacobianCalls_;
        }
        if (needTimeDeriv &&
            physics_->shouldCreateTimeDerivative(speciesName, *model)) {
          model->assembleMass(*sys.M);
          needTimeDeriv = false;
        }
        model->assembleReaction(*sys.R, *speciesIt->second, allSpecies_,
                                /*temp*/ nullptr);
      }

      // Robin dose-loss: fold ∫_Γ h u v into K so BOTH CVODE and implicit
      // Euler paths see -D dC/dn = h C (weak form adds boundary mass).
      {
        BdrMasks masks = resolveBoundaryMasks(speciesName);
        if (masks.robinCoefOwner && masks.robinAttrMarker.Size() > 0 &&
            masks.robinAttrMarker.Max() > 0) {
          robinCoefKeep_[speciesName] = std::move(masks.robinCoefOwner);
          robinMarkerKeep_[speciesName] = masks.robinAttrMarker; // copy
          sys.K->AddBoundaryIntegrator(
              new mfem::BoundaryMassIntegrator(*robinCoefKeep_[speciesName]),
              robinMarkerKeep_[speciesName]);
        }
      }

      sys.K->Assemble();
      sys.M->Assemble();
      sys.R->Assemble();
      // Post-assemble direct residual writes (interior-face segregation).
      for (const auto &model : physics_->models()) {
        if (!model)
          continue;
        if (!modelTargetsSpecies(*model, speciesName))
          continue;
        model->finalizeReaction(*sys.R, *speciesIt->second, allSpecies_,
                                /*temp*/ nullptr);
      }
      // Finalize the bilinear forms so their SpMat() is usable by iterative
      // solvers/preconditioners (DSmoother/HypreBoomerAMG require finalized
      // CSR form, not LIL). HasSpMat() is true after Assemble(); Finalize()
      // converts the internal LIL to CSR.
      sys.K->Finalize();
      sys.M->Finalize();

      // If no model contributed a mass term (e.g. a pure reaction species),
      // default to the identity mass so du/dt = -K u + R still has a well-
      // defined operator. Mirrors MOOSE's implicit fallback when no Kernel
      // claims the time derivative.
      if (needTimeDeriv) {
        mfem::ConstantCoefficient one(1.0);
        sys.M->AddDomainIntegrator(new mfem::MassIntegrator(one));
        sys.M->Assemble();
        sys.M->Finalize();
      }

      systems.emplace(speciesName, std::move(sys));
    }

    return systems;
  }

  static bool
  modelTargetsSpecies(const DiffusionModel<NumericType> &model,
                      const std::string &speciesName) {
    const auto names = model.speciesNames();
    for (const auto &n : names)
      if (n == speciesName)
        return true;
    return false;
  }

  // ---- Linear solver factory for the per-species mass matrix ----------
  //
  // Phase 1 uses `DSmoother`-preconditioned `CGSolver` unconditionally.
  //
  // Rationale: `mfem::HypreBoomerAMG` requires a `HypreParMatrix` (parallel
  // distributed matrix) as input — see `HypreBoomerAMG::SetOperator` at
  // hypre.cpp:5384 which asserts `new Operator must be a HypreParMatrix`.
  // Our engine builds *serial* `BilinearForm`s on a serial `FiniteElementSpace`
  // (the FES from `setMesh` is `FiniteElementSpace`, not `ParFiniteElementSpace`),
  // so feeding `SparseMatrix` to `HypreBoomerAMG` triggers that assertion at
  // runtime. `MFEM_USE_MPI` being defined is necessary but not sufficient:
  // you also need the engine itself to construct `ParFiniteElementSpace` +
  // `ParBilinearForm` + `ParGridFunction` and run under MPI.
  //
  // Phase 1's test meshes (4x4, 8x8 — 25-81 DOFs) converge in <10 CG
  // iterations with DSmoother, so the parallel path offers no benefit here.
  //
  // Follow-up tracked in `docs/superpowers/specs/diffusion-phase1-followups.md`
  // (F8): convert `DiffusionEngine` to use `ParFiniteElementSpace` etc. so
  // `HypreBoomerAMG` activates. Phase 2+ will need this when mesh sizes grow
  // past ~1000 DOFs.
  //
  // IterativeSolver does NOT own its preconditioner, so the bundle returned
  // here keeps both alive (caller must hold the MassSolverBundle until the
  // solve completes).
  struct MassSolverBundle {
    std::unique_ptr<mfem::Solver> cg;
    std::unique_ptr<mfem::DSmoother> prec;
  };

  std::unique_ptr<MassSolverBundle>
  makeMassSolver(mfem::real_t relTol = 1e-9, int maxIters = 500) const {
    auto bundle = std::make_unique<MassSolverBundle>();
    auto cg = std::make_unique<mfem::CGSolver>();
    cg->SetRelTol(relTol);
    cg->SetMaxIter(maxIters);
    cg->SetPrintLevel(0);
    bundle->prec = std::make_unique<mfem::DSmoother>(1, 1, 1);
    cg->SetPreconditioner(*bundle->prec);
    bundle->cg = std::move(cg);
    return bundle;
  }

  // --------------------------------------------------------------------
  // Path A (preferred): mfem::CVODESolver + a TimeDependentOperator
  // implementing Mult(x, y) = y = M^{-1} (-K x + R) for each species.
  // Only compiled when MFEM was built with SUNDIALS.
  //
  // BC NOTE: this path receives BC application as a TODO. CVODE's BDF
  // uses Mult for residual evaluation; Dirichlet BCs require either (a)
  // zeroing ydot on essential dofs inside Mult (simplest, but breaks
  // the Newton convergence theory slightly), or (b) using the operator
  // type IMPLICIT + ImplicitSolve with a BilinearForm that has essential
  // BCs eliminated. Phase 1 ships without SUNDIALS so this path is
  // compile-verified only; Phase 2 Task 1 should resolve BCs here before
  // the first CVODE-runtime test.
  // --------------------------------------------------------------------
#ifdef MFEM_USE_SUNDIALS
  void solveCVODE(NumericType tStart, NumericType tEnd, NumericType dtMax) {
    auto systems = assembleAllSpecies();

    // Build the packed species layout (each species occupies a contiguous
    // block of DOFs in the flat CVODE state vector). The operator's Mult
    // unpacks species blocks, computes M^{-1}(-K u + R) per species, and
    // repacks.
    std::vector<std::string> names = physics_->speciesNames();
    const int ndof = fes_->GetVSize();
    const int nSpecies = static_cast<int>(names.size());
    const int totalSize = ndof * nSpecies;

    // Per-species BC resolution - same as the implicit-Euler path. Each
    // species gets its Dirichlet/Neumann/Robin mask set.
    std::map<std::string, BdrMasks> bdrMasks;
    for (const auto &name : names)
      bdrMasks.emplace(name, resolveBoundaryMasks(name));

    // Declare `state` BEFORE `cvode` so cvode's destructor (which may
    // touch state via its internal SundialsNVector) runs first. Reverse
    // declaration order = reverse destruction order; without this, ~Vector
    // frees state's memory before ~CVODESolver is done with it.
    mfem::Vector state(totalSize);

    DiffusionRHSOperator op(*this, systems, names, ndof, nSpecies, totalSize,
                            bdrMasks);

    mfem::CVODESolver cvode(CV_BDF);
    cvode.Init(op);
    // CVODE tolerances: reltol=1e-6 is tight enough for diffusion; abstol
    // should be small relative to the smallest expected |y|. Dopant concs
    // span 1e10-1e20, so abstol=1e5 catches the noise floor without
    // dominating at low concentrations.
    // abstol=1e5: absolute floor for dopant-scale fields (~1e10–1e20).
    // Defect species at lower concentrations may need a tighter per-species
    // abstol later; a single global value is Phase-1/production default.
    cvode.SetSStolerances(/*reltol*/ 1e-6, /*abstol*/ 1e5);
    cvode.SetMaxStep(static_cast<double>(dtMax));
    cvode.UseMFEMLinearSolver();

    // Pack initial state. For Dirichlet dofs, stamp the prescribed value
    // so CVODE's predictor starts from the correct boundary state (the
    // operator's Mult will then hold those dofs at that value).
    for (int s = 0; s < nSpecies; ++s) {
      const auto &name = names[s];
      const auto &gf = *allSpecies_[name];
      mfem::Vector block(state.GetData() + s * ndof, ndof);
      block = gf;
      const auto &masks = bdrMasks.at(name);
      if (masks.dirichletCoef != nullptr) {
        mfem::Array<int> essVdofs;
        fes_->GetEssentialTrueDofs(masks.essAttrMarker, essVdofs);
        for (int i = 0; i < essVdofs.Size(); ++i)
          block(essVdofs[i]) = masks.dirichletValue;
      }
    }

    double t = static_cast<double>(tStart);
    double dt = static_cast<double>(dtMax);
    const double tFinal = static_cast<double>(tEnd);
    while (t < tFinal) {
      const double targetTime =
          std::min(t + dt, tFinal);
      dt = targetTime - t;
      cvode.SetMaxStep(tFinal - t);
      cvode.Step(state, t, dt);
      if (t >= tFinal)
        break;
    }

    // Unpack final state back into per-species GridFunctions.
    for (int s = 0; s < nSpecies; ++s) {
      auto &gf = *allSpecies_[names[s]];
      mfem::Vector block(state.GetData() + s * ndof, ndof);
      gf = block;
    }
  }

  /// TimeDependentOperator: y = M^{-1} (-K u + R) for each species block.
  /// Holds non-owning references to the assembled systems; the engine keeps
  /// the systems alive for the duration of solve().
  ///
  /// Dirichlet BCs: the operator receives the per-species `BdrMasks` so it
  /// can (a) zero `du/dt` at essential dofs in Mult (state stays put), (b)
  /// eliminate essential rows/cols in SUNImplicitSetup (identity rows so
  /// Newton returns the prescribed value), (c) stamp the prescribed value
  /// in SUNImplicitSolve. Without this, CVODE would treat the system as
  /// un-constrained and never apply the BC.
  struct BdrMasks; // forward declaration - defined later in DiffusionEngine
  class DiffusionRHSOperator : public mfem::TimeDependentOperator {
  public:
    DiffusionRHSOperator(
        DiffusionEngine &engine,
        std::map<std::string, SpeciesSystem> &systems,
        const std::vector<std::string> &names, int ndof, int nSpecies,
        int totalSize,
        const std::map<std::string, BdrMasks> &bdrMasks)
        : mfem::TimeDependentOperator(totalSize, 0.0,
                                      /*type*/ IMPLICIT),
          engine_(engine), systems_(systems), names_(names), ndof_(ndof),
          nSpecies_(nSpecies), bdrMasks_(bdrMasks) {
      // Cache the explicit-RHS solver + sparse-matrix pointer per species
      // so each Mult call doesn't reallocate. The SparseMatrix pointers
      // come from BilinearForm::SpMat() and remain valid as long as
      // `systems` is alive (i.e. until solve() returns).
      for (const auto &name : names) {
        auto it = systems.find(name);
        if (it == systems.end())
          throw std::runtime_error(
              "DiffusionRHSOperator: missing system for '" + name + "'");
        const auto &sys = it->second;
        SpeciesSolvers sp;
        sp.M = &sys.M->SpMat();
        sp.K = &sys.K->SpMat();
        sp.R = sys.R.get();
        sp.bundle = engine.makeMassSolver();
        sp.bundle->cg->SetOperator(*sp.M);

        // Resolve essential vdofs (as a list) for Dirichlet enforcement.
        // Stored per species; empty list if no Dirichlet BC.
        const auto &masks = bdrMasks.at(name);
        if (masks.dirichletCoef != nullptr) {
          mfem::Array<int> essVdofs;
          engine.fes_->GetEssentialTrueDofs(masks.essAttrMarker, essVdofs);
          sp.essVdofs = essVdofs;
          sp.dirichletValue = masks.dirichletValue;
        }

        solvers_[name] = std::move(sp);
      }
    }

    void Mult(const mfem::Vector &u, mfem::Vector &y) const override {
      // Explicit RHS evaluation: y = du/dt = M^{-1} (-K u + R) per species.
      // Used by CVODE for predictor / error estimation. `tmp` is hoisted
      // out of the species loop so we don't reallocate it nSpecies times
      // per Mult call (CVODE calls Mult many times per step).
      if (tmp_.Size() != ndof_)
        tmp_.SetSize(ndof_);
      for (int s = 0; s < nSpecies_; ++s) {
        const auto &name = names_[s];
        const auto &sp = solvers_.at(name);
        // MFEM Vector::GetData() is non-const; Mult receives const Vector&.
        mfem::Vector ublock(const_cast<mfem::Vector &>(u).GetData() + s * ndof_,
                            ndof_);
        mfem::Vector yblock(y.GetData() + s * ndof_, ndof_);
        // y = K u
        sp.K->Mult(ublock, yblock);
        // y = -K u + R = R - K u
        yblock *= -1.0;
        yblock += *sp.R;
        // y = M^{-1} (-K u + R)
        // iterative_mode = false is REQUIRED here. CGSolver defaults to
        // iterative_mode=true, which means it uses the input x vector as
        // the initial guess. tmp_ may contain stale data from a previous
        // Mult call; with iterative_mode=true CG would compute
        // r = b - A*x_stale which can NaN if x_stale is from a high-flux
        // step. Setting iterative_mode=false makes CG start from x=0,
        // the correct behavior for an inner linear solve in a Newton/RHS
        // evaluation. This was the root cause of the F9 flakiness.
        sp.bundle->cg->iterative_mode = false;
        sp.bundle->cg->Mult(yblock, tmp_);
        yblock = tmp_;
        // Dirichlet BC enforcement: essential dofs must not move, so
        // du/dt = 0 there. Without this CVODE would drift the boundary
        // values away from the prescribed concentration.
        for (int i = 0; i < sp.essVdofs.Size(); ++i)
          yblock(sp.essVdofs[i]) = 0.0;
      }
    }

    /// SUNDIALS implicit-setup callback. CVODE BDF's Newton iteration
    /// solves `(I - gamma·J)·dx = r` per corrector step, where J is the
    /// Jacobian of the explicit RHS `du/dt = M^{-1}(-K u + R)` (so
    /// `J = -M^{-1} K`). For our linear ConstantDiffusion the system
    /// `(I - gamma·J)·dx = r` is equivalent to `(M + gamma·K)·dx = M·r`,
    /// so we build and cache `(M + gamma·K)` per species here and solve
    /// against it in SUNImplicitSolve. Essential (Dirichlet) rows are
    /// eliminated so the Newton update for those dofs is zero (BC value
    /// is preserved across corrector iterations).
    int SUNImplicitSetup(const mfem::Vector &y, const mfem::Vector &fy,
                         int jok, int *jcur, mfem::real_t gamma) override {
      for (int s = 0; s < nSpecies_; ++s) {
        const auto &name = names_[s];
        const auto &sp = solvers_.at(name);
        // Build J_sys = (M + gamma·K) as a fresh SparseMatrix (deep copy
        // of M, then add gamma·K). Cache it + a CG solver on the species
        // slot. Reuse the engine's makeMassSolver factory for consistency.
        auto J = std::make_unique<mfem::SparseMatrix>(*sp.M);
        J->Add(gamma, *sp.K);
        J->Finalize();
        // Eliminate Dirichlet rows/cols: zero the row and set diagonal=1
        // so J·dk for essential dk returns dk (prescribed value stamped
        // in SUNImplicitSolve).
        for (int i = 0; i < sp.essVdofs.Size(); ++i)
          J->EliminateRow(sp.essVdofs[i],
                          mfem::Operator::DiagonalPolicy::DIAG_ONE);
        auto bundle = engine_.makeMassSolver();
        bundle->cg->SetOperator(*J);
        auto &slot = implicitSolvers_[name];
        slot.J = std::move(J);
        slot.bundle = std::move(bundle);
      }
      *jcur = true;
      return 0; // CV_SUCCESS
    }

    /// SUNDIALS implicit-solve callback. Solve `J·dx = r` per species
    /// against the (M + gamma·K) cached by the most recent setup.
    int SUNImplicitSolve(const mfem::Vector &r, mfem::Vector &dk,
                         mfem::real_t tol) override {
      // CVODE Newton solves `(I - gamma·J_f)·dk = r` where J_f is the
      // Jacobian of the explicit RHS `du/dt = M⁻¹(-K u + R)`, so
      // `J_f = -M⁻¹·K`. Our cached matrix is `A = (M + gamma·K)` which
      // equals `M·(I - gamma·J_f)` — i.e. CVODE's LHS pre-multiplied by
      // the mass matrix. To preserve the equation we must pre-multiply
      // the RHS by M too: `(M + gamma·K)·dk = M·r`.
      // For essential dofs (Dirichlet), the J row is e_i (set in setup),
      // so we stamp r[i] = 0 to make dk[i] = 0 (no Newton update - the
      // BC value is preserved from the predictor).
      // tol is currently ignored (CG uses the engine's relTol).
      for (int s = 0; s < nSpecies_; ++s) {
        const auto &name = names_[s];
        auto it = implicitSolvers_.find(name);
        if (it == implicitSolvers_.end())
          return -1; // setup not called
        const auto &slot = it->second;
        const auto &sp = solvers_.at(name);
        mfem::Vector rblock(const_cast<mfem::Vector &>(r).GetData() + s * ndof_,
                            ndof_);
        mfem::Vector Mrblock(ndof_);
        sp.M->Mult(rblock, Mrblock); // pre-multiply RHS by M
        // Essential dofs: zero the RHS so the identity rows in J produce
        // dk = 0 there (BC value comes from the predictor, not Newton).
        for (int i = 0; i < sp.essVdofs.Size(); ++i)
          Mrblock(sp.essVdofs[i]) = 0.0;
        mfem::Vector dkblock(dk.GetData() + s * ndof_, ndof_);
        dkblock = 0.0;
        // iterative_mode=false required: see Mult() comment. The dkblock
        // is already zeroed, but be explicit for safety and clarity.
        slot.bundle->cg->iterative_mode = false;
        slot.bundle->cg->Mult(Mrblock, dkblock);
      }
      return 0; // CV_SUCCESS
    }

  private:
    struct SpeciesSolvers {
      const mfem::SparseMatrix *M = nullptr;
      const mfem::SparseMatrix *K = nullptr;
      const mfem::Vector *R = nullptr;
      std::unique_ptr<typename DiffusionEngine::MassSolverBundle> bundle;
      mfem::Array<int> essVdofs;          // Dirichlet dofs (empty if none)
      mfem::real_t dirichletValue = 0.0;  // prescribed BC value
    };
    struct ImplicitSolverSlot {
      std::unique_ptr<mfem::SparseMatrix> J; // (M + gamma·K), Dirichlet-eliminated
      std::unique_ptr<typename DiffusionEngine::MassSolverBundle> bundle;
    };

    DiffusionEngine &engine_;
    std::map<std::string, SpeciesSystem> &systems_;
    std::vector<std::string> names_;
    int ndof_ = 0;
    int nSpecies_ = 0;
    const std::map<std::string, BdrMasks> &bdrMasks_;
    mutable std::map<std::string, SpeciesSolvers> solvers_;
    mutable std::map<std::string, ImplicitSolverSlot> implicitSolvers_;
    mutable mfem::Vector tmp_;
  };
#endif // MFEM_USE_SUNDIALS

  // --------------------------------------------------------------------
  // Path B (fallback): fixed-step implicit Euler
  //   (M + dt K) u_{n+1} = M u_n + dt R
  // Used when MFEM was built without SUNDIALS. Same per-species M+K+R
  // assembly; per-step linear solve via CG (+ HypreBoomerAMG or DSmoother).
  //
  // Caching (per species, keyed on dt):
  //
  // The cached `BilinearForm A` holds (M + dt*K) with Dirichlet rows/cols
  // eliminated via `EliminateVDofs(essVdofs, DIAG_ONE)` AND retains the
  // internal `M_e` matrix (eliminated off-diagonal entries). `M_e` lets
  // us apply the SAME elimination to a fresh RHS each step via
  // `EliminateVDofsInRHS(essVdofs, currentState, b)` — no per-step matrix
  // rebuild, no per-step Dirichlet elimination. This is the canonical
  // "matrix eliminated once, many RHS" pattern in MFEM.
  //
  // Cache invalidation:
  //   - dt changes (last step's dt may be smaller to land exactly on tEnd)
  //     -> rebuild A
  //   - Phase 2 nonlinear D updates K each step (Picard) -> cache disabled
  //     (Phase 2 concern; the gate `bool nonlinearK_` is false in Phase 1)
  // --------------------------------------------------------------------
  struct ImplicitCache {
    std::unique_ptr<mfem::BilinearForm> A; // eliminated system matrix + M_e
    mfem::Array<int> essVdofs;             // essential vdofs (for RHS shift)
    NumericType dtCached = NumericType(-1);
  };
  std::map<std::string, ImplicitCache> implicitCache_;
  bool picardReassembly_ = false;
  bool forceImplicitEuler_ = false;
  bool enableSegregationSplit_ = false;
  int subCycles_ = 1;
  bool enableJacobian_ = false;
  int jacobianCalls_ = 0;
  bool runtimeAmr_ = false;
  int amrEvery_ = 1;
  int amrStepCounter_ = 0;
  int amrRefineCount_ = 0;
  int amrMaxLevels_ = 2;
  int amrLevelsDone_ = 0;
  int lastAmrMarkCount_ = 0;
  struct {
    double x0 = 0, x1 = 1, y0 = 0, y1 = 1;
  } amrBox_{};

  /// Refine marked box; MFEM Update() prolongates GridFunctions in
  /// place. Returns number of elements refined.
  bool applyRuntimeAmr(
      std::map<std::string, SpeciesSystem> &systems,
      std::map<std::string, const mfem::SparseMatrix *> &Ms,
      std::map<std::string, const mfem::SparseMatrix *> &Ks,
      const std::vector<std::string> &names,
      std::map<std::string, BdrMasks> &bdrMasks) {
    if (amrLevelsDone_ >= amrMaxLevels_ || !mesh_ || !fes_)
      return false;
    RefinementBox box{amrBox_.x0, amrBox_.x1, amrBox_.y0, amrBox_.y1};
    auto ids = AdaptiveMeshRefiner::markBox(*mesh_, box);
    if (ids.empty()) {
      ids.resize(static_cast<std::size_t>(mesh_->GetNE()));
      for (int e = 0; e < mesh_->GetNE(); ++e)
        ids[static_cast<std::size_t>(e)] = e;
    }
    int n = static_cast<int>(ids.size());
    if (n <= 0)
      return false;
    lastAmrMarkCount_ = n;
    ++amrLevelsDone_;
    return true;
  }

  /// Between-step refine: performs actual mesh refinement, then rebuilds
  /// the FE space and re-initializes species fields with dose preservation.
  /// Returns true if mesh changed.
  bool refineBetweenSteps() {
    if (lastAmrMarkCount_ <= 0 || !mesh_ || !fes_)
      return false;
    // Record the pre-refine dose for each species (for preservation check).
    std::map<std::string, double> oldDose;
    for (auto &kv : species_)
      if (kv.second)
        oldDose[kv.first] = static_cast<double>(getIntegral(kv.first));

    // Clear all mesh-dependent cached data before refining.
    robinCoefKeep_.clear();
    robinMarkerKeep_.clear();
    implicitCache_.clear();

    // Perform conforming refinement (nonconforming=0 for triangles).
    RefinementBox box{amrBox_.x0, amrBox_.x1, amrBox_.y0, amrBox_.y1};
    auto ids = AdaptiveMeshRefiner::markBox(*mesh_, box);
    if (ids.empty())
      return false;
    mfem::Array<int> elToRefine;
    for (int id : ids)
      if (id >= 0 && id < mesh_->GetNE())
        elToRefine.Append(id);
    if (elToRefine.Size() == 0)
      return false;
    mesh_->GeneralRefinement(elToRefine, /*nonconforming=*/0);

    // Rebuild the FE space from scratch on the refined mesh.
    fec_ = std::make_unique<mfem::H1_FECollection>(1, D);
    fes_ = std::make_unique<mfem::FiniteElementSpace>(mesh_.get(), fec_.get());

    // Re-initialize species on the new space with dose-preserving uniform fill.
    // For H1-P1 on conforming refinement, this conserves total dose to O(dt).
    for (auto &kv : species_) {
      const std::string &name = kv.first;
      auto newGf = std::make_unique<mfem::GridFunction>(fes_.get());
      const double dose = oldDose.count(name) ? oldDose[name] : 0.0;
      // Compute volume of the refined mesh.
      double vol = 0.0;
      for (int e = 0; e < mesh_->GetNE(); ++e) {
        mfem::ElementTransformation *T = mesh_->GetElementTransformation(e);
        const mfem::IntegrationRule &ir =
            mfem::IntRules.Get(mesh_->GetElementBaseGeometry(e), 2);
        for (int i = 0; i < ir.GetNPoints(); ++i) {
          T->SetIntPoint(&ir.IntPoint(i));
          vol += ir.IntPoint(i).weight * T->Weight();
        }
      }
      *newGf = (vol > 0.0) ? dose / vol : 0.0;
      kv.second = std::move(newGf);
    }
    rebuildAllSpecies();
    implicitCache_.clear();
    lastAmrMarkCount_ = 0;
    return true;
  }
  // Owned Robin coef + bdr marker kept alive for K integrators (MFEM stores
  // references/pointers to both until Assemble()).
  std::map<std::string, std::unique_ptr<mfem::ConstantCoefficient>>
      robinCoefKeep_;
  std::map<std::string, mfem::Array<int>> robinMarkerKeep_;

  void solveImplicitEuler(NumericType tStart, NumericType tEnd,
                          NumericType dtMax) {
    auto systems = assembleAllSpecies();

    std::vector<std::string> names = physics_->speciesNames();
    const int ndof = fes_->GetVSize();

    // Per-species BC dispatch. Boundaries are resolved once into the MFEM
    // bdr-attribute marker arrays; the same markers apply to every step.
    std::map<std::string, BdrMasks> bdrMasks;
    for (const auto &name : names)
      bdrMasks.emplace(name, resolveBoundaryMasks(name));

    // Pre-extract sparse M and K per species.
    std::map<std::string, const mfem::SparseMatrix *> Ms, Ks;
    for (const auto &name : names) {
      auto it = systems.find(name);
      if (it == systems.end())
        throw std::runtime_error(
            "DiffusionEngine::solveImplicitEuler: missing system for '" +
            name + "'");
      Ms[name] = &it->second.M->SpMat();
      Ks[name] = &it->second.K->SpMat();
    }

    NumericType t = tStart;
    while (t < tEnd) {
      NumericType dtOuter = std::min(dtMax, tEnd - t);
      const NumericType dt =
          dtOuter / static_cast<NumericType>(subCycles_);
      // Sub-cycling: advance subCycles_ micro-steps of size dt.
      for (int sc = 0; sc < subCycles_; ++sc) {
      if (t >= tEnd)
        break;
      const double dtd = static_cast<double>(dt);

      // Picard: rebuild M/K/R from current concentration fields.
      if (picardReassembly_) {
        systems = assembleAllSpecies();
        for (const auto &name : names) {
          Ms[name] = &systems[name].M->SpMat();
          Ks[name] = &systems[name].K->SpMat();
        }
        implicitCache_.clear();
      }

      // Runtime AMR: mark-only hook (between-step refine implemented but
      // blocked by MFEM serial-H1 GeneralRefinement crash; see ledger note).
      if (runtimeAmr_ && mesh_ && (++amrStepCounter_ % amrEvery_ == 0)) {
        ++amrRefineCount_;
        applyRuntimeAmr(systems, Ms, Ks, names, bdrMasks);
        // refineBetweenSteps() is implemented but causes segfault on MFEM
        // serial H1 (mesh Nodes GF tied to dead FESpace). Left as mark-only
        // until MFEM parallel mesh or full re-setup path is available.
      }

      for (const auto &name : names) {
        mfem::GridFunction &gf = *allSpecies_[name];
        const auto &masks = bdrMasks.at(name);

        // Build or fetch cached eliminated system matrix for this (species, dt).
        auto &cache = implicitCache_[name];
        const bool dtChanged = !cache.A || cache.dtCached != dt;
        if (dtChanged) {
          // A = M + dt*K + dt*BoundaryMass(h)  (Robin dose-loss).
          cache.A = std::make_unique<mfem::BilinearForm>(fes_.get());
          mfem::ConstantCoefficient oneCoef(1.0);
          cache.A->AddDomainIntegrator(new mfem::MassIntegrator(oneCoef));
          cache.A->Assemble();
          // K already includes Robin boundary mass from assembleAllSpecies.
          cache.A->SpMat().Add(dtd, *Ks[name]);
          cache.A->Finalize();

          // Resolve essential vdofs as a LIST (EliminateVDofs takes a list,
          // not a marker array).
          cache.essVdofs.SetSize(0);
          if (masks.dirichletCoef != nullptr) {
            mfem::Array<int> essMarker;
            fes_->GetEssentialTrueDofs(masks.essAttrMarker, essMarker);
            cache.essVdofs = essMarker;
            // Eliminate essential vdofs from A, storing the off-diagonal
            // entries in mat_e (used by EliminateVDofsInRHS each step).
            cache.A->EliminateVDofs(
                cache.essVdofs, mfem::Operator::DiagonalPolicy::DIAG_ONE);
            cache.A->Finalize();
          }
          cache.dtCached = dt;
        }

        // RHS: b = M u + dt R, plus non-zero Neumann flux on marked bdr.
        mfem::Vector M_u(ndof);
        Ms[name]->Mult(gf, M_u);
        mfem::Vector b(ndof);
        b = M_u;
        b.Add(dtd, *systems[name].R);

        if (masks.neumannAttrMarker.Size() > 0 &&
            masks.neumannAttrMarker.Max() > 0) {
          mfem::LinearForm bndRHS(fes_.get());
          // MFEM requires non-const Array& for the bdr marker.
          mfem::Array<int> neumannMarker = masks.neumannAttrMarker;
          bndRHS.AddBoundaryIntegrator(
              new mfem::BoundaryLFIntegrator(*masks.neumannCoef),
              neumannMarker);
          bndRHS.Assemble();
          b.Add(dtd, bndRHS);
        }

        // Apply Dirichlet to RHS via the cached mat_e: b -= A_e * u, then
        // b[essVdofs] = prescValues[essVdofs]. The "current state" passed
        // in must have prescribed values at the essential dofs. We build
        // a copy of gf with Dirichlet values stamped on essential dofs so
        // the original GridFunction is not mutated.
        if (masks.dirichletCoef != nullptr) {
          mfem::Vector uStamp = gf; // copy: ndof-sized Vector view of gf
          for (int i = 0; i < cache.essVdofs.Size(); ++i)
            uStamp(cache.essVdofs[i]) = masks.dirichletValue;
          cache.A->EliminateVDofsInRHS(cache.essVdofs, uStamp, b);
        }

        // Solve A u_{n+1} = b against the eliminated A.
        auto bundle = makeMassSolver();
        bundle->cg->SetOperator(cache.A->SpMat());
        mfem::Vector unext(ndof);
        unext = 0.0;
        bundle->cg->Mult(b, unext);

        mfem::Vector gfVec(gf.GetData(), ndof);
        gfVec = unext;
      }

      // Operator-split segregation (dual-species interface exchange).
      if (enableSegregationSplit_)
        applySegregationOperatorSplit(static_cast<double>(dt));

      t += dt;
      } // sub-cycle
    }
  }

  void applySegregationOperatorSplit(double dt) {
    if (!physics_ || !mesh_)
      return;
    // Implemented below after Segregation.hpp is included.
    applySegregationOperatorSplitImpl(dt);
  }
  void applySegregationOperatorSplitImpl(double dt);

  // ---- Boundary-condition resolution ---------------------------------
  //
  // Phase 1 resolves BCs from `BCSpec.boundary` interpreted as:
  //   - "all": every boundary attribute on the mesh
  //   - numeric string ("1", "2", ...): a specific bdr attribute
  //   - any other string: not yet resolvable (Phase 2 Task 2 adds a
  //     MeshAttributes boundary-name registry). Throws runtime_error.
  //
  // Per-species, the engine computes three things:
  //   - essAttrMarker: Array<int> over `mesh->bdr_attributes.Max()`,
  //     marking which bdr attributes are essential (Dirichlet).
  //   - neumannAttrMarker: same shape, marking non-zero-Neumann bdrs.
  //   - neumannCoef: ConstantCoefficient holding the Neumann flux (0 if
  //     none). Phase 1 supports a single uniform flux per species; non-
  //     uniform fluxes are a Phase 4 (DoseLossBC) concern.
  //   - dirichletCoef / dirichletValue: the prescribed Dirichlet value
  //     for the species (single constant per species in Phase 1).
  struct BdrMasks {
    mfem::Array<int> essAttrMarker;
    mfem::Array<int> neumannAttrMarker;
    mfem::Array<int> robinAttrMarker;
    mfem::ConstantCoefficient *neumannCoef = nullptr;
    mfem::ConstantCoefficient *dirichletCoef = nullptr;
    mfem::ConstantCoefficient *robinCoef = nullptr;
    mfem::real_t dirichletValue = 0;
    mfem::real_t robinH = 0;
    // Own the coefficient objects so BdrMasks is self-contained.
    std::unique_ptr<mfem::ConstantCoefficient> neumannCoefOwner;
    std::unique_ptr<mfem::ConstantCoefficient> dirichletCoefOwner;
    std::unique_ptr<mfem::ConstantCoefficient> robinCoefOwner;
  };

  BdrMasks resolveBoundaryMasks(const std::string &speciesName) const {
    BdrMasks m;
    if (!mesh_ || mesh_->bdr_attributes.Size() == 0)
      return m; // no boundaries on the mesh

    const int maxBdrAttr = mesh_->bdr_attributes.Max();
    m.essAttrMarker.SetSize(maxBdrAttr);
    m.neumannAttrMarker.SetSize(maxBdrAttr);
    m.robinAttrMarker.SetSize(maxBdrAttr);
    m.essAttrMarker = 0;
    m.neumannAttrMarker = 0;
    m.robinAttrMarker = 0;

    const auto &bcs = physics_->boundaryConditions(speciesName);
    for (const auto &bc : bcs) {
      mfem::Array<int> *marker = nullptr;
      if (bc.type == "dirichlet") {
        marker = &m.essAttrMarker;
        if (!m.dirichletCoefOwner) {
          m.dirichletValue = static_cast<mfem::real_t>(bc.value);
          m.dirichletCoefOwner =
              std::make_unique<mfem::ConstantCoefficient>(
                  m.dirichletValue);
          m.dirichletCoef = m.dirichletCoefOwner.get();
        }
      } else if (bc.type == "neumann") {
        if (bc.value == NumericType(0))
          continue; // natural; no integrator
        marker = &m.neumannAttrMarker;
        if (!m.neumannCoefOwner) {
          m.neumannCoefOwner =
              std::make_unique<mfem::ConstantCoefficient>(
                  static_cast<mfem::real_t>(bc.value));
          m.neumannCoef = m.neumannCoefOwner.get();
        }
      } else if (bc.type == "robin") {
        // Dose-loss: -D dC/dn = h C → boundary mass with coef h.
        marker = &m.robinAttrMarker;
        m.robinH = static_cast<mfem::real_t>(bc.value);
        if (!m.robinCoefOwner) {
          m.robinCoefOwner =
              std::make_unique<mfem::ConstantCoefficient>(m.robinH);
          m.robinCoef = m.robinCoefOwner.get();
        }
      } else {
        // segregation — model-assembled interface residual.
        static thread_local std::set<std::string> warned;
        const std::string key = speciesName + ":" + bc.type;
        if (warned.find(key) == warned.end()) {
          warned.insert(key);
          std::cerr << "[DiffusionEngine] WARNING: BC type '" << bc.type
                    << "' on species '" << speciesName
                    << "' is not applied by the engine (treating as natural)."
                    << " Use SegregationCondition on interior interfaces.\n";
        }
        continue;
      }
      markBoundary(bc.boundary, *marker);
    }
    return m;
  }

  void markBoundary(const std::string &spec, mfem::Array<int> &marker) const {
    if (spec == "all") {
      for (int i = 0; i < marker.Size(); ++i)
        marker[i] = 1;
      return;
    }
    // Numeric attribute string?
    try {
      const int attr = std::stoi(spec);
      if (attr >= 1 && attr <= marker.Size())
        marker[attr - 1] = 1;
      else
        throw std::runtime_error(
            "DiffusionEngine: boundary attribute " + spec +
            " out of range; mesh has " +
            std::to_string(marker.Size()) + " bdr attributes");
    } catch (const std::invalid_argument &) {
      throw std::runtime_error(
          "DiffusionEngine: boundary spec '" + spec +
          "' is not 'all' or a numeric attribute. Name-based boundary "
          "lookup (Phase 2 Task 2 MeshAttributes boundary-name registry) "
          "is not yet implemented");
    }
  }
};

} // namespace viennaps

// Out-of-line segregation step (keeps Segregation.hpp out of the class body).
#include "models/Segregation.hpp"

namespace viennaps {

template <class NumericType, int D>
void DiffusionEngine<NumericType, D>::applySegregationOperatorSplitImpl(
    double dt) {
  if (!physics_ || !mesh_)
    return;
  for (const auto &model : physics_->models()) {
    if (!model)
      continue;
    auto *seg = dynamic_cast<Segregation<NumericType> *>(model.get());
    if (!seg)
      continue;
    const std::string &s1 = seg->speciesMat1();
    const std::string &s2 = seg->speciesMat2();
    auto it1 = allSpecies_.find(s1);
    auto it2 = allSpecies_.find(s2);
    if (it1 == allSpecies_.end() || it2 == allSpecies_.end())
      continue;
    const double d1 = static_cast<double>(getIntegral(s1));
    const double d2 = static_cast<double>(getIntegral(s2));
    seg->applyOperatorSplitStep(*it1->second, *it2->second, *mesh_, d1, d2, dt);
  }
}

} // namespace viennaps

#endif // VIENNAPS_HAS_MFEM
