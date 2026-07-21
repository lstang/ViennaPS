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
///     the default preconditioner for diffusion. The engine uses
///     `mfem::HypreBoomerAMG` + CG when MFEM was built with MPI+hypre
///     (MFEM gates HypreBoomerAMG behind `MFEM_USE_MPI`); otherwise it falls
///     back to `DSmoother`-preconditioned CG so Phase 1 still runs on a
///     serial-MFEM build.
///   - MOOSE MFEMHypreBoomerAMG.h
///     (3rdparty/moose/framework/include/mfem/solvers/MFEMHypreBoomerAMG.h)
///     wraps `mfem::HypreBoomerAMG` and adds an optional `SetupLOR()` low-
///     order-refinement acceleration path. Phase 1 does not require LOR
///     (meshes are small); a Phase 2+ task may add it.
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

#include "DiffusionModel.hpp"
#include "DiffusionPhysics.hpp"
#include "MeshAttributes.hpp"

#include <algorithm>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef VIENNAPS_HAS_MFEM
#include <mfem.hpp>

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

    mesh_ = std::move(mesh);
    attrs_ = &attrs;
    fec_ = std::make_unique<mfem::H1_FECollection>(/*order*/ 1, /*dim*/ D);
    fes_ = std::make_unique<mfem::FiniteElementSpace>(
        mesh_.get(), fec_.get());
    species_.clear();
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

#ifdef MFEM_USE_SUNDIALS
    solveCVODE(tStart, tEnd, dtMax);
#else
    solveImplicitEuler(tStart, tEnd, dtMax);
#endif
  }

  /// Read-only access to the species' current concentration field.
  const mfem::GridFunction &getSolution(const std::string &name) const {
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
    mfem::LinearForm ones(fes_.get());
    ones.AddDomainIntegrator(new mfem::DomainLFIntegrator(oneCoef));
    ones.Assemble();
    const double dose = (*it->second) * ones;
    return static_cast<NumericType>(dose);
  }

private:
  // ---- State ----------------------------------------------------------
  std::unique_ptr<mfem::Mesh> mesh_;
  std::unique_ptr<mfem::FiniteElementCollection> fec_;
  std::unique_ptr<mfem::FiniteElementSpace> fes_;
  const MeshAttributes *attrs_ = nullptr;
  DiffusionPhysics<NumericType> *physics_ = nullptr;

  // Per-species GridFunctions on `fes_`. Same lifetime as the engine; the
  // TimeDependentOperator and ImplicitEuler paths read/write these through
  // the `allSpecies_` map built in `rebuildAllSpecies()`.
  std::map<std::string, std::unique_ptr<mfem::GridFunction>> species_;
  std::map<std::string, mfem::GridFunction *> allSpecies_;

  // ---- Helpers --------------------------------------------------------

  mfem::GridFunction *ensureSpecies(const std::string &name) {
    auto it = species_.find(name);
    if (it != species_.end())
      return it->second.get();
    auto gf = std::make_unique<mfem::GridFunction>(fes_.get());
    *gf = mfem::real_t(0);
    mfem::GridFunction *raw = gf.get();
    species_.emplace(name, std::move(gf));
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
    std::unique_ptr<mfem::BilinearForm> K;
    std::unique_ptr<mfem::BilinearForm> M;
    std::unique_ptr<mfem::LinearForm> R;
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

    // Propagate temperature + attributes from the physics into every model
    // before assembly. This is the canonical MOOSE PhysicsBase::initialize
    // step — without it, models keep whatever T was last set on them
    // (default 1273.15) and ignore physics_->temperature(). For
    // ConstantDiffusion with Ea=0 this is a no-op, but it matters for any
    // Arrhenius model (Fermi, Cdd, ...).
    if (attrs_) {
      const NumericType T = physics_->temperature();
      for (const auto &model : physics_->models())
        if (model)
          model->setup(*attrs_, T);
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
      sys.K = std::make_unique<mfem::BilinearForm>(fes_.get());
      sys.M = std::make_unique<mfem::BilinearForm>(fes_.get());
      sys.R = std::make_unique<mfem::LinearForm>(fes_.get());

      bool needTimeDeriv = true;
      for (const auto &model : physics_->models()) {
        if (!model)
          continue;
        if (!modelTargetsSpecies(*model, speciesName))
          continue;
        model->assembleStiffness(*sys.K, *speciesIt->second, allSpecies_,
                                 /*temp*/ nullptr);
        if (needTimeDeriv &&
            physics_->shouldCreateTimeDerivative(speciesName, *model)) {
          model->assembleMass(*sys.M);
          needTimeDeriv = false;
        }
        model->assembleReaction(*sys.R, *speciesIt->second, allSpecies_,
                                /*temp*/ nullptr);
      }

      sys.K->Assemble();
      sys.M->Assemble();
      sys.R->Assemble();

      // If no model contributed a mass term (e.g. a pure reaction species),
      // default to the identity mass so du/dt = -K u + R still has a well-
      // defined operator. Mirrors MOOSE's implicit fallback when no Kernel
      // claims the time derivative.
      if (needTimeDeriv) {
        mfem::ConstantCoefficient one(1.0);
        sys.M->AddDomainIntegrator(new mfem::MassIntegrator(one));
        sys.M->Assemble();
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
  // MFEM gates HypreBoomerAMG behind MFEM_USE_MPI (hypre.hpp:17). When
  // available, use it as the preconditioner for CG (the mass matrix M is
  // symmetric positive definite). When not available (serial MFEM build,
  // the current Phase 1 configuration), fall back to DSmoother-preconditioned
  // CG — same O(n) cost per iteration, slightly more iterations on large
  // meshes but correct for the small 2D tests in Phase 1.
  //
  // IterativeSolver does NOT own its preconditioner, so the bundle returned
  // here keeps both alive (caller must hold the MassSolverBundle until the
  // solve completes).
  struct MassSolverBundle {
    std::unique_ptr<mfem::Solver> cg;
#ifdef MFEM_USE_MPI
    std::unique_ptr<mfem::HypreBoomerAMG> prec;
#else
    std::unique_ptr<mfem::DSmoother> prec;
#endif
  };

  std::unique_ptr<MassSolverBundle>
  makeMassSolver(mfem::real_t relTol = 1e-9, int maxIters = 500) const {
    auto bundle = std::make_unique<MassSolverBundle>();
    auto cg = std::make_unique<mfem::CGSolver>();
    cg->SetRelTol(relTol);
    cg->SetMaxIter(maxIters);
    cg->SetPrintLevel(0);
#ifdef MFEM_USE_MPI
    bundle->prec = std::make_unique<mfem::HypreBoomerAMG>();
    bundle->prec->SetPrintLevel(0);
    cg->SetPreconditioner(*bundle->prec);
#else
    bundle->prec = std::make_unique<mfem::DSmoother>(1, 1, 1);
    cg->SetPreconditioner(*bundle->prec);
#endif
    bundle->cg = std::move(cg);
    return bundle;
  }

  // --------------------------------------------------------------------
  // Path A (preferred): mfem::CVODESolver + a TimeDependentOperator
  // implementing Mult(x, y) = y = M^{-1} (-K x + R) for each species.
  // Only compiled when MFEM was built with SUNDIALS.
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

    DiffusionRHSOperator op(*this, systems, names, ndof, nSpecies, totalSize);

    mfem::CVODESolver cvode(mfem::CV_BDF);
    cvode.Init(op);
    cvode.SetSStolerances(/*reltol*/ 1e-4, /*abstol*/ 1e-9 *
                                                       static_cast<double>(1e18));
    cvode.SetMaxStep(static_cast<double>(dtMax));
    cvode.UseMFEMLinearSolver();

    // Pack initial state.
    mfem::Vector state(totalSize);
    for (int s = 0; s < nSpecies; ++s) {
      const auto &gf = *allSpecies_[names[s]];
      mfem::Vector block(state.GetData() + s * ndof, ndof);
      block = gf;
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
  class DiffusionRHSOperator : public mfem::TimeDependentOperator {
  public:
    DiffusionRHSOperator(
        DiffusionEngine &engine,
        std::map<std::string, SpeciesSystem> &systems,
        const std::vector<std::string> &names, int ndof, int nSpecies,
        int totalSize)
        : mfem::TimeDependentOperator(totalSize, 0.0, /*type*/ EXPLICIT),
          engine_(engine), systems_(systems), names_(names), ndof_(ndof),
          nSpecies_(nSpecies) {
      // Cache a solver + sparse-matrix pointer per species so each Mult
      // call doesn't reallocate. The SparseMatrix pointers come from
      // BilinearForm::SpMat() and remain valid as long as `systems`
      // is alive (i.e. until solve() returns).
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
        solvers_[name] = std::move(sp);
      }
    }

    void Mult(const mfem::Vector &u, mfem::Vector &y) const override {
      // For each species block, compute Ku, subtract from -R, then apply
      // M^{-1}. `tmp` is hoisted out of the species loop so we don't
      // reallocate it nSpecies times per Mult call (CVODE calls Mult
      // many times per step).
      if (tmp_.Size() != ndof_)
        tmp_.SetSize(ndof_);
      for (int s = 0; s < nSpecies_; ++s) {
        const auto &name = names_[s];
        const auto &sp = solvers_.at(name);
        mfem::Vector ublock(const_cast<mfem::Vector &>(u).GetData() + s * ndof_,
                            ndof_);
        mfem::Vector yblock(y.GetData() + s * ndof_, ndof_);
        // y = K u
        sp.K->Mult(ublock, yblock);
        // y = -K u + R = R - K u
        yblock *= -1.0;
        yblock += *sp.R;
        // y = M^{-1} (-K u + R)
        sp.bundle->cg->Mult(yblock, tmp_);
        yblock = tmp_;
      }
    }

  private:
    struct SpeciesSolvers {
      const mfem::SparseMatrix *M = nullptr;
      const mfem::SparseMatrix *K = nullptr;
      const mfem::Vector *R = nullptr;
      std::unique_ptr<typename DiffusionEngine::MassSolverBundle> bundle;
    };

    DiffusionEngine &engine_;
    std::map<std::string, SpeciesSystem> &systems_;
    std::vector<std::string> names_;
    int ndof_ = 0;
    int nSpecies_ = 0;
    mutable std::map<std::string, SpeciesSolvers> solvers_;
    mutable mfem::Vector tmp_;
  };
#endif // MFEM_USE_SUNDIALS

  // --------------------------------------------------------------------
  // Path B (fallback): fixed-step implicit Euler
  //   (M + dt K) u_{n+1} = M u_n + dt R
  // Used when MFEM was built without SUNDIALS. Same per-species M+K+R
  // assembly; per-step linear solve via CG (+ HypreBoomerAMG or DSmoother).
  // --------------------------------------------------------------------
  void solveImplicitEuler(NumericType tStart, NumericType tEnd,
                          NumericType dtMax) {
    auto systems = assembleAllSpecies();

    std::vector<std::string> names = physics_->speciesNames();
    const int ndof = fes_->GetVSize();

    // No essential BC in Phase 1 (zero-flux Neumann is natural). ess_tdof
    // stays empty so FormLinearSystem passes the full system through.
    mfem::Array<int> essTdof;
    // Mass integrator coefficient kept alive as a member-equivalent local —
    // MassIntegrator takes Coefficient& (non-const), so the coefficient
    // object must outlive the BilinearForm assembly.
    mfem::ConstantCoefficient oneCoef(1.0);

    NumericType t = tStart;
    while (t < tEnd) {
      const NumericType dt = std::min(dtMax, tEnd - t);
      const double dtd = static_cast<double>(dt);

      for (const auto &name : names) {
        auto it = systems.find(name);
        if (it == systems.end())
          throw std::runtime_error(
              "DiffusionEngine::solveImplicitEuler: missing system for '" +
              name + "'");
        auto &sys = it->second;
        mfem::GridFunction &gf = *allSpecies_[name];

        // Build A = (M + dt K) by re-assembling M fresh and adding dt*K to
        // its entries. We can recover the assembled matrices from `sys`
        // (per-species M and K were already assembled once and stored).
        mfem::BilinearForm Aform(fes_.get());
        Aform.AddDomainIntegrator(new mfem::MassIntegrator(oneCoef));
        Aform.Assemble();
        mfem::SparseMatrix &Asparse = Aform.SpMat();
        const mfem::SparseMatrix &Ksparse = sys.K->SpMat();
        // A = M + dt * K  (in-place on the freshly assembled M)
        Asparse.Add(dtd, Ksparse);
        Asparse.Finalize();

        // Right-hand side: b = M u + dt R
        const mfem::SparseMatrix &Msparse = sys.M->SpMat();
        mfem::Vector M_u(ndof);
        Msparse.Mult(gf, M_u);
        mfem::Vector b(ndof);
        b = M_u;
        b.Add(dtd, *sys.R);

        // Build the preconditioner + CG solver. Reuse the mass solver
        // factory for consistency (HypreBoomerAMG if available; DSmoother
        // fallback otherwise). Preconditioner choice is robust because the
        // implicit-Euler system (M + dt K) has the same SPD structure as M.
        auto bundle = makeMassSolver();
        bundle->cg->SetOperator(Asparse);
        mfem::Vector unext(ndof);
        unext = 0.0;
        bundle->cg->Mult(b, unext);

        // Push back into the GridFunction for the next species / step.
        mfem::Vector gfVec(gf.GetData(), ndof);
        gfVec = unext;
      }

      t += dt;
    }
  }
};

} // namespace viennaps

#endif // VIENNAPS_HAS_MFEM
