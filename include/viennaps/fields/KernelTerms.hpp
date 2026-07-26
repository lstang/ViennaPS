#pragma once

/// Built-in KernelTerm implementations (MOOSE MatDiffusion / Reaction /
/// CoupledForce / BodyForce equivalents).

#include "KernelTerm.hpp"

#include <cmath>
#include <map>
#include <memory>
#include <string>

namespace viennaps {

#ifdef VIENNAPS_HAS_MFEM

/// ∇·(D ∇C) — MOOSE MatDiffusion.
class DiffusionTerm : public KernelTerm {
public:
  DiffusionTerm(std::string species, double D)
      : species_(std::move(species)), D_(D) {
    setName("DiffusionTerm(" + species_ + ")");
  }

  std::string targetSpecies() const override { return species_; }

  void assembleStiffness(
      mfem::ParBilinearForm &K,
      const std::map<std::string, mfem::ParGridFunction *> & /*species*/,
      const mfem::ParGridFunction * /*temp*/) const override {
    coef_ = std::make_unique<mfem::ConstantCoefficient>(D_);
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*coef_));
  }

  void assembleMass(mfem::ParBilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }

  double diffusivity() const { return D_; }

private:
  std::string species_;
  double D_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> coef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
};

/// Residual contribution λ * C (linear reaction / decay).
/// R += λ * Mass * C  via DomainLFIntegrator with GridFunctionCoefficient,
/// or equivalently stiffness contribution λ * Mass for implicit form.
/// Here: residual load R_i = λ * C_i projected (collocation-friendly).
class ReactionTerm : public KernelTerm {
public:
  ReactionTerm(std::string species, double lambda)
      : species_(std::move(species)), lambda_(lambda) {
    setName("ReactionTerm(" + species_ + ")");
  }

  std::string targetSpecies() const override { return species_; }

  void assembleResidual(
      mfem::ParLinearForm &R,
      const std::map<std::string, mfem::ParGridFunction *> &species,
      const mfem::ParGridFunction * /*temp*/) const override {
    auto it = species.find(species_);
    if (it == species.end() || !it->second)
      return;
    coef_ = std::make_unique<ScaledGFCoef>(*it->second, lambda_);
    R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*coef_));
  }

  /// Pure-math helper for unit tests: r = lambda * C
  static double eval(double lambda, double C) { return lambda * C; }

  double rate() const { return lambda_; }

private:
  class ScaledGFCoef : public mfem::Coefficient {
  public:
    ScaledGFCoef(const mfem::ParGridFunction &gf, double s)
        : gf_(&gf), s_(s) {}
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      return s_ * gf_->GetValue(T, ip);
    }

  private:
    const mfem::ParGridFunction *gf_;
    double s_;
  };

  std::string species_;
  double lambda_;
  mutable std::unique_ptr<ScaledGFCoef> coef_;
};

/// Residual contribution sigma * C_other (MOOSE CoupledForce).
/// For recombination-style coupling: target residual gets sigma * C_coupled.
class CoupledForceTerm : public KernelTerm {
public:
  CoupledForceTerm(std::string target, std::string coupled, double sigma)
      : target_(std::move(target)), coupled_(std::move(coupled)),
        sigma_(sigma) {
    setName("CoupledForceTerm(" + target_ + "," + coupled_ + ")");
  }

  std::string targetSpecies() const override { return target_; }

  void assembleResidual(
      mfem::ParLinearForm &R,
      const std::map<std::string, mfem::ParGridFunction *> &species,
      const mfem::ParGridFunction * /*temp*/) const override {
    auto it = species.find(coupled_);
    if (it == species.end() || !it->second)
      return;
    coef_ = std::make_unique<ScaledGFCoef>(*it->second, sigma_);
    R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*coef_));
  }

  static double eval(double sigma, double C_coupled) {
    return sigma * C_coupled;
  }

  double sigma() const { return sigma_; }

private:
  class ScaledGFCoef : public mfem::Coefficient {
  public:
    ScaledGFCoef(const mfem::ParGridFunction &gf, double s)
        : gf_(&gf), s_(s) {}
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      return s_ * gf_->GetValue(T, ip);
    }

  private:
    const mfem::ParGridFunction *gf_;
    double s_;
  };

  std::string target_;
  std::string coupled_;
  double sigma_;
  mutable std::unique_ptr<ScaledGFCoef> coef_;
};

/// Constant (or callable) body force source f.
class SourceTerm : public KernelTerm {
public:
  SourceTerm(std::string species, double f)
      : species_(std::move(species)), f_(f) {
    setName("SourceTerm(" + species_ + ")");
  }

  std::string targetSpecies() const override { return species_; }

  void assembleResidual(
      mfem::ParLinearForm &R,
      const std::map<std::string, mfem::ParGridFunction *> & /*species*/,
      const mfem::ParGridFunction * /*temp*/) const override {
    coef_ = std::make_unique<mfem::ConstantCoefficient>(f_);
    R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*coef_));
  }

  double source() const { return f_; }

private:
  std::string species_;
  double f_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> coef_;
};

/// Mass-matrix time derivative claim for a species (MOOSE TimeDerivative).
class TimeDerivativeTerm : public KernelTerm {
public:
  explicit TimeDerivativeTerm(std::string species)
      : species_(std::move(species)) {
    setName("TimeDerivativeTerm(" + species_ + ")");
  }

  std::string targetSpecies() const override { return species_; }

  void assembleMass(mfem::ParBilinearForm &M) const override {
    coef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*coef_));
  }

private:
  std::string species_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> coef_;
};

/// SUPG stabilization for advection-like pair diffusion (MOOSE LevelSetAdvectionSUPG).
/// tau = hmin / (2 ||v||); residual contribution tau (v·∇u)(v·∇test).
/// Implemented as ConvectionIntegrator with velocity = v scaled by tau when
/// ||v|| is estimated from a constant drift magnitude.
class SupgAdvectionTerm : public KernelTerm {
public:
  SupgAdvectionTerm(std::string species, double vx, double vy,
                    double hmin = 0.1)
      : species_(std::move(species)), vx_(vx), vy_(vy), hmin_(hmin) {
    setName("SupgAdvectionTerm(" + species_ + ")");
  }

  std::string targetSpecies() const override { return species_; }

  /// Stabilization parameter tau = h/(2|v|).
  double tau() const {
    const double vmag = std::sqrt(vx_ * vx_ + vy_ * vy_);
    if (vmag < 1e-30)
      return 0.0;
    return hmin_ / (2.0 * vmag);
  }

  void assembleStiffness(
      mfem::ParBilinearForm &K,
      const std::map<std::string, mfem::ParGridFunction *> & /*species*/,
      const mfem::ParGridFunction * /*temp*/) const override {
    const double t = tau();
    if (t <= 0.0)
      return;
    // Effective artificial diffusion ~ tau |v|^2 for streamline diffusion.
    const double v2 = vx_ * vx_ + vy_ * vy_;
    artDiff_ = std::make_unique<mfem::ConstantCoefficient>(t * v2);
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*artDiff_));
  }

  double driftX() const { return vx_; }
  double driftY() const { return vy_; }

private:
  std::string species_;
  double vx_, vy_, hmin_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> artDiff_;
};

#else // !VIENNAPS_HAS_MFEM

// Pure-math stubs so unit tests compile without MFEM.
class DiffusionTerm : public KernelTerm {
public:
  DiffusionTerm(std::string species, double D)
      : species_(std::move(species)), D_(D) {}
  std::string targetSpecies() const override { return species_; }
  double diffusivity() const { return D_; }

private:
  std::string species_;
  double D_;
};

class ReactionTerm : public KernelTerm {
public:
  ReactionTerm(std::string species, double lambda)
      : species_(std::move(species)), lambda_(lambda) {}
  std::string targetSpecies() const override { return species_; }
  static double eval(double lambda, double C) { return lambda * C; }
  double rate() const { return lambda_; }

private:
  std::string species_;
  double lambda_;
};

class CoupledForceTerm : public KernelTerm {
public:
  CoupledForceTerm(std::string target, std::string coupled, double sigma)
      : target_(std::move(target)), coupled_(std::move(coupled)),
        sigma_(sigma) {}
  std::string targetSpecies() const override { return target_; }
  static double eval(double sigma, double C) { return sigma * C; }
  double sigma() const { return sigma_; }

private:
  std::string target_, coupled_;
  double sigma_;
};

class SourceTerm : public KernelTerm {
public:
  SourceTerm(std::string species, double f)
      : species_(std::move(species)), f_(f) {}
  std::string targetSpecies() const override { return species_; }
  double source() const { return f_; }

private:
  std::string species_;
  double f_;
};

class TimeDerivativeTerm : public KernelTerm {
public:
  explicit TimeDerivativeTerm(std::string species)
      : species_(std::move(species)) {}
  std::string targetSpecies() const override { return species_; }

private:
  std::string species_;
};

class SupgAdvectionTerm : public KernelTerm {
public:
  SupgAdvectionTerm(std::string species, double vx, double vy,
                    double hmin = 0.1)
      : species_(std::move(species)), vx_(vx), vy_(vy), hmin_(hmin) {}
  std::string targetSpecies() const override { return species_; }
  double tau() const {
    const double vmag = std::sqrt(vx_ * vx_ + vy_ * vy_);
    return (vmag < 1e-30) ? 0.0 : hmin_ / (2.0 * vmag);
  }

private:
  std::string species_;
  double vx_, vy_, hmin_;
};

#endif // VIENNAPS_HAS_MFEM

} // namespace viennaps
