#pragma once

/// Built-in KernelTerm implementations (MOOSE MatDiffusion / Reaction /
/// CoupledForce / BodyForce equivalents).

#include "KernelTerm.hpp"

#include <cmath>
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
      mfem::BilinearForm &K,
      const std::map<std::string, mfem::GridFunction *> & /*species*/,
      const mfem::GridFunction * /*temp*/) const override {
    coef_ = std::make_unique<mfem::ConstantCoefficient>(D_);
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*coef_));
  }

  void assembleMass(mfem::BilinearForm &M) const override {
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
      mfem::LinearForm &R,
      const std::map<std::string, mfem::GridFunction *> &species,
      const mfem::GridFunction * /*temp*/) const override {
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
    ScaledGFCoef(const mfem::GridFunction &gf, double s)
        : gf_(&gf), s_(s) {}
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      return s_ * gf_->GetValue(T, ip);
    }

  private:
    const mfem::GridFunction *gf_;
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
      mfem::LinearForm &R,
      const std::map<std::string, mfem::GridFunction *> &species,
      const mfem::GridFunction * /*temp*/) const override {
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
    ScaledGFCoef(const mfem::GridFunction &gf, double s)
        : gf_(&gf), s_(s) {}
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      return s_ * gf_->GetValue(T, ip);
    }

  private:
    const mfem::GridFunction *gf_;
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
      mfem::LinearForm &R,
      const std::map<std::string, mfem::GridFunction *> & /*species*/,
      const mfem::GridFunction * /*temp*/) const override {
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

  void assembleMass(mfem::BilinearForm &M) const override {
    coef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*coef_));
  }

private:
  std::string species_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> coef_;
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

#endif // VIENNAPS_HAS_MFEM

} // namespace viennaps
