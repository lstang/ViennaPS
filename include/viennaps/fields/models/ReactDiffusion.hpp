#pragma once

/// ReactDiffusion — interstitial + vacancy bulk recombination.
///
/// dC_I/dt -= k * C_I * C_V
/// dC_V/dt -= k * C_I * C_V
/// plus constant-D diffusion for each species.

#include "../DiffusionModel.hpp"
#include "../KernelTerms.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class ReactDiffusion : public DiffusionModel<NumericType> {
public:
  ReactDiffusion(std::string interstitial = "Interstitial",
                  std::string vacancy = "Vacancy")
      : I_(std::move(interstitial)), V_(std::move(vacancy)) {
    this->setName("ReactDiffusion(" + I_ + "," + V_ + ")");
  }

  void setDiffusivities(NumericType D_I, NumericType D_V) {
    D_I_ = D_I;
    D_V_ = D_V;
  }

  void setRecombinationRate(NumericType k) { k_ = k; }

  NumericType recombinationRate() const { return k_; }

  /// One reaction step (host arrays) for unit tests: both decrease by k*I*V*dt.
  void applyReactionStep(std::vector<NumericType> &I,
                         std::vector<NumericType> &V, NumericType dt) const {
    const std::size_t n = std::min(I.size(), V.size());
    for (std::size_t i = 0; i < n; ++i) {
      const NumericType r = k_ * I[i] * V[i] * dt;
      I[i] = std::max(NumericType(0), I[i] - r);
      V[i] = std::max(NumericType(0), V[i] - r);
    }
  }

  int numSpecies() const override { return 2; }
  std::vector<std::string> speciesNames() const override {
    return {I_, V_};
  }

#ifdef VIENNAPS_HAS_MFEM
  void assembleStiffness(
      mfem::BilinearForm &K, const mfem::GridFunction &speciesGF,
      const std::map<std::string, mfem::GridFunction *> &allSpecies,
      const mfem::GridFunction * /*temp*/) const override {
    // Detect which species residual is open by pointer identity.
    double D = static_cast<double>(D_I_);
    auto itI = allSpecies.find(I_);
    auto itV = allSpecies.find(V_);
    if (itV != allSpecies.end() && itV->second == &speciesGF) {
      D = static_cast<double>(D_V_);
    } else if (itI != allSpecies.end() && itI->second == &speciesGF) {
      D = static_cast<double>(D_I_);
    }
    stiffCoef_ = std::make_unique<mfem::ConstantCoefficient>(D);
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*stiffCoef_));
  }

  void assembleReaction(
      mfem::LinearForm &R, const mfem::GridFunction &speciesGF,
      const std::map<std::string, mfem::GridFunction *> &allSpecies,
      const mfem::GridFunction * /*temp*/) const override {
    auto itI = allSpecies.find(I_);
    auto itV = allSpecies.find(V_);
    if (itI == allSpecies.end() || itV == allSpecies.end() || !itI->second ||
        !itV->second)
      return;
    // Residual contribution -k*C_I*C_V for BOTH species (recombination sink).
    recombCoef_ = std::make_unique<ProductCoef>(*itI->second, *itV->second,
                                                -static_cast<double>(k_));
    R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*recombCoef_));
    (void)speciesGF;
  }

  void assembleMass(mfem::BilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
#ifdef VIENNAPS_HAS_MFEM
  class ProductCoef : public mfem::Coefficient {
  public:
    ProductCoef(const mfem::GridFunction &a, const mfem::GridFunction &b,
                double scale)
        : a_(&a), b_(&b), scale_(scale) {}
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      return scale_ * a_->GetValue(T, ip) * b_->GetValue(T, ip);
    }

  private:
    const mfem::GridFunction *a_;
    const mfem::GridFunction *b_;
    double scale_;
  };
  mutable std::unique_ptr<mfem::ConstantCoefficient> stiffCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
  mutable std::unique_ptr<ProductCoef> recombCoef_;
#endif

  std::string I_, V_;
  NumericType D_I_ = NumericType(1e-10);
  NumericType D_V_ = NumericType(1e-12);
  NumericType k_ = NumericType(1e-15); // cm^3/s scale (test-tuned)
};

} // namespace viennaps
