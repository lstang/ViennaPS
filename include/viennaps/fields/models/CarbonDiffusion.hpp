#pragma once

/// CarbonDiffusion — C traps interstitials (C + I → CI), suppressing TED.
/// FEM residual: r = kf * C * I on product/reactants (ImpurityCluster pattern).

#include "../DiffusionModel.hpp"

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class CarbonDiffusion : public DiffusionModel<NumericType> {
public:
  CarbonDiffusion(std::string carbon = "Carbon",
                  std::string interstitial = "Interstitial",
                  std::string complex = "CarbonInterstitial")
      : C_(std::move(carbon)), I_(std::move(interstitial)),
        CI_(std::move(complex)) {
    this->setName("CarbonDiffusion");
  }

  void setTrapRate(NumericType kf) { kf_ = kf; }
  void setDiffusivities(NumericType D_C, NumericType D_I) {
    D_C_ = D_C;
    D_I_ = D_I;
  }

  void applyTrapStep(std::vector<NumericType> &C, std::vector<NumericType> &I,
                     std::vector<NumericType> &CI, NumericType dt) const {
    const std::size_t n = std::min({C.size(), I.size(), CI.size()});
    for (std::size_t i = 0; i < n; ++i) {
      const NumericType r = kf_ * C[i] * I[i] * dt;
      C[i] = std::max(NumericType(0), C[i] - r);
      I[i] = std::max(NumericType(0), I[i] - r);
      CI[i] += r;
    }
  }

  int numSpecies() const override { return 3; }
  std::vector<std::string> speciesNames() const override {
    return {C_, I_, CI_};
  }

#ifdef VIENNAPS_HAS_MFEM
  void assembleStiffness(
      mfem::BilinearForm &K, const mfem::GridFunction &speciesGF,
      const std::map<std::string, mfem::GridFunction *> &allSpecies,
      const mfem::GridFunction * /*temp*/) const override {
    double D = 0.0;
    auto itC = allSpecies.find(C_);
    auto itI = allSpecies.find(I_);
    auto itCI = allSpecies.find(CI_);
    if (itC != allSpecies.end() && itC->second == &speciesGF)
      D = static_cast<double>(D_C_);
    else if (itI != allSpecies.end() && itI->second == &speciesGF)
      D = static_cast<double>(D_I_);
    else if (itCI != allSpecies.end() && itCI->second == &speciesGF)
      D = 0.0; // immobile complex default
    else
      return;
    stiffCoef_ = std::make_unique<mfem::ConstantCoefficient>(D);
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*stiffCoef_));
  }

  void assembleReaction(
      mfem::LinearForm &R, const mfem::GridFunction &speciesGF,
      const std::map<std::string, mfem::GridFunction *> &allSpecies,
      const mfem::GridFunction * /*temp*/) const override {
    auto itC = allSpecies.find(C_);
    auto itI = allSpecies.find(I_);
    auto itCI = allSpecies.find(CI_);
    if (itC == allSpecies.end() || itI == allSpecies.end() ||
        itCI == allSpecies.end() || !itC->second || !itI->second ||
        !itCI->second)
      return;

    double scale = 0.0;
    if (&speciesGF == itCI->second)
      scale = 1.0;
    else if (&speciesGF == itC->second || &speciesGF == itI->second)
      scale = -1.0;
    else
      return;

    coefs_.clear();
    coefs_.push_back(std::make_unique<TrapCoef>(
        *itC->second, *itI->second, static_cast<double>(kf_), scale));
    R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*coefs_.back()));
  }

  void assembleMass(mfem::BilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
#ifdef VIENNAPS_HAS_MFEM
  class TrapCoef : public mfem::Coefficient {
  public:
    TrapCoef(const mfem::GridFunction &C, const mfem::GridFunction &I,
             double kf, double scale)
        : C_(&C), I_(&I), kf_(kf), scale_(scale) {}
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      return scale_ * kf_ * C_->GetValue(T, ip) * I_->GetValue(T, ip);
    }

  private:
    const mfem::GridFunction *C_;
    const mfem::GridFunction *I_;
    double kf_, scale_;
  };
  mutable std::vector<std::unique_ptr<TrapCoef>> coefs_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> stiffCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif

  std::string C_, I_, CI_;
  NumericType kf_ = NumericType(1e-18);
  NumericType D_C_ = NumericType(1e-14);
  NumericType D_I_ = NumericType(1e-12);
};

} // namespace viennaps
