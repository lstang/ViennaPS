#pragma once

/// ImpurityCluster (BIC) — B + I ⇌ BIC
/// dC_BIC/dt = k_f * C_B * C_I - k_r * C_BIC
/// dC_B/dt  -= r; dC_I/dt -= r

#include "../DiffusionModel.hpp"

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class ImpurityCluster : public DiffusionModel<NumericType> {
public:
  ImpurityCluster(std::string bic = "BIC", std::string boron = "Boron",
                  std::string interstitial = "Interstitial")
      : bic_(std::move(bic)), B_(std::move(boron)),
        I_(std::move(interstitial)) {
    this->setName("ImpurityCluster");
  }

  void setRates(NumericType kf, NumericType kr) {
    kf_ = kf;
    kr_ = kr;
  }

  void applyReactionStep(std::vector<NumericType> &B,
                         std::vector<NumericType> &I,
                         std::vector<NumericType> &BIC,
                         NumericType dt) const {
    const std::size_t n = std::min({B.size(), I.size(), BIC.size()});
    for (std::size_t i = 0; i < n; ++i) {
      const NumericType form = kf_ * B[i] * I[i];
      const NumericType diss = kr_ * BIC[i];
      const NumericType dC = (form - diss) * dt;
      BIC[i] = std::max(NumericType(0), BIC[i] + dC);
      B[i] = std::max(NumericType(0), B[i] - dC);
      I[i] = std::max(NumericType(0), I[i] - dC);
    }
  }

  int numSpecies() const override { return 3; }
  std::vector<std::string> speciesNames() const override {
    return {bic_, B_, I_};
  }

#ifdef VIENNAPS_HAS_MFEM
  void assembleReaction(
      mfem::LinearForm &R, const mfem::GridFunction &speciesGF,
      const std::map<std::string, mfem::GridFunction *> &allSpecies,
      const mfem::GridFunction * /*temp*/) const override {
    auto itB = allSpecies.find(B_);
    auto itI = allSpecies.find(I_);
    auto itBIC = allSpecies.find(bic_);
    if (itB == allSpecies.end() || itI == allSpecies.end() ||
        itBIC == allSpecies.end() || !itB->second || !itI->second ||
        !itBIC->second)
      return;

    // r = kf*B*I - kr*BIC; BIC gets +r, B and I get -r
    double scale = 0.0;
    if (&speciesGF == itBIC->second)
      scale = 1.0;
    else if (&speciesGF == itB->second || &speciesGF == itI->second)
      scale = -1.0;
    else
      return;

    coefs_.clear();
    coefs_.push_back(std::make_unique<BicRateCoef>(
        *itB->second, *itI->second, *itBIC->second, static_cast<double>(kf_),
        static_cast<double>(kr_), scale));
    R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*coefs_.back()));
  }

  void assembleMass(mfem::BilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
#ifdef VIENNAPS_HAS_MFEM
  class BicRateCoef : public mfem::Coefficient {
  public:
    BicRateCoef(const mfem::GridFunction &B, const mfem::GridFunction &I,
                const mfem::GridFunction &BIC, double kf, double kr,
                double scale)
        : B_(&B), I_(&I), BIC_(&BIC), kf_(kf), kr_(kr), scale_(scale) {}
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      const double r = kf_ * B_->GetValue(T, ip) * I_->GetValue(T, ip) -
                       kr_ * BIC_->GetValue(T, ip);
      return scale_ * r;
    }

  private:
    const mfem::GridFunction *B_;
    const mfem::GridFunction *I_;
    const mfem::GridFunction *BIC_;
    double kf_, kr_, scale_;
  };
  mutable std::vector<std::unique_ptr<BicRateCoef>> coefs_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif

  std::string bic_, B_, I_;
  NumericType kf_ = NumericType(1e-18);
  NumericType kr_ = NumericType(1e-3);
};

} // namespace viennaps
