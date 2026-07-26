#pragma once

/// CarbonDiffusion — C traps interstitials with detailed balance.
/// Forward: r_f = kf * C * I
/// Reverse: r_r = kr * CI with default kr = kf * C_I_eq (SProcess 4.193)
/// so C_sI_eq = C_s * C_I / C*_I at equilibrium.

#include "../DiffusionModel.hpp"
#include "../PointDefectEquilibrium.hpp"

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
  void setReverseRate(NumericType kr) {
    kr_ = kr;
    krManual_ = true;
  }
  void setDiffusivities(NumericType D_C, NumericType D_I) {
    D_C_ = D_C;
    D_I_ = D_I;
  }

  NumericType forwardRate() const { return kf_; }
  NumericType reverseRate(NumericType T) const {
    if (krManual_)
      return kr_;
    // Detailed balance: kr = kf * C*_I so C_sI / (C_s C_I) = 1/C*_I at eq.
    const NumericType Ceq =
        PointDefectEquilibrium<NumericType>{}.C_I_eq(T, "Si");
    return kf_ * std::max(Ceq, NumericType(0));
  }

  void applyTrapStep(std::vector<NumericType> &C, std::vector<NumericType> &I,
                     std::vector<NumericType> &CI, NumericType dt) const {
    const NumericType kr = reverseRate(this->T_);
    const std::size_t n = std::min({C.size(), I.size(), CI.size()});
    for (std::size_t i = 0; i < n; ++i) {
      const NumericType form = kf_ * C[i] * I[i];
      const NumericType diss = kr * CI[i];
      const NumericType r = (form - diss) * dt;
      C[i] = std::max(NumericType(0), C[i] - r);
      I[i] = std::max(NumericType(0), I[i] - r);
      CI[i] = std::max(NumericType(0), CI[i] + r);
    }
  }

  int numSpecies() const override { return 3; }
  std::vector<std::string> speciesNames() const override {
    return {C_, I_, CI_};
  }

#ifdef VIENNAPS_HAS_MFEM
  void assembleStiffness(
      mfem::ParBilinearForm &K, const mfem::ParGridFunction &speciesGF,
      const std::map<std::string, mfem::ParGridFunction *> &allSpecies,
      const mfem::ParGridFunction * /*temp*/) const override {
    double D = 0.0;
    auto itC = allSpecies.find(C_);
    auto itI = allSpecies.find(I_);
    auto itCI = allSpecies.find(CI_);
    if (itC != allSpecies.end() && itC->second == &speciesGF)
      D = static_cast<double>(D_C_);
    else if (itI != allSpecies.end() && itI->second == &speciesGF)
      D = static_cast<double>(D_I_);
    else if (itCI != allSpecies.end() && itCI->second == &speciesGF)
      D = 0.0;
    else
      return;
    stiffCoef_ = std::make_unique<mfem::ConstantCoefficient>(D);
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*stiffCoef_));
  }

  void assembleReaction(
      mfem::ParLinearForm &R, const mfem::ParGridFunction &speciesGF,
      const std::map<std::string, mfem::ParGridFunction *> &allSpecies,
      const mfem::ParGridFunction * /*temp*/) const override {
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

    const double kr = static_cast<double>(reverseRate(this->T_));
    coefs_.clear();
    coefs_.push_back(std::make_unique<TrapCoef>(
        *itC->second, *itI->second, *itCI->second, static_cast<double>(kf_),
        kr, scale));
    R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*coefs_.back()));
  }

  void assembleMass(mfem::ParBilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
#ifdef VIENNAPS_HAS_MFEM
  class TrapCoef : public mfem::Coefficient {
  public:
    TrapCoef(const mfem::ParGridFunction &C, const mfem::ParGridFunction &I,
             const mfem::ParGridFunction &CI, double kf, double kr, double scale)
        : C_(&C), I_(&I), CI_(&CI), kf_(kf), kr_(kr), scale_(scale) {}
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      const double c = std::max(0.0, C_->GetValue(T, ip));
      const double i = std::max(0.0, I_->GetValue(T, ip));
      const double ci = std::max(0.0, CI_->GetValue(T, ip));
      return scale_ * (kf_ * c * i - kr_ * ci);
    }

  private:
    const mfem::ParGridFunction *C_;
    const mfem::ParGridFunction *I_;
    const mfem::ParGridFunction *CI_;
    double kf_, kr_, scale_;
  };
  mutable std::vector<std::unique_ptr<TrapCoef>> coefs_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> stiffCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif

  std::string C_, I_, CI_;
  NumericType kf_ = NumericType(1e-18);
  NumericType kr_ = NumericType(0);
  bool krManual_ = false;
  NumericType D_C_ = NumericType(1e-14);
  NumericType D_I_ = NumericType(1e-12);
};

} // namespace viennaps
