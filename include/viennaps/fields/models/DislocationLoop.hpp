#pragma once

/// DislocationLoop — loop growth from interstitial supersaturation.
/// dC_loop/dt = k * max(C_I/C_I_eq - 1, 0)^p
/// dC_I/dt   -= same (1:1 atomistic mapping)

#include "../DiffusionModel.hpp"
#include "../PointDefectEquilibrium.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class DislocationLoop : public DiffusionModel<NumericType> {
public:
  explicit DislocationLoop(std::string loop = "DislocationLoop",
                           std::string interstitial = "Interstitial")
      : loop_(std::move(loop)), I_(std::move(interstitial)) {
    this->setName("DislocationLoop");
  }

  void setGrowth(NumericType k, NumericType p = NumericType(1)) {
    k_ = k;
    p_ = p;
  }

  void setCIEq(NumericType Ceq) { Ceq_ = Ceq; }

  void applyReactionStep(std::vector<NumericType> &I,
                         std::vector<NumericType> &loop,
                         NumericType dt) const {
    const std::size_t n = std::min(I.size(), loop.size());
    const NumericType Ceq =
        (Ceq_ > NumericType(0))
            ? Ceq_
            : PointDefectEquilibrium<NumericType>{}.C_I_eq(this->T_, "Si");
    for (std::size_t i = 0; i < n; ++i) {
      const NumericType super =
          std::max(NumericType(0), I[i] / std::max(Ceq, NumericType(1)) -
                                       NumericType(1));
      const NumericType dL = k_ * std::pow(super, p_) * dt;
      loop[i] += dL;
      I[i] = std::max(NumericType(0), I[i] - dL);
    }
  }

  int numSpecies() const override { return 2; }
  std::vector<std::string> speciesNames() const override {
    return {loop_, I_};
  }

#ifdef VIENNAPS_HAS_MFEM
  void assembleReaction(
      mfem::LinearForm &R, const mfem::GridFunction &speciesGF,
      const std::map<std::string, mfem::GridFunction *> &allSpecies,
      const mfem::GridFunction * /*temp*/) const override {
    auto itI = allSpecies.find(I_);
    auto itL = allSpecies.find(loop_);
    if (itI == allSpecies.end() || itL == allSpecies.end() || !itI->second ||
        !itL->second)
      return;

    double scale = 0.0;
    if (&speciesGF == itL->second)
      scale = 1.0;
    else if (&speciesGF == itI->second)
      scale = -1.0;
    else
      return;

    double Ceq = static_cast<double>(Ceq_);
    if (Ceq <= 0.0) {
      Ceq = static_cast<double>(
          PointDefectEquilibrium<NumericType>{}.C_I_eq(this->T_, "Si"));
    }

    coefs_.clear();
    coefs_.push_back(std::make_unique<LoopRateCoef>(
        *itI->second, static_cast<double>(k_), static_cast<double>(p_), Ceq,
        scale));
    R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*coefs_.back()));
  }

  void assembleMass(mfem::BilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
#ifdef VIENNAPS_HAS_MFEM
  class LoopRateCoef : public mfem::Coefficient {
  public:
    LoopRateCoef(const mfem::GridFunction &I, double k, double p, double Ceq,
                 double scale)
        : I_(&I), k_(k), p_(p), Ceq_(std::max(Ceq, 1.0)), scale_(scale) {}
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      const double cI = std::max(0.0, I_->GetValue(T, ip));
      const double super = std::max(0.0, cI / Ceq_ - 1.0);
      return scale_ * k_ * std::pow(super, p_);
    }

  private:
    const mfem::GridFunction *I_;
    double k_, p_, Ceq_, scale_;
  };
  mutable std::vector<std::unique_ptr<LoopRateCoef>> coefs_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif

  std::string loop_, I_;
  NumericType k_ = NumericType(1e12);
  NumericType p_ = NumericType(1);
  NumericType Ceq_ = NumericType(0); // 0 => lookup
};

} // namespace viennaps
