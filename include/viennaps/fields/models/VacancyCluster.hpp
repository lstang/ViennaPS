#pragma once

/// VacancyCluster — vacancy cluster (VC) growth / dissociation.
/// dC_VC/dt = k_f * C_V^m - k_r * C_VC
/// dC_V/dt  -= m * (k_f * C_V^m - k_r * C_VC)

#include "../DiffusionModel.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class VacancyCluster : public DiffusionModel<NumericType> {
public:
  VacancyCluster(std::string cluster = "VC", std::string vacancy = "Vacancy")
      : cluster_(std::move(cluster)), V_(std::move(vacancy)) {
    this->setName("VacancyCluster");
  }

  void setRates(NumericType kf, NumericType kr, int m = 2) {
    kf_ = kf;
    kr_ = kr;
    m_ = m;
  }

  void applyReactionStep(std::vector<NumericType> &V,
                         std::vector<NumericType> &VC,
                         NumericType dt) const {
    const std::size_t n = std::min(V.size(), VC.size());
    for (std::size_t i = 0; i < n; ++i) {
      const NumericType growth =
          kf_ * std::pow(std::max(V[i], NumericType(0)), m_);
      const NumericType dissoc = kr_ * VC[i];
      const NumericType dC = (growth - dissoc) * dt;
      VC[i] = std::max(NumericType(0), VC[i] + dC);
      V[i] = std::max(NumericType(0),
                      V[i] - static_cast<NumericType>(m_) * dC);
    }
  }

  int numSpecies() const override { return 2; }
  std::vector<std::string> speciesNames() const override {
    return {cluster_, V_};
  }

#ifdef VIENNAPS_HAS_MFEM
  void assembleReaction(
      mfem::LinearForm &R, const mfem::GridFunction &speciesGF,
      const std::map<std::string, mfem::GridFunction *> &allSpecies,
      const mfem::GridFunction * /*temp*/) const override {
    auto itV = allSpecies.find(V_);
    auto itC = allSpecies.find(cluster_);
    if (itV == allSpecies.end() || itC == allSpecies.end() || !itV->second ||
        !itC->second)
      return;

    double scale = 0.0;
    if (&speciesGF == itC->second)
      scale = 1.0;
    else if (&speciesGF == itV->second)
      scale = -static_cast<double>(m_);
    else
      return;

    coefs_.clear();
    coefs_.push_back(std::make_unique<RateCoef>(
        *itV->second, *itC->second, static_cast<double>(kf_),
        static_cast<double>(kr_), m_, scale));
    R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*coefs_.back()));
  }

  void assembleMass(mfem::BilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
#ifdef VIENNAPS_HAS_MFEM
  class RateCoef : public mfem::Coefficient {
  public:
    RateCoef(const mfem::GridFunction &V, const mfem::GridFunction &VC,
             double kf, double kr, int m, double scale)
        : V_(&V), VC_(&VC), kf_(kf), kr_(kr), m_(m), scale_(scale) {}
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      const double cV = std::max(0.0, V_->GetValue(T, ip));
      const double cVC = VC_->GetValue(T, ip);
      return scale_ * (kf_ * std::pow(cV, m_) - kr_ * cVC);
    }

  private:
    const mfem::GridFunction *V_;
    const mfem::GridFunction *VC_;
    double kf_, kr_;
    int m_;
    double scale_;
  };
  mutable std::vector<std::unique_ptr<RateCoef>> coefs_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif

  std::string cluster_, V_;
  NumericType kf_ = NumericType(1e-20);
  NumericType kr_ = NumericType(1e-3);
  int m_ = 2;
};

} // namespace viennaps
