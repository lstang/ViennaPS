#pragma once

/// Cluster311 — {311} interstitial cluster growth / dissociation.
///
/// dC_311/dt = k_f * C_I^n - k_r * C_311
/// dC_I/dt  -= n * (k_f * C_I^n - k_r * C_311)   (atom conservation)
///
/// FEM: residual product/power coefficients (MOOSE CoupledForce/Reaction).

#include "../DiffusionModel.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class Cluster311 : public DiffusionModel<NumericType> {
public:
  Cluster311(std::string cluster = "311",
             std::string interstitial = "Interstitial")
      : cluster_(std::move(cluster)), I_(std::move(interstitial)) {
    this->setName("Cluster311");
  }

  void setRates(NumericType kf, NumericType kr, int n = 2) {
    kf_ = kf;
    kr_ = kr;
    n_ = n;
  }

  void applyReactionStep(std::vector<NumericType> &I,
                         std::vector<NumericType> &C311,
                         NumericType dt) const {
    const std::size_t n = std::min(I.size(), C311.size());
    for (std::size_t i = 0; i < n; ++i) {
      const NumericType growth =
          kf_ * std::pow(std::max(I[i], NumericType(0)), n_);
      const NumericType dissoc = kr_ * C311[i];
      const NumericType dC = (growth - dissoc) * dt;
      C311[i] = std::max(NumericType(0), C311[i] + dC);
      I[i] = std::max(NumericType(0),
                      I[i] - static_cast<NumericType>(n_) * dC);
    }
  }

  int numSpecies() const override { return 2; }
  std::vector<std::string> speciesNames() const override {
    return {cluster_, I_};
  }

  NumericType kf() const { return kf_; }
  NumericType kr() const { return kr_; }

#ifdef VIENNAPS_HAS_MFEM
  void assembleReaction(
      mfem::ParLinearForm &R, const mfem::ParGridFunction &speciesGF,
      const std::map<std::string, mfem::ParGridFunction *> &allSpecies,
      const mfem::ParGridFunction * /*temp*/) const override {
    auto itI = allSpecies.find(I_);
    auto itC = allSpecies.find(cluster_);
    if (itI == allSpecies.end() || itC == allSpecies.end() || !itI->second ||
        !itC->second)
      return;

    // net rate r = kf*C_I^n - kr*C_311
    // cluster residual gets +r; interstitial residual gets -n*r
    double scale = 0.0;
    if (&speciesGF == itC->second)
      scale = 1.0;
    else if (&speciesGF == itI->second)
      scale = -static_cast<double>(n_);
    else
      return;

    coefs_.clear();
    coefs_.push_back(std::make_unique<ClusterRateCoef>(
        *itI->second, *itC->second, static_cast<double>(kf_),
        static_cast<double>(kr_), n_, scale));
    R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*coefs_.back()));
  }

  void assembleMass(mfem::ParBilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
#ifdef VIENNAPS_HAS_MFEM
  class ClusterRateCoef : public mfem::Coefficient {
  public:
    ClusterRateCoef(const mfem::ParGridFunction &I, const mfem::ParGridFunction &C311,
                    double kf, double kr, int n, double scale)
        : I_(&I), C311_(&C311), kf_(kf), kr_(kr), n_(n), scale_(scale) {}
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      const double cI = std::max(0.0, I_->GetValue(T, ip));
      const double c311 = C311_->GetValue(T, ip);
      const double r = kf_ * std::pow(cI, n_) - kr_ * c311;
      return scale_ * r;
    }

  private:
    const mfem::ParGridFunction *I_;
    const mfem::ParGridFunction *C311_;
    double kf_, kr_;
    int n_;
    double scale_;
  };
  mutable std::vector<std::unique_ptr<ClusterRateCoef>> coefs_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif

  std::string cluster_, I_;
  NumericType kf_ = NumericType(1e-20);
  NumericType kr_ = NumericType(1e-3);
  int n_ = 2;
};

} // namespace viennaps
