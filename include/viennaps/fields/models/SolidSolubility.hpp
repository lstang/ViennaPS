#pragma once

/// SolidSolubility — caps active dopant concentration; excess → cluster.
///
/// C_ss(T) = C_ss0 * exp(-Ea_ss / (kB * T))
/// For each dof (or element sample): C_excess = max(0, C - C_ss).
/// Reaction contribution (instant deactivation, rate 1/dt when dt provided;
/// otherwise unit rate so the residual encodes excess directly):
///   R[species]  -= C_excess / tau
///   R[cluster]  += C_excess / tau
///
/// This model owns two species names: the active dopant and the cluster.

#include "../DiffusionModel.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class SolidSolubility : public DiffusionModel<NumericType> {
public:
  SolidSolubility(const std::string &species = "Boron",
                  const std::string &cluster = "BoronCluster") {
    this->setName("SolidSolubility(" + species + "->" + cluster + ")");
    species_ = species;
    cluster_ = cluster;
  }

  void setSolidSolubility(NumericType Css0, NumericType Ea_eV = NumericType(0)) {
    Css0_ = Css0;
    Ea_ = Ea_eV;
  }

  /// Optional deactivation time scale. Default 1 so residual ≈ ±excess
  /// when the caller interprets a unit time step.
  void setTimeScale(NumericType tau) {
    tau_ = (tau > NumericType(0)) ? tau : NumericType(1);
  }

  NumericType solidSolubility(NumericType T) const {
    const NumericType kB = static_cast<NumericType>(8.617333262145e-5);
    if (T <= NumericType(0) || Ea_ == NumericType(0)) {
      return Css0_;
    }
    return Css0_ * std::exp(-Ea_ / (kB * T));
  }

  NumericType getCss0() const { return Css0_; }

  /// Apply one deactivation step to packed host arrays (test helper).
  /// active[i] and cluster[i] are updated in place for N dofs.
  void applyReactionStep(std::vector<NumericType> &active,
                         std::vector<NumericType> &cluster,
                         NumericType T) const {
    const NumericType Css = solidSolubility(T);
    const std::size_t n = std::min(active.size(), cluster.size());
    for (std::size_t i = 0; i < n; ++i) {
      const NumericType excess = std::max(NumericType(0), active[i] - Css);
      active[i] -= excess;
      cluster[i] += excess;
    }
  }

  int numSpecies() const override { return 2; }
  std::vector<std::string> speciesNames() const override {
    return {species_, cluster_};
  }
  std::vector<int> applicableAttributes() const override { return attrs_; }
  void setApplicableAttributes(std::vector<int> a) { attrs_ = std::move(a); }

  const std::string &activeSpecies() const { return species_; }
  const std::string &clusterSpecies() const { return cluster_; }

#ifdef VIENNAPS_HAS_MFEM
  /// Per-species residual. Engine invokes this once per owned species.
  /// Excess is always computed from the *active* concentration field
  /// (looked up in allSpecies); sign depends on which residual is open:
  ///   active  residual:  -excess/tau
  ///   cluster residual:  +excess/tau
  void assembleReaction(
      mfem::LinearForm &R, const mfem::GridFunction &speciesGF,
      const std::map<std::string, mfem::GridFunction *> &allSpecies,
      const mfem::GridFunction * /*temp*/) const override {
    const auto activeIt = allSpecies.find(species_);
    const mfem::GridFunction *activeGF =
        (activeIt != allSpecies.end() && activeIt->second)
            ? activeIt->second
            : &speciesGF;

    const NumericType Css = solidSolubility(this->T_);
    const double invTau = 1.0 / static_cast<double>(tau_);

    // Which residual is being assembled? Compare GridFunction identity.
    double scale = 0.0;
    if (&speciesGF == activeGF) {
      scale = -invTau;
    } else {
      // Cluster (or other owned) residual: deposit excess.
      scale = +invTau;
    }

    excessCoef_ = std::make_unique<ExcessCoef>(*activeGF,
                                               static_cast<double>(Css),
                                               scale);
    R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*excessCoef_));
  }
#endif

private:
#ifdef VIENNAPS_HAS_MFEM
  /// Coefficient: value = sign * max(0, C - Css) / tau at qp.
  class ExcessCoef : public mfem::Coefficient {
  public:
    ExcessCoef(const mfem::GridFunction &c, double Css, double scale)
        : conc_(&c), Css_(Css), scale_(scale) {}

    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      const double C = conc_->GetValue(T, ip);
      const double excess = std::max(0.0, C - Css_);
      return scale_ * excess;
    }

  private:
    const mfem::GridFunction *conc_;
    double Css_;
    double scale_;
  };
  mutable std::unique_ptr<ExcessCoef> excessCoef_;
#endif

  std::string species_;
  std::string cluster_;
  NumericType Css0_ = NumericType(1e20);
  NumericType Ea_ = NumericType(0);
  NumericType tau_ = NumericType(1);
  std::vector<int> attrs_;
};

} // namespace viennaps
