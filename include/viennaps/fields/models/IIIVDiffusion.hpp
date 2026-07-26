#pragma once

/// IIIVDiffusion - GaAs/InP donor/acceptor diffusivity on Ga sublattice.
/// Implements ATHENA eq. 3-239 (donor, V_Ga mechanism) and 3-240 (acceptor,
/// I_Ga mechanism) with carrier-concentration-dependent D:
///   D_donor   = D_AV*(n/ni) + D_AV^2*(n/ni)^2        (ATHENA 3-239)
///   D_acceptor = D_AI*(p/ni) + D_AI^2*(p/ni)^2       (ATHENA 3-240)
/// Si/Se are concentration-independent (D_AV^2 ~ 0); Zn is strongly
/// concentration-dependent ((p/ni)^2 dominant). The Arrhenius T-dependence
/// is folded into D_AV/D_AV^2/D_AI/D_AI^2 via the material property DB.

#include "../DiffusionModel.hpp"
#include "../IntrinsicCarrier.hpp"
#include "../MaterialPropertySystem.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class IIIVDiffusion : public DiffusionModel<NumericType> {
public:
  enum class Mechanism { Donor, Acceptor };

  explicit IIIVDiffusion(std::string material = "GaAs",
                         std::string species = "Dopant")
      : material_(std::move(material)), species_(std::move(species)) {
    this->setName("IIIVDiffusion(" + material_ + "," + species_ + ")");
    // Pick donor vs acceptor mechanism from species name.
    // Donors (Si, Se, Ge) diffuse via V_Ga; acceptors (Be, Mg, Zn, C) via I_Ga.
    if (species_ == "Si" || species_ == "Se" || species_ == "Ge")
      mechanism_ = Mechanism::Donor;
    else
      mechanism_ = Mechanism::Acceptor;
    loadParametersFromDB();
  }

  void setMechanism(Mechanism m) {
    mechanism_ = m;
    loadParametersFromDB();
  }
  Mechanism mechanism() const { return mechanism_; }

  /// Load D_AV/D_AV^2/D_AI/D_AI^2 from MaterialPropertySystem (per material
  /// and species). Falls back to generic Arrhenius if keys are absent.
  void loadParametersFromDB() {
    MaterialPropertySystem<NumericType> mps;
    const std::string avKey = species_ + "_D_AV";
    const std::string av2Key = species_ + "_D_AV2";
    const std::string aiKey = species_ + "_D_AI";
    const std::string ai2Key = species_ + "_D_AI2";
    D_AV_ = mps.getProperty(material_, avKey, 0);
    D_AV2_ = mps.getProperty(material_, av2Key, 0);
    D_AI_ = mps.getProperty(material_, aiKey, 0);
    D_AI2_ = mps.getProperty(material_, ai2Key, 0);
    // Generic Arrhenius fallback (for species without DB entries).
    D0_ = mps.getProperty(material_, species_ + "_D", NumericType(1e-15));
    Ea_ = NumericType(2.0);
  }

  void setDiffusivity(NumericType D0, NumericType Ea) {
    D0_ = D0;
    Ea_ = Ea;
  }

  /// Carrier-dependent diffusivity per eq. 3-239/3-240.
  /// C is the local dopant concentration; n, p are electron/hole densities.
  /// For donors (n-type), n ≈ C; for acceptors (p-type), p ≈ C.
  NumericType getDiffusivity(NumericType C, NumericType T) const {
    const NumericType kB = static_cast<NumericType>(8.617333262145e-5);
    if (T <= 0)
      return D0_;
    IntrinsicCarrier<NumericType> ic;
    const NumericType ni = ic.ni(T, material_);
    if (ni <= NumericType(0))
      return D0_ * std::exp(-Ea_ / (kB * T)); // fallback
    // Majority carrier depends on dopant type:
    //   donor  -> n ≈ C (electrons), p = ni^2/n
    //   acceptor -> p ≈ C (holes),    n = ni^2/p
    NumericType n, p;
    if (mechanism_ == Mechanism::Donor) {
      n = std::max(C, ni);
      p = ni * ni / std::max(n, NumericType(1));
    } else {
      p = std::max(C, ni);
      n = ni * ni / std::max(p, NumericType(1));
    }
    const NumericType nRatio = n / ni;
    const NumericType pRatio = p / ni;
    if (mechanism_ == Mechanism::Donor) {
      // ATHENA 3-239: D = D_AV*(n/ni) + D_AV^2*(n/ni)^2
      return D_AV_ * nRatio + D_AV2_ * nRatio * nRatio;
    }
    // ATHENA 3-240: D = D_AI*(p/ni) + D_AI^2*(p/ni)^2
    return D_AI_ * pRatio + D_AI2_ * pRatio * pRatio;
  }

  const std::string &material() const { return material_; }

  /// Compound-sublattice I/V equilibrium: Ga and As sublattices have
  /// distinct interstitial and vacancy species (4 total). These are
  /// registered as engine species when sublattice tracking is enabled.
  void enableSublatticeDefects(bool on = true) {
    sublatticeDefects_ = on;
  }
  bool sublatticeDefectsEnabled() const { return sublatticeDefects_; }

  /// Ga-sublattice interstitial equilibrium (V_Ga mechanism for donors).
  NumericType C_I_Ga_eq(NumericType T) const {
    const NumericType kB = static_cast<NumericType>(8.617333262145e-5);
    return NumericType(1e22) * std::exp(-NumericType(3.0) / (kB * std::max(T, NumericType(1))));
  }
  /// Ga-sublattice vacancy equilibrium.
  NumericType C_V_Ga_eq(NumericType T) const {
    const NumericType kB = static_cast<NumericType>(8.617333262145e-5);
    return NumericType(1e22) * std::exp(-NumericType(2.5) / (kB * std::max(T, NumericType(1))));
  }
  /// As-sublattice interstitial equilibrium (I_Ga mechanism for acceptors).
  NumericType C_I_As_eq(NumericType T) const {
    const NumericType kB = static_cast<NumericType>(8.617333262145e-5);
    return NumericType(1e22) * std::exp(-NumericType(3.5) / (kB * std::max(T, NumericType(1))));
  }
  /// As-sublattice vacancy equilibrium.
  NumericType C_V_As_eq(NumericType T) const {
    const NumericType kB = static_cast<NumericType>(8.617333262145e-5);
    return NumericType(1e22) * std::exp(-NumericType(3.0) / (kB * std::max(T, NumericType(1))));
  }

  // Legacy proxies (kept for backward compatibility).
  NumericType C_I_eq(NumericType T) const { return C_I_Ga_eq(T); }
  NumericType C_V_eq(NumericType T) const { return C_V_Ga_eq(T); }

  int numSpecies() const override {
    return sublatticeDefects_ ? 5 : 1;
  }
  std::vector<std::string> speciesNames() const override {
    if (sublatticeDefects_)
      return {species_, "I_Ga", "V_Ga", "I_As", "V_As"};
    return {species_};
  }

#ifdef VIENNAPS_HAS_MFEM
  /// QP-local coefficient: reads the dopant GridFunction at each QP and
  /// computes D via eq. 3-239/3-240 (mirror of ChargedFermiDCoef pattern).
  class IIIVDCoef : public mfem::Coefficient {
  public:
    IIIVDCoef(const IIIVDiffusion *m, double T)
        : m_(m), T_(T), conc_(nullptr) {}
    void SetField(const mfem::GridFunction *c) { conc_ = c; }
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      double C = 0.0;
      if (conc_)
        C = std::max(0.0, conc_->GetValue(T, ip));
      return static_cast<double>(
          m_->getDiffusivity(static_cast<NumericType>(C),
                             static_cast<NumericType>(T_)));
    }

  private:
    const IIIVDiffusion *m_;
    double T_;
    const mfem::GridFunction *conc_;
  };

  void assembleStiffness(
      mfem::BilinearForm &K, const mfem::GridFunction &speciesGF,
      const std::map<std::string, mfem::GridFunction *> & /*allSpecies*/,
      const mfem::GridFunction * /*temp*/) const override {
    stiffCoef_ = std::make_unique<IIIVDCoef>(this, static_cast<double>(this->T_));
    stiffCoef_->SetField(&speciesGF);
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*stiffCoef_));
  }
  void assembleMass(mfem::BilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
  std::string material_;
  std::string species_;
  Mechanism mechanism_ = Mechanism::Donor;
  bool sublatticeDefects_ = false;
  // Eq. 3-239/3-240 prefactors (loaded from MaterialPropertySystem).
  NumericType D_AV_ = NumericType(1e-15);
  NumericType D_AV2_ = NumericType(0);
  NumericType D_AI_ = NumericType(1e-15);
  NumericType D_AI2_ = NumericType(0);
  // Generic Arrhenius fallback.
  NumericType D0_ = NumericType(1e-15);
  NumericType Ea_ = NumericType(2.0);
#ifdef VIENNAPS_HAS_MFEM
  mutable std::unique_ptr<IIIVDCoef> stiffCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif
};

} // namespace viennaps
