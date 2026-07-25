#pragma once

/// MobileImpurity — general mobile impurity + optional ion-pairing.
/// FEM: QP D with optional coupled dopant field (same pattern as Copper).

#include "../DiffusionModel.hpp"
#include "../IntrinsicCarrier.hpp"

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class MobileImpurity : public DiffusionModel<NumericType> {
public:
  explicit MobileImpurity(std::string species = "Impurity")
      : species_(std::move(species)) {
    this->setName("MobileImpurity(" + species_ + ")");
  }

  void setD0(NumericType D0) { D0_ = D0; }
  void setIonPairing(NumericType beta) { beta_ = beta; }
  void setDopantConcentration(NumericType C) { C_dopant_ = C; }
  void setDopantSpecies(std::string s) { dopantSpecies_ = std::move(s); }

  NumericType getDiffusivity(NumericType C_dopant, NumericType T) const {
    if (beta_ == NumericType(0))
      return D0_;
    IntrinsicCarrier<NumericType> ic;
    const NumericType ni = ic.ni(T, "Si");
    return D0_ * (NumericType(1) +
                  beta_ * C_dopant / std::max(ni, NumericType(1)));
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {species_};
  }

#ifdef VIENNAPS_HAS_MFEM
  class ImpDCoef : public mfem::Coefficient {
  public:
    ImpDCoef(const MobileImpurity *m, double T, const mfem::GridFunction *dop,
             double Cconst)
        : m_(m), T_(T), dop_(dop), Cconst_(Cconst) {}
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      double C = Cconst_;
      if (dop_)
        C = dop_->GetValue(T, ip);
      return static_cast<double>(
          m_->getDiffusivity(static_cast<NumericType>(C),
                             static_cast<NumericType>(T_)));
    }

  private:
    const MobileImpurity *m_;
    double T_;
    const mfem::GridFunction *dop_;
    double Cconst_;
  };

  void assembleStiffness(
      mfem::BilinearForm &K, const mfem::GridFunction & /*speciesGF*/,
      const std::map<std::string, mfem::GridFunction *> &allSpecies,
      const mfem::GridFunction * /*temp*/) const override {
    const mfem::GridFunction *dop = nullptr;
    if (!dopantSpecies_.empty()) {
      auto it = allSpecies.find(dopantSpecies_);
      if (it != allSpecies.end())
        dop = it->second;
    }
    stiffCoef_ = std::make_unique<ImpDCoef>(
        this, static_cast<double>(this->T_), dop,
        static_cast<double>(C_dopant_));
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*stiffCoef_));
  }

  void assembleMass(mfem::BilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
  std::string species_;
  NumericType D0_ = NumericType(1e-10);
  NumericType beta_ = NumericType(0);
  NumericType C_dopant_ = NumericType(0);
  std::string dopantSpecies_;
#ifdef VIENNAPS_HAS_MFEM
  mutable std::unique_ptr<ImpDCoef> stiffCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif
};

} // namespace viennaps
