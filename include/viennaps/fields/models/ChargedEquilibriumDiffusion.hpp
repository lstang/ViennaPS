#pragma once

/// ChargedEquilibriumDiffusion — D = D0 * f_eq(T, n, p) at charge equilibrium.
/// FEM: QP-local D from concentration field (FermiDCoef-style).

#include "../DiffusionModel.hpp"
#include "../IntrinsicCarrier.hpp"

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class ChargedEquilibriumDiffusion : public DiffusionModel<NumericType> {
public:
  explicit ChargedEquilibriumDiffusion(std::string species = "Boron")
      : species_(std::move(species)) {
    this->setName("ChargedEquilibriumDiffusion(" + species_ + ")");
  }

  void setD0(NumericType D0) { D0_ = D0; }
  void setAlpha(NumericType a) { alpha_ = a; }
  void setNi(NumericType ni) { ni_ = ni; }

  NumericType getDiffusivity(NumericType C, NumericType T) const {
    NumericType ni = ni_;
    if (ni <= NumericType(0)) {
      IntrinsicCarrier<NumericType> ic;
      ni = ic.ni(T, "Si");
    }
    const NumericType n = std::max(C, ni);
    return D0_ * (NumericType(1) + alpha_ * n / std::max(ni, NumericType(1)));
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {species_};
  }

#ifdef VIENNAPS_HAS_MFEM
  class EqDCoef : public mfem::Coefficient {
  public:
    EqDCoef(const ChargedEquilibriumDiffusion *m, double T)
        : m_(m), T_(T), conc_(nullptr) {}
    void SetField(const mfem::GridFunction *c) { conc_ = c; }
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      double C = 0.0;
      if (conc_)
        C = conc_->GetValue(T, ip);
      return static_cast<double>(
          m_->getDiffusivity(static_cast<NumericType>(C),
                             static_cast<NumericType>(T_)));
    }

  private:
    const ChargedEquilibriumDiffusion *m_;
    double T_;
    const mfem::GridFunction *conc_;
  };

  void assembleStiffness(
      mfem::BilinearForm &K, const mfem::GridFunction &speciesGF,
      const std::map<std::string, mfem::GridFunction *> & /*allSpecies*/,
      const mfem::GridFunction * /*temp*/) const override {
    stiffCoef_ =
        std::make_unique<EqDCoef>(this, static_cast<double>(this->T_));
    stiffCoef_->SetField(&speciesGF);
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*stiffCoef_));
  }

  void assembleMass(mfem::BilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
  std::string species_;
  NumericType D0_ = NumericType(1e-14);
  NumericType alpha_ = NumericType(1);
  NumericType ni_ = NumericType(0);
#ifdef VIENNAPS_HAS_MFEM
  mutable std::unique_ptr<EqDCoef> stiffCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif
};

} // namespace viennaps
