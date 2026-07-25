#pragma once

/// IIIVDiffusion — GaAs/InP species diffusivities on native sublattices.

#include "../DiffusionModel.hpp"

#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class IIIVDiffusion : public DiffusionModel<NumericType> {
public:
  explicit IIIVDiffusion(std::string material = "GaAs",
                         std::string species = "Dopant")
      : material_(std::move(material)), species_(std::move(species)) {
    this->setName("IIIVDiffusion(" + material_ + "," + species_ + ")");
    // Rough defaults
    if (material_ == "GaAs") {
      D0_ = NumericType(1e-4);
      Ea_ = NumericType(2.5);
    } else if (material_ == "InP") {
      D0_ = NumericType(5e-5);
      Ea_ = NumericType(2.2);
    } else {
      D0_ = NumericType(1e-5);
      Ea_ = NumericType(2.0);
    }
  }

  void setDiffusivity(NumericType D0, NumericType Ea) {
    D0_ = D0;
    Ea_ = Ea;
  }

  NumericType getDiffusivity(NumericType T) const {
    const NumericType kB = static_cast<NumericType>(8.617333262145e-5);
    if (T <= 0)
      return D0_;
    return D0_ * std::exp(-Ea_ / (kB * T));
  }

  const std::string &material() const { return material_; }

  /// Compound-sublattice I/V equilibrium proxies (GaAs-like).
  NumericType C_I_eq(NumericType T) const {
    const NumericType kB = static_cast<NumericType>(8.617333262145e-5);
    // Higher formation energy than Si for III-V native defects.
    return NumericType(1e22) * std::exp(-NumericType(3.0) / (kB * std::max(T, NumericType(1))));
  }
  NumericType C_V_eq(NumericType T) const {
    const NumericType kB = static_cast<NumericType>(8.617333262145e-5);
    return NumericType(1e22) * std::exp(-NumericType(2.5) / (kB * std::max(T, NumericType(1))));
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {species_};
  }

#ifdef VIENNAPS_HAS_MFEM
  void assembleStiffness(
      mfem::BilinearForm &K, const mfem::GridFunction & /*speciesGF*/,
      const std::map<std::string, mfem::GridFunction *> & /*allSpecies*/,
      const mfem::GridFunction * /*temp*/) const override {
    stiffCoef_ = std::make_unique<mfem::ConstantCoefficient>(
        static_cast<double>(getDiffusivity(this->T_)));
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
  NumericType D0_ = NumericType(1e-4);
  NumericType Ea_ = NumericType(2.5);
#ifdef VIENNAPS_HAS_MFEM
  mutable std::unique_ptr<mfem::ConstantCoefficient> stiffCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif
};

} // namespace viennaps
