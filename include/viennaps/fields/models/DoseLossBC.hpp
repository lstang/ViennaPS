#pragma once

/// DoseLossBC — surface evaporation Robin BC: -D dC/dn = h * C
/// (MOOSE ADRobinBC with coef = h/D folded in).

#include "../DiffusionModel.hpp"

#include <algorithm>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class DoseLossBC : public DiffusionModel<NumericType> {
public:
  explicit DoseLossBC(std::string species = "Boron")
      : species_(std::move(species)) {
    this->setName("DoseLossBC(" + species_ + ")");
  }

  void setTransferCoefficient(NumericType h) { h_ = h; }
  void setDiffusivity(NumericType D) { D_ = D; }

  NumericType h() const { return h_; }
  /// Robin coefficient as used by MFEM: coef = h/D so dC/dn = coef * C.
  NumericType robinCoefficient() const {
    return (D_ > NumericType(0)) ? (h_ / D_) : NumericType(0);
  }

  /// 0D host mass-balance step: d(dose)/dt = -h * C_surface * area.
  void applyLossStep(NumericType &C_surface, NumericType &dose,
                     NumericType area, NumericType dt) const {
    const NumericType flux = h_ * C_surface; // out of domain
    const NumericType dDose = flux * area * dt;
    dose = std::max(NumericType(0), dose - dDose);
    C_surface = std::max(NumericType(0), C_surface - flux * dt);
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {species_};
  }

private:
  std::string species_;
  NumericType h_ = NumericType(0);
  NumericType D_ = NumericType(1e-13);
};

/// DissociationLossBC — species lost at surface transforming to another
/// (MOOSE DissociationFluxBC): residual ~ -Kd * C.
template <class NumericType>
class DissociationLossBC : public DiffusionModel<NumericType> {
public:
  DissociationLossBC(std::string species = "Boron", NumericType Kd = 0)
      : species_(std::move(species)), Kd_(Kd) {
    this->setName("DissociationLossBC(" + species_ + ")");
  }

  void setKd(NumericType Kd) { Kd_ = Kd; }
  NumericType Kd() const { return Kd_; }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {species_};
  }

private:
  std::string species_;
  NumericType Kd_ = NumericType(0);
};

} // namespace viennaps
