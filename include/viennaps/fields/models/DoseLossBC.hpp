#pragma once

/// DoseLossBC — surface evaporation Robin BC: -D dC/dn = h * C
/// Register with DiffusionPhysics::addRobinBC(species, boundary, h).

#include "../DiffusionModel.hpp"
#include "../DiffusionPhysics.hpp"

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
  void setBoundary(std::string b) { boundary_ = std::move(b); }

  NumericType h() const { return h_; }
  /// Robin coefficient for forms written as dC/dn = coef * C: coef = h/D.
  NumericType robinCoefficient() const {
    return (D_ > NumericType(0)) ? (h_ / D_) : NumericType(0);
  }

  /// Attach this BC to physics as an engine-applied Robin condition.
  void registerWith(DiffusionPhysics<NumericType> &physics) const {
    physics.addRobinBC(species_, boundary_, h_);
  }

  void applyLossStep(NumericType &C_surface, NumericType &dose,
                     NumericType area, NumericType dt) const {
    const NumericType flux = h_ * C_surface;
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
  std::string boundary_ = "all";
  NumericType h_ = NumericType(0);
  NumericType D_ = NumericType(1e-13);
};

template <class NumericType>
class DissociationLossBC : public DiffusionModel<NumericType> {
public:
  DissociationLossBC(std::string species = "Boron", NumericType Kd = 0)
      : species_(std::move(species)), Kd_(Kd) {
    this->setName("DissociationLossBC(" + species_ + ")");
  }

  void setKd(NumericType Kd) { Kd_ = Kd; }
  NumericType Kd() const { return Kd_; }

  /// Register as Robin with h = Kd (same weak form for first-order loss).
  void registerWith(DiffusionPhysics<NumericType> &physics,
                    const std::string &boundary = "all") const {
    physics.addRobinBC(species_, boundary, Kd_);
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {species_};
  }

private:
  std::string species_;
  NumericType Kd_ = NumericType(0);
};

} // namespace viennaps
