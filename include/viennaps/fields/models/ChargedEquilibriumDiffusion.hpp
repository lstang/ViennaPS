#pragma once

/// ChargedEquilibriumDiffusion — D = D0 * f_eq(T, n, p) at charge equilibrium.

#include "../DiffusionModel.hpp"
#include "../IntrinsicCarrier.hpp"

#include <algorithm>
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

  NumericType getDiffusivity(NumericType C, NumericType T) const {
    IntrinsicCarrier<NumericType> ic;
    const NumericType ni = ic.ni(T, "Si");
    const NumericType n = std::max(C, ni);
    return D0_ * (NumericType(1) + alpha_ * n / std::max(ni, NumericType(1)));
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {species_};
  }

private:
  std::string species_;
  NumericType D0_ = NumericType(1e-14);
  NumericType alpha_ = NumericType(1);
};

} // namespace viennaps
