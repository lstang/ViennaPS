#pragma once

/// MobileImpurity — general mobile impurity + optional ion-pairing.

#include "../DiffusionModel.hpp"
#include "../IntrinsicCarrier.hpp"

#include <algorithm>
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

private:
  std::string species_;
  NumericType D0_ = NumericType(1e-10);
  NumericType beta_ = NumericType(0);
};

} // namespace viennaps
