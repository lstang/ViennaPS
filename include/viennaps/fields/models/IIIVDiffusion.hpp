#pragma once

/// IIIVDiffusion — GaAs/InP species diffusivities on native sublattices.

#include "../DiffusionModel.hpp"

#include <cmath>
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

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {species_};
  }

private:
  std::string material_;
  std::string species_;
  NumericType D0_ = NumericType(1e-4);
  NumericType Ea_ = NumericType(2.5);
};

} // namespace viennaps
