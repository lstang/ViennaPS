#pragma once

/// BandgapModel — E_g(x_Ge, T, strain) for SiGe and related materials.

#include <algorithm>
#include <cmath>

namespace viennaps {

template <class NumericType>
class BandgapModel {
public:
  /// SiGe bandgap with bowing: Eg = (1-x)Eg_Si + x*Eg_Ge - b*x*(1-x)
  NumericType Eg_SiGe(NumericType x_Ge, NumericType T = NumericType(300),
                      NumericType strain = NumericType(0)) const {
    const NumericType x = std::min(NumericType(1), std::max(NumericType(0), x_Ge));
    const NumericType EgSi = Eg_Si(T);
    const NumericType EgGe = Eg_Ge(T);
    const NumericType bow = NumericType(0.21); // eV bowing parameter
    NumericType Eg = (NumericType(1) - x) * EgSi + x * EgGe - bow * x * (NumericType(1) - x);
    // Strain: compressive (negative) increases Eg slightly (simple linear model)
    Eg += strainCoef_ * strain;
    return Eg;
  }

  NumericType Eg_Si(NumericType T = NumericType(300)) const {
    // Varshni-like mild T dependence around 1.12 eV
    return NumericType(1.17) -
           NumericType(4.73e-4) * T * T / (T + NumericType(636));
  }

  NumericType Eg_Ge(NumericType T = NumericType(300)) const {
    return NumericType(0.742) -
           NumericType(4.8e-4) * T * T / (T + NumericType(235));
  }

  void setStrainCoefficient(NumericType c) { strainCoef_ = c; }

  /// Intrinsic carrier concentration scaling: ni ∝ exp(-Eg/(2kT))
  NumericType niRatioToSi(NumericType x_Ge, NumericType T) const {
    const NumericType kB = static_cast<NumericType>(8.617333262145e-5);
    const NumericType dEg = Eg_SiGe(x_Ge, T) - Eg_Si(T);
    return std::exp(-dEg / (NumericType(2) * kB * T));
  }

private:
  NumericType strainCoef_ = NumericType(-1.0); // eV per unit strain
};

} // namespace viennaps
