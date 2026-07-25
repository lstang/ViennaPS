#pragma once

/// SiGeDiffusion — Ge interdiffusion + B diffusivity modified by x_Ge.

#include "../BandgapModel.hpp"
#include "../DiffusionModel.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class SiGeDiffusion : public DiffusionModel<NumericType> {
public:
  SiGeDiffusion() { this->setName("SiGeDiffusion"); }

  void setGeDiffusivity(NumericType D0, NumericType Ea) {
    D0_Ge_ = D0;
    Ea_Ge_ = Ea;
  }

  void setBoronBaseDiffusivity(NumericType D0) { D0_B_ = D0; }

  NumericType geDiffusivity(NumericType T) const {
    const NumericType kB = static_cast<NumericType>(8.617333262145e-5);
    if (T <= 0)
      return D0_Ge_;
    return D0_Ge_ * std::exp(-Ea_Ge_ / (kB * T));
  }

  /// B diffusivity enhancement/suppression vs pure Si via bandgap/ni.
  NumericType boronDiffusivity(NumericType x_Ge, NumericType T) const {
    BandgapModel<NumericType> bg;
    // Smaller Eg → larger ni → typically larger extrinsic enhancement factor.
    const NumericType ratio = bg.niRatioToSi(x_Ge, T);
    return D0_B_ * ratio;
  }

  /// Host 1D intermixing step: smooth Ge profile with constant D.
  void applyIntermixStep(std::vector<NumericType> &xGe, NumericType T,
                         NumericType dx, NumericType dt) const {
    if (xGe.size() < 3 || dx <= 0)
      return;
    const NumericType D = geDiffusivity(T);
    const NumericType alpha = D * dt / (dx * dx);
    std::vector<NumericType> next = xGe;
    for (std::size_t i = 1; i + 1 < xGe.size(); ++i) {
      next[i] = xGe[i] + alpha * (xGe[i - 1] - NumericType(2) * xGe[i] + xGe[i + 1]);
      next[i] = std::min(NumericType(1), std::max(NumericType(0), next[i]));
    }
    xGe.swap(next);
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {"Germanium"};
  }

private:
  NumericType D0_Ge_ = NumericType(1e-3);
  NumericType Ea_Ge_ = NumericType(4.0);
  NumericType D0_B_ = NumericType(1e-13);
};

/// Ge-B pairing: reduces mobile B.
template <class NumericType>
class GeBPairing {
public:
  void setRates(NumericType kf, NumericType kr) {
    kf_ = kf;
    kr_ = kr;
  }

  void applyStep(std::vector<NumericType> &B, std::vector<NumericType> &Ge,
                 std::vector<NumericType> &pair, NumericType dt) const {
    const std::size_t n = std::min({B.size(), Ge.size(), pair.size()});
    for (std::size_t i = 0; i < n; ++i) {
      const NumericType form = kf_ * B[i] * Ge[i];
      const NumericType diss = kr_ * pair[i];
      const NumericType d = (form - diss) * dt;
      pair[i] = std::max(NumericType(0), pair[i] + d);
      B[i] = std::max(NumericType(0), B[i] - d);
      Ge[i] = std::max(NumericType(0), Ge[i] - d);
    }
  }

private:
  NumericType kf_ = NumericType(1e-20);
  NumericType kr_ = NumericType(1e-3);
};

/// Strain-modified diffusivity: D = D0 * exp(-alpha * strain / kT)
template <class NumericType>
class StrainDiffusionModifier {
public:
  void setAlpha(NumericType a) { alpha_ = a; }

  NumericType modifyD(NumericType D0, NumericType strain,
                      NumericType T) const {
    const NumericType kB = static_cast<NumericType>(8.617333262145e-5);
    if (T <= 0)
      return D0;
    return D0 * std::exp(-alpha_ * strain / (kB * T));
  }

private:
  NumericType alpha_ = NumericType(0.5); // eV per unit strain
};

} // namespace viennaps
