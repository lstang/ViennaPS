#pragma once

/// psLithography - Simple aerial-image / threshold lithography for mask generation.
/// Produces a binary (or aerial) intensity field usable as implant/etch mask
/// without full OPC. Phase-2 parity placeholder for structure generation.

#include "fields/PhysicsField.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class LithographyModel {
public:
  LithographyModel() = default;

  void setWavelength(NumericType nm) { wavelength_ = nm; }
  void setNA(NumericType na) { NA_ = na; }
  void setDose(NumericType d) { dose_ = d; }
  void setThreshold(NumericType t) { threshold_ = t; }
  void setLineWidth(NumericType w) { lineWidth_ = w; }
  void setPitch(NumericType p) { pitch_ = p; }
  void setNumLines(int n) { nLines_ = n; }

  /// Generate aerial image along x in [-0.5,0.5], store mask in field.
  void apply(PhysicsField<NumericType>& field, int nX = 64) {
    nX = std::max(8, nX);
    std::vector<NumericType> aerial(nX, 0);
    std::vector<NumericType> mask(nX, 0);

    // Diffraction-limited blur sigma ~ k1 * lambda / NA (normalized)
    NumericType sigma =
        std::max(NumericType(0.01), NumericType(0.6) * wavelength_ / (NA_ * NumericType(1000)));
    // Map pitch/line into normalized domain width 1
    NumericType half = lineWidth_ / NumericType(2);

    for (int i = 0; i < nX; ++i) {
      NumericType x = (static_cast<NumericType>(i) + NumericType(0.5)) / static_cast<NumericType>(nX) -
                      NumericType(0.5);
      NumericType I = 0;
      for (int L = 0; L < nLines_; ++L) {
        NumericType cx = (static_cast<NumericType>(L) - static_cast<NumericType>(nLines_ - 1) / 2) *
                         pitch_;
        // Soft rectangle via erf edges
        NumericType e1 = std::erf(static_cast<double>((x - (cx - half)) / sigma));
        NumericType e2 = std::erf(static_cast<double>((x - (cx + half)) / sigma));
        I += NumericType(0.5) * (e1 - e2);
      }
      I = std::clamp(I * dose_, NumericType(0), NumericType(2));
      aerial[i] = I;
      mask[i] = (I >= threshold_) ? NumericType(1) : NumericType(0);
    }

    field.injectImplantProfile("AerialImage", aerial);
    field.injectImplantProfile("LithoMask", mask);

    NumericType open = 0;
    for (auto m : mask) open += m;
    openFraction_ = open / static_cast<NumericType>(nX);
    std::cout << "[LithographyModel] lambda=" << wavelength_ << "nm NA=" << NA_
              << " threshold=" << threshold_ << " openFraction=" << openFraction_
              << "\n";
  }

  NumericType getOpenFraction() const { return openFraction_; }

private:
  NumericType wavelength_ = 193; // nm
  NumericType NA_ = 0.93;
  NumericType dose_ = 1.0;
  NumericType threshold_ = 0.3;
  NumericType lineWidth_ = 0.08;
  NumericType pitch_ = 0.16;
  int nLines_ = 3;
  NumericType openFraction_ = 0;
};

} // namespace viennaps
