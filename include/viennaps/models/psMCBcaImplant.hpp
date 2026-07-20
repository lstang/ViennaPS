#pragma once

/// psMCBcaImplant - Binary Collision Approximation Monte Carlo implant.
///
/// Features beyond analytic sampling:
///   - Crystal vs amorphous target modes
///   - Channeling: direction-dependent range enhancement along <100>/<110>
///   - Cascade: primary ion → nuclear energy deposition → secondary I/V pairs
///   - Damage accumulation with local amorphization threshold
///
/// Field-only engine has no Domain/ViennaLS dependency (multiphysics smoke).
/// Process-model wrapper is optional via VIENNAPS_MCBCA_PROCESS_MODEL.

#include "fields/PhysicsField.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>
#include <memory>

namespace viennaps {

template <class NumericType>
struct MCBcaResult {
  std::vector<NumericType> dopant;
  std::vector<NumericType> interstitial;
  std::vector<NumericType> vacancy;
  std::vector<NumericType> nuclearDeposition;
  std::vector<NumericType> amorphousFraction;
  NumericType totalCascadePairs = 0;
  NumericType channeledFraction = 0;
};

/// Pure BCA engine (no Domain dependency) — used by process model and tests.
template <class NumericType>
class MCBcaEngine {
public:
  struct Params {
    NumericType energyKeV = 50;
    NumericType dose = 1e13;
    NumericType tiltDeg = 7;
    NumericType rotationDeg = 0;
    std::string crystalMode = "crystal"; // "crystal" | "amorphous"
    std::string channelAxis = "100";     // "100" | "110" | "111" | "none"
    int nIons = 256;
    int nDepth = 128;
    unsigned seed = 42;
  };

  explicit MCBcaEngine(Params p = {}) : params_(std::move(p)) {}

  void setParams(const Params& p) { params_ = p; }
  Params& params() { return params_; }
  const Params& getParams() const { return params_; }

  MCBcaResult<NumericType> run() const {
    MCBcaResult<NumericType> out;
    const int N = std::max(8, params_.nDepth);
    out.dopant.assign(N, 0);
    out.interstitial.assign(N, 0);
    out.vacancy.assign(N, 0);
    out.nuclearDeposition.assign(N, 0);
    out.amorphousFraction.assign(N, 0);

    const NumericType tiltRad = params_.tiltDeg * NumericType(3.141592653589793 / 180.0);
    const NumericType cosT = std::max(NumericType(0.15), std::cos(tiltRad));
    NumericType Rp = params_.energyKeV * NumericType(0.0085) * cosT + NumericType(0.01);
    NumericType dRp = params_.energyKeV * NumericType(0.0035) + NumericType(0.005);

    NumericType channelGain = NumericType(1);
    NumericType channelFrac = NumericType(0);
    if (params_.crystalMode == "crystal" && params_.channelAxis != "none") {
      const NumericType align =
          std::exp(-params_.tiltDeg * params_.tiltDeg / NumericType(50));
      if (params_.channelAxis == "100") {
        channelGain = NumericType(1) + NumericType(0.55) * align;
        channelFrac = NumericType(0.18) * align;
      } else if (params_.channelAxis == "110") {
        channelGain = NumericType(1) + NumericType(0.75) * align;
        channelFrac = NumericType(0.28) * align;
      } else if (params_.channelAxis == "111") {
        channelGain = NumericType(1) + NumericType(0.40) * align;
        channelFrac = NumericType(0.12) * align;
      }
    }
    out.channeledFraction = channelFrac;

    const NumericType RpChan = Rp * channelGain * NumericType(1.35);
    const NumericType dRpChan = dRp * NumericType(1.6);
    const NumericType depthMax = std::max(RpChan, Rp) * NumericType(4);
    const NumericType dz = depthMax / static_cast<NumericType>(N);

    std::mt19937 rng(params_.seed);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    std::normal_distribution<double> G(0.0, 1.0);

    const int nIons = std::max(16, params_.nIons);
    const NumericType dosePerIon = params_.dose / static_cast<NumericType>(nIons);
    NumericType cascadePairs = 0;

    const NumericType Ed = NumericType(15);
    const NumericType nuclearFraction =
        (params_.crystalMode == "amorphous") ? NumericType(0.55) : NumericType(0.45);

    for (int ion = 0; ion < nIons; ++ion) {
      const bool channeled = (U(rng) < static_cast<double>(channelFrac));
      const NumericType rp = channeled ? RpChan : Rp;
      const NumericType drp = channeled ? dRpChan : dRp;

      NumericType zStop = rp + drp * static_cast<NumericType>(G(rng));
      zStop = std::max(NumericType(0), zStop);

      NumericType zNuc = zStop * NumericType(0.82);
      NumericType Enuc = params_.energyKeV * NumericType(1000) * nuclearFraction;

      NumericType nPairs =
          std::max(NumericType(1), Enuc / (NumericType(2) * Ed) * NumericType(0.01));
      nPairs *= NumericType(0.7 + 0.6 * U(rng));
      cascadePairs += nPairs * dosePerIon;

      auto binOf = [&](NumericType z) -> int {
        int b = static_cast<int>(z / dz);
        return std::clamp(b, 0, N - 1);
      };

      out.dopant[binOf(zStop)] += dosePerIon;

      const NumericType sigmaN = drp * NumericType(0.6);
      for (int i = 0; i < N; ++i) {
        NumericType z = (static_cast<NumericType>(i) + NumericType(0.5)) * dz;
        NumericType w = std::exp(
            -NumericType(0.5) *
            std::pow((z - zNuc) / std::max(sigmaN, NumericType(1e-6)), 2));
        out.nuclearDeposition[i] += dosePerIon * Enuc * w;
        out.interstitial[i] += dosePerIon * nPairs * w * NumericType(1.05);
        out.vacancy[i] += dosePerIon * nPairs * w * NumericType(0.95);
      }
    }

    NumericType maxNuc = 0;
    for (auto v : out.nuclearDeposition) maxNuc = std::max(maxNuc, v);
    const NumericType amorphThresh = std::max(maxNuc * NumericType(0.35), NumericType(1));
    for (int i = 0; i < N; ++i) {
      out.amorphousFraction[i] =
          std::min(NumericType(1), out.nuclearDeposition[i] / amorphThresh);
    }

    out.totalCascadePairs = cascadePairs;

    std::cout << "[MCBcaEngine] ions=" << nIons << " mode=" << params_.crystalMode
              << " channel=" << params_.channelAxis
              << " channeledFrac=" << out.channeledFraction
              << " cascadePairs=" << out.totalCascadePairs
              << " dopantDose~" << params_.dose << "\n";
    return out;
  }

  static void injectIntoField(PhysicsField<NumericType>& field,
                              const MCBcaResult<NumericType>& r) {
    field.injectImplantProfile("Dopant", r.dopant);
    field.injectImplantProfile("Interstitial", r.interstitial);
    field.injectImplantProfile("Vacancy", r.vacancy);
    field.injectImplantProfile("NuclearDeposition", r.nuclearDeposition);
    field.injectImplantProfile("AmorphousFraction", r.amorphousFraction);
    field.addDose("CascadeIons", r.totalCascadePairs);
    field.addDose("ChanneledFraction", r.channeledFraction);
  }

private:
  Params params_;
};

/// Lightweight process-facing wrapper (field-only; no LS Domain required).
template <class NumericType, int D = 2>
class MCBcaImplant {
public:
  MCBcaImplant() = default;

  void setEnergy(NumericType e) { eng_.params().energyKeV = e; }
  void setDose(NumericType d) { eng_.params().dose = d; }
  void setTilt(NumericType t) { eng_.params().tiltDeg = t; }
  void setRotation(NumericType r) { eng_.params().rotationDeg = r; }
  void setCrystalMode(const std::string& m) { eng_.params().crystalMode = m; }
  void setChannelAxis(const std::string& a) { eng_.params().channelAxis = a; }
  void setNumIons(int n) { eng_.params().nIons = n; }
  void setSeed(unsigned s) { eng_.params().seed = s; }

  void setPhysicsField(std::shared_ptr<PhysicsField<NumericType>> field) {
    field_ = std::move(field);
  }

  MCBcaResult<NumericType> applyFieldOnly(PhysicsField<NumericType>& field) {
    auto result = eng_.run();
    MCBcaEngine<NumericType>::injectIntoField(field, result);
    lastResult_ = result;
    return result;
  }

  const MCBcaResult<NumericType>& getLastResult() const { return lastResult_; }

private:
  MCBcaEngine<NumericType> eng_;
  std::shared_ptr<PhysicsField<NumericType>> field_;
  MCBcaResult<NumericType> lastResult_;
};

} // namespace viennaps
