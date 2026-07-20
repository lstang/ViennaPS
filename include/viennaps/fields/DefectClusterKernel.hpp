#pragma once

/// DefectClusterKernel - Point defect clustering models for high-fid diffusion.
/// Named mechanisms (selectable via setModel):
///   "recomb"  - I+V recombination only
///   "311"     - excess I aggregation into {311} clusters
///   "bic"     - boron-interstitial clustering (B + I -> BIC)
///   "loop"    - dislocation-loop growth from interstitial supersaturation
/// Default "combined" runs recomb + 311 (backward compatible).

#include "PhysicsKernel.hpp"
#include "PhysicsField.hpp"
#include "MaterialPropertySystem.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

namespace viennaps {

template <class NumericType>
class DefectClusterKernel : public PhysicsKernel<NumericType> {
public:
  explicit DefectClusterKernel(NumericType temperatureK = 1273.15,
                               std::string model = "combined")
      : T_(temperatureK), model_(std::move(model)) {
    this->setName("DefectCluster(" + model_ + ")");
  }

  void setTemperature(NumericType T_K) { T_ = T_K; }
  void setModel(const std::string& m) {
    model_ = m;
    this->setName("DefectCluster(" + model_ + ")");
  }
  const std::string& getModel() const { return model_; }

  void setup() override {
    std::cout << "[DefectClusterKernel] setup model=" << model_ << " at T=" << T_ << "\n";
  }

  void evolve(NumericType dt) override {
    if (!this->field_) return;

    this->field_->addSpecies("Interstitial");
    this->field_->addSpecies("Vacancy");

    if (model_ == "recomb") {
      runRecombination(dt);
    } else if (model_ == "311") {
      run311(dt);
    } else if (model_ == "bic") {
      runBIC(dt);
    } else if (model_ == "loop") {
      runLoop(dt);
    } else {
      // combined: recomb then 311
      runRecombination(dt);
      run311(dt);
    }

    std::cout << "[DefectClusterKernel] model=" << model_ << " evolved dt=" << dt
              << "  I=" << this->field_->getTotalDose("Interstitial")
              << "  V=" << this->field_->getTotalDose("Vacancy")
              << "  Cluster311=" << this->field_->getTotalDose("Cluster311")
              << "  BIC=" << this->field_->getTotalDose("BIC")
              << "  Loop=" << this->field_->getTotalDose("DislocationLoop") << "\n";
  }

  void addToRHS(NumericType /*t*/, const std::vector<NumericType>& y,
                std::vector<NumericType>& ydot, int myIndex) override {
    if (myIndex < 0 || myIndex >= static_cast<int>(ydot.size()) ||
        myIndex >= static_cast<int>(y.size()))
      return;
    NumericType rate = NumericType(0.05);
    if (model_ == "recomb") rate = NumericType(0.08);
    else if (model_ == "311") rate = NumericType(0.05);
    else if (model_ == "bic") rate = NumericType(0.04);
    else if (model_ == "loop") rate = NumericType(0.03);
    ydot[myIndex] += -rate * y[myIndex];
  }

  void addToFieldRHS(NumericType /*t*/, const std::vector<NumericType>& y,
                     std::vector<NumericType>& ydot) override {
    if (!this->field_) return;
    auto offI = this->field_->getSpeciesOffset("Interstitial");
    auto offV = this->field_->getSpeciesOffset("Vacancy");
    auto n = this->field_->getProfileSize();
    if (offI == static_cast<std::size_t>(-1) || offI + n > y.size()) return;

    NumericType rate = NumericType(0.05);
    if (model_ == "recomb") rate = NumericType(0.08);
    else if (model_ == "311") rate = NumericType(0.05);
    else if (model_ == "bic") rate = NumericType(0.04);
    else if (model_ == "loop") rate = NumericType(0.03);

    for (std::size_t i = 0; i < n; ++i) {
      NumericType Ii = std::max(NumericType(0), y[offI + i]);
      ydot[offI + i] += -rate * Ii;
      if (offV != static_cast<std::size_t>(-1) && offV + n <= y.size() &&
          (model_ == "recomb" || model_ == "combined")) {
        NumericType Vi = std::max(NumericType(0), y[offV + i]);
        NumericType r = std::min(Ii, Vi) * rate;
        ydot[offI + i] += -r;
        ydot[offV + i] += -r;
      }
    }
  }

private:
  void runRecombination(NumericType dt) {
    NumericType I = std::max(NumericType(0), this->field_->getTotalDose("Interstitial"));
    NumericType V = std::max(NumericType(0), this->field_->getTotalDose("Vacancy"));
    NumericType recomb = std::min(I, V) * NumericType(0.15) * std::min(dt, NumericType(1));
    recomb = std::min(recomb, std::min(I, V) * NumericType(0.9));
    if (recomb <= 0) return;
    if (I > 0) this->field_->scaleProfile("Interstitial", (I - recomb) / I);
    if (V > 0) this->field_->scaleProfile("Vacancy", (V - recomb) / V);
    this->field_->addSpecies("RecombinedIV");
    this->field_->injectImplantProfile("RecombinedIV",
                                       std::vector<NumericType>(1, recomb));
  }

  void run311(NumericType dt) {
    this->field_->addSpecies("Cluster311");
    NumericType I = this->field_->getTotalDose("Interstitial");
    NumericType V = this->field_->getTotalDose("Vacancy");
    NumericType excessI = I - V;
    if (excessI <= 0) return;
    NumericType cluster = excessI * NumericType(0.08) * dt;
    this->field_->injectImplantProfile("Cluster311",
                                       std::vector<NumericType>(1, cluster));
    if (I > 0)
      this->field_->scaleProfile("Interstitial",
                                 std::max(NumericType(0.5), (I - cluster) / I));
  }

  void runBIC(NumericType dt) {
    // Boron-interstitial clustering: consumes Dopant + Interstitial -> BIC
    this->field_->addSpecies("BIC");
    this->field_->addSpecies("Dopant");
    this->field_->addSpecies("Interstitial");
    NumericType B = this->field_->getTotalDose("Dopant");
    NumericType I = this->field_->getTotalDose("Interstitial");
    NumericType form = std::min(B, I) * NumericType(0.06) * dt;
    if (form <= 0) return;
    this->field_->injectImplantProfile("BIC", std::vector<NumericType>(1, form));
    if (B > 0) this->field_->scaleProfile("Dopant", (B - form) / B);
    if (I > 0) this->field_->scaleProfile("Interstitial", (I - form) / I);
  }

  void runLoop(NumericType dt) {
    // Dislocation loop growth from interstitial supersaturation
    this->field_->addSpecies("DislocationLoop");
    this->field_->addSpecies("Interstitial");
    NumericType I = this->field_->getTotalDose("Interstitial");
    if (I <= 0) return;
    NumericType loop = I * NumericType(0.04) * dt;
    this->field_->injectImplantProfile("DislocationLoop",
                                       std::vector<NumericType>(1, loop));
    this->field_->scaleProfile("Interstitial",
                               std::max(NumericType(0.6), (I - loop) / I));
  }

  NumericType T_ = 1273.15;
  std::string model_ = "combined";
};

} // namespace viennaps
