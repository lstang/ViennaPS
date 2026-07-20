#pragma once

/// DefectClusterKernel - Basic point defect clustering (start of BIC/311/Loop models).
/// Simple: I + V -> recombination, or I aggregation to clusters (311 like).
/// Reduces mobile I/V, produces "Cluster" species.
/// This is the beginning of high-fid clustering per manuals.

#include "PhysicsKernel.hpp"
#include "PhysicsField.hpp"
#include "MaterialPropertySystem.hpp"

#include <iostream>

namespace viennaps {

template <class NumericType>
class DefectClusterKernel : public PhysicsKernel<NumericType> {
public:
  DefectClusterKernel(NumericType temperatureK = 1273.15)
      : T_(temperatureK) {
    this->setName("DefectCluster");
  }

  void setTemperature(NumericType T_K) { T_ = T_K; }

  void setup() override {
    std::cout << "[DefectClusterKernel] setup at T=" << T_ << " (recomb + simple clustering)\n";
  }

  void evolve(NumericType dt) override {
    if (!this->field_) return;

    this->field_->addSpecies("Interstitial");
    this->field_->addSpecies("Vacancy");
    this->field_->addSpecies("Cluster311"); // example cluster

    NumericType I = this->field_->getTotalDose("Interstitial");
    NumericType V = this->field_->getTotalDose("Vacancy");

    // Simple recombination I + V -> 0
    NumericType recomb = std::min(I, V) * 0.1 * dt; // rate
    if (recomb > 0) {
      // reduce by scaling
      if (I > 0) this->field_->scaleProfile("Interstitial", (I - recomb) / I );
      if (V > 0) this->field_->scaleProfile("Vacancy", (V - recomb) / V );
    }

    // Simple clustering: excess I -> Cluster311
    NumericType excessI = this->field_->getTotalDose("Interstitial") - this->field_->getTotalDose("Vacancy");
    if (excessI > 0) {
      NumericType cluster = excessI * 0.05 * dt;
      std::vector<NumericType> cl(1, cluster);
      this->field_->injectImplantProfile("Cluster311", cl);
      // remove from I
      this->field_->scaleProfile("Interstitial", 0.95);
    }

    std::cout << "[DefectClusterKernel] I/V recomb + clustering dt=" << dt << "\n";
  }

  // Contribute recombination/clustering rate
  virtual void addToRHS(NumericType t, const std::vector<NumericType>& y, std::vector<NumericType>& ydot, int myIndex) override {
    if (myIndex >= 0 && myIndex < static_cast<int>(ydot.size()) && myIndex < static_cast<int>(y.size())) {
      // Simple: negative rate for clustering/recomb
      ydot[myIndex] += -0.05 * y[myIndex];
    }
  }

private:
  NumericType T_ = 1273.15;
};

} // namespace viennaps
