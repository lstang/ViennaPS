#pragma once

/// FermiDiffusionKernel - Fermi (extrinsic) diffusion model.
/// In high concentration regime, diffusivity increases with dopant concentration.
/// Simple model: D_eff = D_intrinsic * (1 + (C / C_ref))
/// (Real Fermi models use more terms with ni, charged fractions etc. from manuals.)
/// Uses MaterialPropertySystem for base D, and PhysicsField for local concentration.

#include "PhysicsKernel.hpp"
#include "PhysicsField.hpp"
#include "MaterialPropertySystem.hpp"

#include <iostream>

namespace viennaps {

template <class NumericType>
class FermiDiffusionKernel : public PhysicsKernel<NumericType> {
public:
  FermiDiffusionKernel(const std::string& species, NumericType temperatureK = 1273.15,
                       NumericType concRef = 1e18)
      : species_(species), T_(temperatureK), Cref_(concRef) {
    this->setName("FermiDiffusion(" + species + ")");
  }

  void setTemperature(NumericType T_K) { T_ = T_K; }

  void setup() override {
    if (this->material_) {
      D0_ = this->material_->getDiffusivity(species_, "Si", T_);
    } else {
      D0_ = NumericType(1e-14);
    }
    std::cout << "[FermiDiffusionKernel] " << species_ << " base D0=" << D0_
              << " at T=" << T_ << "K, Cref=" << Cref_ << "\n";
  }

  void evolve(NumericType dt) override {
    if (!this->field_) return;

    // Get approximate local conc (use total dose / volume proxy or profile avg)
    NumericType dose = this->field_->getTotalDose(species_);
    NumericType avgConc = dose / 1e-8; // very rough normalization for demo
    if (avgConc < 1) avgConc = 1;

    // Effective D
    NumericType Deff = D0_ * (NumericType(1) + (avgConc / Cref_));

    // Couple to stress (from ViscoelasticStressKernel)
    NumericType stress = this->field_->getTotalDose("HydrostaticStress");
    if (this->material_) {
      NumericType sFactor = this->material_->getProperty("Si", "stressFactor", T_);
      Deff *= (1.0 + (stress / 1e9) * 0.05 * sFactor); // stress in Pa -> GPa, small effect
    }
    NumericType totalBefore = this->field_->getTotalDose(species_);

    // Simple explicit relaxation scaled by Deff*dt (toy, like before)
    NumericType alpha = Deff * dt * NumericType(5e9);
    NumericType factor = NumericType(1) / (NumericType(1) + alpha * 0.05);
    if (factor < NumericType(0.6)) factor = NumericType(0.6);

    this->field_->scaleProfile(species_, factor);

    std::cout << "[FermiDiffusionKernel] " << species_ << " Deff=" << Deff
              << "  evolved dt=" << dt << "  total " << totalBefore << " -> "
              << this->field_->getTotalDose(species_) << "\n";
  }

  // Contribute to RHS: concentration-enhanced decay
  virtual void addToRHS(NumericType t, const std::vector<NumericType>& y, std::vector<NumericType>& ydot, int myIndex) override {
    if (myIndex >= 0 && myIndex < static_cast<int>(ydot.size()) && myIndex < static_cast<int>(y.size())) {
      NumericType yval = y[myIndex];
      NumericType Deff = D0_ * (1.0 + (yval / Cref_));
      NumericType alpha = 0.01 * (Deff / (D0_ > 0 ? D0_ : 1.0));  // scale
      ydot[myIndex] += -alpha * yval;
    }
  }

private:
  std::string species_;
  NumericType T_ = 1273.15;
  NumericType Cref_ = 1e18;
  NumericType D0_ = 0;
};

} // namespace viennaps
