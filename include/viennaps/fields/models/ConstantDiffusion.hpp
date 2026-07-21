#pragma once

#include "../DiffusionModel.hpp"

#include <cmath>
#include <string>
#include <vector>

namespace viennaps {

/// Concrete diffusion model with constant (concentration-independent)
/// diffusivity given by the Arrhenius relation
///   D = D0 * exp(-Ea / (kB * T))
/// where D0 is the pre-exponential factor, Ea is the activation energy in eV,
/// kB = 8.617333262145e-5 eV/K is the Boltzmann constant, and T is the
/// absolute temperature in Kelvin.
///
/// Design reference: MOOSE `MatDiffusion` Kernel — the canonical "register a
/// diffusion integrator with a scalar coefficient" pattern. The MFEM analog is
/// `mfem::DiffusionIntegrator` driven by a `mfem::ConstantCoefficient`, which
/// plays the role of MOOSE's `D` material property. See ADR-0001 for the
/// rationale behind the standalone-MFEM-engine architecture of which this
/// model is the first concrete instance.
template <class NumericType>
class ConstantDiffusion : public DiffusionModel<NumericType> {
public:
  explicit ConstantDiffusion(const std::string& species = "Dopant") {
    this->setName("ConstantDiffusion(" + species + ")");
    species_ = species;
  }

  void setSpecies(const std::string& s) { species_ = s; }
  void setDiffusivity(NumericType D0, NumericType Ea_eV) {
    D0_ = D0;
    Ea_ = Ea_eV;
  }

  NumericType getDiffusivity() const {
    const NumericType kB = static_cast<NumericType>(8.617333262145e-5);
    if (this->T_ <= 0)
      return D0_;
    return D0_ * std::exp(-Ea_ / (kB * this->T_));
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {species_};
  }
  std::vector<int> applicableAttributes() const override { return attrs_; }
  void setApplicableAttributes(std::vector<int> a) { attrs_ = std::move(a); }

#ifdef VIENNAPS_HAS_MFEM
  void assembleStiffness(mfem::BilinearForm& K, const mfem::GridFunction& speciesGF,
                         const std::map<std::string, mfem::GridFunction*>& allSpecies,
                         const mfem::GridFunction* temp) const override {
    double D = static_cast<double>(getDiffusivity());
    mfem::ConstantCoefficient Dcoef(D);
    // DiffusionIntegrator adds D * grad(phi_i) . grad(phi_j)
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(Dcoef));
  }

  void assembleMass(mfem::BilinearForm& M) const override {
    mfem::ConstantCoefficient one(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(one));
  }
#endif

private:
  std::string species_;
  NumericType D0_ = NumericType(1e-14);
  NumericType Ea_ = NumericType(3.0);
  std::vector<int> attrs_;
};

} // namespace viennaps
