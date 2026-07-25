#pragma once

/// DiffusivityMaterial — concentration-dependent diffusivity coefficients
/// for MFEM assembly (MOOSE MatDiffusionBase + DerivativeMaterialInterface).
///
/// Provides:
///   - FermiDCoef   : D(C) = D_i * (1 + alpha * n/ni)  at quadrature points
///   - FermiDdCCoef : dD/dC for Newton Jacobian path (b1)
///
/// Coefficients must outlive the integrators that reference them (MFEM stores
/// Coefficient& and evaluates lazily during Assemble()). Models that own these
/// coefficients keep them as members.

#include <algorithm>
#include <cmath>

#ifdef VIENNAPS_HAS_MFEM
#include <mfem.hpp>
#endif

namespace viennaps {

/// Pure-math helpers shared by Fermi diffusion (no MFEM required).
template <class NumericType>
struct FermiDiffusivity {
  static NumericType D(NumericType C, NumericType D_i, NumericType alpha,
                       NumericType ni) {
    const NumericType n = std::max(C, ni);
    const NumericType niSafe = std::max(ni, NumericType(1));
    return D_i * (NumericType(1) + alpha * n / niSafe);
  }

  /// d/dC [D_i*(1 + alpha*n/ni)] with n = max(C, ni).
  /// Extrinsic (C > ni): dn/dC = 1 → dD/dC = D_i * alpha / ni.
  /// Intrinsic (C ≤ ni): dn/dC = 0 → dD/dC = 0.
  static NumericType dDdC(NumericType C, NumericType D_i, NumericType alpha,
                          NumericType ni) {
    if (C > ni) {
      const NumericType niSafe = std::max(ni, NumericType(1));
      return D_i * alpha / niSafe;
    }
    return NumericType(0);
  }
};

#ifdef VIENNAPS_HAS_MFEM
/// D(C) coefficient: D = D_i * (1 + alpha * n/ni), n = max(C, ni).
class FermiDCoef : public mfem::Coefficient {
public:
  FermiDCoef(double D_i, double alpha, double ni, double /*T*/ = 0.0)
      : D_i_(D_i), alpha_(alpha), ni_(ni), conc_(nullptr) {}

  void SetConcentrationField(const mfem::GridFunction *c) { conc_ = c; }

  double Eval(mfem::ElementTransformation &T,
              const mfem::IntegrationPoint &ip) override {
    double C = 0.0;
    if (conc_) {
      C = conc_->GetValue(T, ip);
    }
    const double n = std::max(C, ni_);
    const double niSafe = std::max(ni_, 1.0);
    return D_i_ * (1.0 + alpha_ * n / niSafe);
  }

private:
  double D_i_, alpha_, ni_;
  const mfem::GridFunction *conc_;
};

/// dD/dC coefficient for Jacobian chain-rule term (strategy b1).
class FermiDdCCoef : public mfem::Coefficient {
public:
  FermiDdCCoef(double D_i, double alpha, double ni)
      : D_i_(D_i), alpha_(alpha), ni_(ni), conc_(nullptr) {}

  void SetConcentrationField(const mfem::GridFunction *c) { conc_ = c; }

  double Eval(mfem::ElementTransformation &T,
              const mfem::IntegrationPoint &ip) override {
    double C = 0.0;
    if (conc_) {
      C = conc_->GetValue(T, ip);
    }
    if (C > ni_) {
      const double niSafe = std::max(ni_, 1.0);
      return D_i_ * alpha_ / niSafe;
    }
    return 0.0;
  }

private:
  double D_i_, alpha_, ni_;
  const mfem::GridFunction *conc_;
};
#endif // VIENNAPS_HAS_MFEM

} // namespace viennaps
