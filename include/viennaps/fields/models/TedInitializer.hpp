#pragma once

/// TedInitializer — project implant damage onto continuum fields with
/// integral-preserving rescaling (MOOSE IntegralPreservingFunctionIC).

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

#ifdef VIENNAPS_HAS_MFEM
#include <mfem.hpp>
#endif

namespace viennaps {

template <class NumericType>
class TedInitializer {
public:
  struct DamageProfile {
    std::vector<NumericType> values;
  };

  static void projectIntegralPreserving(const DamageProfile &source,
                                        NumericType dxSource,
                                        std::vector<NumericType> &target,
                                        NumericType dxTarget) {
    if (source.values.empty() || target.empty())
      return;

    const std::size_t ns = source.values.size();
    const std::size_t nt = target.size();
    for (std::size_t i = 0; i < nt; ++i) {
      const double xi =
          (nt == 1) ? 0.0
                    : static_cast<double>(i) / static_cast<double>(nt - 1);
      const double sj = xi * static_cast<double>(ns - 1);
      const std::size_t j0 = static_cast<std::size_t>(sj);
      const std::size_t j1 = std::min(j0 + 1, ns - 1);
      const double f = sj - static_cast<double>(j0);
      target[i] = static_cast<NumericType>(
          (1.0 - f) * static_cast<double>(source.values[j0]) +
          f * static_cast<double>(source.values[j1]));
    }

    const NumericType doseSrc =
        std::accumulate(source.values.begin(), source.values.end(),
                        NumericType(0)) *
        dxSource;
    const NumericType doseTgt =
        std::accumulate(target.begin(), target.end(), NumericType(0)) *
        dxTarget;
    if (doseTgt > NumericType(0) && doseSrc > NumericType(0)) {
      const NumericType scale = doseSrc / doseTgt;
      for (auto &v : target)
        v *= scale;
    }
  }

  static NumericType totalDose(const std::vector<NumericType> &C,
                               NumericType dx) {
    return std::accumulate(C.begin(), C.end(), NumericType(0)) * dx;
  }

#ifdef VIENNAPS_HAS_MFEM
  /// Project damage onto an MFEM GridFunction and rescale to targetDose = ∫C.
  static void projectToGridFunction(const DamageProfile &source,
                                    mfem::GridFunction &gf,
                                    double targetDose) {
    if (source.values.empty() || gf.Size() == 0)
      return;
    const std::size_t ns = source.values.size();
    for (int i = 0; i < gf.Size(); ++i) {
      const double xi =
          (gf.Size() == 1)
              ? 0.0
              : static_cast<double>(i) / static_cast<double>(gf.Size() - 1);
      const double sj = xi * static_cast<double>(ns - 1);
      const std::size_t j0 = static_cast<std::size_t>(sj);
      const std::size_t j1 = std::min(j0 + 1, ns - 1);
      const double f = sj - static_cast<double>(j0);
      gf(i) = static_cast<mfem::real_t>(
          (1.0 - f) * static_cast<double>(source.values[j0]) +
          f * static_cast<double>(source.values[j1]));
    }
    mfem::ConstantCoefficient one(1.0);
    mfem::LinearForm mass(gf.FESpace());
    mass.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
    mass.Assemble();
    const double doseNow = gf * mass;
    if (doseNow > 0.0 && targetDose > 0.0)
      gf *= static_cast<mfem::real_t>(targetDose / doseNow);
  }
#endif
};

} // namespace viennaps
