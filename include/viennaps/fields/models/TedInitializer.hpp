#pragma once

/// TedInitializer — project implant damage (I/V) onto FEM dofs with
/// integral-preserving rescaling (MOOSE IntegralPreservingFunctionIC).

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class TedInitializer {
public:
  /// Analytic / mock damage profile samples (positions optional).
  struct DamageProfile {
    std::vector<NumericType> values; ///< C_I or C_V samples
  };

  /// Project source onto target mesh samples and rescale so
  /// sum(target)*dx_target == sum(source)*dx_source (dose conservation).
  static void projectIntegralPreserving(const DamageProfile &source,
                                        NumericType dxSource,
                                        std::vector<NumericType> &target,
                                        NumericType dxTarget) {
    if (source.values.empty() || target.empty())
      return;

    // Nearest-neighbor / linear sample from source onto target length.
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
};

} // namespace viennaps
