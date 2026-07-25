#pragma once

/// NeutralReactDiffusion — I+V recombination without charge enhancement.
/// Thin alias of ReactDiffusion for API completeness (Phase 3 Task 11).

#include "ReactDiffusion.hpp"

#include <string>

namespace viennaps {

template <class NumericType>
class NeutralReactDiffusion : public ReactDiffusion<NumericType> {
public:
  NeutralReactDiffusion(std::string interstitial = "Interstitial",
                        std::string vacancy = "Vacancy")
      : ReactDiffusion<NumericType>(std::move(interstitial),
                                    std::move(vacancy)) {
    this->setName("NeutralReactDiffusion");
  }
};

} // namespace viennaps
