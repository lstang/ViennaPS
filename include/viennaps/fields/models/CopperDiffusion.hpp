#pragma once

/// CopperDiffusion — MobileImpurity specialization for Cu (Phase 4 Task 7).
///
/// D = D0 * (1 + beta * C_dopant / ni) with optional Nernst–Planck drift
/// J = −D (∇C + (q/kT) z C E) and Cu + acceptor ⇌ CuA pairing.
/// Uses CopperTag for compile-time defaults (D0=1e-5, charge=+1).

#include "MobileImpurityTags.hpp"
#include "MobileImpurity.hpp"

#include <string>

namespace viennaps {

template <class NumericType>
using CopperDiffusion = MobileImpurity<NumericType, CopperTag>;

} // namespace viennaps