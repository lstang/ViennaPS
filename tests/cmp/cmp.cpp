// cmp.cpp — CMP Preston planarization model tests (GAP_ANALYSIS §4.4).
#include <lsTestAsserts.hpp>
#include <models/psCMP.hpp>
#include <psDomain.hpp>
#include <vcTestAsserts.hpp>

namespace viennacore {

using namespace viennaps;

template <typename NumericType, int D>
Vec3D<NumericType> makeCoord(NumericType x, NumericType h) {
  Vec3D<NumericType> c{NumericType(0), NumericType(0), NumericType(0)};
  c[D - 1] = h;
  c[0] = x;
  return c;
}

template <class NumericType, int D> void RunTest() {
  using impl::CmpVelocityField;
  constexpr NumericType eps = NumericType(1e-9);

  // --- Test 1: Preston rate law + pattern-density modulation + hard stops.
  {
    // base rate = K_p * P * v_rel = 0.1 um/s; α = 0.5; L_p = 1.0 um;
    // reference height = 5.0 um.
    const std::vector<std::pair<Material, NumericType>> polish = {
        {Material::SiO2, NumericType(1.0)}, {Material::Mask, NumericType(2.0)}};
    const std::vector<Material> hardStop = {Material::Si};
    CmpVelocityField<NumericType, D> vf(/*prestonRate=*/NumericType(0.1),
                                        /*alpha=*/NumericType(0.5),
                                        /*planarizationLength=*/NumericType(1.0),
                                        /*refHeight=*/NumericType(5.0), polish,
                                        hardStop);

    // Protrusion (h=6) polishes faster than depression (h=4):
    // V(6) = 0.1*(1+0.5*1) = 0.15 ; V(4) = 0.1*(1-0.5*1) = 0.05.
    const auto vProtrusion = vf.getScalarVelocity(
        makeCoord<NumericType, D>(NumericType(0), NumericType(6)),
        static_cast<int>(Material::SiO2), {NumericType(0), NumericType(1), NumericType(0)}, 0);
    const auto vDepression = vf.getScalarVelocity(
        makeCoord<NumericType, D>(NumericType(0), NumericType(4)),
        static_cast<int>(Material::SiO2), {NumericType(0), NumericType(1), NumericType(0)}, 0);
    VC_TEST_ASSERT(std::abs(vProtrusion - NumericType(0.15)) < eps);
    VC_TEST_ASSERT(std::abs(vDepression - NumericType(0.05)) < eps);
    VC_TEST_ASSERT(vProtrusion > vDepression);

    // Selectivity: Mask (sel 2.0) at the reference height → 0.2.
    const auto vMask = vf.getScalarVelocity(
        makeCoord<NumericType, D>(NumericType(0), NumericType(5)),
        static_cast<int>(Material::Mask), {NumericType(0), NumericType(1), NumericType(0)}, 0);
    VC_TEST_ASSERT(std::abs(vMask - NumericType(0.2)) < eps);

    // Hard stop: Si never polishes.
    const auto vSi = vf.getScalarVelocity(
        makeCoord<NumericType, D>(NumericType(0), NumericType(6)),
        static_cast<int>(Material::Si), {NumericType(0), NumericType(1), NumericType(0)}, 0);
    VC_TEST_ASSERT(std::abs(vSi) < eps);

    // Non-polish material (e.g. Undefined / metal not in list) → 0.
    const auto vOther = vf.getScalarVelocity(
        makeCoord<NumericType, D>(NumericType(0), NumericType(6)),
        static_cast<int>(Material::W), {NumericType(0), NumericType(1), NumericType(0)}, 0);
    VC_TEST_ASSERT(std::abs(vOther) < eps);

    // Clamp: h=20 → f = 1+0.5*15 = 8.5 → clamped to 2.0 → V = 0.2.
    const auto vHigh = vf.getScalarVelocity(
        makeCoord<NumericType, D>(NumericType(0), NumericType(20)),
        static_cast<int>(Material::SiO2), {NumericType(0), NumericType(1), NumericType(0)}, 0);
    VC_TEST_ASSERT(std::abs(vHigh - NumericType(0.2)) < eps);

    // Clamp: h=-10 → f = 1-0.5*15 = -6.5 → clamped to 0.1 → V = 0.01.
    const auto vLow = vf.getScalarVelocity(
        makeCoord<NumericType, D>(NumericType(0), NumericType(-10)),
        static_cast<int>(Material::SiO2), {NumericType(0), NumericType(1), NumericType(0)}, 0);
    VC_TEST_ASSERT(std::abs(vLow - NumericType(0.01)) < eps);
  }

  // --- Test 2 (added in Task 3): integration run on a Mask-bump stack.
}

} // namespace viennacore

int main() { VC_RUN_ALL_TESTS }
