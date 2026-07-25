#pragma once

/// MaterialConverter — re-tag mesh materials (e.g. Si → GaAs).

#include "MeshAttributes.hpp"

#include <string>

#ifdef VIENNAPS_HAS_MFEM
#include <mfem.hpp>
#endif

namespace viennaps {

class MaterialConverter {
public:
  /// Rename attribute mapping and optionally retag all elements with oldAttr.
  static void convert(MeshAttributes &attrs, int attrId,
                      const std::string &newMaterial) {
    attrs.setAttributeName(attrId, newMaterial);
  }

#ifdef VIENNAPS_HAS_MFEM
  static void convertMesh(mfem::Mesh &mesh, MeshAttributes &attrs, int fromAttr,
                          int toAttr, const std::string &newMaterial) {
    attrs.setAttributeName(toAttr, newMaterial);
    for (int e = 0; e < mesh.GetNE(); ++e) {
      if (mesh.GetAttribute(e) == fromAttr) {
        mesh.SetAttribute(e, toAttr);
      }
    }
  }
#endif
};

} // namespace viennaps
