#pragma once

#include <string>
#include <map>

namespace viennaps {

class MeshAttributes {
public:
  void setAttributeName(int attr, const std::string& name) {
    attrToName_[attr] = name;
    nameToAttr_[name] = attr;
  }

  const std::string& materialName(int attr) const {
    auto it = attrToName_.find(attr);
    if (it == attrToName_.end()) {
      static const std::string empty;
      return empty;
    }
    return it->second;
  }

  bool isMaterial(int attr, const std::string& name) const {
    return materialName(attr) == name;
  }

  int numMaterials() const { return static_cast<int>(attrToName_.size()); }

  bool hasMaterial(const std::string& name) const {
    return nameToAttr_.find(name) != nameToAttr_.end();
  }

  int attributeOf(const std::string& name) const {
    auto it = nameToAttr_.find(name);
    return it == nameToAttr_.end() ? -1 : it->second;
  }

private:
  std::map<int, std::string> attrToName_;
  std::map<std::string, int> nameToAttr_;
};

} // namespace viennaps
