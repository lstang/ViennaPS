"""Type stubs for ViennaPS multiphysics (header-only C++ APIs).

C++ headers under include/viennaps/fields and models expose:
  PhysicsField, MaterialPropertySystem, ParameterDatabase
  Diffusion/Fermi/Pair/ChargedReact/DefectCluster/SPER kernels
  Viscoelastic/Elastic/MfemElasticity stress kernels
  SundialsTimeIntegrator, AmgclSolver, BandLimitedSolver
  GeometryFieldCoupler, LocosDopingValidator, PhysicsFieldAdapter
  AnalyticImplant, MCBcaImplant, SilicidationModel, LithographyModel
  ProcessOrchestrator

Full pybind11 bindings for every kernel are Phase-2 polish; C++ is the
primary API. Import paths match viennaps.hpp umbrella includes.
"""

from typing import Any

class PhysicsField:
    def add_species(self, name: str) -> None: ...
    def get_total_dose(self, species: str) -> float: ...
    def inject_implant_profile(self, species: str, profile: list[float]) -> None: ...
    def add_dose(self, species: str, amount: float) -> None: ...
    def set_species_dose(self, species: str, dose: float) -> None: ...

class ParameterDatabase:
    def get(self, material: str, key: str, T_K: float = ...) -> float: ...
    def blend(self, mat_a: str, mat_b: str, key: str, weight_b: float, T_K: float = ...) -> float: ...

class MCBcaEngine:
    def run(self) -> Any: ...

class SilicidationModel:
    def evolve(self, field: PhysicsField) -> float: ...

class LithographyModel:
    def apply(self, field: PhysicsField, n_x: int = ...) -> None: ...
    def get_open_fraction(self) -> float: ...
