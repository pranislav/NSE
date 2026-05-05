#ifndef STEP57_CASE_CONFIG_H
#define STEP57_CASE_CONFIG_H

#include <deal.II/grid/tria.h>

#include <map>
#include <set>
#include <string>
#include <vector>

namespace Cht
{
  struct VelocityBoundaryCondition
  {
    enum class Type
    {
      constant,
      parabolic
    };

    dealii::types::boundary_id boundary_id;
    Type                       type;
    unsigned int               component;
    unsigned int               coordinate;
    double                     value;
  };

  struct TemperatureBoundaryCondition
  {
    dealii::types::boundary_id boundary_id;
    double                     value;
  };

  struct MaterialData
  {
    enum class Kind
    {
      fluid,
      solid
    };

    Kind   kind;
    double thermal_diffusivity;
  };

  struct CaseConfig
  {
    std::string                                           mesh_file;
    double                                                reynolds;
    double                                                gamma;
    unsigned int                                          degree;
    unsigned int                                          adaptive_refinement_cycles;
    std::map<dealii::types::material_id, MaterialData>    materials;
    std::set<dealii::types::boundary_id>                  no_slip_boundary_ids;
    std::vector<VelocityBoundaryCondition>                velocity_dirichlet_boundaries;
    std::vector<TemperatureBoundaryCondition>             temperature_dirichlet_boundaries;
  };

  CaseConfig read_case_config(const std::string &case_file);
}

#endif
