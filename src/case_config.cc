#include "case_config.h"

#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>

#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>

namespace Cht
{
  using namespace dealii;

  namespace
  {
    std::vector<std::string> split_string(const std::string &input,
                                          const char         delimiter)
    {
      std::vector<std::string> parts;
      std::stringstream        stream(input);
      std::string              item;

      while (std::getline(stream, item, delimiter))
        {
          item = Utilities::trim(item);
          if (!item.empty())
            parts.push_back(item);
        }

      return parts;
    }

    std::filesystem::path resolve_relative_path(const std::string &base_file,
                                                const std::string &path)
    {
      const std::filesystem::path file_path(path);
      if (file_path.is_absolute())
        return file_path;

      return std::filesystem::path(base_file).parent_path() / file_path;
    }

    std::set<types::material_id>
    discover_material_ids_from_case_file(const std::string &case_file)
    {
      std::ifstream input(case_file);
      AssertThrow(input, ExcFileNotOpen(case_file));

      std::set<types::material_id> material_ids;
      const std::string            prefix = "subsection Material ";
      std::string                  line;

      while (std::getline(input, line))
        {
          const std::string trimmed = Utilities::trim(line);
          if (trimmed.rfind(prefix, 0) != 0)
            continue;

          const std::string suffix =
            Utilities::trim(trimmed.substr(prefix.size()));
          AssertThrow(!suffix.empty(),
                      ExcMessage("Missing material id in subsection declaration."));

          material_ids.insert(static_cast<types::material_id>(std::stoi(suffix)));
        }

      return material_ids;
    }

    std::set<types::material_id> parse_material_id_set(const std::string &text)
    {
      std::set<types::material_id> ids;
      for (const auto &token : split_string(text, ','))
        ids.insert(static_cast<types::material_id>(std::stoi(token)));
      return ids;
    }

    MaterialData::Kind parse_material_kind(const std::string &text)
    {
      if (text == "fluid")
        return MaterialData::Kind::fluid;
      if (text == "solid")
        return MaterialData::Kind::solid;

      AssertThrow(false, ExcMessage("Unsupported material kind: " + text));
      return MaterialData::Kind::solid;
    }

    std::vector<VelocityBoundaryCondition>
    parse_velocity_boundary_conditions(const std::string &text)
    {
      std::vector<VelocityBoundaryCondition> conditions;

      for (const auto &entry : split_string(text, ';'))
        {
          const auto parts = split_string(entry, ':');
          AssertThrow(parts.size() == 5,
                      ExcMessage("Velocity boundary entries must be "
                                 "boundary_id:type:component:coordinate:value."));

          VelocityBoundaryCondition condition;
          condition.boundary_id =
            static_cast<types::boundary_id>(std::stoi(parts[0]));
          if (parts[1] == "constant")
            condition.type = VelocityBoundaryCondition::Type::constant;
          else if (parts[1] == "parabolic")
            condition.type = VelocityBoundaryCondition::Type::parabolic;
          else
            AssertThrow(false,
                        ExcMessage("Unsupported velocity boundary type: " +
                                   parts[1]));
          condition.component  = std::stoi(parts[2]);
          condition.coordinate = std::stoi(parts[3]);
          condition.value      = std::stod(parts[4]);
          conditions.push_back(condition);
        }

      return conditions;
    }

    std::vector<TemperatureBoundaryCondition>
    parse_temperature_boundary_conditions(const std::string &text)
    {
      std::vector<TemperatureBoundaryCondition> conditions;

      for (const auto &entry : split_string(text, ';'))
        {
          const auto parts = split_string(entry, ':');
          AssertThrow(parts.size() == 2,
                      ExcMessage(
                        "Temperature boundary entries must be boundary_id:value."));

          conditions.push_back(
            {static_cast<types::boundary_id>(std::stoi(parts[0])),
             std::stod(parts[1])});
        }

      return conditions;
    }
  }

  CaseConfig read_case_config(const std::string &case_file)
  {
    ParameterHandler prm;

    prm.enter_subsection("Mesh");
    prm.declare_entry("File", "", Patterns::Anything());
    prm.leave_subsection();

    prm.enter_subsection("Solver");
    prm.declare_entry("Reynolds number", "7500", Patterns::Double(0.0));
    prm.declare_entry("Gamma", "1.0", Patterns::Double(0.0));
    prm.declare_entry("Polynomial degree", "1", Patterns::Integer(1));
    prm.declare_entry("Adaptive refinement cycles", "4", Patterns::Integer(0));
    prm.declare_entry("Use MMS", "false", Patterns::Bool());
    prm.leave_subsection();

    prm.enter_subsection("Materials");
    prm.declare_entry("Ids", "0", Patterns::Anything());
    prm.leave_subsection();

    prm.enter_subsection("Boundary conditions");
    prm.declare_entry("Velocity Dirichlet", "", Patterns::Anything());
    prm.declare_entry("Temperature Dirichlet", "", Patterns::Anything());
    prm.leave_subsection();

    const auto material_ids = discover_material_ids_from_case_file(case_file);

    for (const auto material_id : material_ids)
      {
        prm.enter_subsection("Materials");
        prm.enter_subsection("Material " + std::to_string(material_id));
        prm.declare_entry("Kind", "solid", Patterns::Selection("fluid|solid"));
        prm.declare_entry("Thermal diffusivity", "1.0", Patterns::Double(0.0));
        prm.leave_subsection();
        prm.leave_subsection();
      }

    prm.parse_input(case_file);

    CaseConfig config;

    prm.enter_subsection("Mesh");
    config.mesh_file = resolve_relative_path(case_file, prm.get("File")).string();
    prm.leave_subsection();

    prm.enter_subsection("Solver");
    config.reynolds                   = prm.get_double("Reynolds number");
    config.gamma                      = prm.get_double("Gamma");
    config.degree                     = prm.get_integer("Polynomial degree");
    config.adaptive_refinement_cycles =
      prm.get_integer("Adaptive refinement cycles");
    config.use_mms = prm.get_bool("Use MMS");
    prm.leave_subsection();

    prm.enter_subsection("Materials");
    const auto declared_material_ids = parse_material_id_set(prm.get("Ids"));
    prm.leave_subsection();

    AssertThrow(declared_material_ids == material_ids,
                ExcMessage("Mismatch between 'Materials/Ids' and "
                           "'subsection Material <id>' declarations."));

    for (const auto material_id : material_ids)
      {
        prm.enter_subsection("Materials");
        prm.enter_subsection("Material " + std::to_string(material_id));

        MaterialData material;
        material.kind                = parse_material_kind(prm.get("Kind"));
        material.thermal_diffusivity = prm.get_double("Thermal diffusivity");

        config.materials.emplace(material_id, material);

        prm.leave_subsection();
        prm.leave_subsection();
      }

    prm.enter_subsection("Boundary conditions");
    config.velocity_dirichlet_boundaries =
      parse_velocity_boundary_conditions(prm.get("Velocity Dirichlet"));
    config.temperature_dirichlet_boundaries =
      parse_temperature_boundary_conditions(prm.get("Temperature Dirichlet"));
    prm.leave_subsection();

    AssertThrow(!config.mesh_file.empty(),
                ExcMessage("Case config must define a mesh file."));
    AssertThrow(!config.materials.empty(),
                ExcMessage("Case config must define at least one material."));

    return config;
  }
}
