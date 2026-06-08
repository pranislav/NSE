#include "conjugate_heat_transfer_solver.h"

#include "boundary_conditions.h"

#include "mms.h"

#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <sstream>

namespace Cht
{
  using namespace dealii;

  namespace
  {
    std::string join_with_underscore(const std::vector<std::string> &parts)
    {
      std::ostringstream out;
      for (std::size_t i = 0; i < parts.size(); ++i)
        {
          if (i > 0)
            out << "_";
          out << parts[i];
        }
      return out.str();
    }

    template <typename IdType>
    std::string join_ids(const std::set<IdType> &ids)
    {
      std::ostringstream out;
      bool               first = true;

      for (const auto id : ids)
        {
          if (!first)
            out << ", ";
          out << static_cast<unsigned int>(id);
          first = false;
        }

      return out.str();
    }

    template <class PreconditionerMp>
    class BlockSchurPreconditioner : public EnableObserverPointer
    {
    public:
      BlockSchurPreconditioner(double                           gamma,
                               double                           viscosity,
                               const BlockSparseMatrix<double> &S,
                               const SparseMatrix<double>      &P,
                               const PreconditionerMp          &Mppreconditioner)
        : gamma(gamma)
        , viscosity(viscosity)
        , stokes_matrix(S)
        , pressure_mass_matrix(P)
        , mp_preconditioner(Mppreconditioner)
      {
        A_inverse.initialize(stokes_matrix.block(0, 0));
      }

      void vmult(BlockVector<double>       &dst,
                 const BlockVector<double> &src) const
      {
        Vector<double> utmp(src.block(0));

        {
          SolverControl solver_control(1000, 1e-6 * src.block(1).l2_norm());
          SolverCG<Vector<double>> cg(solver_control);

          dst.block(1) = 0.0;
          cg.solve(pressure_mass_matrix,
                   dst.block(1),
                   src.block(1),
                   mp_preconditioner);
          dst.block(1) *= -(viscosity + gamma);
        }

        {
          stokes_matrix.block(0, 1).vmult(utmp, dst.block(1));
          utmp *= -1.0;
          utmp += src.block(0);
        }

        A_inverse.vmult(dst.block(0), utmp);
      }

    private:
      const double                     gamma;
      const double                     viscosity;
      const BlockSparseMatrix<double> &stokes_matrix;
      const SparseMatrix<double>      &pressure_mass_matrix;
      const PreconditionerMp          &mp_preconditioner;
      SparseDirectUMFPACK              A_inverse;
    };
  }

  template <int dim>
  ConjugateHeatTransferSolver<dim>::BoundaryExtent::BoundaryExtent()
    : initialized(false)
  {
    min.fill(std::numeric_limits<double>::max());
    max.fill(std::numeric_limits<double>::lowest());
  }

  template <int dim>
  ConjugateHeatTransferSolver<dim>::ConjugateHeatTransferSolver(
    const CaseConfig  &config,
    const std::string &output_directory,
    const bool         save_mesh_output,
    const bool         output_partial_solutions,
    const RefinementMode refinement_mode)
    : config(config)
    , refinement_mode(refinement_mode)
    , viscosity(1.0 / config.reynolds)
    , gamma(config.gamma)
    , degree(config.degree)
    , output_directory(output_directory)
    , save_mesh_output(save_mesh_output)
    , output_partial_solutions(output_partial_solutions)
    , triangulation(Triangulation<dim>::maximum_smoothing)
    , fe_fluid(FE_Q<dim>(degree + 1) ^ dim, FE_Q<dim>(degree))
    , fe_solid(FE_Nothing<dim>(), dim, FE_Nothing<dim>(), 1)
    , temperature_fe(degree + 1)
    , dof_handler(triangulation)
    , temperature_dof_handler(triangulation)
  {
    fe_collection.push_back(fe_fluid);
    fe_collection.push_back(fe_solid);
  }

  template <int dim>
  std::string ConjugateHeatTransferSolver<dim>::format_compact_double(
    const double value)
  {
    std::ostringstream out;
    out << std::fixed << std::setprecision(6) << value;
    std::string text = out.str();
    while (!text.empty() && text.back() == '0')
      text.pop_back();
    if (!text.empty() && text.back() == '.')
      text.pop_back();
    if (text.empty())
      text = "0";
    return text;
  }

  template <int dim>
  std::string ConjugateHeatTransferSolver<dim>::case_tag() const
  {
    std::vector<std::string> diffusivity_parts;
    for (const auto &[material_id, material] : config.materials)
      {
        (void)material_id;
        diffusivity_parts.push_back(
          format_compact_double(material.thermal_diffusivity));
      }

    std::vector<std::string> velocity_parts;
    for (const auto &boundary : config.velocity_dirichlet_boundaries)
      velocity_parts.push_back(format_compact_double(boundary.value));

    std::vector<std::string> temperature_parts;
    for (const auto &boundary : config.temperature_dirichlet_boundaries)
      temperature_parts.push_back(std::to_string(boundary.boundary_id) + "_" +
                                  format_compact_double(boundary.value));

    const std::string temperature_tag =
      temperature_parts.empty() ? "none" : join_with_underscore(temperature_parts);

    return "kappa_" + join_with_underscore(diffusivity_parts) + "-v_" +
           join_with_underscore(velocity_parts) + "-Re_" +
           format_compact_double(config.reynolds) + "-T_" + temperature_tag;
  }

  template <int dim>
  const MaterialData &
  ConjugateHeatTransferSolver<dim>::material_data(
    const types::material_id material_id) const
  {
    const auto it = config.materials.find(material_id);
    AssertThrow(it != config.materials.end(),
                ExcMessage("Missing material definition for material id " +
                           std::to_string(material_id) + "."));
    return it->second;
  }

  template <int dim>
  bool ConjugateHeatTransferSolver<dim>::cell_is_in_fluid_domain(
    const typename DoFHandler<dim>::cell_iterator &cell) const
  {
    return material_data(cell->material_id()).kind == MaterialData::Kind::fluid;
  }

  template <int dim>
  bool ConjugateHeatTransferSolver<dim>::cell_is_in_solid_domain(
    const typename DoFHandler<dim>::cell_iterator &cell) const
  {
    return material_data(cell->material_id()).kind == MaterialData::Kind::solid;
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::set_active_fe_indices()
  {
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell_is_in_fluid_domain(cell))
          cell->set_active_fe_index(0);
        else if (cell_is_in_solid_domain(cell))
          cell->set_active_fe_index(1);
        else
          DEAL_II_NOT_IMPLEMENTED();
      }
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::add_fluid_solid_interface_constraints(
    AffineConstraints<double> &constraints) const
  {
    std::vector<types::global_dof_index> local_face_dof_indices(
      fe_fluid.n_dofs_per_face());

    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell_is_in_fluid_domain(cell))
        for (const auto face_no : cell->face_indices())
          if (!cell->face(face_no)->at_boundary())
            {
              bool face_is_on_interface = false;

              if (!cell->neighbor(face_no)->has_children() &&
                  cell_is_in_solid_domain(cell->neighbor(face_no)))
                face_is_on_interface = true;
              else if (cell->neighbor(face_no)->has_children())
                {
                  for (unsigned int sf = 0;
                       sf < cell->face(face_no)->n_children();
                       ++sf)
                    if (cell_is_in_solid_domain(
                          cell->neighbor_child_on_subface(face_no, sf)))
                      {
                        face_is_on_interface = true;
                        break;
                      }
                }

              if (face_is_on_interface)
                {
                  cell->face(face_no)->get_dof_indices(local_face_dof_indices,
                                                       0);

                  for (unsigned int i = 0; i < local_face_dof_indices.size(); ++i)
                    if (fe_fluid.face_system_to_component_index(i).first < dim)
                      constraints.constrain_dof_to_zero(local_face_dof_indices[i]);
                }
            }
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::collect_boundary_extents()
  {
    boundary_extents.clear();

    for (const auto &cell : triangulation.active_cell_iterators())
      for (const auto face_no : cell->face_indices())
        if (cell->face(face_no)->at_boundary())
          {
            auto &extent = boundary_extents[cell->face(face_no)->boundary_id()];

            for (const auto vertex_no : cell->face(face_no)->vertex_indices())
              {
                const Point<dim> &vertex = cell->face(face_no)->vertex(vertex_no);
                extent.initialized       = true;

                for (unsigned int d = 0; d < dim; ++d)
                  {
                    extent.min[d] = std::min(extent.min[d], vertex[d]);
                    extent.max[d] = std::max(extent.max[d], vertex[d]);
                  }
              }
          }
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::validate_case_against_mesh() const
  {
    std::set<types::material_id> mesh_material_ids;
    std::set<types::material_id> case_material_ids;
    std::set<types::material_id> missing_in_case;
    std::set<types::material_id> unused_in_case;

    for (const auto &entry : config.materials)
      case_material_ids.insert(entry.first);

    for (const auto &cell : triangulation.active_cell_iterators())
      mesh_material_ids.insert(cell->material_id());

    std::set_difference(mesh_material_ids.begin(),
                        mesh_material_ids.end(),
                        case_material_ids.begin(),
                        case_material_ids.end(),
                        std::inserter(missing_in_case, missing_in_case.end()));

    std::set_difference(case_material_ids.begin(),
                        case_material_ids.end(),
                        mesh_material_ids.begin(),
                        mesh_material_ids.end(),
                        std::inserter(unused_in_case, unused_in_case.end()));

    AssertThrow(missing_in_case.empty(),
                ExcMessage("Mesh uses material ids {" +
                           join_ids(mesh_material_ids) +
                           "}, but case file defines {" +
                           join_ids(case_material_ids) +
                           "}. Missing material definitions for: {" +
                           join_ids(missing_in_case) + "}."));

    AssertThrow(unused_in_case.empty(),
                ExcMessage("Case file defines material ids {" +
                           join_ids(case_material_ids) + "}, but mesh uses {" +
                           join_ids(mesh_material_ids) +
                           "}. Unused material definitions: {" +
                           join_ids(unused_in_case) + "}."));
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::setup_temperature_dofs()
  {
    temperature_matrix.clear();

    temperature_dof_handler.distribute_dofs(temperature_fe);

    temperature_constraints.clear();
    DoFTools::make_hanging_node_constraints(temperature_dof_handler,
                                            temperature_constraints);
    for (const auto &boundary : config.temperature_dirichlet_boundaries)
      VectorTools::interpolate_boundary_values(
        temperature_dof_handler,
        boundary.boundary_id,
        TemperatureBoundaryValues<dim>(boundary.value),
        temperature_constraints);

    temperature_constraints.close();

    DynamicSparsityPattern dsp(temperature_dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(temperature_dof_handler,
                                    dsp,
                                    temperature_constraints);
    temperature_sparsity_pattern.copy_from(dsp);

    temperature_matrix.reinit(temperature_sparsity_pattern);
    temperature_solution.reinit(temperature_dof_handler.n_dofs());
    temperature_rhs.reinit(temperature_dof_handler.n_dofs());
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::setup_dofs()
  {
    setup_flow_dofs();
    setup_temperature_dofs();
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::setup_flow_dofs()
  {
    system_matrix.clear();
    pressure_mass_matrix.clear();
    collect_boundary_extents();

    set_active_fe_indices();

    dof_handler.distribute_dofs(fe_collection);

    std::vector<unsigned int> block_component(dim + 1, 0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise(dof_handler, block_component);

    dofs_per_block =
      DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
    const unsigned int dof_u = dofs_per_block[0];
    const unsigned int dof_p = dofs_per_block[1];

    const FEValuesExtractors::Vector velocities(0);
    {
      nonzero_constraints.clear();

      DoFTools::make_hanging_node_constraints(dof_handler, nonzero_constraints);

      if (config.use_mms)
        {
          for (const auto &entry : boundary_extents)
            {
              const auto boundary_id = entry.first;
              VectorTools::interpolate_boundary_values(
                dof_handler,
                boundary_id,
                MMS::VelocityBoundaryValues<dim>(),
                nonzero_constraints,
                fe_collection.component_mask(velocities));
            }
        }
      else
        {
          for (const auto &boundary : config.velocity_dirichlet_boundaries)
            {
              const auto extent_it = boundary_extents.find(boundary.boundary_id);
              AssertThrow(extent_it != boundary_extents.end() &&
                            extent_it->second.initialized,
                          ExcMessage("Missing boundary extent for configured "
                                    "velocity boundary."));

              AssertIndexRange(boundary.component, dim);
              AssertIndexRange(boundary.coordinate, dim);

              const auto &extent = extent_it->second;

              VectorTools::interpolate_boundary_values(
                dof_handler,
                boundary.boundary_id,
                VelocityBoundaryValues<dim>(boundary,
                                            extent.min[boundary.coordinate],
                                            extent.max[boundary.coordinate]),
                nonzero_constraints,
                fe_collection.component_mask(velocities));
            }
          }

      if (not config.use_mms)
        add_fluid_solid_interface_constraints(nonzero_constraints);
    }
    nonzero_constraints.close();

    {
      zero_constraints.clear();

      DoFTools::make_hanging_node_constraints(dof_handler, zero_constraints);

      if (config.use_mms)
        {
          for (const auto &entry : boundary_extents)
            {
              const auto boundary_id = entry.first;
              VectorTools::interpolate_boundary_values(
                dof_handler,
                boundary_id,
                Functions::ZeroFunction<dim>(dim + 1),
                zero_constraints,
                fe_collection.component_mask(velocities));
            }
        }
      else
        {
          for (const auto &boundary : config.velocity_dirichlet_boundaries)
            VectorTools::interpolate_boundary_values(
              dof_handler,
              boundary.boundary_id,
              Functions::ZeroFunction<dim>(dim + 1),
              zero_constraints,
              fe_collection.component_mask(velocities));
        }

      add_fluid_solid_interface_constraints(zero_constraints);
    }
    zero_constraints.close();

    std::cout << "Number of active cells: " << triangulation.n_active_cells()
              << std::endl
              << "Number of degrees of freedom: " << dof_handler.n_dofs()
              << " (" << dof_u << " + " << dof_p << ')' << std::endl;
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::initialize_system()
  {
    {
      BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);
      DoFTools::make_sparsity_pattern(dof_handler, dsp, nonzero_constraints);
      sparsity_pattern.copy_from(dsp);
    }

    system_matrix.reinit(sparsity_pattern);

    flow_solution.reinit(dofs_per_block);
    newton_update.reinit(dofs_per_block);
    flow_rhs.reinit(dofs_per_block);
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::assemble_nse(const bool initial_step,
                                                  const bool assemble_matrix)
  {
    if (assemble_matrix)
      system_matrix = 0;

    flow_rhs = 0;

    const UpdateFlags update_flags =
      update_values | update_quadrature_points | update_JxW_values |
      update_gradients;

    hp::QCollection<dim> q_collection;
    q_collection.push_back(QGauss<dim>(degree + 2));
    q_collection.push_back(QGauss<dim>(degree + 2));

    hp::FEValues<dim> hp_fe_values(fe_collection, q_collection, update_flags);

    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar pressure(dim);

    Cht::MMS::RightHandSide<dim> manufactured_rhs(viscosity);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        hp_fe_values.reinit(cell);
        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
        const unsigned int   n_q_points = fe_values.n_quadrature_points;

        std::vector<Tensor<1, dim>> present_velocity_values(n_q_points);
        std::vector<Tensor<2, dim>> present_velocity_gradients(n_q_points);
        std::vector<double>         present_pressure_values(n_q_points);

        const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();

        FullMatrix<double>                 local_matrix(dofs_per_cell, dofs_per_cell);
        Vector<double>                     local_rhs(dofs_per_cell);
        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

        std::vector<double>         div_phi_u(dofs_per_cell);
        std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
        std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell);
        std::vector<double>         phi_p(dofs_per_cell);

        local_matrix = 0;
        local_rhs    = 0;

        fe_values[velocities].get_function_values(evaluation_point,
                                                  present_velocity_values);

        fe_values[velocities].get_function_gradients(evaluation_point,
                                                     present_velocity_gradients);

        fe_values[pressure].get_function_values(evaluation_point,
                                                present_pressure_values);

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            for (unsigned int k = 0; k < dofs_per_cell; ++k)
              {
                div_phi_u[k]  = fe_values[velocities].divergence(k, q);
                grad_phi_u[k] = fe_values[velocities].gradient(k, q);
                phi_u[k]      = fe_values[velocities].value(k, q);
                phi_p[k]      = fe_values[pressure].value(k, q);
              }
              const Tensor<1, dim> f = 
                config.use_mms ?
                  manufactured_rhs.value(fe_values.quadrature_point(q)) :
                  Tensor<1, dim>();

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                if (assemble_matrix)
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    local_matrix(i, j) +=
                      (viscosity * scalar_product(grad_phi_u[i], grad_phi_u[j]) +
                       phi_u[i] *
                         (present_velocity_gradients[q] * phi_u[j]) +
                       phi_u[i] *
                         (grad_phi_u[j] * present_velocity_values[q]) -
                       div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j] +
                       gamma * div_phi_u[i] * div_phi_u[j] +
                       phi_p[i] * phi_p[j]) *
                      fe_values.JxW(q);

                const double present_velocity_divergence =
                  trace(present_velocity_gradients[q]);
                local_rhs(i) +=
                  (-viscosity *
                     scalar_product(grad_phi_u[i],
                                    present_velocity_gradients[q]) -
                   phi_u[i] *
                     (present_velocity_gradients[q] * present_velocity_values[q]) +
                   div_phi_u[i] * present_pressure_values[q] +
                   phi_p[i] * present_velocity_divergence -
                   gamma * div_phi_u[i] * present_velocity_divergence) *
                  fe_values.JxW(q);
                if (config.use_mms)
                  {
                    local_rhs(i) += (phi_u[i] * f) * fe_values.JxW(q);
                  }
              }
          }

        cell->get_dof_indices(local_dof_indices);

        const AffineConstraints<double> &constraints_used =
          initial_step ? nonzero_constraints : zero_constraints;

        if (assemble_matrix)
          constraints_used.distribute_local_to_global(local_matrix,
                                                      local_rhs,
                                                      local_dof_indices,
                                                      system_matrix,
                                                      flow_rhs);
        else
          constraints_used.distribute_local_to_global(local_rhs,
                                                      local_dof_indices,
                                                      flow_rhs);
      }

    if (assemble_matrix)
      {
        pressure_mass_matrix.reinit(sparsity_pattern.block(1, 1));
        pressure_mass_matrix.copy_from(system_matrix.block(1, 1));
        system_matrix.block(1, 1) = 0;
      }
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::assemble_flow_system(const bool initial_step)
  {
    assemble_nse(initial_step, true);
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::assemble_flow_rhs(const bool initial_step)
  {
    assemble_nse(initial_step, false);
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::solve_flow(const bool initial_step)
  {
    const AffineConstraints<double> &constraints_used =
      initial_step ? nonzero_constraints : zero_constraints;

    SolverControl solver_control(system_matrix.m(),
                                 1e-4 * flow_rhs.l2_norm(),
                                 true);

    SolverFGMRES<BlockVector<double>> gmres(solver_control);
    SparseILU<double>                 pmass_preconditioner;
    pmass_preconditioner.initialize(pressure_mass_matrix,
                                    SparseILU<double>::AdditionalData());

    const BlockSchurPreconditioner<SparseILU<double>> preconditioner(
      gamma,
      viscosity,
      system_matrix,
      pressure_mass_matrix,
      pmass_preconditioner);

    gmres.solve(system_matrix, newton_update, flow_rhs, preconditioner);
    std::cout << "FGMRES steps: " << solver_control.last_step() << std::endl;

    constraints_used.distribute(newton_update);
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::assemble_temperature_system()
  {
    temperature_matrix = 0;
    temperature_rhs    = 0;

    const QGauss<dim> quadrature_formula(degree + 2);
    hp::QCollection<dim> flow_quadratures;
    flow_quadratures.push_back(quadrature_formula);
    flow_quadratures.push_back(quadrature_formula);
    hp::FEValues<dim> flow_fe_values(fe_collection,
                                     flow_quadratures,
                                     update_values);
    FEValues<dim> temperature_fe_values(temperature_fe,
                                        quadrature_formula,
                                        update_values | update_gradients |
                                          update_JxW_values);
    const FEValuesExtractors::Vector velocities(0);

    const unsigned int dofs_per_cell = temperature_fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double>                 local_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>                     local_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    std::vector<Tensor<1, dim>>           grad_phi_T(dofs_per_cell);
    std::vector<double>                   phi_T(dofs_per_cell);
    std::vector<Tensor<1, dim>>           velocity_values(n_q_points);
    std::vector<typename DoFHandler<dim>::active_cell_iterator> flow_cells(
      triangulation.n_active_cells(), dof_handler.end());

    for (const auto &flow_cell : dof_handler.active_cell_iterators())
      flow_cells[flow_cell->active_cell_index()] = flow_cell;

    for (const auto &cell : temperature_dof_handler.active_cell_iterators())
      {
        local_matrix = 0;
        local_rhs    = 0;
        temperature_fe_values.reinit(cell);

        const auto flow_cell = flow_cells[cell->active_cell_index()];
        Assert(flow_cell != dof_handler.end(), ExcInternalError());

        const bool   in_fluid_domain = cell_is_in_fluid_domain(flow_cell);
        const double thermal_diffusivity =
          material_data(cell->material_id()).thermal_diffusivity;

        if (in_fluid_domain)
          {
            flow_fe_values.reinit(flow_cell);
            flow_fe_values.get_present_fe_values()[velocities].get_function_values(
              flow_solution, velocity_values);
          }

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            for (unsigned int k = 0; k < dofs_per_cell; ++k)
              {
                phi_T[k]      = temperature_fe_values.shape_value(k, q);
                grad_phi_T[k] = temperature_fe_values.shape_grad(k, q);
              }

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  local_matrix(i, j) +=
                    (thermal_diffusivity * grad_phi_T[i] * grad_phi_T[j] *
                     temperature_fe_values.JxW(q));

                  if (in_fluid_domain)
                    local_matrix(i, j) +=
                      ((velocity_values[q] * grad_phi_T[j]) * phi_T[i] *
                       temperature_fe_values.JxW(q));
                }
          }

        cell->get_dof_indices(local_dof_indices);
        temperature_constraints.distribute_local_to_global(local_matrix,
                                                           local_rhs,
                                                           local_dof_indices,
                                                           temperature_matrix,
                                                           temperature_rhs);
      }
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::solve_temperature()
  {
    SolverControl solver_control(temperature_matrix.m(),
                                 1e-6 * temperature_rhs.l2_norm() + 1e-12);
    SolverGMRES<Vector<double>> gmres(solver_control);
    SparseILU<double>           preconditioner;
    preconditioner.initialize(temperature_matrix,
                              SparseILU<double>::AdditionalData());
    gmres.solve(temperature_matrix,
                temperature_solution,
                temperature_rhs,
                preconditioner);
    temperature_constraints.distribute(temperature_solution);
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::update_temperature_field()
  {
    assemble_temperature_system();
    solve_temperature();
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::mark_cells_for_adaptive_refinement()
  {
    Vector<float> flow_error_per_cell(triangulation.n_active_cells());
    Vector<float> temperature_error_per_cell(triangulation.n_active_cells());
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
    const FEValuesExtractors::Vector velocity(0);

    hp::QCollection<dim - 1> face_quadratures;
    for (unsigned int fe_index = 0; fe_index < fe_collection.size(); ++fe_index)
      face_quadratures.push_back(QGauss<dim - 1>(degree + 1));

    KellyErrorEstimator<dim>::estimate(
      dof_handler,
      face_quadratures,
      std::map<types::boundary_id, const Function<dim> *>(),
      flow_solution,
      flow_error_per_cell,
      fe_collection.component_mask(velocity));

    if (temperature_solution.size() == temperature_dof_handler.n_dofs())
      KellyErrorEstimator<dim>::estimate(
        temperature_dof_handler,
        QGauss<dim - 1>(degree + 1),
        std::map<types::boundary_id, const Function<dim> *>(),
        temperature_solution,
        temperature_error_per_cell);

    const float max_flow_error        = flow_error_per_cell.linfty_norm();
    const float max_temperature_error = temperature_error_per_cell.linfty_norm();

    for (unsigned int cell = 0; cell < estimated_error_per_cell.size(); ++cell)
      {
        const float normalized_flow_error =
          max_flow_error > 0.0f ? flow_error_per_cell[cell] / max_flow_error :
                                  0.0f;
        const float normalized_temperature_error =
          max_temperature_error > 0.0f ?
            temperature_error_per_cell[cell] / max_temperature_error :
            0.0f;

        estimated_error_per_cell[cell] =
          std::max(normalized_flow_error, normalized_temperature_error);
      }

    GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                    estimated_error_per_cell,
                                                    0.3,
                                                    0.0);
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::mark_cells_for_global_refinement()
  {
    for (const auto &cell : triangulation.active_cell_iterators())
      cell->set_refine_flag();
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::execute_refinement()
  {
    triangulation.execute_coarsening_and_refinement();
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::refine_mesh(
    const unsigned int refinement_cycle)
  {
    SolutionTransfer<dim, BlockVector<double>> solution_transfer(dof_handler);

    switch (refinement_mode)
      {
        case global_refinement:
          mark_cells_for_global_refinement();
          break;
        case adaptive_refinement:
          mark_cells_for_adaptive_refinement();
          break;
        default:
          AssertThrow(false, ExcMessage("Unhandled refinement_mode in switch"));
      }

    triangulation.prepare_coarsening_and_refinement();
    solution_transfer.prepare_for_coarsening_and_refinement(flow_solution);
    execute_refinement();

    setup_dofs();

    BlockVector<double> tmp(dofs_per_block);
    solution_transfer.interpolate(tmp);
    nonzero_constraints.distribute(tmp);

    initialize_system();
    flow_solution = tmp;

    if (save_mesh_output)
      output_mesh(refinement_cycle);
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::newton_iteration(
    const double       tolerance,
    const unsigned int max_n_line_searches,
    const unsigned int max_n_refinements,
    const bool         is_initial_step,
    const bool         output_result)
  {
    bool first_step = is_initial_step;

    for (unsigned int refinement_n = 0; refinement_n < max_n_refinements + 1;
         ++refinement_n)
      {
        unsigned int line_search_n = 0;
        double       last_res      = 1.0;
        double       current_res   = 1.0;
        std::cout << "grid refinements: " << refinement_n << std::endl
                  << "viscosity: " << viscosity << std::endl;

        while ((first_step || (current_res > tolerance)) &&
               line_search_n < max_n_line_searches)
          {
            if (first_step)
              {
                setup_dofs();
                initialize_system();
                evaluation_point = flow_solution;
                assemble_flow_system(first_step);
                solve_flow(first_step);
                flow_solution = newton_update;
                nonzero_constraints.distribute(flow_solution);
                first_step       = false;
                evaluation_point = flow_solution;
                assemble_flow_rhs(first_step);
                current_res = flow_rhs.l2_norm();
                std::cout << "The residual of initial guess is " << current_res
                          << std::endl;
                last_res = current_res;
              }
            else
              {
                evaluation_point = flow_solution;
                assemble_flow_system(first_step);
                solve_flow(first_step);

                for (double alpha = 1.0; alpha > 1e-5; alpha *= 0.5)
                  {
                    evaluation_point = flow_solution;
                    evaluation_point.add(alpha, newton_update);
                    nonzero_constraints.distribute(evaluation_point);
                    assemble_flow_rhs(first_step);
                    current_res = flow_rhs.l2_norm();
                    std::cout << "  alpha: " << std::setw(10) << alpha
                              << std::setw(0)
                              << "  residual: " << current_res << std::endl;
                    if (current_res < last_res)
                      break;
                  }

                flow_solution = evaluation_point;
                std::cout << "  number of line searches: " << line_search_n
                          << "  residual: " << current_res << std::endl;
                last_res = current_res;
                ++line_search_n;
              }

            if (output_result &&
                (output_partial_solutions ||
                 (refinement_n == max_n_refinements && current_res <= tolerance)))
              {
                update_temperature_field();
                output_results(refinement_n, line_search_n);
              }
          }
        
        if (config.use_mms)
          compute_errors(refinement_n);
        
        if (refinement_n < max_n_refinements)
          refine_mesh(refinement_n + 1);
      }
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::compute_initial_guess(double step_size)
  {
    const double target_Re = 1.0 / viscosity;
    bool         is_initial_step = true;

    for (double Re = 1000.0; Re < target_Re;
         Re        = std::min(Re + step_size, target_Re))
      {
        viscosity = 1.0 / Re;
        std::cout << "Searching for initial guess with Re = " << Re
                  << std::endl;
        newton_iteration(1e-12, 50, 0, is_initial_step, false);
        is_initial_step = false;
      }
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::output_results(
    const unsigned int refinement_cycle,
    const unsigned int newton_step) const
  {
    std::filesystem::create_directories(output_directory);

    std::vector<std::string> solution_names(dim, "velocity");
    solution_names.emplace_back("pressure");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(flow_solution,
                             solution_names,
                             DataOut<dim>::type_dof_data,
                             data_component_interpretation);

    BlockVector<double> flow_error;
    if (config.use_mms)
      {
        flow_error.reinit(flow_solution);
        VectorTools::interpolate(dof_handler,
                                 Cht::MMS::Solution<dim>(),
                                 flow_error);

        BlockVector<double> normalized_flow_solution(flow_solution);
        const double        mean_pressure = VectorTools::compute_mean_value(
          dof_handler,
          QGauss<dim>(degree + 2),
          normalized_flow_solution,
          dim);
        normalized_flow_solution.block(1).add(-mean_pressure);

        flow_error -= normalized_flow_solution;

        std::vector<std::string> error_names(dim, "velocity_error");
        error_names.emplace_back("pressure_error");
        data_out.add_data_vector(flow_error,
                                 error_names,
                                 DataOut<dim>::type_dof_data,
                                 data_component_interpretation);
      }

    Vector<double> temperature_error;
    if (temperature_solution.size() == temperature_dof_handler.n_dofs() &&
        temperature_dof_handler.n_dofs() > 0)
      {
        data_out.add_data_vector(temperature_dof_handler,
                                 temperature_solution,
                                 "temperature");

        if (config.use_mms)
          {
            temperature_error.reinit(temperature_solution.size());
            VectorTools::interpolate(temperature_dof_handler,
                                     Cht::MMS::TemperatureSolution<dim>(),
                                     temperature_error);
            temperature_error -= temperature_solution;
            data_out.add_data_vector(temperature_dof_handler,
                                     temperature_error,
                                     "temperature_error");
          }
      }
    data_out.build_patches();

    std::ofstream output(output_directory + "/" + case_tag() + "_ref" +
                         std::to_string(refinement_cycle) + "_newt" +
                         std::to_string(newton_step) + ".vtk");
    data_out.write_vtk(output);

  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::output_mesh(
    const unsigned int refinement_cycle) const
  {
    GridOut        grid_out;
    const std::string mesh_stem =
      std::filesystem::path(config.mesh_file).stem().string();
    std::ofstream output(output_directory + "/mesh_" + mesh_stem + "_ref" +
                         std::to_string(refinement_cycle) + ".vtu");
    grid_out.write_vtu(triangulation, output);
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::compute_errors(const unsigned int cycle)
  {
    // Compute the mean pressure $\frac{1}{\Omega} \int_{\Omega} p(x) dx $
    // and then subtract it from each pressure coefficient. This will result
    // in a pressure with mean value zero. Here we make use of the fact that
    // the pressure is component $dim$ and that the finite element space
    // is nodal.
    const double mean_pressure = VectorTools::compute_mean_value(
      dof_handler, QGauss<dim>(degree + 2), flow_solution, dim); // TODO former pressure_degree (right change?)
    flow_solution.block(1).add(-mean_pressure);
    std::cout << "   Note: The mean value was adjusted by " << -mean_pressure
              << std::endl;

    const ComponentSelectFunction<dim> pressure_mask(dim, dim + 1);
    const ComponentSelectFunction<dim> velocity_mask(std::make_pair(0, dim),
                                                     dim + 1);

    Vector<float> difference_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler,
                                      flow_solution,
                                      Cht::MMS::Solution<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(degree + 2),
                                      VectorTools::L2_norm,
                                      &velocity_mask);

    const double Velocity_L2_error =
      VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::L2_norm);

    VectorTools::integrate_difference(dof_handler,
                                      flow_solution,
                                      Cht::MMS::Solution<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(degree + 2),
                                      VectorTools::L2_norm,
                                      &pressure_mask);

    const double Pressure_L2_error =
      VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::L2_norm);

    VectorTools::integrate_difference(dof_handler,
                                      flow_solution,
                                      Cht::MMS::Solution<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(degree + 2),
                                      VectorTools::H1_norm,
                                      &velocity_mask);

    const double Velocity_H1_error =
      VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::H1_norm);

    std::cout << std::endl
              << "   Velocity L2 Error: " << Velocity_L2_error << std::endl
              << "   Pressure L2 Error: " << Pressure_L2_error << std::endl
              << "   Velocity H1 Error: " << Velocity_H1_error << std::endl;
    
    const unsigned int n_active_cells = triangulation.n_active_cells();
    const unsigned int n_dofs         = dof_handler.n_dofs();
              
    convergence_table.add_value("cycle", cycle);
    convergence_table.add_value("cells", n_active_cells);
    convergence_table.add_value("dofs", n_dofs);
    convergence_table.add_value("L2_velocity", Velocity_L2_error);
    convergence_table.add_value("L2_pressure", Pressure_L2_error);
    convergence_table.add_value("H1_velocity", Velocity_H1_error);
    // convergence_table.add_value("Linfty", Linfty_error);
  }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::make_error_table()
    {
        
      convergence_table.set_scientific("L2_velocity", true);
      convergence_table.set_scientific("L2_pressure", true);
      convergence_table.set_scientific("H1_velocity", true);
      
      std::string error_filename = make_error_filename();
      std::ofstream org_mode_table(error_filename + ".org");
      convergence_table.write_text(org_mode_table, TableHandler::TextOutputFormat::org_mode_table);


      convergence_table.set_precision("L2_velocity", 3);
      convergence_table.set_precision("L2_pressure", 3);
      convergence_table.set_precision("H1_velocity", 3);
  
      convergence_table.set_tex_caption("cells", "\\# cells");
      convergence_table.set_tex_caption("dofs", "\\# dofs");
      convergence_table.set_tex_caption("L2_velocity", "$L^2$ velocity");
      convergence_table.set_tex_caption("L2_pressure", "$L^2$ pressure");
      convergence_table.set_tex_caption("H1_velocity", "$H^1$ velocity");
  
      // convergence_table.set_tex_format("cells", "r");
      // convergence_table.set_tex_format("dofs", "r");
  
      std::cout << std::endl;
      convergence_table.write_text(std::cout);

      std::ofstream error_table_file(error_filename + ".tex");
  
      convergence_table.write_tex(error_table_file);
  
    }

    template <int dim>
    std::string ConjugateHeatTransferSolver<dim>::make_error_filename()
    {
      std::string error_filename = output_directory + "/error";
      switch (refinement_mode)
        {
          case global_refinement:
            error_filename += "-global";
            break;
          case adaptive_refinement:
            error_filename += "-adaptive";
            break;
          default:
            DEAL_II_ASSERT_UNREACHABLE();
        }
  
      error_filename += "-q" + std::to_string(degree);

      return error_filename;
    }

  template <int dim>
  void ConjugateHeatTransferSolver<dim>::run(const unsigned int refinement)
  {
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(triangulation);
    std::ifstream input_file(config.mesh_file);
    AssertThrow(input_file, ExcFileNotOpen(config.mesh_file));
    Assert(dim == 2, ExcNotImplemented());

    grid_in.read_msh(input_file);
    validate_case_against_mesh();

    const double Re = 1.0 / viscosity;

    if (Re > 1000.0)
      {
        std::cout << "Searching for initial guess ..." << std::endl;
        const double step_size = 2000.0;
        compute_initial_guess(step_size);
        std::cout << "Found initial guess." << std::endl;
        std::cout << "Computing solution with target Re = " << Re << std::endl;
        viscosity = 1.0 / Re;
        newton_iteration(1e-12, 50, refinement, false, true);
      }
    else
      newton_iteration(1e-12, 50, refinement, true, true);

    if (config.use_mms)
      make_error_table();

  }

  template class ConjugateHeatTransferSolver<2>;
}
