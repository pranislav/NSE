/*
 * Simulation of a viscose incompressible fluid flow in a tube
 * based on step-57 of dealii tutorial
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/exceptions.h>
 
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
 
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
 
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
 
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
 
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
 
#include <deal.II/numerics/solution_transfer.h>
 
#include <deal.II/lac/sparse_direct.h>
 
#include <deal.II/lac/sparse_ilu.h>
 
#include <deal.II/grid/grid_in.h>
 
#include <filesystem>
#include <fstream>
#include <regex>
#include <iostream>
 
namespace Step57
{
  using namespace dealii;

  struct InflowsParameters
  {
    double pipe_len;
    double pipe_diameter;
    double lower_inflow_height;
    double higher_inflow_height;
  };

  InflowsParameters read_inflows_parameters_from_geo(const std::string &geo_file)
  {
    std::ifstream in(geo_file);
    AssertThrow(in, ExcFileNotOpen(geo_file));

    const std::regex parameter_regex(
      R"(^\s*(length|h_fluid|h_insul|h_membrane)\s*=\s*([0-9]+(?:\.[0-9]*)?(?:[eE][-+]?[0-9]+)?)\s*;\s*$)");

    double length = 0.0;
    double h_fluid = 0.0;
    double h_insul = 0.0;
    double h_membrane = 0.0;

    std::string line;
    while (std::getline(in, line))
      {
        std::smatch match;
        if (!std::regex_match(line, match, parameter_regex))
          continue;

        const double value = std::stod(match[2].str());
        const std::string name = match[1].str();
        if (name == "length")
          length = value;
        else if (name == "h_fluid")
          h_fluid = value;
        else if (name == "h_insul")
          h_insul = value;
        else if (name == "h_membrane")
          h_membrane = value;
      }

    AssertThrow(length > 0.0, ExcMessage("Missing length in geo file."));
    AssertThrow(h_fluid > 0.0, ExcMessage("Missing h_fluid in geo file."));
    AssertThrow(h_insul > 0.0, ExcMessage("Missing h_insul in geo file."));
    AssertThrow(h_membrane >= 0.0,
                ExcMessage("Missing h_membrane in geo file."));

    return {length, h_fluid, h_insul, h_insul + h_fluid + h_membrane};
  }
 
 
  template <int dim>
  class StationaryNavierStokes
  {
  public:
    StationaryNavierStokes(const unsigned int degree,
                           const std::string &output_directory = ".",
                           const bool         save_mesh_output = false);
    void run(const unsigned int refinement);
 
  private:

    enum
    {
      fluid_material_id,
      insulator_material_id,
      conductor_material_id
    };

    static bool cell_is_in_fluid_domain(
      const typename DoFHandler<dim>::cell_iterator &cell);

    static bool cell_is_in_solid_domain(
      const typename DoFHandler<dim>::cell_iterator &cell);

    void set_active_fe_indices();
    void add_fluid_solid_interface_constraints(
      AffineConstraints<double> &constraints) const;
    void setup_temperature_dofs();
    void assemble_temperature_system();
    void solve_temperature();
    void update_temperature_field();

    void setup_dofs();
 
    void initialize_system();
 
    void assemble(const bool initial_step, const bool assemble_matrix);
 
    void assemble_system(const bool initial_step);
 
    void assemble_rhs(const bool initial_step);
 
    void solve(const bool initial_step);
 
    void refine_mesh();
 
    void process_solution(unsigned int refinement);
 
    void output_results(const unsigned int refinement_cycle) const;
    void output_mesh(const unsigned int output_index) const;
 
    void newton_iteration(const double       tolerance,
                          const unsigned int max_n_line_searches,
                          const unsigned int max_n_refinements,
                          const bool         is_initial_step,
                          const bool         output_result);
 
    void compute_initial_guess(double step_size);

    double solid_thermal_conductivity(dealii::types::material_id material_id);
 
    double                               viscosity;
    double                               gamma;
    const double                         fluid_thermal_conductivity;
    const double                         conductor_thermal_conductivity;
    const double                         insulator_thermal_conductivity;
    const InflowsParameters              inflows_params;
    const unsigned int                   degree;
    const std::string                    output_directory;
    const bool                           save_mesh_output;
    std::vector<types::global_dof_index> dofs_per_block;
 
    Triangulation<dim>    triangulation;
    const FESystem<dim>   fe_fluid;
    const FESystem<dim>   fe_solid;
    const FE_Q<dim>       temperature_fe;
    hp::FECollection<dim> fe_collection;
    DoFHandler<dim>       dof_handler;
    DoFHandler<dim>       temperature_dof_handler;
 
    AffineConstraints<double> zero_constraints;
    AffineConstraints<double> nonzero_constraints;
    AffineConstraints<double> temperature_constraints;
 
    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;
    SparseMatrix<double>      pressure_mass_matrix;
    SparsityPattern           temperature_sparsity_pattern;
    SparseMatrix<double>      temperature_matrix;
 
    BlockVector<double> present_solution;
    BlockVector<double> newton_update;
    BlockVector<double> system_rhs;
    BlockVector<double> evaluation_point;
    Vector<double>      temperature_solution;
    Vector<double>      temperature_rhs;
  };
 
 
  template <int dim>
  class VelocityBoundaryValues : public Function<dim>
  {
  public:
    VelocityBoundaryValues(const InflowsParameters &inflows_params)
      : Function<dim>(dim + 1)
      , inflows_params(inflows_params)
    {}
    virtual double value(const Point<dim>  &p,
                         const unsigned int component) const override;

  private:
    const InflowsParameters &inflows_params;
  };

    template <int dim>
    double VelocityBoundaryValues<dim>::value(const Point<dim> &p,
                                      const unsigned int component) const
      {
        const double y = p[1];
        if (std::abs(p[0]) < 1e-12) // lower pipe
          {
          const double z0 = inflows_params.lower_inflow_height;
          const double z1 = z0 + inflows_params.pipe_diameter;
          if (component == 0 && y >= z0 && y <= z1)
            return 10 * (y - z0) * (z1 - y);
          }
        else if (std::abs(p[0] - inflows_params.pipe_len) < 1e-12)
          {
          const double z0 = inflows_params.higher_inflow_height;
          const double z1 = z0 + inflows_params.pipe_diameter;
          if (component == 0 && y >= z0 && y <= z1)
            return -10 * (y - z0) * (z1 - y);
          }
        return 0;
      }

  template <int dim>
  class TemperatureBoundaryValues : public Function<dim>
  {
  public:
    TemperatureBoundaryValues(const double temperature)
      : Function<dim>(1)
      , temperature(temperature)
    {}

    virtual double value(const Point<dim> & /*p*/,
                         const unsigned int component) const override
    {
      AssertIndexRange(component, 1);
      return temperature;
    }

  private:
    const double temperature;
  };

  template <int dim>
  class LocalizedInletTemperatureBoundaryValues : public Function<dim>
  {
  public:
    LocalizedInletTemperatureBoundaryValues()
      : Function<dim>(1)
    {}

    virtual double value(const Point<dim>  &p,
                         const unsigned int component) const override
    {
      AssertIndexRange(component, 1);

      const double center = 1.5;
      const double half_width = 0.15;
      const double distance = std::abs(p[1] - center);

      if (distance >= half_width)
        return 0.0;

      const double phase = numbers::PI * distance / half_width;
      return 0.5 * (1.0 + std::cos(phase));
    }
  };
 
  // {
  //   Assert(component < this->n_components,
  //          ExcIndexRange(component, 0, this->n_components));
  //   if (component == 0) // x-velocity
  //   {
  //     if (std::abs(p[0]) < 1e-12) // inlet
  //     {
  //       const double y = p[1];
  //       return y * (1.0 - y); // Poiseuille profile
  //     }
  //   }

  //   return 0.0;
  // }
  // {
  //   Assert(component < this->n_components,
  //          ExcIndexRange(component, 0, this->n_components));
  //   if (component == 0 && std::abs(p[dim - 1] - 1.0) < 1e-10)
  //     return 1.0;

  //   return 0;
  // }
 
  template <class PreconditionerMp>
  class BlockSchurPreconditioner : public EnableObserverPointer
  {
  public:
    BlockSchurPreconditioner(double                           gamma,
                             double                           viscosity,
                             const BlockSparseMatrix<double> &S,
                             const SparseMatrix<double>      &P,
                             const PreconditionerMp          &Mppreconditioner);
 
    void vmult(BlockVector<double> &dst, const BlockVector<double> &src) const;
 
  private:
    const double                     gamma;
    const double                     viscosity;
    const BlockSparseMatrix<double> &stokes_matrix;
    const SparseMatrix<double>      &pressure_mass_matrix;
    const PreconditionerMp          &mp_preconditioner;
    SparseDirectUMFPACK              A_inverse;
  };
 
 
  template <class PreconditionerMp>
  BlockSchurPreconditioner<PreconditionerMp>::BlockSchurPreconditioner(
    double                           gamma,
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
 
  template <class PreconditionerMp>
  void BlockSchurPreconditioner<PreconditionerMp>::vmult(
    BlockVector<double>       &dst,
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
 
  template <int dim>
  StationaryNavierStokes<dim>::StationaryNavierStokes(
    const unsigned int degree,
    const std::string &output_directory,
    const bool         save_mesh_output)
    : viscosity(1.0 / 7500.0)
    , gamma(1.0)
    , fluid_thermal_conductivity(1.0)
    , conductor_thermal_conductivity(0.5)
    , insulator_thermal_conductivity(5.0)
    , inflows_params(
        read_inflows_parameters_from_geo("../meshes/thermal_exchanger.geo"))
    , degree(degree)
    , output_directory(output_directory)
    , save_mesh_output(save_mesh_output)
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
  bool StationaryNavierStokes<dim>::cell_is_in_fluid_domain(
    const typename DoFHandler<dim>::cell_iterator &cell)
  {
    return (cell->material_id() == fluid_material_id);
  }


  template <int dim>
  bool StationaryNavierStokes<dim>::cell_is_in_solid_domain(
    const typename DoFHandler<dim>::cell_iterator &cell)
  {
    return (cell->material_id() == insulator_material_id ||
            cell->material_id() == conductor_material_id);
  }

  template <int dim>
  void StationaryNavierStokes<dim>::set_active_fe_indices()
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
  void StationaryNavierStokes<dim>::add_fluid_solid_interface_constraints(
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

                  for (unsigned int i = 0; i < local_face_dof_indices.size();
                       ++i)
                    if (fe_fluid.face_system_to_component_index(i).first < dim)
                      constraints.constrain_dof_to_zero(
                        local_face_dof_indices[i]);
                }
            }
  }

  template <int dim>
  void StationaryNavierStokes<dim>::setup_temperature_dofs()
  {
    temperature_matrix.clear();

    temperature_dof_handler.distribute_dofs(temperature_fe);

    temperature_constraints.clear();
    DoFTools::make_hanging_node_constraints(temperature_dof_handler,
                                            temperature_constraints);
    VectorTools::interpolate_boundary_values(temperature_dof_handler,
                                             30,
                                             TemperatureBoundaryValues<dim>(1.0),
                                             temperature_constraints);
    VectorTools::interpolate_boundary_values(temperature_dof_handler,
                                             40,
                                             TemperatureBoundaryValues<dim>(0.0),
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
  void StationaryNavierStokes<dim>::setup_dofs()
  {
    system_matrix.clear();
    pressure_mass_matrix.clear();

    set_active_fe_indices();
 
    dof_handler.distribute_dofs(fe_collection);
 
    std::vector<unsigned int> block_component(dim + 1, 0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise(dof_handler, block_component);
 
    dofs_per_block =
      DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
    unsigned int dof_u = dofs_per_block[0];
    unsigned int dof_p = dofs_per_block[1];
 
    const FEValuesExtractors::Vector velocities(0);
    {
      nonzero_constraints.clear();
 
      DoFTools::make_hanging_node_constraints(dof_handler, nonzero_constraints);

      for (int id: {30, 40}) {
        VectorTools::interpolate_boundary_values(
          dof_handler,
          id,
          VelocityBoundaryValues<dim>(inflows_params),
          nonzero_constraints,
          fe_collection.component_mask(velocities));
        }

      add_fluid_solid_interface_constraints(nonzero_constraints);

    }
    nonzero_constraints.close();
 
    {
      zero_constraints.clear();
 
      DoFTools::make_hanging_node_constraints(dof_handler, zero_constraints);

      for (int id: {30, 40}) {
        VectorTools::interpolate_boundary_values(
          dof_handler,
          id,
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

    setup_temperature_dofs();
  }
 
  template <int dim>
  void StationaryNavierStokes<dim>::initialize_system()
  {
    {
      BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);
      DoFTools::make_sparsity_pattern(dof_handler, dsp, nonzero_constraints);
      sparsity_pattern.copy_from(dsp);
    }
 
    system_matrix.reinit(sparsity_pattern);
 
    present_solution.reinit(dofs_per_block);
    newton_update.reinit(dofs_per_block);
    system_rhs.reinit(dofs_per_block);
  }
 
  template <int dim>
  void StationaryNavierStokes<dim>::assemble(const bool initial_step,
                                             const bool assemble_matrix)
  {
    if (assemble_matrix)
      system_matrix = 0;
 
    system_rhs = 0;
 
    const UpdateFlags update_flags =
      update_values | update_quadrature_points |
      update_JxW_values | update_gradients;

    hp::QCollection<dim> q_collection;
    q_collection.push_back(QGauss<dim>(degree + 2));
    q_collection.push_back(QGauss<dim>(degree + 2));

    hp::FEValues<dim> hp_fe_values(fe_collection,
                                  q_collection,
                                  update_flags);
  
    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar pressure(dim);  

 
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        hp_fe_values.reinit(cell);
        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
        const unsigned int n_q_points = fe_values.n_quadrature_points;

        std::vector<Tensor<1, dim>> present_velocity_values(n_q_points);
        std::vector<Tensor<2, dim>> present_velocity_gradients(n_q_points);
        std::vector<double>         present_pressure_values(n_q_points);

        const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();

        FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
        Vector<double>     local_rhs(dofs_per_cell);
        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

        std::vector<double>         div_phi_u(dofs_per_cell);
        std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
        std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell);
        std::vector<double>         phi_p(dofs_per_cell);
 
        local_matrix = 0;
        local_rhs    = 0;
 
        fe_values[velocities].get_function_values(evaluation_point,
                                                  present_velocity_values);
 
        fe_values[velocities].get_function_gradients(
          evaluation_point, present_velocity_gradients);
 
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
 
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                if (assemble_matrix)
                  {
                    for (unsigned int j = 0; j < dofs_per_cell; ++j)
                      {
                        local_matrix(i, j) +=
                          (viscosity *
                             scalar_product(grad_phi_u[i], grad_phi_u[j]) +
                           phi_u[i] *
                             (present_velocity_gradients[q] * phi_u[j]) +
                           phi_u[i] *
                             (grad_phi_u[j] * present_velocity_values[q]) -
                           div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j] +
                           gamma * div_phi_u[i] * div_phi_u[j] +
                           phi_p[i] * phi_p[j]) *
                          fe_values.JxW(q);
                      }
                  }
 
                double present_velocity_divergence =
                  trace(present_velocity_gradients[q]);
                local_rhs(i) +=
                  (-viscosity * scalar_product(grad_phi_u[i],
                                               present_velocity_gradients[q]) -
                   phi_u[i] * (present_velocity_gradients[q] *
                               present_velocity_values[q]) +
                   div_phi_u[i] * present_pressure_values[q] +
                   phi_p[i] * present_velocity_divergence -
                   gamma * div_phi_u[i] * present_velocity_divergence) *
                  fe_values.JxW(q);
              }
          }
 
        cell->get_dof_indices(local_dof_indices);
 
        const AffineConstraints<double> &constraints_used =
          initial_step ? nonzero_constraints : zero_constraints;
 
        if (assemble_matrix)
          {
            constraints_used.distribute_local_to_global(local_matrix,
                                                        local_rhs,
                                                        local_dof_indices,
                                                        system_matrix,
                                                        system_rhs);
          }
        else
          {
            constraints_used.distribute_local_to_global(local_rhs,
                                                        local_dof_indices,
                                                        system_rhs);
          }
      }
 
    if (assemble_matrix)
      {
        pressure_mass_matrix.reinit(sparsity_pattern.block(1, 1));
        pressure_mass_matrix.copy_from(system_matrix.block(1, 1));
 
        system_matrix.block(1, 1) = 0;
      }
  }
 
  template <int dim>
  void StationaryNavierStokes<dim>::assemble_system(const bool initial_step)
  {
    assemble(initial_step, true);
  }
 
  template <int dim>
  void StationaryNavierStokes<dim>::assemble_rhs(const bool initial_step)
  {
    assemble(initial_step, false);
  }
 
  template <int dim>
  void StationaryNavierStokes<dim>::solve(const bool initial_step)
  {
    const AffineConstraints<double> &constraints_used =
      initial_step ? nonzero_constraints : zero_constraints;
 
    SolverControl solver_control(system_matrix.m(),
                                 1e-4 * system_rhs.l2_norm(),
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
 
    gmres.solve(system_matrix, newton_update, system_rhs, preconditioner);
    std::cout << "FGMRES steps: " << solver_control.last_step() << std::endl;
 
    constraints_used.distribute(newton_update);
  }

  template <int dim>
  double StationaryNavierStokes<dim>::
    solid_thermal_conductivity(dealii::types::material_id material_id)
  {
    if (material_id == 1)
      return insulator_thermal_conductivity;
    if(material_id == 2)
      return conductor_thermal_conductivity;
    throw std::runtime_error("unexpected material id in solid domain");
  }

  template <int dim>
  void StationaryNavierStokes<dim>::assemble_temperature_system()
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

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     local_rhs(dofs_per_cell);
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
        const double thermal_conductivity =
          in_fluid_domain ? fluid_thermal_conductivity :
                            solid_thermal_conductivity(cell->material_id());

        if (in_fluid_domain)
          {
            flow_fe_values.reinit(flow_cell);
            flow_fe_values.get_present_fe_values()[velocities].get_function_values(
              present_solution, velocity_values);
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
                    (thermal_conductivity * grad_phi_T[i] * grad_phi_T[j] *
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
  void StationaryNavierStokes<dim>::solve_temperature()
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
  void StationaryNavierStokes<dim>::update_temperature_field()
  {
    assemble_temperature_system();
    solve_temperature();
  }
 
  template <int dim>
  void StationaryNavierStokes<dim>::refine_mesh()
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
      present_solution,
      flow_error_per_cell,
      fe_collection.component_mask(velocity));

    if (temperature_solution.size() == temperature_dof_handler.n_dofs())
      KellyErrorEstimator<dim>::estimate(
        temperature_dof_handler,
        QGauss<dim - 1>(degree + 1),
        std::map<types::boundary_id, const Function<dim> *>(),
        temperature_solution,
        temperature_error_per_cell);

    const float max_flow_error = flow_error_per_cell.linfty_norm();
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
 
    triangulation.prepare_coarsening_and_refinement();
    SolutionTransfer<dim, BlockVector<double>> solution_transfer(dof_handler);
    solution_transfer.prepare_for_coarsening_and_refinement(present_solution);
    triangulation.execute_coarsening_and_refinement();
 
    setup_dofs();
 
    BlockVector<double> tmp(dofs_per_block);
 
    solution_transfer.interpolate(tmp);
    nonzero_constraints.distribute(tmp);
 
    initialize_system();
    present_solution = tmp;
  }
 
  template <int dim>
  void StationaryNavierStokes<dim>::newton_iteration(
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
                evaluation_point = present_solution;
                assemble_system(first_step);
                solve(first_step);
                present_solution = newton_update;
                nonzero_constraints.distribute(present_solution);
                first_step       = false;
                evaluation_point = present_solution;
                assemble_rhs(first_step);
                current_res = system_rhs.l2_norm();
                std::cout << "The residual of initial guess is " << current_res
                          << std::endl;
                last_res = current_res;
              }
            else
              {
                evaluation_point = present_solution;
                assemble_system(first_step);
                solve(first_step);
 
                for (double alpha = 1.0; alpha > 1e-5; alpha *= 0.5)
                  {
                    evaluation_point = present_solution;
                    evaluation_point.add(alpha, newton_update);
                    nonzero_constraints.distribute(evaluation_point);
                    assemble_rhs(first_step);
                    current_res = system_rhs.l2_norm();
                    std::cout << "  alpha: " << std::setw(10) << alpha
                              << std::setw(0) << "  residual: " << current_res
                              << std::endl;
                    if (current_res < last_res)
                      break;
                  }
                {
                  present_solution = evaluation_point;
                  std::cout << "  number of line searches: " << line_search_n
                            << "  residual: " << current_res << std::endl;
                  last_res = current_res;
                }
                ++line_search_n;
              }
 
            if (output_result)
              {
                update_temperature_field();
                output_results(max_n_line_searches * refinement_n +
                               line_search_n);
 
                if (current_res <= tolerance)
                  process_solution(refinement_n);
              }
          }
 
        if (refinement_n < max_n_refinements)
          {
            refine_mesh();
          }
      }
  }
 
  template <int dim>
  void StationaryNavierStokes<dim>::compute_initial_guess(double step_size)
  {
    const double target_Re = 1.0 / viscosity;
 
    bool is_initial_step = true;
 
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
  void StationaryNavierStokes<dim>::output_results(
    const unsigned int output_index) const
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
    data_out.add_data_vector(present_solution,
                             solution_names,
                             DataOut<dim>::type_dof_data,
                             data_component_interpretation);
    if (temperature_solution.size() == temperature_dof_handler.n_dofs() &&
        temperature_dof_handler.n_dofs() > 0)
      data_out.add_data_vector(temperature_dof_handler,
                               temperature_solution,
                               "temperature");
    data_out.build_patches();
 
    std::ofstream output(output_directory + "/" +
                         std::to_string(1.0 / viscosity) + "-solution-" +
                         Utilities::int_to_string(output_index, 4) + ".vtk");
    data_out.write_vtk(output);

    if (save_mesh_output)
      output_mesh(output_index);
  }

  template <int dim>
  void StationaryNavierStokes<dim>::output_mesh(
    const unsigned int output_index) const
  {
    GridOut grid_out;
    std::ofstream output(output_directory + "/" +
                         std::to_string(1.0 / viscosity) + "-mesh-" +
                         Utilities::int_to_string(output_index, 4) + ".vtu");
    grid_out.write_vtu(triangulation, output);
  }
 
  template <int dim>
  void StationaryNavierStokes<dim>::process_solution(unsigned int refinement)
  {
    std::filesystem::create_directories(output_directory);

    std::ofstream f(output_directory + "/" + std::to_string(1.0 / viscosity) +
                    "-line-" + std::to_string(refinement) + ".txt");
    f << "# y u_x u_y" << std::endl;
 
    Point<dim> p;
    p[0] = 0.5;
    p[1] = 0.5;
 
    f << std::scientific;
 
    for (unsigned int i = 0; i <= 100; ++i)
      {
        p[dim - 1] = i / 100.0;
 
        Vector<double> tmp_vector(dim + 1);
        VectorTools::point_value(dof_handler, present_solution, p, tmp_vector);
        f << p[dim - 1];
 
        for (int j = 0; j < dim; ++j)
          f << ' ' << tmp_vector(j);
        f << std::endl;
      }
  }

  template <int dim>
  void set_pipe_boundary_ids(Triangulation<dim> &triangulation)
  {
    for (auto &cell : triangulation.active_cell_iterators())
      for (auto f : cell->face_iterators())
        if (f->at_boundary())
        {
          const auto center = f->center();

          if (std::abs(center[0]) < 1e-12)
            f->set_boundary_id(40); // inlet
          else if (std::abs(center[0] - 1.0) < 1e-12)
            f->set_boundary_id(20); // outlet
          else
            f->set_boundary_id(10); // walls
        }
  }

  template <int dim>
  void StationaryNavierStokes<dim>::run(const unsigned int refinement)
  {
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(triangulation);
    std::ifstream input_file("../meshes/thermal_exchanger.msh");
    Assert(dim == 2, ExcNotImplemented());

    grid_in.read_msh(input_file);

    // GridGenerator::hyper_cube(triangulation);
    // triangulation.refine_global(5);
    // set_pipe_boundary_ids(triangulation);
  
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
      {
 
        newton_iteration(1e-12, 50, refinement, true, true);
      }
  }
} // namespace Step57
 
int main(int argc, char **argv)
{
  try
    {
      using namespace Step57;

      std::string output_directory = ".";
      bool        save_mesh_output = false;

      for (int i = 1; i < argc; ++i)
        {
          const std::string argument = argv[i];

          if (argument == "--output-dir")
            {
              AssertThrow(i + 1 < argc,
                          ExcMessage("Missing value for --output-dir."));
              output_directory = argv[++i];
            }
          else if (argument == "--save-mesh")
            save_mesh_output = true;
          else
            AssertThrow(false,
                        ExcMessage("Unknown command line argument: " +
                                   argument));
        }
 
      StationaryNavierStokes<2> flow(/* degree = */ 1,
                                     output_directory,
                                     save_mesh_output);
      flow.run(4);
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}
