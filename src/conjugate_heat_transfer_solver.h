#ifndef STEP57_CONJUGATE_HEAT_TRANSFER_SOLVER_H
#define STEP57_CONJUGATE_HEAT_TRANSFER_SOLVER_H

#include "case_config.h"

#include <deal.II/base/function.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/base/convergence_table.h>

#include <array>
#include <map>
#include <string>
#include <vector>

namespace Cht
{
  template <int dim>
  class ConjugateHeatTransferSolver
  {
  public:
    enum RefinementMode
    {
      global_refinement,
      adaptive_refinement
    };

    ConjugateHeatTransferSolver(const CaseConfig  &config,
                                const std::string &output_directory,
                                const bool         save_mesh_output,
                                const RefinementMode refinement_mode);

    void run(const unsigned int refinement);

  private:
    struct BoundaryExtent
    {
      BoundaryExtent();

      std::array<double, dim> min;
      std::array<double, dim> max;
      bool                    initialized;
    };

    const MaterialData &material_data(const dealii::types::material_id material_id) const;
    bool cell_is_in_fluid_domain(
      const typename dealii::DoFHandler<dim>::cell_iterator &cell) const;
    bool cell_is_in_solid_domain(
      const typename dealii::DoFHandler<dim>::cell_iterator &cell) const;

    void set_active_fe_indices();
    void add_fluid_solid_interface_constraints(
      dealii::AffineConstraints<double> &constraints) const;
    void validate_case_against_mesh() const;
    void collect_boundary_extents();
    void setup_temperature_dofs();
    void setup_flow_dofs();
    void setup_dofs();
    void assemble_temperature_system();
    void solve_temperature();
    void update_temperature_field();
    void initialize_system();
    void assemble_nse(const bool initial_step, const bool assemble_matrix);
    void assemble_flow_system(const bool initial_step);
    void assemble_flow_rhs(const bool initial_step);
    void solve_flow(const bool initial_step);
    void mark_cells_for_adaptive_refinement();
    void mark_cells_for_global_refinement();
    void execute_refinement();
    void refine_mesh(unsigned int refinement_cycle);
    void process_solution(unsigned int refinement);
    void compute_errors(const unsigned int cycle);
    void output_results(const unsigned int refinement_cycle,
                        const unsigned int newton_step) const;
    void output_mesh(const unsigned int refinement_cycle) const;
    void newton_iteration(const double       tolerance,
                          const unsigned int max_n_line_searches,
                          const unsigned int max_n_refinements,
                          const bool         is_initial_step,
                          const bool         output_result);
    void compute_initial_guess(double step_size);
    std::string case_tag() const;
    static std::string format_compact_double(double value);

    const CaseConfig                     config;
    const RefinementMode                 refinement_mode;
    double                               viscosity;
    double                               gamma;
    const unsigned int                   degree;
    const std::string                    output_directory;
    const bool                           save_mesh_output;
    std::vector<dealii::types::global_dof_index> dofs_per_block;
    std::map<dealii::types::boundary_id, BoundaryExtent> boundary_extents;

    dealii::Triangulation<dim>    triangulation;
    const dealii::FESystem<dim>   fe_fluid;
    const dealii::FESystem<dim>   fe_solid;
    const dealii::FE_Q<dim>       temperature_fe;
    dealii::hp::FECollection<dim> fe_collection;
    dealii::DoFHandler<dim>       dof_handler;
    dealii::DoFHandler<dim>       temperature_dof_handler;

    dealii::AffineConstraints<double> zero_constraints;
    dealii::AffineConstraints<double> nonzero_constraints;
    dealii::AffineConstraints<double> temperature_constraints;

    dealii::BlockSparsityPattern      sparsity_pattern;
    dealii::BlockSparseMatrix<double> system_matrix;
    dealii::SparseMatrix<double>      pressure_mass_matrix;
    dealii::SparsityPattern           temperature_sparsity_pattern;
    dealii::SparseMatrix<double>      temperature_matrix;

    dealii::BlockVector<double> flow_solution;
    dealii::BlockVector<double> newton_update;
    dealii::BlockVector<double> flow_rhs;
    dealii::BlockVector<double> evaluation_point;
    dealii::Vector<double>      temperature_solution;
    dealii::Vector<double>      temperature_rhs;

    dealii::ConvergenceTable convergence_table;
  };
}

#endif
