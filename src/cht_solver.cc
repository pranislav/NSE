#include "case_config.h"
#include "conjugate_heat_transfer_solver.h"

#include <deal.II/base/exceptions.h>

#include <iostream>
#include <string>

int main(int argc, char **argv)
{
  try
    {
      using namespace Cht;

      std::string output_directory = ".";
      bool        save_mesh_output = false;
      std::string case_file        = "../cases/heat_exchanger.prm";

      const auto print_help = [&]() {
        std::cout << "Usage: cht_solver [--case FILE] [--output-dir DIR] "
                     "[--save-mesh] [--help]\n"
                  << "  --case FILE        Case file to read\n"
                  << "  --output-dir DIR   Output directory (default: .)\n"
                  << "  --save-mesh        Save mesh after each refinement\n"
                  << "  --help             Show this message\n";
      };

      for (int i = 1; i < argc; ++i)
        {
          const std::string argument = argv[i];

          if (argument == "--output-dir")
            {
              AssertThrow(i + 1 < argc,
                          dealii::ExcMessage("Missing value for --output-dir."));
              output_directory = argv[++i];
            }
          else if (argument == "--case")
            {
              AssertThrow(i + 1 < argc,
                          dealii::ExcMessage("Missing value for --case."));
              case_file = argv[++i];
            }
          else if (argument == "--save-mesh")
            save_mesh_output = true;
          else if (argument == "--help")
            {
              print_help();
              return 0;
            }
          else
            AssertThrow(false,
                        dealii::ExcMessage("Unknown command line argument: " +
                                           argument));
        }

      const CaseConfig case_config = read_case_config(case_file);

      ConjugateHeatTransferSolver<2> solver(case_config,
                                            output_directory,
                                            save_mesh_output);
      solver.run(case_config.adaptive_refinement_cycles);
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
