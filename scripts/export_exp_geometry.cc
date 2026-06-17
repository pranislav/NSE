#include <deal.II/base/point.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>

#include <fstream>
#include <map>
#include <set>

using namespace dealii;

int main()
{
    constexpr int dim = 2;

    Triangulation<dim> tria;

    {
        GridIn<dim> grid_in;
        grid_in.attach_triangulation(tria);

        std::ifstream in("../meshes/experiment_layout.msh");
        grid_in.read_msh(in);
    }

    std::ofstream out("material_interfaces.vtk");

    std::vector<std::array<Point<dim>, 2>> segments;
    std::set<std::pair<unsigned int, unsigned int>> visited;

    for (const auto &cell : tria.active_cell_iterators())
    {
        for (unsigned int f = 0;
             f < GeometryInfo<dim>::faces_per_cell;
             ++f)
        {
            if (cell->at_boundary(f))
                continue;

            const auto neighbor = cell->neighbor(f);

            if (cell->material_id() == neighbor->material_id())
                continue;

            const auto key =
                std::minmax(cell->active_cell_index(),
                            neighbor->active_cell_index());

            if (!visited.insert(key).second)
                continue;

            auto face = cell->face(f);

            segments.push_back(
                {face->vertex(0), face->vertex(1)});
        }
    }

    const unsigned int n_points = 2 * segments.size();

    out << "# vtk DataFile Version 2.0\n";
    out << "Material interfaces\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";

    out << "POINTS " << n_points << " double\n";

    for (const auto &s : segments)
    {
        out << s[0][0] << " " << s[0][1] << " 0\n";
        out << s[1][0] << " " << s[1][1] << " 0\n";
    }

    out << "\nLINES "
        << segments.size() << " "
        << 3 * segments.size() << "\n";

    for (unsigned int i = 0; i < segments.size(); ++i)
        out << "2 " << 2 * i << " " << 2 * i + 1 << "\n";

    std::cout << "Wrote "
              << segments.size()
              << " interface segments to material_interfaces.vtk\n";
}