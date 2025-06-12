#include <mtetcol/contour.h>
#include <mtetcol/io.h>
#include <mtetcol/logger.h>

#include <mtet/grid.h>
#include <mtet/io.h>

#include <stf/stf.h>

#include <iostream>
#include <vector>

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <sample_file> <coeff_file>" << std::endl;
        return 1;
    }
    mtetcol::logger().set_level(spdlog::level::debug);

    constexpr size_t resolution = 64;
    auto tet_mesh = mtet::generate_tet_grid(
        {resolution, resolution, resolution},
        {0, 0, 0},
        {1, 1, 1},
        mtet::TET5);

    stf::Duchon fn(argv[1], argv[2], {0.5, 0.5, 0.5}, 0.4, true);

    std::vector<stf::Scalar> values;
    values.reserve(tet_mesh.get_num_vertices());

    tet_mesh.seq_foreach_vertex([&](mtet::VertexId vid, std::span<const mtet::Scalar, 3> data) {
        values.push_back(fn.value({data[0], data[1], data[2]}));
    });

    mtet::save_mesh("grid.msh", tet_mesh, "value", {values.data(), values.size()});

    return 0;
}
