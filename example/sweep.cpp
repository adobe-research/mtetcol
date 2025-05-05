#include <mtetcol/contour.h>
#include <mtetcol/simplicial_column.h>
#include <mtetcol/io.h>
#include <mtetcol/logger.h>

#include <ankerl/unordered_dense.h>
#include <mtet/io.h>

#include "grid.h"
#include "implicit_function.h"
#include "sweep_function.h"
#include "transform.h"

template <typename Scalar, typename Index>
std::tuple<std::vector<Scalar>, std::vector<Index>> generate_simpicial_column(
    const mtet::MTetMesh& tet_mesh)
{
    using VertexId = mtet::VertexId;
    using TetId = mtet::TetId;
    using IndexMap = ankerl::unordered_dense::map<uint64_t, size_t>;
    IndexMap vertex_tag_map;
    vertex_tag_map.reserve(tet_mesh.get_num_vertices());
    std::vector<Scalar> vertices;
    vertices.reserve(tet_mesh.get_num_vertices() * 3);
    std::vector<Index> tets;
    tets.reserve(tet_mesh.get_num_tets() * 4);

    tet_mesh.seq_foreach_vertex([&](VertexId vid, std::span<const mtet::Scalar, 3> data) {
        size_t vertex_tag = vertex_tag_map.size() + 1;
        vertex_tag_map[value_of(vid)] = vertex_tag;
        vertices.push_back(static_cast<Scalar>(data[0]));
        vertices.push_back(static_cast<Scalar>(data[1]));
        vertices.push_back(static_cast<Scalar>(data[2]));
    });

    tet_mesh.seq_foreach_tet([&](TetId, std::span<const VertexId, 4> data) {
        tets.push_back(static_cast<Index>(value_of(data[0])));
        tets.push_back(static_cast<Index>(value_of(data[1])));
        tets.push_back(static_cast<Index>(value_of(data[2])));
        tets.push_back(static_cast<Index>(value_of(data[3])));
    });

    return {vertices, tets};
}

void sample_time_derivative(
    mtetcol::SimplicialColumn<4>& column,
    const mtetcol::SweepFunction<3>& sweep_function,
    size_t num_time_samples)
{
    using Index = mtetcol::Index;
    using Scalar = mtetcol::Scalar;

    std::vector<Scalar> time_samples;
    for (size_t i = 0; i < num_time_samples; ++i) {
        time_samples.push_back(static_cast<Scalar>(i) / (num_time_samples - 1));
    }

    std::vector<Scalar> function_values;
    function_values.reserve(time_samples.size());

    column.set_time_samples(
        [&](Index) { return std::span<Scalar>(time_samples.data(), time_samples.size()); },
        [&](Index vi) {
            function_values.clear();
            auto pos = column.get_spatial_vertex(vi);
            for (const auto& t : time_samples) {
                function_values.push_back(
                    sweep_function.time_derivative({pos[0], pos[1], pos[2]}, t));
            }
            return std::span<Scalar>(function_values.data(), function_values.size());
        });
}

int main(int argc, char** argv)
{
    mtetcol::logger().set_level(spdlog::level::debug);

    auto tet_mesh = grid::generate_tet_mesh({10, 10, 10}, {0, 0, 0}, {1, 1, 1}, grid::TET5);
    mtet::save_mesh("grid.msh", tet_mesh);

    auto [vertices, tets] = generate_simpicial_column<mtetcol::Scalar, mtetcol::Index>(tet_mesh);
    mtetcol::SimplicialColumn<4> column;
    column.set_vertices(vertices);
    column.set_simplices(tets);

    constexpr size_t num_time_samples_per_vertex = 10;

    mtetcol::ImplicitSphere base_shape(0.2, {0.25, 0.25, 0.25});
    mtetcol::Translation<3> translation({0.5, 0.5, 0.6});
    mtetcol::SweepFunction<3> sweep_function(base_shape, translation);
    sample_time_derivative(column, sweep_function, num_time_samples_per_vertex);

    auto contour = column.extract_contour(0.0);
    contour.triangulate_cycles();

    size_t num_contour_vertices = contour.get_num_vertices();
    std::vector<mtetcol::Scalar> function_values(num_contour_vertices);
    for (size_t i=0; i< num_contour_vertices; ++i) {
        auto pos = contour.get_vertex(i);
        function_values[i] = sweep_function.value({pos[0], pos[1], pos[2]}, pos[3]);
    }
    auto isocontour = contour.isocontour(function_values);
    mtetcol::save_contour("tmp.obj", isocontour);

    isocontour.triangulate_cycles();

    return 0;
}
