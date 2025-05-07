#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>

#include <mtetcol/contour.h>
#include <mtetcol/implicit_function.h>
#include <mtetcol/logger.h>
#include <mtetcol/simplicial_column.h>
#include <mtetcol/sweep_function.h>
#include <mtetcol/transform.h>

#include <ankerl/unordered_dense.h>
#include <mtet/grid.h>

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

auto sample_time_derivative(
    mtetcol::SimplicialColumn<4>& columns,
    const mtetcol::SweepFunction<3>& sweep_function,
    size_t num_time_samples) -> std::
    tuple<std::vector<mtetcol::Scalar>, std::vector<mtetcol::Scalar>, std::vector<mtetcol::Index>>
{
    using Index = mtetcol::Index;
    using Scalar = mtetcol::Scalar;

    std::vector<Scalar> time_samples;
    std::vector<Scalar> function_values;
    std::vector<Index> vertex_start_indices;

    const size_t num_vertices = columns.get_num_spatial_vertices();
    time_samples.reserve(num_time_samples * num_vertices);
    function_values.reserve(num_time_samples * num_vertices);
    vertex_start_indices.reserve(num_vertices + 1);
    vertex_start_indices.push_back(0);

    for (size_t vi = 0; vi < num_vertices; vi++) {
        auto pos = columns.get_spatial_vertex(vi);
        for (size_t i = 0; i < num_time_samples; ++i) {
            Scalar t = static_cast<Scalar>(i) / (num_time_samples - 1);
            time_samples.push_back(t);
            function_values.push_back(sweep_function.time_derivative({pos[0], pos[1], pos[2]}, t));
        }
        vertex_start_indices.push_back(static_cast<Index>(time_samples.size()));
    }

    return {time_samples, function_values, vertex_start_indices};
}


TEST_CASE("benchmark", "[mtetcol][.benchmark]")
{
    constexpr size_t resolution = 64;
    constexpr size_t num_time_samples_per_vertex = 64;

    mtetcol::logger().set_level(spdlog::level::warn);

    auto grid = mtet::generate_tet_grid(
        {resolution, resolution, resolution},
        {0, 0, 0},
        {1, 1, 1},
        mtet::TET5);
    auto [vertices, tets] = generate_simpicial_column<mtetcol::Scalar, mtetcol::Index>(grid);
    mtetcol::SimplicialColumn<4> columns;
    columns.set_vertices(vertices);
    columns.set_simplices(tets);

    mtetcol::ImplicitTorus base_shape(0.2, 0.04, {0.25, 0.5, 0.5});
    mtetcol::Rotation<3> rotation({0.25, 0.5, 0.5}, {1, 0, 0});
    mtetcol::Translation<3> translation({-0.5, 0, 0});
    mtetcol::Compose<3> flip(translation, rotation);
    mtetcol::SweepFunction<3> sweep_function(base_shape, flip);

    auto [time_samples, function_values, vertex_start_indices] =
        sample_time_derivative(columns, sweep_function, num_time_samples_per_vertex);

    columns.set_time_samples(time_samples, function_values, vertex_start_indices);

    BENCHMARK("Extract time derivative contour")
    {
        return columns.extract_contour(0.0, false);
    };
}
