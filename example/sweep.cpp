#include <mtetcol/contour.h>
#include <mtetcol/implicit_function.h>
#include <mtetcol/io.h>
#include <mtetcol/logger.h>
#include <mtetcol/simplicial_column.h>
#include <mtetcol/sweep_function.h>
#include <mtetcol/transform.h>

#include <ankerl/unordered_dense.h>
#include <mtet/grid.h>
#include <mtet/io.h>

#include <chrono>

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
    size_t num_time_samples,
    bool cyclic = false) -> std::
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

    if (cyclic) {
        // For cyclic time samples, we need to set the last time sample to the first one
        for (size_t vi = 0; vi < num_vertices; vi++) {
            function_values[vertex_start_indices[vi + 1] - 1] =
                function_values[vertex_start_indices[vi]];
        }
    }

    return {time_samples, function_values, vertex_start_indices};
}

mtetcol::Contour<4> sphere_translation(mtetcol::SimplicialColumn<4>& columns)
{
    using Scalar = mtetcol::Scalar;
    using Index = mtetcol::Index;

    constexpr size_t num_time_samples_per_vertex = 10;

    mtetcol::ImplicitSphere base_shape(0.2, {0.25, 0.5, 0.5});
    mtetcol::Translation<3> translation({-0.5, 0, 0});
    mtetcol::SweepFunction<3> sweep_function(base_shape, translation);

    auto [time_samples, function_values, vertex_start_indices] =
        sample_time_derivative(columns, sweep_function, num_time_samples_per_vertex);

    columns.set_time_samples(time_samples, function_values, vertex_start_indices);

    auto contour = columns.extract_contour(0.0, false);
    contour.triangulate_cycles();

    size_t num_contour_vertices = contour.get_num_vertices();

    function_values.clear();
    function_values.resize(num_contour_vertices);
    for (size_t i = 0; i < num_contour_vertices; ++i) {
        auto pos = contour.get_vertex(i);
        function_values[i] = sweep_function.value({pos[0], pos[1], pos[2]}, pos[3]);
    }

    return contour.isocontour(function_values);
}

mtetcol::Contour<4> sphere_rotation(mtetcol::SimplicialColumn<4>& columns)
{
    using Scalar = mtetcol::Scalar;
    using Index = mtetcol::Index;

    constexpr size_t num_time_samples_per_vertex = 10;

    mtetcol::ImplicitSphere base_shape(0.35, {0.25, 0.5, 0.5});
    mtetcol::Rotation<3> rotation({0.5, 0.5, 0.5}, {0, 0, 1});
    mtetcol::SweepFunction<3> sweep_function(base_shape, rotation);

    auto [time_samples, function_values, vertex_start_indices] =
        sample_time_derivative(columns, sweep_function, num_time_samples_per_vertex, true);

    columns.set_time_samples(time_samples, function_values, vertex_start_indices);

    auto contour = columns.extract_contour(0.0, true);
    contour.triangulate_cycles();

    size_t num_contour_vertices = contour.get_num_vertices();
    function_values.clear();
    function_values.resize(num_contour_vertices);
    for (size_t i = 0; i < num_contour_vertices; ++i) {
        auto pos = contour.get_vertex(i);
        function_values[i] = sweep_function.value({pos[0], pos[1], pos[2]}, pos[3]);
    }

    return contour.isocontour(function_values);
}

mtetcol::Contour<4> torus_flip(mtetcol::SimplicialColumn<4>& columns)
{
    using Scalar = mtetcol::Scalar;
    using Index = mtetcol::Index;

    constexpr size_t num_time_samples_per_vertex = 64;

    mtetcol::ImplicitTorus base_shape(0.2, 0.04, {0.25, 0.5, 0.5});
    mtetcol::Rotation<3> rotation({0.25, 0.5, 0.5}, {1, 0, 0});
    mtetcol::Translation<3> translation({-0.5, 0, 0});
    mtetcol::Compose<3> flip(translation, rotation);
    mtetcol::SweepFunction<3> sweep_function(base_shape, flip);

    auto t0 = std::chrono::high_resolution_clock::now();
    auto [time_samples, function_values, vertex_start_indices] =
        sample_time_derivative(columns, sweep_function, num_time_samples_per_vertex);

    auto t1 = std::chrono::high_resolution_clock::now();

    columns.set_time_samples(time_samples, function_values, vertex_start_indices);

    auto t2 = std::chrono::high_resolution_clock::now();

    auto contour = columns.extract_contour(0.0, false);

    auto t3 = std::chrono::high_resolution_clock::now();

    contour.triangulate_cycles();

    auto t4 = std::chrono::high_resolution_clock::now();

    size_t num_contour_vertices = contour.get_num_vertices();
    function_values.clear();
    function_values.resize(num_contour_vertices);
    for (size_t i = 0; i < num_contour_vertices; ++i) {
        auto pos = contour.get_vertex(i);
        function_values[i] = sweep_function.value({pos[0], pos[1], pos[2]}, pos[3]);
    }

    auto t5 = std::chrono::high_resolution_clock::now();

    auto result = contour.isocontour(function_values);
    auto t6 = std::chrono::high_resolution_clock::now();

    mtetcol::logger().info(
        "sample_time_derivative: {} ms",
        std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count());
    mtetcol::logger().info(
        "set_time_samples: {} ms",
        std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count());
    mtetcol::logger().info(
        "extract_contour: {} ms",
        std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count());
    mtetcol::logger().info(
        "triangulate_cycles: {} ms",
        std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count());
    mtetcol::logger().info(
        "initialize function value: {} ms",
        std::chrono::duration_cast<std::chrono::milliseconds>(t5 - t4).count());
    mtetcol::logger().info(
        "isocontour: {} ms",
        std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count());

    return result;
}

int main(int argc, char** argv)
{
    mtetcol::logger().set_level(spdlog::level::debug);

    constexpr size_t resolution = 64;
    auto tet_mesh = mtet::generate_tet_grid(
        {resolution, resolution, resolution},
        {0, 0, 0},
        {1, 1, 1},
        mtet::TET5);
    mtet::save_mesh("grid.msh", tet_mesh);

    auto [vertices, tets] = generate_simpicial_column<mtetcol::Scalar, mtetcol::Index>(tet_mesh);
    mtetcol::SimplicialColumn<4> columns;
    columns.set_vertices(vertices);
    columns.set_simplices(tets);

    // auto isocontour = sphere_translation(columns);
    auto isocontour = sphere_rotation(columns);
    // auto isocontour = torus_flip(columns);

    isocontour.triangulate_cycles();
    mtetcol::save_contour("contour.msh", isocontour);

    return 0;
}
