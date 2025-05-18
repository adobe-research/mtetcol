#include <mtetcol/contour.h>
#include <mtetcol/io.h>
#include <mtetcol/logger.h>
#include <mtetcol/simplicial_column.h>

#include "flipping_donut.h"

#include <ankerl/unordered_dense.h>
#include <mtet/grid.h>
#include <mtet/io.h>
#include <spdlog/fmt/ranges.h>
#include <stf/stf.h>

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
    const stf::SpaceTimeFunction<3>& sweep_function,
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

mtetcol::Contour<4> generate_contour(
    mtetcol::SimplicialColumn<4>& columns,
    const stf::SpaceTimeFunction<3>& sweep_function)
{
    using Scalar = mtetcol::Scalar;
    using Index = mtetcol::Index;

    constexpr size_t num_time_samples_per_vertex = 128;

    auto t0 = std::chrono::high_resolution_clock::now();
    auto [time_samples, function_values, vertex_start_indices] =
        sample_time_derivative(columns, sweep_function, num_time_samples_per_vertex);

    auto t1 = std::chrono::high_resolution_clock::now();

    columns.set_time_samples(time_samples, function_values, vertex_start_indices);

    auto t2 = std::chrono::high_resolution_clock::now();

    auto contour = columns.extract_contour(0.0, false);

    auto t3 = std::chrono::high_resolution_clock::now();

    // contour.triangulate_cycles();
    // mtetcol::save_contour("contour_grid.msh", contour);

    auto t4 = std::chrono::high_resolution_clock::now();

    size_t num_contour_vertices = contour.get_num_vertices();
    function_values.clear();
    function_values.resize(num_contour_vertices);
    for (size_t i = 0; i < num_contour_vertices; ++i) {
        auto pos = contour.get_vertex(i);
        function_values[i] = sweep_function.value({pos[0], pos[1], pos[2]}, pos[3]);
    }

    std::vector<Scalar> function_gradients;
    function_gradients.reserve(num_contour_vertices * 4);
    for (size_t i = 0; i < num_contour_vertices; ++i) {
        auto pos = contour.get_vertex(i);
        std::array<Scalar, 4> gradient = sweep_function.gradient({pos[0], pos[1], pos[2]}, pos[3]);
        for (int j = 0; j < 4; ++j) {
            function_gradients.push_back(gradient[j]);
        }
    }

    auto t5 = std::chrono::high_resolution_clock::now();

    auto result = contour.isocontour(function_values, function_gradients, true);
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

mtetcol::Contour<4> sphere_translation(mtetcol::SimplicialColumn<4>& columns)
{
    stf::ImplicitSphere base_shape(0.2, {0.25, 0.5, 0.5});
    stf::Translation<3> translation({-0.5, 0, 0});
    stf::SweepFunction<3> sweep_function(base_shape, translation);

    return generate_contour(columns, sweep_function);
}

mtetcol::Contour<4> sphere_rotation(mtetcol::SimplicialColumn<4>& columns)
{
    stf::ImplicitSphere base_shape(0.3, {0.4, 0.5, 0.5});
    stf::Rotation<3> rotation({0.5, 0.5, 0.5}, {0, 0, 1});
    stf::SweepFunction<3> sweep_function(base_shape, rotation);

    return generate_contour(columns, sweep_function);
}

mtetcol::Contour<4> torus_rotation(mtetcol::SimplicialColumn<4>& columns)
{
    stf::ImplicitTorus base_shape(0.2, 0.04, {0.25, 0.5, 0.5});
    stf::Rotation<3> rotation({0.25, 0.5, 0.5}, {1, 0, 0});
    stf::Translation<3> translation({-0.5, 0, 0});
    stf::Compose<3> flip(translation, rotation);
    stf::SweepFunction<3> sweep_function(base_shape, flip);

    return generate_contour(columns, sweep_function);
}

mtetcol::Contour<4> torus_flip(mtetcol::SimplicialColumn<4>& columns)
{
    using Scalar = mtetcol::Scalar;
    stf::ExplicitForm<3> sweep_function(
        [](std::array<Scalar, 3> pos, Scalar t) {
            return mtetcol::raw_flipping_donut({pos[0], pos[1], pos[2], t}).first;
        },
        [](std::array<Scalar, 3> pos, Scalar t) {
            return mtetcol::raw_flipping_donut({pos[0], pos[1], pos[2], t}).second[3];
        });

    return generate_contour(columns, sweep_function);
}

mtetcol::Contour<4> elbow(mtetcol::SimplicialColumn<4>& columns)
{
    // stf::ImplicitTorus base_shape(0.2, 0.05, {0.0, 0.0, 0.0});
    stf::ImplicitSphere base_shape(0.2, {0.0, 0.0, 0.0});
    stf::Polyline<3> polyline({{0.3, 0.3, 0.3}, {0.7, 0.3, 0.3}, {0.7, 0.7, 0.3}});
    stf::SweepFunction<3> sweep_function(base_shape, polyline);

    return generate_contour(columns, sweep_function);
}

mtetcol::Contour<4> bezier(mtetcol::SimplicialColumn<4>& columns)
{
    stf::ImplicitTorus base_shape(0.07, 0.03, {0.0, 0.0, 0.0});
    // stf::ImplicitSphere base_shape(0.05, {0.0, 0.0, 0.0});
    stf::PolyBezier<3> bezier(
        {{0.2, 0.2, 0.3}, {1.4, 0.8, 0.3}, {-0.4, 0.8, 0.3}, {0.8, 0.2, 0.3}});
    stf::SweepFunction<3> sweep_function(base_shape, bezier);

    return generate_contour(columns, sweep_function);
}

mtetcol::Contour<4> blending(mtetcol::SimplicialColumn<4>& columns)
{
    stf::ImplicitSphere sphere(0.08, {0.0, 0.0, 0.0});
    stf::ImplicitTorus torus(0.1, 0.05, {0.0, 0.0, 0.0});
    stf::Polyline<3> polyline({{0.2, 0.5, 0.5}, {0.8, 0.5, 0.5}});
    // stf::ImplicitSphere base_shape(0.05, {0.0, 0.0, 0.0});
    stf::SweepFunction<3> sphere_sweep(sphere, polyline);
    stf::SweepFunction<3> torus_sweep(torus, polyline);
    stf::InterpolateFunction<3> blend(sphere_sweep, torus_sweep);

    return generate_contour(columns, blend);
}

mtetcol::Contour<4> blending_spheres(mtetcol::SimplicialColumn<4>& columns)
{
    stf::ImplicitSphere sphere(0.2, {0.0, 0.0, 0.0}, 2);
    stf::ImplicitSphere sphere2(0.2, {0.0, 0.3, 0.0}, 2);
    stf::ImplicitSphere sphere3(0.2, {0.0, -0.3, 0.0}, 2);
    stf::ImplicitUnion two_spheres(sphere2, sphere3, 0.03);

    stf::Polyline<3> polyline({{0.2, 0.5, 0.5}, {0.8, 0.5, 0.5}});
    stf::SweepFunction<3> sphere_sweep(sphere, polyline);
    stf::SweepFunction<3> two_sphere_sweep(two_spheres, polyline);
    stf::InterpolateFunction<3> blend(sphere_sweep, two_sphere_sweep);

    return generate_contour(columns, blend);
}

mtetcol::Contour<4> union_of_sweeps(mtetcol::SimplicialColumn<4>& columns)
{
    stf::ImplicitSphere sphere(0.1, {0.0, 0.0, 0.0});

    stf::Polyline<3> path1({{0.31, 0.5, 0.5}, {0.7, 0.75, 0.5}});
    stf::Polyline<3> path2({{0.3, 0.5, 0.5}, {0.7, 0.25, 0.5}});
    stf::SweepFunction<3> sphere_sweep_1(sphere, path1);
    stf::SweepFunction<3> sphere_sweep_2(sphere, path2);
    stf::UnionFunction<3> merged_sweep(sphere_sweep_1, sphere_sweep_2, 0.05);

    return generate_contour(columns, merged_sweep);
}

mtetcol::Contour<4> sphere_spiral(mtetcol::SimplicialColumn<4>& columns)
{
    stf::ImplicitSphere base_shape(0.05, {0.0, 0.0, 0.0});
    // stf::ImplicitCapsule<3> base_shape(0.01, {0.0, 0.0, 0.0}, {0.1, 0, 0});
    // clang-format off
    std::vector<std::array<stf::Scalar, 3>> samples {
        { 0.500000, 0.050000, 0.500000},
        { 0.522888, 0.056137, 0.570442},
        { 0.381791, 0.074382, 0.585884},
        { 0.326728, 0.104237, 0.374110},
        { 0.585411, 0.144887, 0.237132},
        { 0.831076, 0.195223, 0.500000},
        { 0.616414, 0.253873, 0.858287},
        { 0.166606, 0.319237, 0.742225},
        { 0.147082, 0.389532, 0.243590},
        { 0.638583, 0.462839, 0.073486},
        { 0.948463, 0.537161, 0.500000},
        { 0.634803, 0.610468, 0.914880},
        { 0.166606, 0.680763, 0.742225},
        { 0.195223, 0.746127, 0.278567},
        { 0.602308, 0.804777, 0.185128},
        { 0.776396, 0.855113, 0.500000},
        { 0.566184, 0.895763, 0.703694},
        { 0.381791, 0.925618, 0.585884},
        { 0.440078, 0.943863, 0.456464},
        { 0.500000, 0.950000, 0.500000},
        { 0.425932, 0.943863, 0.500000},
    };
    auto curve = stf::PolyBezier<3>::from_samples(samples);
    //auto curve = stf::Polyline<3>(samples);
    // clang-format on
    stf::SweepFunction<3> sweep_function(base_shape, curve);

    return generate_contour(columns, sweep_function);
}

mtetcol::Contour<4> knot(mtetcol::SimplicialColumn<4>& columns)
{
    // stf::ImplicitSphere base_shape(0.025, {0.0, 0.0, 0.0});
    stf::ImplicitCapsule<3> base_shape(0.01, {-0.1, 0.0, 0.0}, {0.1, 0, 0});
    std::vector<std::array<stf::Scalar, 3>> samples{
        {0.4000, 0.0000, 0.0000},    {0.4000, 0.0832, 0.1464},    {0.1708, 0.0915, 0.2424},
        {0.0174, 0.0985, 0.1732},    {-0.1140, 0.1045, 0.1139},   {-0.0335, 0.1510, -0.1139},
        {-0.0940, 0.0342, -0.1732},  {-0.1646, -0.1021, -0.2424}, {-0.2721, -0.3048, -0.1464},
        {-0.2000, -0.3464, -0.0000}, {-0.1279, -0.3880, 0.1464},  {-0.0062, -0.1936, 0.2424},
        {0.0766, -0.0643, 0.1732},   {0.1475, 0.0465, 0.1139},    {0.1475, -0.0465, -0.1139},
        {0.0766, 0.0643, -0.1732},   {-0.0062, 0.1936, -0.2424},  {-0.1279, 0.3880, -0.1464},
        {-0.2000, 0.3464, -0.0000},  {-0.2721, 0.3048, 0.1464},   {-0.1646, 0.1021, 0.2424},
        {-0.0940, -0.0342, 0.1732},  {-0.0335, -0.1510, 0.1139},  {-0.1140, -0.1045, -0.1139},
        {0.0174, -0.0985, -0.1732},  {0.1708, -0.0915, -0.2424},  {0.4000, -0.0832, -0.1464},
        {0.4000, 0.0000, 0.0000},
    };
    for (auto& p : samples) {
        p[0] += 0.5;
        p[1] += 0.5;
        p[2] += 0.5;
    }
    auto curve = stf::PolyBezier<3>(samples);
    //auto curve = stf::Polyline<3>(samples);
    stf::SweepFunction<3> sweep_function(base_shape, curve);

    return generate_contour(columns, sweep_function);
}


int main(int argc, char** argv)
{
    mtetcol::logger().set_level(spdlog::level::debug);

    constexpr size_t resolution = 128;
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
    // auto isocontour = sphere_rotation(columns);
    // auto isocontour = torus_rotation(columns);
    // auto isocontour = torus_flip(columns);
    // auto isocontour = elbow(columns);
    // auto isocontour = bezier(columns);
    // auto isocontour = blending(columns);
    // auto isocontour = blending_spheres(columns);
    // auto isocontour = sphere_spiral(columns);
    // auto isocontour = union_of_sweeps(columns);
    auto isocontour = knot(columns);

    isocontour.triangulate_cycles();
    mtetcol::save_contour("contour.msh", isocontour);

    return 0;
}
