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

mtetcol::Contour<4> blending_sphere_torus(mtetcol::SimplicialColumn<4>& columns)
{
    stf::ImplicitSphere sphere(0.08, {0.0, 0.0, 0.0});
    stf::ImplicitTorus torus(0.1, 0.05, {0.0, 0.0, 0.0});
    stf::Polyline<3> polyline({{0.2, 0.5, 0.5}, {0.8, 0.5, 0.5}});
    // stf::ImplicitSphere base_shape(0.05, {0.0, 0.0, 0.0});
    stf::SweepFunction<3> sphere_sweep(sphere, polyline);
    stf::SweepFunction<3> torus_sweep(torus, polyline);
    stf::InterpolateFunction<3> blend(
        sphere_sweep,
        torus_sweep,
        [](stf::Scalar t) { return (std::sin(t * 2 * 2 * M_PI - M_PI / 2) + 1) / 2; },
        [](stf::Scalar t) { return 2 * M_PI * std::cos(t * 2 * 2 * M_PI - M_PI / 2); });

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
    stf::ImplicitSphere base_shape(0.05, {0.0, 0.0, 0.0});
    // stf::ImplicitCapsule<3> base_shape(0.01, {-0.1, 0.0, 0.0}, {0.1, 0, 0});
    // clang-format off
    std::vector<std::array<stf::Scalar, 3>> samples{
        { 0.7000, 0.5000, 0.5000 },
        { 0.7000, 0.5416, 0.5732 },
        { 0.5799, 0.5460, 0.6187 },
        { 0.5087, 0.5493, 0.5866 },
        { 0.4375, 0.5525, 0.5545 },
        { 0.4858, 0.5804, 0.4455 },
        { 0.4530, 0.5171, 0.4134 },
        { 0.4202, 0.4538, 0.3813 },
        { 0.3639, 0.3476, 0.4268 },
        { 0.4000, 0.3268, 0.5000 },
        { 0.4360, 0.3060, 0.5732 },
        { 0.4999, 0.4078, 0.6187 },
        { 0.5383, 0.4678, 0.5866 },
        { 0.5767, 0.5279, 0.5545 },
        { 0.5767, 0.4721, 0.4455 },
        { 0.5383, 0.5322, 0.4134 },
        { 0.4999, 0.5922, 0.3813 },
        { 0.4360, 0.6940, 0.4268 },
        { 0.4000, 0.6732, 0.5000 },
        { 0.3639, 0.6524, 0.5732 },
        { 0.4202, 0.5462, 0.6187 },
        { 0.4530, 0.4829, 0.5866 },
        { 0.4858, 0.4196, 0.5545 },
        { 0.4375, 0.4475, 0.4455 },
        { 0.5087, 0.4507, 0.4134 },
        { 0.5799, 0.4540, 0.3813 },
        { 0.7000, 0.4584, 0.4268 },
        { 0.7000, 0.5000, 0.5000 },
    };
    // clang-format on
    auto curve = stf::PolyBezier<3>(samples);
    // auto curve = stf::Polyline<3>(samples);
    stf::SweepFunction<3> sweep_function(base_shape, curve);

    return generate_contour(columns, sweep_function);
}

mtetcol::Contour<4> concentric_rings(mtetcol::SimplicialColumn<4>& columns)
{
    stf::ImplicitSphere sphere1(0.05, {0.6, 0.5, 0.5});
    stf::ImplicitSphere sphere2(0.05, {0.71, 0.5, 0.5});
    stf::ImplicitUnion base_shape(sphere1, sphere2);

    stf::Rotation<3> rotation({0.5, 0.5, 0.5}, {0, 0, 1});
    stf::SweepFunction<3> sweep_function(base_shape, rotation);

    return generate_contour(columns, sweep_function);
}

mtetcol::Contour<4> spinning_rings(mtetcol::SimplicialColumn<4>& columns)
{
    stf::ImplicitSphere sphere1(0.05, {0.6, 0.5, 0.5});
    stf::ImplicitSphere sphere2(0.05, {0.75, 0.5, 0.5});
    stf::ImplicitUnion base_shape(sphere1, sphere2);

    stf::Rotation<3> spin({0.675, 0.5, 0.5}, {0, 1, 0}, 360 * 3);
    stf::Rotation<3> rotation({0.5, 0.5, 0.5}, {0, 0, 1});
    stf::Compose<3> transform(rotation, spin);
    stf::SweepFunction<3> sweep_function(base_shape, transform);

    return generate_contour(columns, sweep_function);
}

mtetcol::Contour<4> spinning_rods(mtetcol::SimplicialColumn<4>& columns)
{
    stf::ImplicitCapsule<3> rod(0.02, {0.5, 0.5, 0.5}, {0.7, 0.5, 0.5});

    stf::Rotation<3> spin({0.6, 0.5, 0.5}, {0, 1, 0}, 360 * 3);
    stf::Rotation<3> rotation({0.5, 0.5, 0.5}, {0, 0, 1});
    stf::Compose<3> transform(rotation, spin);
    stf::SweepFunction<3> sweep_function(rod, transform);

    return generate_contour(columns, sweep_function);
}

mtetcol::Contour<4> brush_stroke(mtetcol::SimplicialColumn<4>& columns)
{
    stf::ImplicitSphere base_shape(0.045, {0.0, 0.0, 0.0});

    // clang-format off
    std::vector<std::array<stf::Scalar, 3>> samples{
        {0.3090, 0.6504, 0.5},
        {0.3074, 0.6294, 0.5},
        {0.3000, 0.5009, 0.5},
        {0.3934, 0.4126, 0.5},
        {0.4897, 0.3216, 0.5},
        {0.6325, 0.3306, 0.5},
        {0.6360, 0.3432, 0.5},
        {0.6389, 0.3537, 0.5},
        {0.5428, 0.3618, 0.5},
        {0.4755, 0.4415, 0.5},
        {0.3973, 0.5340, 0.5},
        {0.4045, 0.6679, 0.5},
        {0.4223, 0.6732, 0.5},
        {0.4402, 0.6784, 0.5},
        {0.4594, 0.5506, 0.5},
        {0.5693, 0.4865, 0.5},
        {0.6152, 0.4597, 0.5},
        {0.6628, 0.4525, 0.5},
        {0.7000, 0.4514, 0.5},
    };
    // clang-format on
    stf::PolyBezier<3> transform(samples, false);
    stf::SweepFunction<3> sweep_function(base_shape, transform);

    return generate_contour(columns, sweep_function);
}

mtetcol::Contour<4> torus_double_rotation(mtetcol::SimplicialColumn<4>& columns)
{
    stf::ImplicitTorus base_shape(0.25, 0.05, {0.5, 0.5, 0.25});
    stf::Rotation<3> rotation({0.5, 0.5, 0.25}, {1, 0, 0});
    stf::Translation<3> translation({0, 0, -0.5});
    stf::Compose<3> flip(translation, rotation);
    stf::SweepFunction<3> sweep_function(base_shape, flip);

    return generate_contour(columns, sweep_function);
}

mtetcol::Contour<4> rotating_rods(mtetcol::SimplicialColumn<4>& columns)
{
    stf::ImplicitCapsule<3> rod1(0.02, {-0.1, 0.0, 0.0}, {0.1, 0.0, 0.0});
    stf::ImplicitCapsule<3> rod2(0.02, {0.0, -0.1, 0.0}, {0.0, 0.1, 0.0});
    stf::Rotation<3> rotation1({0.0, 0.0, 0.0}, {0, 0, 1});
    stf::Rotation<3> rotation2({0.0, 0.0, 0.0}, {0, 0, -1});
    stf::Polyline<3> polyline({{0.25, 0.5, 0.5}, {0.75, 0.5, 0.5}});
    stf::Compose<3> transform1(polyline, rotation1);
    stf::Compose<3> transform2(polyline, rotation2);

    stf::SweepFunction<3> sweep_function1(rod1, transform1);
    stf::SweepFunction<3> sweep_function2(rod2, transform2);
    stf::UnionFunction<3> union_sweeps(sweep_function1, sweep_function2, 0.01);

    return generate_contour(columns, union_sweeps);
}

mtetcol::Contour<4> blending_spheres_nonlinear(mtetcol::SimplicialColumn<4>& columns)
{
    stf::ImplicitSphere sphere(0.07, {0.0, 0.0, 0.0}, 2);
    stf::ImplicitSphere sphere2(0.07, {0.0, 0.1, 0.0}, 2);
    stf::ImplicitSphere sphere3(0.07, {0.0, -0.1, 0.0}, 2);
    stf::ImplicitUnion two_spheres(sphere2, sphere3, 0.004);

    stf::Polyline<3> polyline({{0.2, 0.5, 0.5}, {0.8, 0.5, 0.5}});
    stf::SweepFunction<3> sphere_sweep(sphere, polyline);
    stf::SweepFunction<3> two_sphere_sweep(two_spheres, polyline);
    stf::InterpolateFunction<3> blend(
        sphere_sweep,
        two_sphere_sweep,
        [](stf::Scalar t) { return (std::sin(t * 2 * 2 * M_PI - M_PI / 2) + 1) / 2; },
        [](stf::Scalar t) { return 2 * M_PI * std::cos(t * 2 * 2 * M_PI - M_PI / 2); });

    return generate_contour(columns, blend);
}

stf::PolyBezier<3> letter_L_curve(bool use_frame = false)
{
    stf::PolyBezier<3> curve(
        {
            {0.6941, 0.4189, 0.45},
            {0.6457, 0.3864, 0.4532520646014161},
            {0.5952, 0.3504, 0.4567053138602828},
            {0.5448, 0.3504, 0.45951681021477053},
            {0.5076, 0.3504, 0.461587292026215},
            {0.4752, 0.3696, 0.46368795156949405},
            {0.4447, 0.3899, 0.46573106669305236},
            {0.4110, 0.4126, 0.4679918240358986},
            {0.3755, 0.4392, 0.4704676882285434},
            {0.3422, 0.4392, 0.47232022458615164},
            {0.3180, 0.4392, 0.47367148639993645},
            {0.3000, 0.4204, 0.4751204549163221},
            {0.3000, 0.3993, 0.4762973603670379},
            {0.3000, 0.3782, 0.4774742658177537},
            {0.3192, 0.3555, 0.4791290732784213},
            {0.3567, 0.3555, 0.4812213496352494},
            {0.3893, 0.3555, 0.48304119417478214},
            {0.4193, 0.3729, 0.48497004222424356},
            {0.4439, 0.3966, 0.48687369938850344},
            {0.4670, 0.4196, 0.48869220574687294},
            {0.5081, 0.4804, 0.49278055529388515},
            {0.5323, 0.5170, 0.49522667680735677},
            {0.5946, 0.6112, 0.5015202760795618},
            {0.6257, 0.6496, 0.5042756679944657},
            {0.6683, 0.6496, 0.5066512734412809},
            {0.6855, 0.6496, 0.5076102334381605},
            {0.7000, 0.6363, 0.508705395789782},
            {0.7000, 0.6159, 0.5098387121497305},
            {0.7000, 0.5741, 0.5121707285057785},
            {0.6157, 0.4908, 0.5187744847481601},
            {0.5025, 0.4873, 0.5250870538782416},
            {0.4029, 0.4842, 0.5306437666433996},
            {0.3395, 0.5440, 0.5355028764354673},
            {0.3395, 0.5933, 0.5382489891538043},
            {0.3395, 0.6245, 0.5399925527844943},
            {0.3649, 0.6488, 0.541950304765777},
            {0.4017, 0.6488, 0.5439989920318379},
            {0.4498, 0.6488, 0.5466797211140239},
            {0.4799, 0.6022, 0.5497688624782774},
            {0.4834, 0.6002, 0.55},
        },
        use_frame);
    return curve;
}

mtetcol::Contour<4> letter_L(mtetcol::SimplicialColumn<4>& columns)
{
    stf::ImplicitSphere sphere(0.02, {0.0, 0.0, 0.0});
    stf::PolyBezier<3> curve = letter_L_curve();

    stf::SweepFunction<3> sweep_function(sphere, curve);
    return generate_contour(columns, sweep_function);
}

mtetcol::Contour<4> letter_L_blend(mtetcol::SimplicialColumn<4>& columns)
{
    static stf::ImplicitTorus torus(0.03, 0.015, {0.0, 0.0, 0.0});
    stf::PolyBezier<3> curve = letter_L_curve();

    stf::Rotation<3> torus_rotation({0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 360 * 3);
    stf::Compose<3> torus_curve(curve, torus_rotation);
    stf::SweepFunction<3> torus_sweep(torus, torus_curve);
    stf::OffsetFunction<3> sweep_function(
        torus_sweep,
        [](stf::Scalar t) { return - 0.01 * std::sin(t * 6 * M_PI) - 0.01; },
        [](stf::Scalar t) { return - 0.01 * std::cos(t * 6 * M_PI) * 6 * M_PI; });

    return generate_contour(columns, sweep_function);
}

mtetcol::Contour<4> loopDloop_with_offset(mtetcol::SimplicialColumn<4>& columns)
{
    stf::ImplicitTorus base_shape(0.07, 0.03, {0.0, 0.0, 0.0});
    // stf::ImplicitSphere base_shape(0.05, {0.0, 0.0, 0.0});
    stf::PolyBezier<3> bezier(
        {{0.2, 0.2, 0.3}, {1.4, 0.8, 0.3}, {-0.4, 0.8, 0.3}, {0.8, 0.2, 0.3}});
    stf::SweepFunction<3> sweep_function(base_shape, bezier);
    stf::OffsetFunction<3> offset_function(
        sweep_function,
        [](stf::Scalar t) { return -0.02 * std::cos(t * 2 * M_PI) - 0.02; },
        [](stf::Scalar t) { return 0.02 * std::sin(t * 2 * M_PI) * 2 * M_PI; });

    return generate_contour(columns, offset_function);
}

mtetcol::Contour<4> doghead(mtetcol::SimplicialColumn<4>& columns)
{
    stf::Duchon base_shape(
        "vipss_data/doghead_800_shifted.xyz",
        "vipss_data/doghead_800_shifted_coeff",
        {0.0, 0.0, 0.0},
        0.2,
        true);
    stf::Translation<3> translation({-0.5, 0.0, 0.0});
    stf::Rotation<3> rotation_Y({0.5, 0.5, 0.5}, {0.0, 1.0, 0.0}, 180);
    stf::Rotation<3> rotation_X({0.25, 0.5, 0.5}, {1.0, 0.0, 0.0}, 360);
    stf::Rotation<3> rotation_Z({0.5, 0.5, 0.5}, {0.0, 0.0, 1.0}, 180);
    stf::Compose<3> transform(translation, rotation_Y);
    stf::PolyBezier<3> bezier(
        {{0.2, 0.2, 0.3}, {1.4, 0.8, 0.4}, {-0.4, 0.8, 0.5}, {0.8, 0.2, 0.6}});
    stf::SweepFunction<3> sweep_function(base_shape, bezier);

    return generate_contour(columns, sweep_function);
}

mtetcol::Contour<4> rotating_rod(mtetcol::SimplicialColumn<4>& columns)
{
    stf::ImplicitCapsule<3> rod(0.1, {0.5, 0.5, 0.5}, {0.85, 0.5, 0.5});
    stf::Rotation<3> rotation({0.5, 0.5, 0.5}, {0.0, 0.0, 1.0}, 180);
    stf::SweepFunction<3> sweep_function(rod, rotation);

    return generate_contour(columns, sweep_function);
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
    // auto isocontour = sphere_rotation(columns);
    // auto isocontour = torus_rotation(columns);
    // auto isocontour = torus_flip(columns);
    // auto isocontour = elbow(columns);
    // auto isocontour = bezier(columns);
    // auto isocontour = blending_sphere_torus(columns);
    // auto isocontour = blending_spheres(columns);
    // auto isocontour = blending_spheres_nonlinear(columns);
    // auto isocontour = sphere_spiral(columns);
    // auto isocontour = union_of_sweeps(columns);
    // auto isocontour = knot(columns);
    // auto isocontour = spinning_rings(columns);
    // auto isocontour = concentric_rings(columns);
    // auto isocontour = brush_stroke(columns);
    // auto isocontour = torus_double_rotation(columns);
    // auto isocontour = rotating_rods(columns);
    // auto isocontour = spinning_rods(columns);
    // auto isocontour = letter_L(columns);
    // auto isocontour = loopDloop_with_offset(columns);
    // auto isocontour = doghead(columns);
    // auto isocontour = letter_L_blend(columns);
    auto isocontour = rotating_rod(columns);

    isocontour.triangulate_cycles();
    mtetcol::save_contour("contour.msh", isocontour);

    return 0;
}
