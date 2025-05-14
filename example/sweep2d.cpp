#include <mtetcol/contour.h>
#include <mtetcol/implicit_function.h>
#include <mtetcol/io.h>
#include <mtetcol/logger.h>
#include <mtetcol/simplicial_column.h>
#include <mtetcol/sweep_function.h>
#include <mtetcol/transform.h>

using Scalar = mtetcol::Scalar;
using Index = mtetcol::Index;

auto generate_grid(
    size_t resolution,
    const std::array<Scalar, 2>& min,
    const std::array<Scalar, 2>& max) -> std::tuple<std::vector<Scalar>, std::vector<Index>>
{
    std::vector<Scalar> vertices;
    std::vector<Index> triangles;
    vertices.reserve((resolution + 1) * (resolution + 1) * 2);
    triangles.reserve(resolution * resolution * 2);

    for (size_t i = 0; i <= resolution; ++i) {
        for (size_t j = 0; j <= resolution; ++j) {
            vertices.push_back(min[0] + (max[0] - min[0]) * i / resolution);
            vertices.push_back(min[1] + (max[1] - min[1]) * j / resolution);
        }
    }

    for (size_t i = 0; i < resolution; ++i) {
        for (size_t j = 0; j < resolution; ++j) {
            triangles.push_back(i * (resolution + 1) + j);
            triangles.push_back((i + 1) * (resolution + 1) + j);
            triangles.push_back(i * (resolution + 1) + (j + 1));

            triangles.push_back((i + 1) * (resolution + 1) + j);
            triangles.push_back((i + 1) * (resolution + 1) + (j + 1));
            triangles.push_back(i * (resolution + 1) + (j + 1));
        }
    }

    return std::make_tuple(vertices, triangles);
}

auto sample_time_derivative(
    mtetcol::SimplicialColumn<3>& columns,
    const mtetcol::SpaceTimeFunction<2>& sweep_function,
    size_t num_time_samples,
    bool cyclic) -> std::
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
            function_values.push_back(sweep_function.time_derivative({pos[0], pos[1]}, t));
        }
        vertex_start_indices.push_back(static_cast<Index>(time_samples.size()));
    }

    if (cyclic) {
        // Ensure df(0) = df(1) for cyclic time series
        for (size_t vi = 0; vi < num_vertices; vi++) {
            function_values[vertex_start_indices[vi]] =
                function_values[vertex_start_indices[vi + 1] - 1];
        }
    }

    return {time_samples, function_values, vertex_start_indices};
}

mtetcol::Contour<3> circle_rotation(mtetcol::SimplicialColumn<3>& columns)
{
    using Scalar = mtetcol::Scalar;
    using Index = mtetcol::Index;

    constexpr size_t num_time_samples_per_vertex = 64;

    std::array<Scalar, 2> center = {0.50, 0.50};
    mtetcol::ImplicitCircle base_shape(0.3, {0.4, 0.5});
    mtetcol::Rotation<2> rotation(center, {0, 0});
    mtetcol::SweepFunction<2> sweep_function(base_shape, rotation);

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
        function_values[i] = sweep_function.value({pos[0], pos[1]}, pos[2]);
    }

    std::vector<Scalar> function_gradients;
    function_gradients.reserve(num_contour_vertices * 3);
    for (size_t i = 0; i < num_contour_vertices; ++i) {
        auto pos = contour.get_vertex(i);
        std::array<Scalar, 3> gradient = sweep_function.gradient({pos[0], pos[1]}, pos[2]);
        for (int j = 0; j < 3; ++j) {
            function_gradients.push_back(gradient[j]);
        }
    }

    return contour.isocontour(function_values, function_gradients);
}

mtetcol::Contour<3> lemniscate_of_bernoulli(mtetcol::SimplicialColumn<3>& columns)
{
    using Scalar = mtetcol::Scalar;
    using Index = mtetcol::Index;

    constexpr size_t num_time_samples_per_vertex = 128; 

    mtetcol::ExplicitForm<2> sweep_function(
        [](std::array<Scalar, 2> pos, Scalar t) {
        Scalar x = std::abs(pos[0] - 0.5);
        Scalar y = pos[1] - 0.5;
        t = (t - 0.5) * 20;
        Scalar a = 0.1;
        return (x - a * std::cosh(t)) * (x - a * std::cosh(t)) +
               (y - a * std::sinh(t)) * (y - a * std::sinh(t)) - a * a * std::cosh(2 * t);
        });

    auto [time_samples, function_values, vertex_start_indices] =
        sample_time_derivative(columns, sweep_function, num_time_samples_per_vertex, true);

    columns.set_time_samples(time_samples, function_values, vertex_start_indices);

    auto contour = columns.extract_contour(0.0, true);
    contour.triangulate_cycles();
    mtetcol::save_contour("contour_grid2d.msh", contour);

    size_t num_contour_vertices = contour.get_num_vertices();
    function_values.clear();
    function_values.resize(num_contour_vertices);
    for (size_t i = 0; i < num_contour_vertices; ++i) {
        auto pos = contour.get_vertex(i);
        function_values[i] = sweep_function.value({pos[0], pos[1]}, pos[2]);
    }

    std::vector<Scalar> function_gradients;
    function_gradients.reserve(num_contour_vertices * 3);
    for (size_t i = 0; i < num_contour_vertices; ++i) {
        auto pos = contour.get_vertex(i);
        std::array<Scalar, 3> gradient = sweep_function.gradient({pos[0], pos[1]}, pos[2]);
        for (int j = 0; j < 3; ++j) {
            function_gradients.push_back(gradient[j]);
        }
    }

    return contour.isocontour(function_values, function_gradients);
}

int main(int argc, char** argv)
{
    mtetcol::logger().set_level(spdlog::level::debug);

    constexpr size_t resolution = 128;

    auto [vertices, triangles] = generate_grid(resolution, {0, 0}, {1, 1});
    mtetcol::SimplicialColumn<3> columns;
    columns.set_vertices(vertices);
    columns.set_simplices(triangles);

    //auto isocontour = circle_rotation(columns);
    auto isocontour = lemniscate_of_bernoulli(columns);
    mtetcol::save_contour("contour2d.msh", isocontour);

    return 0;
}
