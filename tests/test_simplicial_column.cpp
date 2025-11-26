#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mtetcol/logger.h>
#include <mtetcol/simplicial_column.h>

#include <stf/stf.h>

#include <array>
#include <cmath>
#include <vector>

template <int dim>
void check_contour(mtetcol::Contour<dim>& contour)
{
    using namespace mtetcol;
    size_t num_vertices = contour.get_num_vertices();
    size_t num_segments = contour.get_num_segments();
    size_t num_cycles = contour.get_num_cycles();

    for (size_t si = 0; si < num_segments; si++) {
        auto segment = contour.get_segment(si);
        REQUIRE(segment.size() == 2);
        REQUIRE(segment[0] < num_vertices);
        REQUIRE(segment[1] < num_vertices);
    }

    for (size_t ci = 0; ci < num_cycles; ci++) {
        auto cycle = contour.get_cycle(ci);
        size_t cycle_size = cycle.size();
        for (size_t i = 0; i < cycle_size; i++) {
            SignedIndex curr_signed_si = cycle[i];
            SignedIndex next_signed_si = cycle[(i + 1) % cycle_size];

            Index curr_si = index(curr_signed_si);
            Index next_si = index(next_signed_si);

            bool curr_ori = orientation(curr_signed_si);
            bool next_ori = orientation(next_signed_si);

            auto curr_segment = contour.get_segment(curr_si);
            auto next_segment = contour.get_segment(next_si);

            Index v0 = curr_ori ? curr_segment[1] : curr_segment[0];
            Index v1 = next_ori ? next_segment[0] : next_segment[1];

            REQUIRE(v0 == v1);
        }
    }
}

template <int dim, typename Scalar>
void finite_difference(
    std::array<Scalar, dim> pos,
    Scalar t,
    stf::SweepFunction<dim>& sweep_function)
{
    constexpr Scalar delta = 1e-3;
    Scalar val = sweep_function.value(pos, t);
    Scalar val_prev = sweep_function.value(pos, t - delta);
    Scalar val_next = sweep_function.value(pos, t + delta);
    Scalar fd = (val_next - val_prev) / (2 * delta);
    if (std::abs(fd) < 1e-6) {
        REQUIRE_THAT(sweep_function.time_derivative(pos, t), Catch::Matchers::WithinAbs(fd, 1e-9));
    } else {
        REQUIRE_THAT(sweep_function.time_derivative(pos, t), Catch::Matchers::WithinRel(fd, 1e-3));
    }
}

template <int dim>
void populate(
    mtetcol::SimplicialColumn<dim + 1>& columns,
    stf::SweepFunction<dim>& sweep_function,
    size_t num_time_samples_per_vertex)
{
    using namespace mtetcol;
    size_t num_vertices = columns.get_num_spatial_vertices();

    std::vector<Scalar> time_samples;
    std::vector<Scalar> function_values;
    std::vector<Index> start_indices;
    time_samples.reserve(num_time_samples_per_vertex * num_vertices);
    function_values.reserve(num_time_samples_per_vertex * num_vertices);
    start_indices.reserve(num_vertices + 1);
    start_indices.push_back(0);

    for (size_t i = 0; i < num_vertices; i++) {
        std::span<const Scalar, dim> pos = columns.get_spatial_vertex(i);
        std::array<Scalar, dim> p;
        std::copy(pos.begin(), pos.end(), p.begin());

        for (size_t j = 0; j < num_time_samples_per_vertex; j++) {
            Scalar t = static_cast<Scalar>(j) / (num_time_samples_per_vertex - 1);
            time_samples.push_back(t);
            function_values.push_back(sweep_function.time_derivative(p, t));
            finite_difference<dim>(p, t, sweep_function);
        }
        start_indices.push_back(static_cast<Index>(time_samples.size()));
    }

    columns.set_time_samples(
        std::span<Scalar>(time_samples.data(), time_samples.size()),
        std::span<Scalar>(function_values.data(), function_values.size()),
        std::span<Index>(start_indices.data(), start_indices.size()));
}

template <int dim>
void check_translation(mtetcol::SimplicialColumn<dim + 1>& columns)
{
    if constexpr (dim == 2) {
        stf::ImplicitCircle base_shape(0.1, {0.3, 0.5});
        stf::Translation<2> translation({0.4, 0});
        stf::SweepFunction<2> sweep_function(base_shape, translation);

        populate(columns, sweep_function, 10);
    } else if (dim == 3) {
        stf::ImplicitSphere base_shape(0.1, {0.3, 0.5, 0.5});
        stf::Translation<3> translation({0.4, 0.0, 0.0});
        stf::SweepFunction<3> sweep_function(base_shape, translation);

        populate(columns, sweep_function, 10);
    }

    auto contour = columns.extract_contour(0, false);
    // mtetcol::logger().set_level(spdlog::level::warn);

    REQUIRE(contour.get_num_vertices() == columns.get_num_spatial_vertices());
    REQUIRE(contour.get_num_segments() == columns.get_num_spatial_edges());
    REQUIRE(contour.get_num_cycles() == columns.get_num_spatial_triangles());

    check_contour(contour);
}

template <int dim>
void check_rotation(mtetcol::SimplicialColumn<dim + 1>& columns)
{
    if constexpr (dim == 2) {
        stf::ImplicitCircle base_shape(0.1, {0.3, 0.5});
        stf::Rotation<2> rotation({0.5, 0.5}, {0, 0});
        stf::SweepFunction<2> sweep_function(base_shape, rotation);

        populate(columns, sweep_function, 10);
    } else if (dim == 3) {
        stf::ImplicitSphere base_shape(0.1, {0.3, 0.5, 0.5});
        stf::Rotation<3> rotation({0.5, 0.5, 0.5}, {0, 0, 1});
        stf::SweepFunction<3> sweep_function(base_shape, rotation);

        populate(columns, sweep_function, 10);
    }

    auto contour = columns.extract_contour(0, false);
    // mtetcol::logger().set_level(spdlog::level::warn);

    check_contour(contour);
}

template <int dim>
void check_overlapped_rotation(mtetcol::SimplicialColumn<dim + 1>& columns)
{
    if constexpr (dim == 2) {
        stf::ImplicitCircle base_shape(0.35, {0.25, 0.5});
        stf::Rotation<2> rotation({0.5, 0.5}, {0, 0});
        stf::SweepFunction<2> sweep_function(base_shape, rotation);

        populate(columns, sweep_function, 10);
    } else if (dim == 3) {
        stf::ImplicitSphere base_shape(0.35, {0.25, 0.5, 0.5});
        stf::Rotation<3> rotation({0.5, 0.5, 0.5}, {0, 0, 1});
        stf::SweepFunction<3> sweep_function(base_shape, rotation);

        populate(columns, sweep_function, 10);
    }

    auto contour = columns.extract_contour(0, false);
    // mtetcol::logger().set_level(spdlog::level::warn);

    check_contour(contour);
}

TEST_CASE("simplicial_column", "[mtetcol]")
{
    using Index = mtetcol::Index;
    using Scalar = mtetcol::Scalar;
    using SignedIndex = mtetcol::SignedIndex;

    SECTION("Single tet")
    {
        // clang-format off
        Scalar vertices[] = {
            0, 0, 0,
            1, 0, 0,
            0, 1, 0,
            0, 0, 1,
        };
        Index tets[] = {
            0, 1, 2, 3
        };
        // clang-format on

        mtetcol::SimplicialColumn<4> columns;
        columns.set_vertices(std::span<Scalar>(vertices, sizeof(vertices) / sizeof(Scalar)));
        columns.set_simplices(std::span<Index>(tets, sizeof(tets) / sizeof(Index)));

        REQUIRE(columns.get_num_spatial_vertices() == 4);
        REQUIRE(columns.get_num_spatial_edges() == 6);
        REQUIRE(columns.get_num_spatial_triangles() == 4);
        REQUIRE(columns.get_num_spatial_tetrahedra() == 1);

        check_translation<3>(columns);
        // mtetcol::logger().set_level(spdlog::level::debug);
        check_rotation<3>(columns);
        // mtetcol::logger().set_level(spdlog::level::warn);
    }

    SECTION("Two tets sharing a face")
    {
        // clang-format off
        Scalar vertices[] = {
            0, 0, 0,
            1, 0, 0,
            0, 1, 0,
            0, 0, 1,
            0, 0, -1,
        };
        Index tets[] = {
            0, 1, 2, 3,
            0, 2, 1, 4,
        };
        // clang-format on

        mtetcol::SimplicialColumn<4> columns;
        columns.set_vertices(std::span<Scalar>(vertices, sizeof(vertices) / sizeof(Scalar)));
        columns.set_simplices(std::span<Index>(tets, sizeof(tets) / sizeof(Index)));

        REQUIRE(columns.get_num_spatial_vertices() == 5);
        REQUIRE(columns.get_num_spatial_edges() == 9);
        REQUIRE(columns.get_num_spatial_triangles() == 7);
        REQUIRE(columns.get_num_spatial_tetrahedra() == 2);

        check_translation<3>(columns);
        check_rotation<3>(columns);
    }

    SECTION("Two tets sharing an edge")
    {
        // clang-format off
        Scalar vertices[] = {
            0, 0, 0,
            1, 0, 0,
            0, 1, 0,
            0, 0, 1,
            0, -1, 0,
            0, 0, -1,
        };
        Index tets[] = {
            0, 1, 2, 3,
            1, 0, 4, 5,
        };
        // clang-format on

        mtetcol::SimplicialColumn<4> columns;
        columns.set_vertices(std::span<Scalar>(vertices, sizeof(vertices) / sizeof(Scalar)));
        columns.set_simplices(std::span<Index>(tets, sizeof(tets) / sizeof(Index)));

        REQUIRE(columns.get_num_spatial_vertices() == 6);
        REQUIRE(columns.get_num_spatial_edges() == 11);
        REQUIRE(columns.get_num_spatial_triangles() == 8);
        REQUIRE(columns.get_num_spatial_tetrahedra() == 2);

        check_translation<3>(columns);
        check_rotation<3>(columns);
    }

    SECTION("Two tets sharing a vertex")
    {
        // clang-format off
        Scalar vertices[] = {
            0, 0, 0,
            1, 0, 0,
            0, 1, 0,
            0, 0, 1,
            -1, 0, 0,
            0, -1, 0,
            0, 0, -1,
        };
        Index tets[] = {
            0, 1, 2, 3,
            1, 4, 5, 6,
        };
        // clang-format on

        mtetcol::SimplicialColumn<4> columns;
        columns.set_vertices(std::span<Scalar>(vertices, sizeof(vertices) / sizeof(Scalar)));
        columns.set_simplices(std::span<Index>(tets, sizeof(tets) / sizeof(Index)));

        REQUIRE(columns.get_num_spatial_vertices() == 7);
        REQUIRE(columns.get_num_spatial_edges() == 12);
        REQUIRE(columns.get_num_spatial_triangles() == 8);
        REQUIRE(columns.get_num_spatial_tetrahedra() == 2);

        check_translation<3>(columns);
        check_rotation<3>(columns);
    }

    SECTION("Single triangle")
    {
        // clang-format off
        Scalar vertices[] = {
            0, 0,
            1, 0,
            0, 1,
        };
        Index tris[] = {
            0, 1, 2,
        };
        // clang-format on

        mtetcol::SimplicialColumn<3> columns;
        columns.set_vertices(std::span<Scalar>(vertices, sizeof(vertices) / sizeof(Scalar)));
        columns.set_simplices(std::span<Index>(tris, sizeof(tris) / sizeof(Index)));

        REQUIRE(columns.get_num_spatial_vertices() == 3);
        REQUIRE(columns.get_num_spatial_edges() == 3);
        REQUIRE(columns.get_num_spatial_triangles() == 1);

        check_translation<2>(columns);
        check_rotation<2>(columns);
    }

    SECTION("Two triangle sharing a vertex")
    {
        // clang-format off
        Scalar vertices[] = {
            0, 0,
            1, 0,
            0, 1,
            -1, 0,
            0, -1,
        };
        Index tris[] = {
            0, 1, 2,
            0, 3, 4,
        };
        // clang-format on

        mtetcol::SimplicialColumn<3> columns;
        columns.set_vertices(std::span<Scalar>(vertices, sizeof(vertices) / sizeof(Scalar)));
        columns.set_simplices(std::span<Index>(tris, sizeof(tris) / sizeof(Index)));

        REQUIRE(columns.get_num_spatial_vertices() == 5);
        REQUIRE(columns.get_num_spatial_edges() == 6);
        REQUIRE(columns.get_num_spatial_triangles() == 2);

        check_translation<2>(columns);
        check_rotation<2>(columns);
    }

    SECTION("Two triangle sharing an edge")
    {
        // clang-format off
        Scalar vertices[] = {
            0, 0,
            1, 0,
            0, 1,
            0, -1,
        };
        Index tris[] = {
            0, 1, 2,
            1, 0, 3,
        };
        // clang-format on

        mtetcol::SimplicialColumn<3> columns;
        columns.set_vertices(std::span<Scalar>(vertices, sizeof(vertices) / sizeof(Scalar)));
        columns.set_simplices(std::span<Index>(tris, sizeof(tris) / sizeof(Index)));

        REQUIRE(columns.get_num_spatial_vertices() == 4);
        REQUIRE(columns.get_num_spatial_edges() == 5);
        REQUIRE(columns.get_num_spatial_triangles() == 2);

        check_translation<2>(columns);
        check_rotation<2>(columns);
    }

    SECTION("Bubble")
    {
        // clang-format off
        Scalar vertices[] = {
            0, 0,
            1, 0,
            0, 1,
        };
        Index tris[] = {
            0, 1, 2,
        };
        // clang-format on

        mtetcol::SimplicialColumn<3> columns;
        columns.set_vertices(std::span<Scalar>(vertices, sizeof(vertices) / sizeof(Scalar)));
        columns.set_simplices(std::span<Index>(tris, sizeof(tris) / sizeof(Index)));

        REQUIRE(columns.get_num_spatial_vertices() == 3);
        REQUIRE(columns.get_num_spatial_edges() == 3);
        REQUIRE(columns.get_num_spatial_triangles() == 1);

        std::vector<Scalar> time_samples = {0, 0.4, 0.5, 0.6, 1, 0, 1, 0, 1};
        std::vector<Scalar> function_values = {1, 1, -1, 1, 1, 1, 1, 1, 1};
        std::vector<Index> vertex_start_indices = {0, 5, 7, 9};
        columns.set_time_samples(
            std::span<Scalar>(time_samples.data(), time_samples.size()),
            std::span<Scalar>(function_values.data(), function_values.size()),
            std::span<Index>(vertex_start_indices.data(), vertex_start_indices.size()));

        mtetcol::Contour<3> contour = columns.extract_contour(0, false);
        REQUIRE(contour.get_num_vertices() == 5);
        REQUIRE(contour.get_num_segments() == 4);

        check_translation<2>(columns);
        check_rotation<2>(columns);
    }

    SECTION("Singularity 2D")
    {
        // clang-format off
        Scalar vertices[] = {
            0, 0.5,
            0.25, 0.5,
            0, 0.75,
        };
        Index tris[] = {
            0, 1, 2,
        };
        // clang-format on

        mtetcol::SimplicialColumn<3> columns;
        columns.set_vertices(std::span<Scalar>(vertices, 6));
        columns.set_simplices(std::span<Index>(tris, 3));

        REQUIRE(columns.get_num_spatial_vertices() == 3);
        REQUIRE(columns.get_num_spatial_edges() == 3);
        REQUIRE(columns.get_num_spatial_triangles() == 1);

        check_translation<2>(columns);
        check_rotation<2>(columns);
    }

    SECTION("cusp 3D")
    {
        // clang-format off
        Scalar vertices[] = {
            0.25, 0.5, 0,
            0.5, 0.5, 0,
            0.5, 0.75, 0,
            0.5, 0.5, 0.25,
        };
        Index tets[] = {
            0, 1, 2, 3
        };
        // clang-format on

        mtetcol::SimplicialColumn<4> columns;
        columns.set_vertices(std::span<Scalar>(vertices, 12));
        columns.set_simplices(std::span<Index>(tets, 4));

        check_overlapped_rotation<3>(columns);
    }

    SECTION("Non-manifold mesh: triangle shared by 3 tets")
    {
        // clang-format off
        // Create a configuration where 3 tetrahedra share the same triangle [0, 1, 2]
        // This is non-manifold and should throw an exception
        Scalar vertices[] = {
            0, 0, 0,    // vertex 0
            1, 0, 0,    // vertex 1
            0, 1, 0,    // vertex 2
            0, 0, 1,    // vertex 3
            0, 0, -1,   // vertex 4
            0.5, 0.5, 2 // vertex 5
        };
        Index tets[] = {
            0, 1, 2, 3,  // first tet with triangle [0, 2, 1]
            0, 2, 1, 4,  // second tet with triangle [0, 1, 2] (same as above but reversed)
            0, 1, 2, 5,  // third tet with triangle [0, 2, 1] (same orientation as first)
        };
        // clang-format on

        mtetcol::SimplicialColumn<4> columns;
        columns.set_vertices(std::span<Scalar>(vertices, sizeof(vertices) / sizeof(Scalar)));

        // This should throw an exception because triangle [0, 2, 1] appears 3 times
        REQUIRE_THROWS_AS(
            columns.set_simplices(std::span<Index>(tets, sizeof(tets) / sizeof(Index))),
            std::invalid_argument);
    }
}
