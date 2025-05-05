#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <vector>
#include <array>
#include <cmath>

#include <mtetcol/logger.h>
#include <mtetcol/simplicial_column.h>

TEST_CASE("simplicial_column", "[mtetcol]")
{
    using Index = mtetcol::Index;
    using Scalar = mtetcol::Scalar;
    using SignedIndex = mtetcol::SignedIndex;

    auto finite_difference = [](auto& f, auto& df, Scalar x, Scalar y, Scalar z, Scalar t) {
        Scalar val = f(x, y, z, t);
        Scalar val_prev = f(x, y, z, t - 1e-3);
        Scalar val_next = f(x, y, z, t + 1e-3);
        Scalar fd = (val_next - val_prev) / (2 * 1e-3);
        if (std::abs(fd) < 1e-6) {
            REQUIRE_THAT(df(x, y, z, t), Catch::Matchers::WithinAbs(fd, 1e-9));
        } else {
            REQUIRE_THAT(df(x, y, z, t), Catch::Matchers::WithinRel(fd, 1e-3));
        }
    };

    // Sphere of radius 0.5 translating in the x direction by 1 unit.
    auto sphere_translation = [](Scalar x, Scalar y, Scalar z, Scalar t) -> Scalar {
        const Scalar radius = 0.5;
        Scalar center[3] = {t - 0.5, 0, 0};
        Scalar d = std::sqrt(
            (x - center[0]) * (x - center[0]) + (y - center[1]) * (y - center[1]) +
            (z - center[2]) * (z - center[2]));
        return d - radius;
    };

    auto sphere_translation_time_derivative = [](Scalar x, Scalar y, Scalar z, Scalar t) -> Scalar {
        const Scalar radius = 0.5;
        Scalar center[3] = {t - 0.5, 0, 0};
        Scalar d = std::sqrt(
            (x - center[0]) * (x - center[0]) + (y - center[1]) * (y - center[1]) +
            (z - center[2]) * (z - center[2]));
        if (d == 0) return 0;

        Scalar spatial_grad[3] = {(x - center[0]) / d, (y - center[1]) / d, (z - center[2]) / d};
        return -spatial_grad[0];
    };

    auto check_sphere_translation_time_derivative = [&](Scalar x, Scalar y, Scalar z, Scalar t) {
        return finite_difference(
            sphere_translation,
            sphere_translation_time_derivative,
            x,
            y,
            z,
            t);
    };

    auto sphere_rotation = [](Scalar x, Scalar y, Scalar z, Scalar t) -> Scalar {
        constexpr Scalar offset = 0.1;
        const Scalar radius = 0.5;
        const Scalar theta = t * 2 * M_PI;
        Scalar center[3] = {std::cos(theta) + offset, std::sin(theta) + offset, offset};
        Scalar d = std::sqrt(
            (x - center[0]) * (x - center[0]) + (y - center[1]) * (y - center[1]) +
            (z - center[2]) * (z - center[2]));
        return d - radius;
    };

    auto sphere_rotation_time_derivative = [](Scalar x, Scalar y, Scalar z, Scalar t) -> Scalar {
        constexpr Scalar offset = 0.1;
        const Scalar radius = 0.5;
        const Scalar theta = t * 2 * M_PI;
        Scalar center[3] = {std::cos(theta) + offset, std::sin(theta) + offset, offset};
        Scalar d = std::sqrt(
            (x - center[0]) * (x - center[0]) + (y - center[1]) * (y - center[1]) +
            (z - center[2]) * (z - center[2]));
        if (d < 1e-12) return 0;

        Scalar spatial_grad[3] = {(x - center[0]) / d, (y - center[1]) / d, (z - center[2]) / d};
        return (spatial_grad[0] * std::sin(theta) - spatial_grad[1] * std::cos(theta)) * 2 * M_PI;
    };

    auto check_sphere_rotation_time_derivative = [&](Scalar x, Scalar y, Scalar z, Scalar t) {
        finite_difference(sphere_rotation, sphere_rotation_time_derivative, x, y, z, t);
    };

    auto populate_columns = [&](mtetcol::SimplicialColumn<4>& columns,
                                size_t num_time_samples_per_vertex,
                                auto& f,
                                auto& df) {
        auto vertices = columns.get_spatial_vertices();
        size_t num_vertices = vertices.size() / 3;

        std::vector<Scalar> time_samples;
        std::vector<Scalar> function_values;
        std::vector<Index> start_indices;
        time_samples.reserve(num_time_samples_per_vertex * num_vertices);
        function_values.reserve(num_time_samples_per_vertex * num_vertices);
        start_indices.reserve(num_vertices + 1);
        start_indices.push_back(0);

        for (size_t i = 0; i < num_vertices; i++) {
            Scalar x = vertices[i * 3];
            Scalar y = vertices[i * 3 + 1];
            Scalar z = vertices[i * 3 + 2];
            for (size_t j = 0; j < num_time_samples_per_vertex; j++) {
                Scalar t = static_cast<Scalar>(j) / (num_time_samples_per_vertex - 1);
                time_samples.push_back(t);
                function_values.push_back(df(x, y, z, t));
                finite_difference(f, df, x, y, z, t);
            }
            start_indices.push_back(static_cast<Index>(time_samples.size()));
        }

        columns.set_time_samples(
            std::span<Scalar>(time_samples.data(), time_samples.size()),
            std::span<Scalar>(function_values.data(), function_values.size()),
            std::span<Index>(start_indices.data(), start_indices.size()));
    };

    auto check_contour = [&](mtetcol::Contour<4>& contour) {
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
    };

    auto check_translation = [&](mtetcol::SimplicialColumn<4>& columns,
                                 size_t num_time_samples_per_vertex) {
        populate_columns(
            columns,
            num_time_samples_per_vertex,
            sphere_translation,
            sphere_translation_time_derivative);
        //mtetcol::logger().set_level(spdlog::level::debug);
        auto contour = columns.extract_contour(0, false);
        //mtetcol::logger().set_level(spdlog::level::warn);

        REQUIRE(contour.get_num_vertices() == columns.get_num_spatial_vertices());
        REQUIRE(contour.get_num_segments() == columns.get_num_spatial_edges());
        REQUIRE(contour.get_num_cycles() == columns.get_num_spatial_triangles());

        check_contour(contour);
    };

    auto check_rotation = [&](mtetcol::SimplicialColumn<4>& columns,
                              size_t num_time_samples_per_vertex) {
        populate_columns(
            columns,
            num_time_samples_per_vertex,
            sphere_rotation,
            sphere_rotation_time_derivative);
        auto contour = columns.extract_contour(0, false);
        check_contour(contour);

        // mtetcol::logger().set_level(spdlog::level::debug);
        auto cyclic_contour = columns.extract_contour(0, true);
        check_contour(cyclic_contour);
        // mtetcol::logger().set_level(spdlog::level::warn);
    };

    SECTION("Derivative test")
    {
        // finite_difference(
        //     sphere_rotation, sphere_rotation_time_derivative, 1, 0, 0, 0);
        // Scalar d0 = sphere_rotation_time_derivative(1, 0, 0, 0);
        // Scalar d1 = sphere_rotation_time_derivative(1, 0, 0, 0.001);
        // Scalar d2 = sphere_rotation_time_derivative(1, 0, 0, -0.001);
        // Scalar d3 = sphere_rotation_time_derivative(1, 0, 0, 0.999);
        // REQUIRE_THAT(d0, Catch::Matchers::WithinRel(d1, 1e-3));
        // REQUIRE_THAT(sphere_rotation_time_derivative(1, 0, 0, 1),
        //     Catch::Matchers::WithinAbs(0, 1e-3));
    }

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
        columns.set_vertices(std::span<Scalar>(vertices, sizeof(vertices)/sizeof(Scalar)));
        columns.set_simplices(std::span<Index>(tets, sizeof(tets)/sizeof(Index)));

        REQUIRE(columns.get_num_spatial_vertices() == 4);
        REQUIRE(columns.get_num_spatial_edges() == 6);
        REQUIRE(columns.get_num_spatial_triangles() == 4);
        REQUIRE(columns.get_num_spatial_tetrahedra() == 1);

        check_translation(columns, 10);
        check_rotation(columns, 10);
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
        columns.set_vertices(std::span<Scalar>(vertices, sizeof(vertices)/sizeof(Scalar)));
        columns.set_simplices(std::span<Index>(tets, sizeof(tets)/sizeof(Index)));

        REQUIRE(columns.get_num_spatial_vertices() == 5);
        REQUIRE(columns.get_num_spatial_edges() == 9);
        REQUIRE(columns.get_num_spatial_triangles() == 7);
        REQUIRE(columns.get_num_spatial_tetrahedra() == 2);

        check_translation(columns, 10);
        check_rotation(columns, 10);
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
        columns.set_vertices(std::span<Scalar>(vertices, sizeof(vertices)/sizeof(Scalar)));
        columns.set_simplices(std::span<Index>(tets, sizeof(tets)/sizeof(Index)));

        REQUIRE(columns.get_num_spatial_vertices() == 6);
        REQUIRE(columns.get_num_spatial_edges() == 11);
        REQUIRE(columns.get_num_spatial_triangles() == 8);
        REQUIRE(columns.get_num_spatial_tetrahedra() == 2);

        check_translation(columns, 10);
        check_rotation(columns, 10);
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
        columns.set_vertices(std::span<Scalar>(vertices, sizeof(vertices)/sizeof(Scalar)));
        columns.set_simplices(std::span<Index>(tets, sizeof(tets)/sizeof(Index)));

        REQUIRE(columns.get_num_spatial_vertices() == 7);
        REQUIRE(columns.get_num_spatial_edges() == 12);
        REQUIRE(columns.get_num_spatial_triangles() == 8);
        REQUIRE(columns.get_num_spatial_tetrahedra() == 2);

        check_translation(columns, 10);
        check_rotation(columns, 10);
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
        columns.set_vertices(std::span<Scalar>(vertices, sizeof(vertices)/sizeof(Scalar)));
        columns.set_simplices(std::span<Index>(tris, sizeof(tris)/sizeof(Index)));

        REQUIRE(columns.get_num_spatial_vertices() == 3);
        REQUIRE(columns.get_num_spatial_edges() == 3);
        REQUIRE(columns.get_num_spatial_triangles() == 1);
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
        columns.set_vertices(std::span<Scalar>(vertices, sizeof(vertices)/sizeof(Scalar)));
        columns.set_simplices(std::span<Index>(tris, sizeof(tris)/sizeof(Index)));

        REQUIRE(columns.get_num_spatial_vertices() == 5);
        REQUIRE(columns.get_num_spatial_edges() == 6);
        REQUIRE(columns.get_num_spatial_triangles() == 2);
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
        columns.set_vertices(std::span<Scalar>(vertices, sizeof(vertices)/sizeof(Scalar)));
        columns.set_simplices(std::span<Index>(tris, sizeof(tris)/sizeof(Index)));

        REQUIRE(columns.get_num_spatial_vertices() == 4);
        REQUIRE(columns.get_num_spatial_edges() == 5);
        REQUIRE(columns.get_num_spatial_triangles() == 2);
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
        columns.set_vertices(std::span<Scalar>(vertices, sizeof(vertices)/sizeof(Scalar)));
        columns.set_simplices(std::span<Index>(tris, sizeof(tris)/sizeof(Index)));

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
        {
            mtetcol::logger().set_level(spdlog::level::debug);

            size_t num_vertices = contour.get_num_vertices();
            for (size_t i=0; i<num_vertices; i++) {
                auto pos = contour.get_vertex(i);
                mtetcol::logger().debug("Vertex {}: ({}, {})  t: {}", i, pos[0], pos[1], pos[2]);
            }

            size_t num_segments = contour.get_num_segments();
            for (size_t i=0; i<num_segments; i++) {
                auto segment = contour.get_segment(i);
                mtetcol::logger().debug("Segment {}: ({}, {})", i, segment[0], segment[1]);
            }

            size_t num_cycles = contour.get_num_cycles();
            for (size_t i=0; i<num_cycles; i++) {
                auto cycle = contour.get_cycle(i);
                std::string str = fmt::format("Cycle {}: ", i);
                for (auto si : cycle) {
                    str += fmt::format("{} ", value_of(si));
                }
                mtetcol::logger().debug("{}", str);
            }

            mtetcol::logger().set_level(spdlog::level::warn);
        }
        REQUIRE(contour.get_num_vertices() == 5);
        REQUIRE(contour.get_num_segments() == 4);
    }
}
