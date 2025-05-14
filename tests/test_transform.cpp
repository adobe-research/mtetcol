#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mtetcol/transform.h>

template <int dim>
void check_velocity(
    const mtetcol::Transform<dim>& transform,
    const std::array<mtetcol::Scalar, dim>& pos,
    mtetcol::Scalar t)
{
    auto v = transform.velocity(pos, t);
    auto v_fd = transform.finite_difference_velocity(pos, t);
    for (int i = 0; i < dim; ++i) {
        REQUIRE_THAT(v[i], Catch::Matchers::WithinAbs(v_fd[i], 1e-6));
    }
}

template <int dim>
void check_jacobian(
    const mtetcol::Transform<dim>& transform,
    const std::array<mtetcol::Scalar, dim>& pos,
    mtetcol::Scalar t)
{
    auto J = transform.position_Jacobian(pos, t);
    auto J_fd = transform.finite_difference_Jacobian(pos, t);
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j) {
            REQUIRE_THAT(J[i][j], Catch::Matchers::WithinAbs(J_fd[i][j], 1e-6));
        }
}

TEST_CASE("transform", "[mtetcol]")
{
    SECTION("Rotation 2D") {
        mtetcol::Rotation<2> rotation({0.0, 0.0}, {0, 0});

        auto p0 = rotation.transform({1, 0}, 0);
        REQUIRE_THAT(p0[0], Catch::Matchers::WithinAbs(1, 1e-6));
        REQUIRE_THAT(p0[1], Catch::Matchers::WithinAbs(0, 1e-6));
        check_velocity(rotation, {1, 0}, 0);
        check_jacobian(rotation, {1, 0}, 0);

        auto p1 = rotation.transform({1, 0}, 0.5);
        REQUIRE_THAT(p1[0], Catch::Matchers::WithinAbs(-1, 1e-6));
        REQUIRE_THAT(p1[1], Catch::Matchers::WithinAbs(0, 1e-6));
        check_velocity(rotation, {1, 0}, 0.5);
        check_jacobian(rotation, {1, 0}, 0.5);

        auto p2 = rotation.transform({1, 0}, 0.25);
        REQUIRE_THAT(p2[0], Catch::Matchers::WithinAbs(0, 1e-6));
        REQUIRE_THAT(p2[1], Catch::Matchers::WithinAbs(1, 1e-6));
        check_velocity(rotation, {1, 0}, 0.25);
        check_jacobian(rotation, {1, 0}, 0.25);

        check_velocity(rotation, {1, 1}, 0.75);
        check_jacobian(rotation, {1, 1}, 0.75);
    }

    SECTION("Compose") {
        mtetcol::Translation<3> translation({1, 0, 0});
        mtetcol::Rotation<3> rotation({0, 0, 0}, {0, 0, 1});
        mtetcol::Compose<3> compose(rotation, translation);

        SECTION("Origin at t=0") {
            auto p0 = compose.transform({0, 0, 0}, 0);
            REQUIRE_THAT(p0[0], Catch::Matchers::WithinAbs(0, 1e-6));
            REQUIRE_THAT(p0[1], Catch::Matchers::WithinAbs(0, 1e-6));
            REQUIRE_THAT(p0[2], Catch::Matchers::WithinAbs(0, 1e-6));
            check_velocity(compose, {0, 0, 0}, 0);
            check_jacobian(compose, {0, 0, 0}, 0);
        }

        SECTION("Origin at t=0.5") {
            auto p0 = compose.transform({0, 0, 0}, 0.5);
            REQUIRE_THAT(p0[0], Catch::Matchers::WithinAbs(0.5, 1e-6));
            REQUIRE_THAT(p0[1], Catch::Matchers::WithinAbs(0, 1e-6));
            REQUIRE_THAT(p0[2], Catch::Matchers::WithinAbs(0, 1e-6));
            check_velocity(compose, {0, 0, 0}, 0.5);
            check_jacobian(compose, {0, 0, 0}, 0.5);
        }

        SECTION("Origin at t=1") {
            auto p0 = compose.transform({0, 0, 0}, 1);
            REQUIRE_THAT(p0[0], Catch::Matchers::WithinAbs(1, 1e-6));
            REQUIRE_THAT(p0[1], Catch::Matchers::WithinAbs(0, 1e-6));
            REQUIRE_THAT(p0[2], Catch::Matchers::WithinAbs(0, 1e-6));
            check_velocity(compose, {0, 0, 0}, 1);
            check_jacobian(compose, {0, 0, 0}, 1);
        }

        SECTION("[1, 0, 0] at t=0") {
            auto p0 = compose.transform({1, 0, 0}, 0);
            REQUIRE_THAT(p0[0], Catch::Matchers::WithinAbs(1, 1e-6));
            REQUIRE_THAT(p0[1], Catch::Matchers::WithinAbs(0, 1e-6));
            REQUIRE_THAT(p0[2], Catch::Matchers::WithinAbs(0, 1e-6));
            check_velocity(compose, {1, 0, 0}, 0);
            check_jacobian(compose, {1, 0, 0}, 0);
        }

        SECTION("[1, 0, 0] at t=0.5") {
            auto p0 = compose.transform({1, 0, 0}, 0.5);
            REQUIRE_THAT(p0[0], Catch::Matchers::WithinAbs(-0.5, 1e-6));
            REQUIRE_THAT(p0[1], Catch::Matchers::WithinAbs(0, 1e-6));
            REQUIRE_THAT(p0[2], Catch::Matchers::WithinAbs(0, 1e-6));
            check_velocity(compose, {1, 0, 0}, 0.5);
            check_jacobian(compose, {1, 0, 0}, 0.5);
        }

        SECTION("[1, 0, 0] at t=1.0") {
            auto p0 = compose.transform({1, 0, 0}, 1.0);
            REQUIRE_THAT(p0[0], Catch::Matchers::WithinAbs(2, 1e-6));
            REQUIRE_THAT(p0[1], Catch::Matchers::WithinAbs(0, 1e-6));
            REQUIRE_THAT(p0[2], Catch::Matchers::WithinAbs(0, 1e-6));
            check_velocity(compose, {1, 0, 0}, 1.0);
            check_jacobian(compose, {1, 0, 0}, 1.0);
        }
    }
}
