#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mtetcol/transform.h>

TEST_CASE("transform", "[mtetcol]")
{
    SECTION("Rotation 2D") {
        mtetcol::Rotation<2> rotation({0.0, 0.0}, {0, 0});

        auto p0 = rotation.transform({1, 0}, 0);
        auto v0 = rotation.velocity({1, 0}, 0);
        REQUIRE_THAT(p0[0], Catch::Matchers::WithinAbs(1, 1e-6));
        REQUIRE_THAT(p0[1], Catch::Matchers::WithinAbs(0, 1e-6));
        REQUIRE_THAT(v0[0], Catch::Matchers::WithinAbs(0, 1e-6));
        REQUIRE_THAT(v0[1], Catch::Matchers::WithinAbs(2 * M_PI, 1e-6));

        auto p1 = rotation.transform({1, 0}, 0.5);
        auto v1 = rotation.velocity({1, 0}, 0.5);
        REQUIRE_THAT(p1[0], Catch::Matchers::WithinAbs(-1, 1e-6));
        REQUIRE_THAT(p1[1], Catch::Matchers::WithinAbs(0, 1e-6));
        REQUIRE_THAT(v1[0], Catch::Matchers::WithinAbs(0, 1e-6));
        REQUIRE_THAT(v1[1], Catch::Matchers::WithinAbs(-2 * M_PI, 1e-6));

        auto p2 = rotation.transform({1, 0}, 0.25);
        auto v2 = rotation.velocity({1, 0}, 0.25);
        REQUIRE_THAT(p2[0], Catch::Matchers::WithinAbs(0, 1e-6));
        REQUIRE_THAT(p2[1], Catch::Matchers::WithinAbs(1, 1e-6));
        REQUIRE_THAT(v2[0], Catch::Matchers::WithinAbs(-2 * M_PI, 1e-6));
        REQUIRE_THAT(v2[1], Catch::Matchers::WithinAbs(0, 1e-6));

        auto v3 = rotation.velocity({1, 1}, 0.75);
        auto v3_fd = rotation.finite_difference({1, 1}, 0.75);
        REQUIRE_THAT(v3[0], Catch::Matchers::WithinAbs(v3_fd[0], 1e-6));
        REQUIRE_THAT(v3[1], Catch::Matchers::WithinAbs(v3_fd[1], 1e-6));
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
        }

        SECTION("Origin at t=0.5") {
            auto p0 = compose.transform({0, 0, 0}, 0.5);
            REQUIRE_THAT(p0[0], Catch::Matchers::WithinAbs(0.5, 1e-6));
            REQUIRE_THAT(p0[1], Catch::Matchers::WithinAbs(0, 1e-6));
            REQUIRE_THAT(p0[2], Catch::Matchers::WithinAbs(0, 1e-6));
        }

        SECTION("Origin at t=1") {
            auto p0 = compose.transform({0, 0, 0}, 1);
            REQUIRE_THAT(p0[0], Catch::Matchers::WithinAbs(1, 1e-6));
            REQUIRE_THAT(p0[1], Catch::Matchers::WithinAbs(0, 1e-6));
            REQUIRE_THAT(p0[2], Catch::Matchers::WithinAbs(0, 1e-6));
        }

        SECTION("[1, 0, 0] at t=0") {
            auto p0 = compose.transform({1, 0, 0}, 0);
            REQUIRE_THAT(p0[0], Catch::Matchers::WithinAbs(1, 1e-6));
            REQUIRE_THAT(p0[1], Catch::Matchers::WithinAbs(0, 1e-6));
            REQUIRE_THAT(p0[2], Catch::Matchers::WithinAbs(0, 1e-6));
        }

        SECTION("[1, 0, 0] at t=0.5") {
            auto p0 = compose.transform({1, 0, 0}, 0.5);
            REQUIRE_THAT(p0[0], Catch::Matchers::WithinAbs(-0.5, 1e-6));
            REQUIRE_THAT(p0[1], Catch::Matchers::WithinAbs(0, 1e-6));
            REQUIRE_THAT(p0[2], Catch::Matchers::WithinAbs(0, 1e-6));
        }

        SECTION("[1, 0, 0] at t=1.0") {
            auto p0 = compose.transform({1, 0, 0}, 1.0);
            REQUIRE_THAT(p0[0], Catch::Matchers::WithinAbs(2, 1e-6));
            REQUIRE_THAT(p0[1], Catch::Matchers::WithinAbs(0, 1e-6));
            REQUIRE_THAT(p0[2], Catch::Matchers::WithinAbs(0, 1e-6));
        }
    }
}
