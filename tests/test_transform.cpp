#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mtetcol/transform.h>

TEST_CASE("transform", "[mtetcol]")
{

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
