#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mtetcol/implicit_function.h>

#include <array>


TEST_CASE("implicit_function", "[mtetcol]")
{
    using Scalar = mtetcol::Scalar;
    auto finite_difference = [](auto& f, std::array<Scalar, 3> pos) {
        constexpr Scalar eps = 1e-6;
        auto value_prev_x = f.value({pos[0] - eps, pos[1], pos[2]});
        auto value_prev_y = f.value({pos[0], pos[1] - eps, pos[2]});
        auto value_prev_z = f.value({pos[0], pos[1], pos[2] - eps});
        auto value_next_x = f.value({pos[0] + eps, pos[1], pos[2]});
        auto value_next_y = f.value({pos[0], pos[1] + eps, pos[2]});
        auto value_next_z = f.value({pos[0], pos[1], pos[2] + eps});

        std::array<Scalar, 3> gradient;
        gradient[0] = (value_next_x - value_prev_x) / (2 * eps);
        gradient[1] = (value_next_y - value_prev_y) / (2 * eps);
        gradient[2] = (value_next_z - value_prev_z) / (2 * eps);

        return gradient;
    };

    SECTION("Torus") {
        mtetcol::ImplicitTorus torus(1, 0.2, {1, 0, 0});
        REQUIRE_THAT(torus.value({0, 0, 0}), Catch::Matchers::WithinRel(-0.2, 1e-6));
        REQUIRE_THAT(torus.value({1, 0, 0}), Catch::Matchers::WithinRel(0.8, 1e-6));

        SECTION("Gradient check 1") {
            std::array<Scalar, 3> pos = {0.2, 0, 0};
            auto gradient = torus.gradient(pos);
            auto gradient_fd = finite_difference(torus, pos);
            REQUIRE_THAT(gradient[0], Catch::Matchers::WithinAbs(gradient_fd[0], 1e-6));
            REQUIRE_THAT(gradient[1], Catch::Matchers::WithinAbs(gradient_fd[1], 1e-6));
            REQUIRE_THAT(gradient[2], Catch::Matchers::WithinAbs(gradient_fd[2], 1e-6));
        }

        SECTION("Gradient check 2") {
            std::array<Scalar, 3> pos = {0.2, 0.2, 0.2};
            auto gradient = torus.gradient(pos);
            auto gradient_fd = finite_difference(torus, pos);
            REQUIRE_THAT(gradient[0], Catch::Matchers::WithinAbs(gradient_fd[0], 1e-6));
            REQUIRE_THAT(gradient[1], Catch::Matchers::WithinAbs(gradient_fd[1], 1e-6));
            REQUIRE_THAT(gradient[2], Catch::Matchers::WithinAbs(gradient_fd[2], 1e-6));
        }
    }
}
