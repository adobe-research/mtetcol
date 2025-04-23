#include <catch2/catch_test_macros.hpp>

#include <mtetcol/contour.h>

TEST_CASE("contour", "[mtetcol]")
{
    SECTION("4D") {
        mtetcol::Contour<4> contour;
        REQUIRE(contour.get_num_vertices() == 0);
        REQUIRE(contour.get_num_segments() == 0);
        REQUIRE(contour.get_num_cycles() == 0);
        REQUIRE(contour.get_num_polyhedra() == 0);
    }
}
