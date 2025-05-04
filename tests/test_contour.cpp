#include <catch2/catch_test_macros.hpp>

#include <mtetcol/contour.h>

TEST_CASE("contour", "[mtetcol]")
{
    SECTION("4D tet")
    {
        mtetcol::Contour<4> contour;
        REQUIRE(contour.get_num_vertices() == 0);
        REQUIRE(contour.get_num_segments() == 0);
        REQUIRE(contour.get_num_cycles() == 0);
        REQUIRE(contour.get_num_polyhedra() == 0);

        contour.add_vertex({0, 0, 0, 0});
        contour.add_vertex({1, 0, 0, 0});
        contour.add_vertex({0, 1, 0, 0});
        contour.add_vertex({0, 0, 1, 0});

        REQUIRE(contour.get_num_vertices() == 4);
        REQUIRE(contour.get_vertex(1)[0] == 1);
        REQUIRE(contour.get_vertex(0)[3] == 0);

        contour.add_segment(0, 1); // 0
        contour.add_segment(0, 2); // 1
        contour.add_segment(0, 3); // 2
        contour.add_segment(1, 2); // 3
        contour.add_segment(1, 3); // 4
        contour.add_segment(2, 3); // 5

        REQUIRE(contour.get_num_segments() == 6);

        contour.add_cycle(
            {mtetcol::signed_index(1, true),
             mtetcol::signed_index(3, false),
             mtetcol::signed_index(0, false)});
        contour.add_cycle(
            {mtetcol::signed_index(3, true),
             mtetcol::signed_index(5, true),
             mtetcol::signed_index(4, false)});
        contour.add_cycle(
            {mtetcol::signed_index(2, true),
             mtetcol::signed_index(5, false),
             mtetcol::signed_index(1, false)});
        contour.add_cycle(
            {mtetcol::signed_index(0, true),
             mtetcol::signed_index(4, true),
             mtetcol::signed_index(2, false)});
        REQUIRE(contour.get_num_cycles() == 4);

        contour.add_polyhedron(
            {mtetcol::signed_index(0, true),
             mtetcol::signed_index(1, true),
             mtetcol::signed_index(2, true),
             mtetcol::signed_index(3, true)});
        REQUIRE(contour.get_num_polyhedra() == 1);
    }

    SECTION("4D cube")
    {
        mtetcol::Contour<4> contour;
        REQUIRE(contour.get_num_vertices() == 0);
        REQUIRE(contour.get_num_segments() == 0);
        REQUIRE(contour.get_num_cycles() == 0);
        REQUIRE(contour.get_num_polyhedra() == 0);

        contour.add_vertex({0, 0, 0, 0});
        contour.add_vertex({1, 0, 0, 0});
        contour.add_vertex({1, 1, 0, 0});
        contour.add_vertex({0, 1, 0, 0});
        contour.add_vertex({0, 0, 1, 0});
        contour.add_vertex({1, 0, 1, 0});
        contour.add_vertex({1, 1, 1, 0});
        contour.add_vertex({0, 1, 1, 0});

        contour.add_segment(0, 1);
        contour.add_segment(1, 2);
        contour.add_segment(2, 3);
        contour.add_segment(3, 0);

        contour.add_segment(4, 5);
        contour.add_segment(5, 6);
        contour.add_segment(6, 7);
        contour.add_segment(7, 4);

        contour.add_segment(1, 5);
        contour.add_segment(2, 6);
        contour.add_segment(3, 7);
        contour.add_segment(0, 4);

        contour.add_cycle(
            {mtetcol::signed_index(0, false),
             mtetcol::signed_index(3, false),
             mtetcol::signed_index(2, false),
             mtetcol::signed_index(1, false)});
        contour.add_cycle(
            {mtetcol::signed_index(4, true),
             mtetcol::signed_index(5, true),
             mtetcol::signed_index(6, true),
             mtetcol::signed_index(7, true)});
        contour.add_cycle(
            {mtetcol::signed_index(1, true),
             mtetcol::signed_index(9, true),
             mtetcol::signed_index(5, false),
             mtetcol::signed_index(8, false)});
        contour.add_cycle(
            {mtetcol::signed_index(3, true),
             mtetcol::signed_index(11, true),
             mtetcol::signed_index(7, false),
             mtetcol::signed_index(10, false)});
        contour.add_cycle(
            {mtetcol::signed_index(0, true),
             mtetcol::signed_index(8, true),
             mtetcol::signed_index(4, false),
             mtetcol::signed_index(11, false)});
        contour.add_cycle(
            {mtetcol::signed_index(2, true),
             mtetcol::signed_index(10, true),
             mtetcol::signed_index(6, false),
             mtetcol::signed_index(9, false)});

        contour.add_polyhedron(
            {mtetcol::signed_index(0, true),
             mtetcol::signed_index(1, true),
             mtetcol::signed_index(2, true),
             mtetcol::signed_index(3, true),
             mtetcol::signed_index(4, true),
             mtetcol::signed_index(5, true)});

        REQUIRE(contour.get_num_vertices() == 8);
        REQUIRE(contour.get_num_segments() == 12);
        REQUIRE(contour.get_num_cycles() == 6);
        REQUIRE(contour.get_num_polyhedra() == 1);

        contour.triangulate_cycles();

        REQUIRE(contour.get_num_vertices() == 8);
        REQUIRE(contour.get_num_segments() == 18);
        REQUIRE(contour.get_num_cycles() == 12);
        REQUIRE(contour.get_num_polyhedra() == 1);
    }

    SECTION("3D")
    {
        mtetcol::Contour<3> contour;
        REQUIRE(contour.get_num_vertices() == 0);
        REQUIRE(contour.get_num_segments() == 0);
        REQUIRE(contour.get_num_cycles() == 0);

        contour.add_vertex({0, 0, 0});
        contour.add_vertex({1, 0, 0});
        contour.add_vertex({0, 1, 0});
        contour.add_vertex({1, 1, 0});

        contour.add_segment(0, 1);
        contour.add_segment(1, 2);
        contour.add_segment(2, 3);
        contour.add_segment(3, 0);

        REQUIRE(contour.get_num_segments() == 4);

        contour.add_cycle(
            {mtetcol::signed_index(0, true),
             mtetcol::signed_index(1, true),
             mtetcol::signed_index(2, true),
             mtetcol::signed_index(3, true)});
        REQUIRE(contour.get_num_cycles() == 1);

        contour.triangulate_cycles();
        REQUIRE(contour.get_num_cycles() == 2);
    }
}
