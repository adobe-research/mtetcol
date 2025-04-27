#include <catch2/catch_test_macros.hpp>

#include <mtetcol/simplicial_column.h>

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
        columns.set_vertices(vertices);
        columns.set_simplices(tets);

        REQUIRE(columns.get_num_spatial_vertices() == 4);
        REQUIRE(columns.get_num_spatial_edges() == 6);
        REQUIRE(columns.get_num_spatial_triangles() == 4);
        REQUIRE(columns.get_num_spatial_tetrahedra() == 1);
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
        columns.set_vertices(vertices);
        columns.set_simplices(tets);

        REQUIRE(columns.get_num_spatial_vertices() == 5);
        REQUIRE(columns.get_num_spatial_edges() == 9);
        REQUIRE(columns.get_num_spatial_triangles() == 7);
        REQUIRE(columns.get_num_spatial_tetrahedra() == 2);
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
        columns.set_vertices(vertices);
        columns.set_simplices(tets);

        REQUIRE(columns.get_num_spatial_vertices() == 6);
        REQUIRE(columns.get_num_spatial_edges() == 11);
        REQUIRE(columns.get_num_spatial_triangles() == 8);
        REQUIRE(columns.get_num_spatial_tetrahedra() == 2);
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
        columns.set_vertices(vertices);
        columns.set_simplices(tets);

        REQUIRE(columns.get_num_spatial_vertices() == 7);
        REQUIRE(columns.get_num_spatial_edges() == 12);
        REQUIRE(columns.get_num_spatial_triangles() == 8);
        REQUIRE(columns.get_num_spatial_tetrahedra() == 2);
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
        columns.set_vertices(vertices);
        columns.set_simplices(tris);

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
        columns.set_vertices(vertices);
        columns.set_simplices(tris);

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
        columns.set_vertices(vertices);
        columns.set_simplices(tris);

        REQUIRE(columns.get_num_spatial_vertices() == 4);
        REQUIRE(columns.get_num_spatial_edges() == 5);
        REQUIRE(columns.get_num_spatial_triangles() == 2);
    }
}
