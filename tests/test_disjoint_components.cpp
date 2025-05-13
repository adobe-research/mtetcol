#include <catch2/catch_test_macros.hpp>

#include <mtetcol/disjoint_components.h>

#include <vector>

TEST_CASE("disjoint_components", "[mtetcol]")
{
    using Index = mtetcol::Index;
    using Scalar = mtetcol::Scalar;
    using SignedIndex = mtetcol::SignedIndex;

    std::vector<SignedIndex> cycles = {
        mtetcol::signed_index(0, true),
        mtetcol::signed_index(1, true),
        mtetcol::signed_index(2, true),

        mtetcol::signed_index(2, false),
        mtetcol::signed_index(3, true),
        mtetcol::signed_index(4, true),

        mtetcol::signed_index(0, true),
        mtetcol::signed_index(4, true),
        mtetcol::signed_index(3, true),
        mtetcol::signed_index(1, true),

        mtetcol::signed_index(5, true),
        mtetcol::signed_index(6, true),
        mtetcol::signed_index(7, true),
        mtetcol::signed_index(8, true),

        mtetcol::signed_index(7, false),
        mtetcol::signed_index(9, true),
        mtetcol::signed_index(10, true),

        mtetcol::signed_index(8, false),
        mtetcol::signed_index(10, false),
        mtetcol::signed_index(9, false),
        mtetcol::signed_index(6, false),
        mtetcol::signed_index(5, false),
    };
    std::vector<Index> cycle_indices = {0, 3, 6, 10, 14, 17, 22};

    mtetcol::DisjointComponents components(11, cycles, cycle_indices);
    std::vector<SignedIndex> polyhedra;
    std::vector<Index> polyhedron_indices;
    polyhedron_indices.push_back(0);

    SECTION("First comp")
    {
        components.register_cycle(mtetcol::signed_index(0, true));
        components.register_cycle(mtetcol::signed_index(1, true));
        components.register_cycle(mtetcol::signed_index(2, false));

        components.extract_components(polyhedra, polyhedron_indices);

        size_t num_polyhedra = polyhedron_indices.size() - 1;
        REQUIRE(num_polyhedra == 1);
        REQUIRE(polyhedron_indices[1] == 3);
    }

    SECTION("Second comp")
    {
        components.register_cycle(mtetcol::signed_index(3, true));
        components.register_cycle(mtetcol::signed_index(4, true));
        components.register_cycle(mtetcol::signed_index(5, true));

        components.extract_components(polyhedra, polyhedron_indices);

        size_t num_polyhedra = polyhedron_indices.size() - 1;
        REQUIRE(num_polyhedra == 1);
        REQUIRE(polyhedron_indices[1] == 3);
    }

    SECTION("Both comps")
    {
        components.register_cycle(mtetcol::signed_index(0, true));
        components.register_cycle(mtetcol::signed_index(1, true));
        components.register_cycle(mtetcol::signed_index(2, false));
        components.register_cycle(mtetcol::signed_index(3, true));
        components.register_cycle(mtetcol::signed_index(4, true));
        components.register_cycle(mtetcol::signed_index(5, true));

        components.extract_components(polyhedra, polyhedron_indices);

        size_t num_polyhedra = polyhedron_indices.size() - 1;
        REQUIRE(num_polyhedra == 2);
    }
}
