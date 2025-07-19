#include <catch2/catch_test_macros.hpp>

#include <mtetcol/nondisjoint_cycles.h>

#include <vector>

TEST_CASE("nondisjoint_cycles", "[mtetcol]")
{
    using Index = mtetcol::Index;
    using Scalar = mtetcol::Scalar;
    using SignedIndex = mtetcol::SignedIndex;

    // clang-format off
    std::vector<Index> segments = {
        0, 1, // 0
        1, 2, // 1
        2, 0, // 2
        2, 3, // 3
        3, 4, // 4
        4, 5, // 5
        5, 3, // 6
    };
    // clang-format on

    mtetcol::NonDisjointCycles cycles(segments);
    std::vector<SignedIndex> cycle_segments;
    std::vector<Index> cycle_indices;
    cycle_indices.push_back(0);

    SECTION("Single triangle")
    {
        cycles.register_segment(mtetcol::signed_index(0, true));
        cycles.register_segment(mtetcol::signed_index(1, true));
        cycles.register_segment(mtetcol::signed_index(2, true));

        cycles.extract_cycles(cycle_segments, cycle_indices);
        size_t num_cycles = cycle_indices.size() - 1;

        REQUIRE(num_cycles == 1);
        REQUIRE(cycle_indices[1] == 3);
    }

    SECTION("Two disjoint triangles")
    {
        cycles.register_segment(mtetcol::signed_index(0, true));
        cycles.register_segment(mtetcol::signed_index(1, true));
        cycles.register_segment(mtetcol::signed_index(2, true));

        cycles.register_segment(mtetcol::signed_index(4, true));
        cycles.register_segment(mtetcol::signed_index(5, true));
        cycles.register_segment(mtetcol::signed_index(6, true));

        cycles.extract_cycles(cycle_segments, cycle_indices);
        size_t num_cycles = cycle_indices.size() - 1;

        REQUIRE(num_cycles == 2);
        REQUIRE(cycle_indices[1] == 3);
        REQUIRE(cycle_indices[2] == 6);
    }

    SECTION("Two connected cycles")
    {
        cycles.register_segment(mtetcol::signed_index(0, true));
        cycles.register_segment(mtetcol::signed_index(1, true));
        cycles.register_segment(mtetcol::signed_index(2, true));

        cycles.register_segment(mtetcol::signed_index(3, true));
        cycles.register_segment(mtetcol::signed_index(3, false));
        cycles.extract_cycles(cycle_segments, cycle_indices);

        size_t num_cycles = cycle_indices.size() - 1;

        REQUIRE(num_cycles == 2);
        REQUIRE(cycle_indices[1] == 3);
        REQUIRE(cycle_indices[2] == 5);
    }

    SECTION("Three connected cycles")
    {
        cycles.register_segment(mtetcol::signed_index(0, true));
        cycles.register_segment(mtetcol::signed_index(1, true));
        cycles.register_segment(mtetcol::signed_index(2, true));

        cycles.register_segment(mtetcol::signed_index(3, true));
        cycles.register_segment(mtetcol::signed_index(3, false));

        cycles.register_segment(mtetcol::signed_index(4, true));
        cycles.register_segment(mtetcol::signed_index(5, true));
        cycles.register_segment(mtetcol::signed_index(6, true));

        cycles.extract_cycles(cycle_segments, cycle_indices);

        size_t num_cycles = cycle_indices.size() - 1;

        REQUIRE(num_cycles == 3);
        REQUIRE(cycle_indices[1] == 3);
        REQUIRE(cycle_indices[2] == 5);
        REQUIRE(cycle_indices[3] == 8);
    }
}
