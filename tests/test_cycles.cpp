#include <catch2/catch_test_macros.hpp>

#include <chain_cycles.h>

#include <vector>

TEST_CASE("cycles", "[mtetcol]")
{
    using Index = mtetcol::Index;
    using Scalar = mtetcol::Scalar;
    using SignedIndex = mtetcol::SignedIndex;

    // clang-format off
    std::vector<Index> segments = {
        0, 1, // 0
        1, 2, // 1
        2, 3, // 2
        0, 3, // 3
        4, 5, // 4
        5, 6, // 5
        4, 6, // 6
    };
    // clang-format on

    mtetcol::Cycles cycles(7, segments);
    std::vector<SignedIndex> cycle_segments;
    std::vector<Index> cycle_indices;
    cycle_indices.push_back(0);

    SECTION("quad")
    {
        cycles.register_segment(mtetcol::signed_index(0, true));
        cycles.register_segment(mtetcol::signed_index(1, true));
        cycles.register_segment(mtetcol::signed_index(2, true));
        cycles.register_segment(mtetcol::signed_index(3, false));

        cycles.extract_cycles(cycle_segments, cycle_indices);
        size_t num_cycles = cycle_indices.size() - 1;

        REQUIRE(num_cycles == 1);
        REQUIRE(cycle_indices[1] == 4);
    }

    SECTION("Quad and triangle")
    {
        cycles.register_segment(mtetcol::signed_index(0, true));
        cycles.register_segment(mtetcol::signed_index(1, true));
        cycles.register_segment(mtetcol::signed_index(2, true));
        cycles.register_segment(mtetcol::signed_index(3, false));
        cycles.register_segment(mtetcol::signed_index(4, true));
        cycles.register_segment(mtetcol::signed_index(5, true));
        cycles.register_segment(mtetcol::signed_index(6, false));

        cycles.extract_cycles(cycle_segments, cycle_indices);
        size_t num_cycles = cycle_indices.size() - 1;

        REQUIRE(num_cycles == 2);
        for (size_t ci = 0; ci < num_cycles; ci++) {
            Index cycle_size = cycle_indices[ci + 1] - cycle_indices[ci];
            if (cycle_size == 4) {
                for (size_t i = cycle_indices[ci]; i < cycle_indices[ci + 1]; i++) {
                    mtetcol::SignedIndex si = cycle_segments[i];
                    Index seg_id = mtetcol::index(si);
                    REQUIRE(seg_id < 4);
                }
            } else {
                REQUIRE(cycle_size == 3);
                for (size_t i = cycle_indices[ci]; i < cycle_indices[ci + 1]; i++) {
                    mtetcol::SignedIndex si = cycle_segments[i];
                    Index seg_id = mtetcol::index(si);
                    REQUIRE(seg_id >= 4);
                }
            }
        }
    }
}
