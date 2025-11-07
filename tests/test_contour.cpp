#include <catch2/catch_test_macros.hpp>

#include <mtetcol/contour.h>

#include <algorithm>
#include <set>

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

        SECTION("Isocontour: single corner")
        {
            std::vector<mtetcol::Scalar> function_values = {-1, 1, 1, 1, 1, 1, 1, 1};
            auto isocontour = contour.isocontour(function_values);
            REQUIRE(isocontour.get_num_cycles() == 1);
        }

        SECTION("Isocontour: two corners")
        {
            std::vector<mtetcol::Scalar> function_values = {-1, 1, 1, 1, 1, 1, -1, 1};
            auto isocontour = contour.isocontour(function_values);
            REQUIRE(isocontour.get_num_cycles() == 2);
        }

        SECTION("Isocontour: top")
        {
            std::vector<mtetcol::Scalar> function_values = {-1, -1, -1, -1, 1, 1, 1, 1};
            auto isocontour = contour.isocontour(function_values);
            REQUIRE(isocontour.get_num_cycles() == 1);
            REQUIRE(isocontour.get_num_segments() == 8);
            REQUIRE(isocontour.get_num_vertices() == 8);
        }
    }

    SECTION("3D square")
    {
        mtetcol::Contour<3> contour;
        REQUIRE(contour.get_num_vertices() == 0);
        REQUIRE(contour.get_num_segments() == 0);
        REQUIRE(contour.get_num_cycles() == 0);

        contour.add_vertex({0, 0, 0});
        contour.add_vertex({1, 0, 0});
        contour.add_vertex({1, 1, 0});
        contour.add_vertex({0, 1, 0});

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


        SECTION("No snapping")
        {
            std::vector<mtetcol::Scalar> function_values = {-1, -1, 1, 1};
            auto isocontour = contour.isocontour(function_values);
            REQUIRE(isocontour.get_num_vertices() == 3);
            REQUIRE(isocontour.get_num_segments() == 2);
        }

        SECTION("Snapping")
        {
            std::vector<mtetcol::Scalar> function_values = {-0.0001, -0.01, 1, 1};
            auto isocontour = contour.isocontour(function_values, {}, true);
            REQUIRE(isocontour.get_num_vertices() == 2);
            REQUIRE(isocontour.get_num_segments() == 1);
        }
        SECTION("Snapping case 2")
        {
            std::vector<mtetcol::Scalar> function_values = {0.0001, -0.01, 1, 1};
            auto isocontour = contour.isocontour(function_values, {}, true);
            REQUIRE(isocontour.get_num_vertices() == 2);
            REQUIRE(isocontour.get_num_segments() == 1);
        }
        SECTION("Snapping case 3")
        {
            std::vector<mtetcol::Scalar> function_values = {1, -0.01, 1, 1};
            auto isocontour = contour.isocontour(function_values, {}, true);
            REQUIRE(isocontour.get_num_vertices() == 1);
            REQUIRE(isocontour.get_num_segments() == 0);
        }
    }

    SECTION("3D double square")
    {
        mtetcol::Contour<3> contour;
        REQUIRE(contour.get_num_vertices() == 0);
        REQUIRE(contour.get_num_segments() == 0);
        REQUIRE(contour.get_num_cycles() == 0);

        contour.add_vertex({0, 0, 0});
        contour.add_vertex({1, 0, 0});
        contour.add_vertex({0, 1, 0});
        contour.add_vertex({1, 1, 0});
        contour.add_vertex({0, 2, 0});
        contour.add_vertex({1, 2, 0});

        contour.add_segment(0, 1); // 0
        contour.add_segment(1, 3); // 1
        contour.add_segment(3, 2); // 2
        contour.add_segment(2, 0); // 3
        contour.add_segment(3, 5); // 4
        contour.add_segment(5, 4); // 5
        contour.add_segment(2, 4); // 6

        REQUIRE(contour.get_num_segments() == 7);

        contour.add_cycle(
            {mtetcol::signed_index(0, true),
             mtetcol::signed_index(1, true),
             mtetcol::signed_index(2, true),
             mtetcol::signed_index(3, true)});
        contour.add_cycle(
            {mtetcol::signed_index(2, false),
             mtetcol::signed_index(4, true),
             mtetcol::signed_index(5, true),
             mtetcol::signed_index(6, false)});
        REQUIRE(contour.get_num_cycles() == 2);

        SECTION("No snapping")
        {
            std::vector<mtetcol::Scalar> function_values = {1, 1, 1, -1, 1, 1};
            auto isocontour = contour.isocontour(function_values);
            REQUIRE(isocontour.get_num_vertices() == 3);
            REQUIRE(isocontour.get_num_segments() == 2);
        }

        SECTION("Snapping")
        {
            std::vector<mtetcol::Scalar> function_values = {1, 1, -0.0001, -0.001, 1, 1};
            auto isocontour = contour.isocontour(function_values, {}, true);
            REQUIRE(isocontour.get_num_vertices() == 2);
            REQUIRE(isocontour.get_num_segments() == 2);
        }

        SECTION("Snapping case 2")
        {
            std::vector<mtetcol::Scalar> function_values = {1, 1, 0.0001, -0.01, 1, 1};
            auto isocontour = contour.isocontour(function_values, {}, true);
            REQUIRE(isocontour.get_num_vertices() == 2);
            REQUIRE(isocontour.get_num_segments() == 2);
        }

        SECTION("Snapping case 3")
        {
            std::vector<mtetcol::Scalar> function_values = {1, 1, 1, -0.01, 1, 1};
            auto isocontour = contour.isocontour(function_values, {}, true);
            REQUIRE(isocontour.get_num_vertices() == 1);
            REQUIRE(isocontour.get_num_segments() == 0);
        }

        SECTION("Snapping case 4")
        {
            std::vector<mtetcol::Scalar> function_values = {1, 1, 1, -1, 1, 1000};
            auto isocontour = contour.isocontour(function_values, {}, true);
            REQUIRE(isocontour.get_num_vertices() == 3);
            REQUIRE(isocontour.get_num_segments() == 2);
        }
    }
}

TEST_CASE("extract_simple_loops", "[mtetcol]")
{
    using namespace mtetcol;

    // Helper function to check if a subcycle has duplicate vertices
    auto has_duplicate_vertices = [](const auto& subcycle, auto get_segment) {
        std::vector<Index> vertices;
        for (auto sid : subcycle) {
            Index seg_idx = index(sid);
            bool seg_ori = orientation(sid);
            auto seg = get_segment(seg_idx);
            Index v0 = seg[seg_ori ? 0 : 1];
            vertices.push_back(v0);
        }

        // Check for duplicates
        std::sort(vertices.begin(), vertices.end());
        auto it = std::adjacent_find(vertices.begin(), vertices.end());
        return it != vertices.end();
    };

    SECTION("Simple cycle without duplicate vertices")
    {
        // Create a simple triangle cycle: 0->1->2->0
        // Segments: 0:(0,1), 1:(1,2), 2:(2,0)
        std::vector<SignedIndex> cycle = {
            signed_index(0, true), // segment 0: 0->1
            signed_index(1, true), // segment 1: 1->2
            signed_index(2, true) // segment 2: 2->0
        };

        auto get_segment = [](Index seg_idx) -> std::array<Index, 2> {
            if (seg_idx == 0) return {0, 1};
            if (seg_idx == 1) return {1, 2};
            if (seg_idx == 2) return {2, 0};
            return {0, 0};
        };

        auto subcycles = extract_simple_loops(cycle, get_segment);

        REQUIRE(subcycles.size() == 1);
        REQUIRE(subcycles[0].size() == 3);
        REQUIRE_FALSE(has_duplicate_vertices(subcycles[0], get_segment));
    }

    SECTION("Cycle with one duplicate vertex (figure-8)")
    {
        // Create a figure-8 cycle with junction at vertex 0:
        // Loop 1: 0->1->2->0
        // Loop 2: 0->3->4->0
        // Segments: 0:(0,1), 1:(1,2), 2:(2,0), 3:(0,3), 4:(3,4), 5:(4,0)
        std::vector<SignedIndex> cycle = {
            signed_index(0, true), // segment 0: 0->1
            signed_index(1, true), // segment 1: 1->2
            signed_index(2, true), // segment 2: 2->0
            signed_index(3, true), // segment 3: 0->3
            signed_index(4, true), // segment 4: 3->4
            signed_index(5, true) // segment 5: 4->0
        };

        auto get_segment = [](Index seg_idx) -> std::array<Index, 2> {
            if (seg_idx == 0) return {0, 1};
            if (seg_idx == 1) return {1, 2};
            if (seg_idx == 2) return {2, 0};
            if (seg_idx == 3) return {0, 3};
            if (seg_idx == 4) return {3, 4};
            if (seg_idx == 5) return {4, 0};
            return {0, 0};
        };

        auto subcycles = extract_simple_loops(cycle, get_segment);

        REQUIRE(subcycles.size() == 2);

        // Check that each subcycle has no duplicate vertices
        for (const auto& subcycle : subcycles) {
            REQUIRE(subcycle.size() >= 3);
            REQUIRE_FALSE(has_duplicate_vertices(subcycle, get_segment));
        }
    }

    SECTION("Cycle with multiple duplicate vertices")
    {
        // Create a more complex cycle with multiple junctions:
        // Vertex 1 and 3 are duplicates.
        std::vector<SignedIndex> cycle = {
            signed_index(0, false), // segment 0: 0->1
            signed_index(1, true), // segment 1: 1->2
            signed_index(2, true), // segment 2: 2->3
            signed_index(3, true), // segment 3: 3->1
            signed_index(4, true), // segment 4: 1->4
            signed_index(5, true), // segment 5: 4->3
            signed_index(6, true), // segment 6: 3->5
            signed_index(7, true), // segment 7: 5->0
        };

        auto get_segment = [](Index seg_idx) -> std::array<Index, 2> {
            if (seg_idx == 0) return {1, 0};
            if (seg_idx == 1) return {1, 2};
            if (seg_idx == 2) return {2, 3};
            if (seg_idx == 3) return {3, 1};
            if (seg_idx == 4) return {1, 4};
            if (seg_idx == 5) return {4, 3};
            if (seg_idx == 6) return {3, 5};
            if (seg_idx == 7) return {5, 0};
            return {0, 0};
        };

        auto subcycles = extract_simple_loops(cycle, get_segment);

        REQUIRE(subcycles.size() == 2);

        // Check that each subcycle has no duplicate vertices
        for (const auto& subcycle : subcycles) {
            REQUIRE(subcycle.size() >= 3);
            REQUIRE_FALSE(has_duplicate_vertices(subcycle, get_segment));
        }

        // Verify that all segments are used
        std::set<Index> used_segments;
        for (const auto& subcycle : subcycles) {
            for (auto sid : subcycle) {
                used_segments.insert(index(sid));
            }
        }
        REQUIRE(used_segments.size() == cycle.size());
    }

    SECTION("Cycle with duplicate but opposite segments")
    {
        std::vector<SignedIndex> cycle = {
            signed_index(0, true), // segment 0: 22864 -> 22866
            signed_index(1, true), // segment 1: 22866 -> 109833
            signed_index(2, true), // segment 2: 109833 -> 109834
            signed_index(3, true), // segment 3: 109834 -> 22867
            signed_index(4, true), // segment 4: 22867 -> 22865
            signed_index(5, true), // segment 5: 22865 -> 109834
            signed_index(6, true), // segment 6: 109834 -> 109833
            signed_index(7, true), // segment 7: 109833 -> 22864
        };

        auto get_segment = [](Index seg_idx) -> std::array<Index, 2> {
            if (seg_idx == 0) return {22864, 22866};
            if (seg_idx == 1) return {22866, 109833};
            if (seg_idx == 2) return {109833, 109834};
            if (seg_idx == 3) return {109834, 22867};
            if (seg_idx == 4) return {22867, 22865};
            if (seg_idx == 5) return {22865, 109834};
            if (seg_idx == 6) return {109834, 109833};
            if (seg_idx == 7) return {109833, 22864};
            return {0, 0};
        };

        auto subcycles = extract_simple_loops(cycle, get_segment);

        REQUIRE(subcycles.size() == 2);

        // Check that each subcycle has no duplicate vertices
        for (const auto& subcycle : subcycles) {
            REQUIRE(subcycle.size() >= 3);
            REQUIRE_FALSE(has_duplicate_vertices(subcycle, get_segment));
        }

        // Verify that all non-canceled segments are used
        // Segments 2 and 6 cancel (opposite directions: 109833↔109834)
        std::set<Index> used_segments;
        for (const auto& subcycle : subcycles) {
            for (auto sid : subcycle) {
                used_segments.insert(index(sid));
            }
        }
        // 8 original segments - 2 canceled = 6 segments in output
        REQUIRE(used_segments.size() == 6);
        // Verify canceled segments are not in output
        REQUIRE(used_segments.count(2) == 0);
        REQUIRE(used_segments.count(6) == 0);
    }

    SECTION("Complex cycle with degenerate loops")
    {
        std::vector<SignedIndex> cycle = {
            signed_index(0, false),  // segment 0: 101217 -> 101215
            signed_index(1, false),  // segment 1: 41702 -> 101217
            signed_index(2, true),   // segment 2: 41702 -> 20246
            signed_index(3, false),  // segment 3: 101218 -> 20246
            signed_index(4, false),  // segment 4: 20246 -> 101218
            signed_index(5, true),   // segment 5: 20246 -> 41702
            signed_index(6, false),  // segment 6: 101215 -> 41702
        };

        auto get_segment = [](Index seg_idx) -> std::array<Index, 2> {
            if (seg_idx == 0) return {101217, 101215};
            if (seg_idx == 1) return {41702, 101217};
            if (seg_idx == 2) return {41702, 20246};
            if (seg_idx == 3) return {101218, 20246};
            if (seg_idx == 4) return {20246, 101218};
            if (seg_idx == 5) return {20246, 41702};
            if (seg_idx == 6) return {101215, 41702};
            return {0, 0};
        };

        auto subcycles = extract_simple_loops(cycle, get_segment);

        REQUIRE(subcycles.size() == 1);

        // Check that each subcycle has no duplicate vertices
        for (const auto& subcycle : subcycles) {
            REQUIRE(subcycle.size() >= 3);
            REQUIRE_FALSE(has_duplicate_vertices(subcycle, get_segment));
        }
    }

    SECTION("Empty cycle")
    {
        std::vector<SignedIndex> cycle;
        auto get_segment = [](Index) -> std::array<Index, 2> { return {0, 0}; };

        auto subcycles = extract_simple_loops(cycle, get_segment);

        REQUIRE(subcycles.empty());
    }

    SECTION("Cycle with reversed segments")
    {
        // Create a triangle with reversed segments: 0<-1<-2<-0
        // Segments: 0:(0,1), 1:(1,2), 2:(2,0) but traversed backwards
        std::vector<SignedIndex> cycle = {
            signed_index(0, false), // segment 0: 1->0
            signed_index(1, false), // segment 1: 2->1
            signed_index(2, false) // segment 2: 0->2
        };

        auto get_segment = [](Index seg_idx) -> std::array<Index, 2> {
            if (seg_idx == 0) return {0, 1};
            if (seg_idx == 1) return {1, 2};
            if (seg_idx == 2) return {2, 0};
            return {0, 0};
        };

        auto subcycles = extract_simple_loops(cycle, get_segment);

        REQUIRE(subcycles.size() == 1);
        REQUIRE(subcycles[0].size() == 3);
        REQUIRE_FALSE(has_duplicate_vertices(subcycles[0], get_segment));
    }
}
