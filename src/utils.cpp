#include "utils.h"
#include "hashmap.h"

#include <mtetcol/disjoint_components.h>
#include <mtetcol/disjoint_cycles.h>
#include <mtetcol/logger.h>

#include <SmallVector.h>

#include <array>
#include <functional>
#include <stdexcept>

namespace mtetcol {

bool tet_mesh_is_manifold(std::span<Index> tets)
{
    assert(tets.size() % 4 == 0);
    size_t num_tets = tets.size() / 4;
    TriangleMap triangle_count;

    for (size_t i = 0; i < num_tets; i++) {
        auto tet = tets.subspan(i * 4, 4);
        // List all 4 faces of the tetrahedron
        std::array<std::array<Index, 3>, 4> faces = {{
            {tet[0], tet[2], tet[1]},
            {tet[1], tet[2], tet[3]},
            {tet[0], tet[1], tet[3]},
            {tet[0], tet[3], tet[2]}
        }};
        for (const auto& face : faces) {
            auto [it, _] = triangle_count.try_emplace(face, 0);
            it->second += 1;
        }
    }

    for (auto [triangle, count] : triangle_count) {
        if (count > 2) {
            logger().error(
                "Non-manifold triangle found: [{}, {}, {}] appears {} times",
                triangle[0],
                triangle[1],
                triangle[2],
                count);
            return false;
        }
    }
    return true;
}

void extract_vertex_zero_crossing(
    std::span<const Scalar> time_samples,
    std::span<const Scalar> function_values,
    Scalar value,
    bool cyclic,
    std::vector<Scalar>& zero_crossing_times)
{
    if (time_samples.empty() || function_values.empty()) {
        throw std::invalid_argument("Time samples and function values must not be empty");
    }
    if (time_samples.size() != function_values.size()) {
        throw std::invalid_argument("Time samples and function values must have the same size");
    }
    if (time_samples.size() < 2) {
        throw std::invalid_argument("At least 2 samples are required");
    }
    if (time_samples.front() != 0 || time_samples.back() != 1) {
        throw std::invalid_argument("Time samples must start at 0 and end at 1");
    }

    if (cyclic) {
        // Do some sanity check for cyclic time samples
        if (function_values.front() != function_values.back()) {
            Scalar f0 = function_values.front();
            Scalar f1 = function_values.back();

            if ((f0 < value && f1 >= value) || (f0 >= value && f1 < value)) {
                logger().error(
                    "Cyclic function have different signs at the t=0 and t=1: {}, {}",
                    f0,
                    f1);
                throw std::invalid_argument(
                    "Cyclic time samples with different function signs at the start and end");
            } else if (f0 != f1) {
                // Not desirable, but not a blocker.
                logger().warn(
                    "Cyclic function have different values at t=0 and t=1: {}, {}",
                    f0,
                    f1);
            }
        }
    }

    zero_crossing_times.reserve(zero_crossing_times.size() + time_samples.size());

    if (!cyclic && function_values[0] >= value) {
        zero_crossing_times.push_back(time_samples[0]);
    }

    for (size_t i = 0; i < time_samples.size(); i++) {
        size_t next_i = (i + 1) % time_samples.size();
        if (!cyclic && next_i == 0) break;

        if ((function_values[i] < value && function_values[next_i] < value) ||
            (function_values[i] >= value && function_values[next_i] >= value)) {
            continue;
        } else {
            const Scalar t_curr = time_samples[i];
            const Scalar t_next = time_samples[next_i];
            const Scalar val_curr = function_values[i] - value;
            const Scalar val_next = function_values[next_i] - value;
            assert(cyclic || t_next > t_curr);
            assert(!cyclic || next_i == 0 || t_next > t_curr);

            zero_crossing_times.push_back(
                t_curr + (t_next - t_curr) * (-val_curr) / (val_next - val_curr));
        }
    }

    if (!cyclic && function_values[time_samples.size() - 1] < value) {
        zero_crossing_times.push_back(time_samples[time_samples.size() - 1]);
    }
}

std::tuple<std::vector<Scalar>, std::vector<size_t>, std::vector<bool>> extract_contour_vertices(
    const std::vector<Scalar>& time_samples,
    const std::vector<Scalar>& function_values,
    const std::vector<Index>& vertex_start_indices,
    Scalar value,
    bool cyclic)
{
    if (vertex_start_indices.empty()) {
        throw std::invalid_argument("Vertex start indices must not be empty");
    }
    if (vertex_start_indices.back() > time_samples.size()) {
        throw std::invalid_argument("Vertex start indices exceed time samples size");
    }

    size_t num_vertices = vertex_start_indices.size() - 1;

    std::vector<Scalar> zero_crossing_times;
    std::vector<size_t> zero_crossing_indices;
    std::vector<bool> initial_signs;
    zero_crossing_times.reserve(time_samples.size());
    zero_crossing_indices.reserve(num_vertices + 1);
    zero_crossing_indices.push_back(0);
    initial_signs.reserve(num_vertices);

    // Process each vertex
    for (size_t vi = 0; vi < num_vertices; vi++) {
        const Index idx_begin = vertex_start_indices[vi];
        const Index idx_end = vertex_start_indices[vi + 1];

        if (idx_begin >= idx_end) {
            throw std::invalid_argument("Invalid vertex start indices");
        }

        std::span<const Scalar> time_samples_i(
            time_samples.data() + idx_begin,
            idx_end - idx_begin);
        std::span<const Scalar> function_values_i(
            function_values.data() + idx_begin,
            idx_end - idx_begin);

        initial_signs.push_back(function_values_i[0] >= value);

        extract_vertex_zero_crossing(
            time_samples_i,
            function_values_i,
            value,
            cyclic,
            zero_crossing_times);
        zero_crossing_indices.push_back(zero_crossing_times.size());
    }

    return {
        std::move(zero_crossing_times),
        std::move(zero_crossing_indices),
        std::move(initial_signs)};
}

std::tuple<
    std::vector<Index>,
    std::vector<SignedIndex>,
    std::vector<Index>,
    std::vector<bool>,
    std::vector<int8_t>>
extract_contour_segments(
    const std::vector<Scalar>& contour_times,
    const std::vector<size_t>& contour_time_indices,
    const std::vector<bool>& initial_signs,
    const std::vector<Index>& edges,
    bool cyclic)
{
    size_t num_edges = edges.size() / 2;

    OrientedEdgeMap vertical_edge_map;
    std::vector<Index> segments;
    std::vector<SignedIndex> segment_on_edges;
    std::vector<Index> segment_on_edges_indices;
    std::vector<bool> is_simple(num_edges, true);
    std::vector<int8_t> shift(num_edges, 0);
    segments.reserve(contour_times.size() * 2);
    segment_on_edges.reserve(contour_times.size()); // Just a guess
    segment_on_edges_indices.reserve(num_edges + 1);
    segment_on_edges_indices.push_back(0);

    struct LocalIndex
    {
        Index vertex_index;
        Index time_index;
    };

    auto merge_time_samples = [](Index v0,
                                 Index v1,
                                 std::span<const Scalar> times_0,
                                 std::span<const Scalar> times_1,
                                 auto& local_indices) {
        size_t num_times_0 = times_0.size();
        size_t num_times_1 = times_1.size();
        assert((num_times_0 + num_times_1) % 2 == 0);

        local_indices.clear();
        local_indices.reserve(num_times_0 + num_times_1);

        size_t i0 = 0, i1 = 0;
        while (i0 + i1 < num_times_0 + num_times_1) {
            if (i0 < num_times_0 && i1 < num_times_1) {
                if (times_0[i0] < times_1[i1]) {
                    local_indices.push_back({v0, static_cast<Index>(i0)});
                    i0++;
                } else {
                    local_indices.push_back({v1, static_cast<Index>(i1)});
                    i1++;
                }
            } else if (i0 < num_times_0) {
                local_indices.push_back({v0, static_cast<Index>(i0)});
                i0++;
            } else if (i1 < num_times_1) {
                local_indices.push_back({v1, static_cast<Index>(i1)});
                i1++;
            }
        }
    };


    llvm_vecsmall::SmallVector<LocalIndex, 1024> local_indices;
    for (size_t ei = 0; ei < num_edges; ei++) {
        Index v0 = edges[ei * 2];
        Index v1 = edges[ei * 2 + 1];

        std::span<const Scalar> times_0(
            contour_times.data() + contour_time_indices[v0],
            contour_time_indices[v0 + 1] - contour_time_indices[v0]);
        std::span<const Scalar> times_1(
            contour_times.data() + contour_time_indices[v1],
            contour_time_indices[v1 + 1] - contour_time_indices[v1]);

        merge_time_samples(v0, v1, times_0, times_1, local_indices);
        assert(local_indices.size() % 2 == 0);

        size_t offset_v0 = 0;
        size_t offset_v1 = 0;

        // Note:
        //
        // Parity represents the function sign change for the first point of merged `local_indices`
        // __after__ the offset.
        //
        // * parity == 0 means function is changing from negative to positive
        // * parity == 1 means function is changing from positive to negative
        //
        // In fact, after offsetting, the initial sign for v0 should be the same as the initial sign
        // for v1. So, parity provides information about the initial sign for both v0 and v1 after
        // offset.
        size_t parity = 0;

        if (cyclic) {
            if (initial_signs[v0] != initial_signs[v1]) {
                // Different sign at the beginning
                if (local_indices[0].vertex_index == v0) {
                    // First merged time sample is on v0.
                    // The offset will shift the samples, and the first time sample available to
                    // match on v0 will be the second one, which have the opposite initial sign.
                    if (!initial_signs[v0]) {
                        parity = 1;
                    }
                    shift[ei] = -1;
                    offset_v0 = 1;
                } else {
                    // First merged time sample is on v1.
                    // Shifting does not change the initial sign of v0's first sample.
                    if (initial_signs[v0]) {
                        parity = 1;
                    }
                    shift[ei] = 1;
                    offset_v1 = 1;
                }
            } else {
                if (initial_signs[v0]) {
                    parity = 1;
                }
            }
        }

        auto add_segment = [&](Index v0, Index v1) -> Index {
            segments.push_back(v0);
            segments.push_back(v1);
            return static_cast<Index>(segments.size() / 2 - 1);
        };

        size_t num_segments = local_indices.size() / 2;
        for (size_t si = 0; si < num_segments; si++) {
            const auto& lid0 = local_indices[(si * 2 + offset_v0 + offset_v1) % local_indices.size()];
            const auto& lid1 = local_indices[(si * 2 + 1 + offset_v0 + offset_v1) % local_indices.size()];

            Index p0 = contour_time_indices[lid0.vertex_index] + lid0.time_index;
            Index p1 = contour_time_indices[lid1.vertex_index] + lid1.time_index;
            Index seg_id = invalid_index;
            bool seg_ori = true;

            // Add segment
            if (lid0.vertex_index == lid1.vertex_index) {
                // Vertical segment on the same vertex
                if (vertical_edge_map.contains({p0, p1})) {
                    seg_id = vertical_edge_map[{p0, p1}];
                } else if (vertical_edge_map.contains({p1, p0})) {
                    seg_id = vertical_edge_map[{p1, p0}];
                    seg_ori = false;
                } else {
                    vertical_edge_map[{p0, p1}] = segments.size() / 2;
                    seg_id = add_segment(p0, p1);
                }
                is_simple[ei] = false;
            } else {
                // cross edge
                seg_id = add_segment(p0, p1);
            }
            assert(seg_id != invalid_index);

            // The seg_id and seg_ori should now reflect the segment [p0, p1]
            assert(segments[seg_id * 2] == (seg_ori ? p0: p1));
            assert(segments[seg_id * 2 + 1] == (seg_ori ? p1: p0));

            // Compute segment orientation
            // Pair local indices up based on even/odd parity.
            if ((lid0.vertex_index == v0 && (lid0.time_index + offset_v0) % 2 == parity) ||
                (lid0.vertex_index == v1 && (lid0.time_index + offset_v1) % 2 != parity)) {
                // lid0 -> lid1, Nothing to do.
            } else {
                // lid1 -> lid0, flipping orientation
                seg_ori = !seg_ori; //p1 < p0;
            }

            segment_on_edges.push_back(signed_index(seg_id, seg_ori));
        }
        segment_on_edges_indices.push_back(static_cast<Index>(segment_on_edges.size()));
    }

    return {
        std::move(segments),
        std::move(segment_on_edges),
        std::move(segment_on_edges_indices),
        std::move(is_simple),
        std::move(shift),
    };
}

auto extract_contour_cycles(
    const size_t num_contour_vertices,
    const std::vector<Index>& contour_segments,
    const std::vector<SignedIndex>& contour_segment_on_edges,
    const std::vector<Index>& contour_segment_on_edges_indices,
    const std::vector<bool>& edge_is_simple,
    const std::vector<int8_t>& edge_shift,
    const std::vector<Index>& edges,
    const std::vector<SignedIndex>& triangles) -> std::
    tuple<std::vector<SignedIndex>, std::vector<Index>, std::vector<Index>, std::vector<bool>>
{
    std::vector<SignedIndex> contour_cycles;
    std::vector<Index> contour_cycle_indices;
    std::vector<Index> contour_cycle_triangle_indices;
    std::vector<bool> triangle_is_simple(triangles.size() / 3, true);

    // Reserve space based on number of triangles and expected cycles per triangle
    contour_cycles.reserve(triangles.size() / 3 * 2); // Rough estimate: 2 cycles per triangle
    contour_cycle_indices.reserve(triangles.size() / 3 + 1);
    contour_cycle_triangle_indices.reserve(triangles.size() / 3 + 1);

    contour_cycle_indices.push_back(0);
    contour_cycle_triangle_indices.push_back(0);

    DisjointCycles cycle_engine(num_contour_vertices, contour_segments);

    auto register_edge = [&](Index e_id, bool ori) {
        Index seg_begin = contour_segment_on_edges_indices[e_id];
        Index seg_end = contour_segment_on_edges_indices[e_id + 1];
        assert(seg_begin <= seg_end);

        for (Index i = seg_begin; i < seg_end; i++) {
            SignedIndex si = contour_segment_on_edges[i];
            if (!ori) si = -si;
            cycle_engine.register_segment(si);
        }
    };

    auto is_connected = [&](SignedIndex seg_01, SignedIndex seg_12) {
        Index seg_01_id = index(seg_01);
        Index seg_12_id = index(seg_12);
        bool seg_01_ori = orientation(seg_01);
        bool seg_12_ori = orientation(seg_12);

        Index seg_01_v1 =
            seg_01_ori ? contour_segments[seg_01_id * 2 + 1] : contour_segments[seg_01_id * 2];
        Index seg_12_v1 =
            seg_12_ori ? contour_segments[seg_12_id * 2] : contour_segments[seg_12_id * 2 + 1];
        return seg_01_v1 == seg_12_v1;
    };

    auto extract_simple_cycle =
        [&](Index e01_id, Index e12_id, Index e20_id, bool e01_ori, bool e12_ori, bool e20_ori) {
            Index e01_seg_begin = contour_segment_on_edges_indices[e01_id];
            Index e01_seg_end = contour_segment_on_edges_indices[e01_id + 1];
            Index e12_seg_begin = contour_segment_on_edges_indices[e12_id];
            Index e12_seg_end = contour_segment_on_edges_indices[e12_id + 1];
            Index e20_seg_begin = contour_segment_on_edges_indices[e20_id];
            Index e20_seg_end = contour_segment_on_edges_indices[e20_id + 1];

            Index num_segments_over_edge = e01_seg_end - e01_seg_begin;
            assert(num_segments_over_edge == e12_seg_end - e12_seg_begin);
            assert(num_segments_over_edge == e20_seg_end - e20_seg_begin);

            for (Index i = 0; i < num_segments_over_edge; i++) {
                SignedIndex seg_01 = contour_segment_on_edges[e01_seg_begin + i];
                SignedIndex seg_12 = contour_segment_on_edges[e12_seg_begin + i];
                SignedIndex seg_20 = contour_segment_on_edges[e20_seg_begin + i];

                if (!e01_ori) seg_01 = -seg_01;
                if (!e12_ori) seg_12 = -seg_12;
                if (!e20_ori) seg_20 = -seg_20;

                if (is_connected(seg_01, seg_12)) {
                    assert(is_connected(seg_12, seg_20));
                    assert(is_connected(seg_20, seg_01));

                    contour_cycles.push_back(seg_01);
                    contour_cycles.push_back(seg_12);
                    contour_cycles.push_back(seg_20);
                } else {
                    assert(is_connected(seg_01, seg_20));
                    assert(is_connected(seg_12, seg_01));
                    assert(is_connected(seg_20, seg_12));

                    contour_cycles.push_back(seg_01);
                    contour_cycles.push_back(seg_20);
                    contour_cycles.push_back(seg_12);
                }

                contour_cycle_indices.push_back(static_cast<Index>(contour_cycles.size()));
            }
        };

    assert(triangles.size() % 3 == 0);
    const size_t num_triangles = triangles.size() / 3;
    for (size_t ti = 0; ti < num_triangles; ti++) {
        cycle_engine.clear();

        SignedIndex e01 = triangles[ti * 3];
        SignedIndex e12 = triangles[ti * 3 + 1];
        SignedIndex e20 = triangles[ti * 3 + 2];

        Index e01_id = index(e01);
        Index e12_id = index(e12);
        Index e20_id = index(e20);

        bool e01_ori = orientation(e01);
        bool e12_ori = orientation(e12);
        bool e20_ori = orientation(e20);

        // More precise check.
        // triangle_is_simple[ti] =
        //     edge_is_simple[e01_id] && edge_is_simple[e12_id] && edge_is_simple[e20_id] &&
        //     ((edge_shift[e01_id] * (e01_ori ? 1 : -1) + edge_shift[e12_id] * (e12_ori ? 1 : -1) +
        //       edge_shift[e20_id] * (e20_ori ? 1 : -1)) == 0);

        // More conservative check.
        triangle_is_simple[ti] = edge_is_simple[e01_id] && edge_is_simple[e12_id] &&
                                 edge_is_simple[e20_id] && edge_shift[e01_id] == 0 &&
                                 edge_shift[e12_id] == 0 && edge_shift[e20_id] == 0;

        if (triangle_is_simple[ti]) {
            extract_simple_cycle(e01_id, e12_id, e20_id, e01_ori, e12_ori, e20_ori);
        } else {
            register_edge(e01_id, e01_ori);
            register_edge(e12_id, e12_ori);
            register_edge(e20_id, e20_ori);

            cycle_engine.extract_cycles(contour_cycles, contour_cycle_indices);
        }

        contour_cycle_triangle_indices.push_back(contour_cycle_indices.size() - 1);
    }

    return {
        contour_cycles,
        contour_cycle_indices,
        contour_cycle_triangle_indices,
        triangle_is_simple};
}

std::tuple<std::vector<SignedIndex>, std::vector<Index>, std::vector<Index>>
extract_contour_polyhedra(
    size_t num_contour_segments,
    const std::vector<SignedIndex>& contour_cycles,
    const std::vector<Index>& contour_cycle_indices,
    const std::vector<Index>& contour_cycle_triangle_indices,
    const std::vector<bool>& triangle_is_simple,
    const std::vector<SignedIndex>& tetrahedra)
{
    size_t num_cycles = contour_cycle_indices.size() - 1;
    size_t num_triangles = contour_cycle_triangle_indices.size() - 1;
    size_t num_tets = tetrahedra.size() / 4;

    std::vector<SignedIndex> polyhedra;
    std::vector<Index> polyhedron_indices;
    std::vector<Index> polyhedron_tet_indices;

    polyhedra.reserve(num_cycles);
    polyhedron_indices.reserve(num_cycles + 1);
    polyhedron_tet_indices.reserve(num_tets + 1);

    polyhedron_indices.push_back(0);
    polyhedron_tet_indices.push_back(0);

    DisjointComponents component_engine(
        num_contour_segments,
        contour_cycles,
        contour_cycle_indices);

    auto register_cycles = [&](Index cycles_begin, Index cycles_end, bool ori) {
        assert(cycles_begin < cycles_end);

        for (Index ci = cycles_begin; ci < cycles_end; ci++) {
            component_engine.register_cycle(signed_index(ci, ori));
        }
    };

    auto extract_simple_components = [&](Index t021_id,
                                         Index t123_id,
                                         Index t013_id,
                                         Index t032_id,
                                         bool t021_ori,
                                         bool t123_ori,
                                         bool t013_ori,
                                         bool t032_ori) {
        Index cycles_021_begin = contour_cycle_triangle_indices[t021_id];
        Index cycles_021_end = contour_cycle_triangle_indices[t021_id + 1];
        Index cycles_123_begin = contour_cycle_triangle_indices[t123_id];
        Index cycles_123_end = contour_cycle_triangle_indices[t123_id + 1];
        Index cycles_013_begin = contour_cycle_triangle_indices[t013_id];
        Index cycles_013_end = contour_cycle_triangle_indices[t013_id + 1];
        Index cycles_032_begin = contour_cycle_triangle_indices[t032_id];
        Index cycles_032_end = contour_cycle_triangle_indices[t032_id + 1];

        Index num_cycles_per_triangle = cycles_021_end - cycles_021_begin;
        assert(num_cycles_per_triangle == cycles_123_end - cycles_123_begin);
        assert(num_cycles_per_triangle == cycles_013_end - cycles_013_begin);
        assert(num_cycles_per_triangle == cycles_032_end - cycles_032_begin);

        for (Index i = 0; i < num_cycles_per_triangle; i++) {
            Index cycle_021_id = cycles_021_begin + i;
            Index cycle_123_id = cycles_123_begin + i;
            Index cycle_013_id = cycles_013_begin + i;
            Index cycle_032_id = cycles_032_begin + i;

            SignedIndex cycle_021 = signed_index(cycle_021_id, t021_ori);
            SignedIndex cycle_123 = signed_index(cycle_123_id, t123_ori);
            SignedIndex cycle_013 = signed_index(cycle_013_id, t013_ori);
            SignedIndex cycle_032 = signed_index(cycle_032_id, t032_ori);

            polyhedra.push_back(cycle_021);
            polyhedra.push_back(cycle_123);
            polyhedra.push_back(cycle_013);
            polyhedra.push_back(cycle_032);

            polyhedron_indices.push_back(static_cast<Index>(polyhedra.size()));
        }
    };

    auto compute_components =
        [&](SignedIndex t021, SignedIndex t123, SignedIndex t013, SignedIndex t032) {
            component_engine.clear();

            Index t021_id = index(t021);
            Index t123_id = index(t123);
            Index t013_id = index(t013);
            Index t032_id = index(t032);

            assert(t021_id < num_triangles);
            assert(t123_id < num_triangles);
            assert(t013_id < num_triangles);
            assert(t032_id < num_triangles);

            bool t021_ori = orientation(t021);
            bool t123_ori = orientation(t123);
            bool t013_ori = orientation(t013);
            bool t032_ori = orientation(t032);

            bool polyhedron_is_simple = triangle_is_simple[t021_id] &&
                                        triangle_is_simple[t123_id] &&
                                        triangle_is_simple[t013_id] && triangle_is_simple[t032_id];

            if (polyhedron_is_simple) {
                extract_simple_components(
                    t021_id,
                    t123_id,
                    t013_id,
                    t032_id,
                    t021_ori,
                    t123_ori,
                    t013_ori,
                    t032_ori);
            } else {
                Index cycles_021_begin = contour_cycle_triangle_indices[t021_id];
                Index cycles_021_end = contour_cycle_triangle_indices[t021_id + 1];
                Index cycles_123_begin = contour_cycle_triangle_indices[t123_id];
                Index cycles_123_end = contour_cycle_triangle_indices[t123_id + 1];
                Index cycles_013_begin = contour_cycle_triangle_indices[t013_id];
                Index cycles_013_end = contour_cycle_triangle_indices[t013_id + 1];
                Index cycles_032_begin = contour_cycle_triangle_indices[t032_id];
                Index cycles_032_end = contour_cycle_triangle_indices[t032_id + 1];

                register_cycles(cycles_021_begin, cycles_021_end, t021_ori);
                register_cycles(cycles_123_begin, cycles_123_end, t123_ori);
                register_cycles(cycles_013_begin, cycles_013_end, t013_ori);
                register_cycles(cycles_032_begin, cycles_032_end, t032_ori);

                component_engine.extract_components(polyhedra, polyhedron_indices);
            }

            polyhedron_tet_indices.push_back(static_cast<Index>(polyhedron_indices.size() - 1));
        };

    for (size_t ti = 0; ti < num_tets; ti++) {
        SignedIndex t021 = tetrahedra[ti * 4];
        SignedIndex t123 = tetrahedra[ti * 4 + 1];
        SignedIndex t013 = tetrahedra[ti * 4 + 2];
        SignedIndex t032 = tetrahedra[ti * 4 + 3];

        compute_components(t021, t123, t013, t032);
    }

    return {std::move(polyhedra), std::move(polyhedron_indices), std::move(polyhedron_tet_indices)};
}

} // namespace mtetcol
