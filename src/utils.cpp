#include "utils.h"

namespace mtetcol {

void extract_vertex_zero_crossing(
    std::span<const Scalar> time_samples,
    std::span<const Scalar> function_values,
    Scalar value,
    bool cyclic,
    std::vector<Scalar>& zero_crossing_times)
{
    assert(time_samples.front() == 0);
    assert(time_samples.back() == 1);
    size_t num_samples = time_samples.size();
    assert(num_samples >= 2);
    assert(num_samples == function_values.size());
    zero_crossing_times.reserve(zero_crossing_times.size() + num_samples);

    if (!cyclic && function_values[0] >= value) {
        zero_crossing_times.push_back(time_samples[0]);
    }

    for (size_t i = 0; i < num_samples; i++) {
        size_t next_i = (i + 1) % num_samples;
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

    if (!cyclic && function_values[num_samples - 1] < value) {
        zero_crossing_times.push_back(time_samples[num_samples - 1]);
    }
}

std::tuple<std::vector<Scalar>, std::vector<size_t>, std::vector<bool>> extract_contour_vertices(
    const std::vector<Scalar>& time_samples,
    const std::vector<Scalar>& function_values,
    const std::vector<Index>& vertex_start_indices,
    Scalar value,
    bool cyclic)
{
    size_t num_vertices = vertex_start_indices.size() - 1;

    std::vector<Scalar> zero_crossing_times;
    std::vector<size_t> zero_crossing_indices;
    std::vector<bool> initial_signs;
    zero_crossing_times.reserve(time_samples.size());
    zero_crossing_indices.reserve(num_vertices + 1);
    zero_crossing_indices.push_back(0);
    initial_signs.reserve(num_vertices);
    for (size_t vi = 0; vi < num_vertices; vi++) {
        const Index idx_begein = vertex_start_indices[vi];
        const Index idx_end = vertex_start_indices[vi + 1];
        std::span<const Scalar> time_samples_i =
            std::span<const Scalar>(time_samples.data() + idx_begein, idx_end - idx_begein);
        std::span<const Scalar> function_values_i =
            std::span<const Scalar>(function_values.data() + idx_begein, idx_end - idx_begein);
        initial_signs.push_back(function_values_i[0] >= value ? true : false);

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

std::tuple<std::vector<Index>, std::vector<Index>> extract_contour_segments(
    const std::vector<Scalar>& contour_times,
    const std::vector<size_t>& contour_time_indices,
    const std::vector<bool>& initial_signs,
    const std::vector<Index>& edges,
    bool cyclic)
{
    size_t num_edges = edges.size() / 2;

    std::vector<Index> segments;
    std::vector<Index> segment_indices;
    segment_indices.reserve(num_edges + 1);
    segment_indices.push_back(0);

    struct LocalIndex
    {
        Index vertex_index;
        Index time_index;
    };

    auto merge_time_samples = [](Index v0,
                                 Index v1,
                                 std::span<const Scalar> times_0,
                                 std::span<const Scalar> times_1,
                                 std::vector<LocalIndex>& local_indices) {
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


    std::vector<LocalIndex> local_indices; // TODO: this dynamic allocation can be avoided.
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

        size_t offset = 0;
        size_t parity = 0;

        if (cyclic) {
            if (initial_signs[v0] != initial_signs[v1]) {
                // Different sign at the beginning
                offset = 1;

                if (local_indices[0].vertex_index == v0) {
                    // First merged time sample is on v0.
                    // The offset will shift the samples, and the first time sample available to
                    // match on v0 will be the second one, which have the opposite initial sign.
                    if (!initial_signs[v0]) {
                        parity = 1;
                    }
                } else {
                    // First merged time sample is on v1.
                    // Shifting does not change the initial sign of v0's first sample.
                    if (initial_signs[v0]) {
                        parity = 1;
                    }
                }
            } else {
                if (initial_signs[v0]) {
                    parity = 1;
                }
            }
        }

        size_t num_segments = local_indices.size() / 2;
        for (size_t si = 0; si < num_segments; si++) {
            const auto& lid0 = local_indices[(si * 2 + offset) % local_indices.size()];
            const auto& lid1 = local_indices[(si * 2 + 1 + offset) % local_indices.size()];

            // Pair local indices up based on even/odd parity.
            if ((lid0.vertex_index == v0 && lid0.time_index % 2 == parity) ||
                (lid0.vertex_index == v1 && lid0.time_index % 2 != parity)) {
                // lid0 -> lid1
                segments.push_back(contour_time_indices[lid0.vertex_index] + lid0.time_index);
                segments.push_back(contour_time_indices[lid1.vertex_index] + lid1.time_index);
            } else {
                // lid1 -> lid0
                segments.push_back(contour_time_indices[lid1.vertex_index] + lid1.time_index);
                segments.push_back(contour_time_indices[lid0.vertex_index] + lid0.time_index);
            }
        }
        segment_indices.push_back(static_cast<Index>(segments.size()));
    }

    return {std::move(segments), std::move(segment_indices)};
}

std::tuple<std::vector<SignedIndex>, std::vector<Index>, std::vector<Index>> extract_contour_cycles(
    const size_t num_contour_vertices,
    const std::vector<Index>& contour_segments,
    const std::vector<Index>& contour_segment_indices,
    const std::vector<Index>& edges,
    const std::vector<SignedIndex>& triangles)
{
    // Map from contour vertex index to contour segement index.
    std::vector<Index> next_index(num_contour_vertices, invalid_index);
    std::vector<SignedIndex> contour_cycles;
    std::vector<Index> contour_cycle_indices;
    std::vector<Index> contour_cycle_triangle_indices;

    contour_cycle_indices.push_back(0);
    contour_cycle_triangle_indices.push_back(0);

    assert(contour_segments.size() % 2 == 0);
    const size_t num_segments = contour_segments.size() / 2;

    auto clear_next_index = [&](Index eid) {
        assert(eid >= 0 && eid + 1 < contour_segment_indices.size());
        Index start = contour_segment_indices[eid];
        Index end = contour_segment_indices[eid + 1];
        assert((end - start) % 2 == 0);
        for (Index i = start; i < end; i++) {
            next_index[contour_segments[i]] = invalid_index;
        }
    };

    auto print_edge = [&](Index eid, bool ori) {
        Index v0 = edges[eid * 2];
        Index v1 = edges[eid * 2 + 1];
        if (!ori) {
            std::swap(v0, v1);
        }
        logger().debug("edge: ({}, {})", v0, v1);
    };

    auto print_segment = [&](Index eid, bool ori) {
        Index start = contour_segment_indices[eid];
        Index end = contour_segment_indices[eid + 1];

        Index num_segments_over_edge = (end - start) / 2;
        for (size_t si = 0; si < num_segments_over_edge; si++) {
            Index v0 = contour_segments[start + si * 2];
            Index v1 = contour_segments[start + si * 2 + 1];
            if (!ori) std::swap(v0, v1);

            logger().debug("segment: ({}, {})", v0, v1);
        }
    };

    auto register_next_index = [&](Index eid, bool ori) {
        logger().debug("eid: {}, ori: {}", eid, ori);
        assert(eid >= 0 && eid + 1 < contour_segment_indices.size());
        Index start = contour_segment_indices[eid];
        Index end = contour_segment_indices[eid + 1];
        assert((end - start) % 2 == 0);
        Index num_segments_over_edge = (end - start) / 2;
        for (size_t si = 0; si < num_segments_over_edge; si++) {
            Index v0 = contour_segments[start + si * 2];
            Index v1 = contour_segments[start + si * 2 + 1];
            if (!ori) std::swap(v0, v1);
            assert(v0 < num_contour_vertices);
            assert(next_index[v0] == invalid_index);
            next_index[v0] = start / 2 + si;
            logger().debug("registering next index: ({}, {}) -> {}", v0, v1, next_index[v0]);
        }
    };

    auto grow_cycle = [&](Index vid) {
        assert(next_index[vid] != invalid_index);
        Index sid = next_index[vid];
        next_index[vid] = invalid_index;

        Index v0 = contour_segments[sid * 2];
        Index v1 = contour_segments[sid * 2 + 1];

        if (vid == v0) {
            contour_cycles.push_back(signed_index(sid, true));
            return v1;
        } else {
            contour_cycles.push_back(signed_index(sid, false));
            return v0;
        }
    };

    auto chain_cycles = [&](Index eid, bool ori) {
        assert(eid >= 0 && eid + 1 < contour_segment_indices.size());
        Index start = contour_segment_indices[eid];
        Index end = contour_segment_indices[eid + 1];
        assert((end - start) % 2 == 0);
        Index num_segments_over_edge = (end - start) / 2;

        for (size_t si = 0; si < num_segments_over_edge; si++) {
            Index v0 = contour_segments[start + si * 2];
            Index v1 = contour_segments[start + si * 2 + 1];
            if (!ori) std::swap(v0, v1);
            if (next_index[v0] == invalid_index) continue;

            while (next_index[v0] != invalid_index) {
                v0 = grow_cycle(v0);
            }
            contour_cycle_indices.push_back(static_cast<Index>(contour_cycles.size()));
        }
    };

    assert(triangles.size() % 3 == 0);
    const size_t num_triangles = triangles.size() / 3;
    for (size_t ti = 0; ti < num_triangles; ti++) {
        logger().debug("Triangle {}", ti);
        SignedIndex e01 = triangles[ti * 3];
        SignedIndex e12 = triangles[ti * 3 + 1];
        SignedIndex e20 = triangles[ti * 3 + 2];

        Index e01_id = index(e01);
        Index e12_id = index(e12);
        Index e20_id = index(e20);

        bool e01_ori = orientation(e01);
        bool e12_ori = orientation(e12);
        bool e20_ori = orientation(e20);

        clear_next_index(e01_id);
        clear_next_index(e12_id);
        clear_next_index(e20_id);

        print_edge(e01_id, e01_ori);
        print_segment(e01_id, e01_ori);
        print_edge(e12_id, e12_ori);
        print_segment(e12_id, e12_ori);
        print_edge(e20_id, e20_ori);
        print_segment(e20_id, e20_ori);

        register_next_index(e01_id, e01_ori);
        register_next_index(e12_id, e12_ori);
        register_next_index(e20_id, e20_ori);

        chain_cycles(e01_id, e01_ori);
        chain_cycles(e12_id, e12_ori);
        chain_cycles(e20_id, e20_ori);

        contour_cycle_triangle_indices.push_back(contour_cycle_indices.size());
    }

    return {contour_cycles, contour_cycle_indices, contour_cycle_triangle_indices};
}

std::tuple<std::vector<SignedIndex>, std::vector<Index>, std::vector<Index>>
extract_contour_polyhedra(
    size_t num_contour_segments,
    const std::vector<SignedIndex>& contour_cycles,
    const std::vector<Index>& contour_cycle_indices,
    const std::vector<Index>& contour_cycle_triangle_indices,
    const std::vector<SignedIndex>& tetrahedra)
{
    size_t num_cycles = contour_cycle_indices.size() - 1;
    size_t num_triangles = contour_cycle_triangle_indices.size() - 1;
    size_t num_tets = tetrahedra.size() / 4;

    std::vector<SignedIndex> polyhedra;
    std::vector<Index> polyhedron_indices;
    std::vector<Index> polyhedron_tet_indices;
    polyhedra.reserve(num_cycles * 2);
    polyhedron_indices.reserve(num_cycles + 1);
    polyhedron_tet_indices.reserve(num_tets + 1);
    polyhedra_indices.push_back(0);
    polyhedron_tet_indices.push_back(0);

    // Map from segment index to the cycle containing the it and is conistently oriented.
    std::vector<SignedIndex> positive_segment_map(num_contour_segments, invalid_signed_index);
    // Map from segment index to the cycle containing the it and has opposite orientation.
    std::vector<SignedIndex> negative_segment_map(num_contour_segments, invalid_signed_index);
    // Map from cycle index to its relative orientation to the tet column.
    std::vector<bool> cycle_orientations(num_cycles, false);
    // Tracking whether a cycle has been processed.
    std::vector<bool> involved(num_cycles, false);

    auto clear_register_cycle = [&](Index cycle_id) {
        for (Index i = cycle_begin; i < cycle_end; i++) {
            SignedIndex si = contour_cycles[i];
            Index seg_id = index(si);
            positive_segment_map[seg_id] = invalid_signed_index;
            negative_segment_map[seg_id] = invalid_signed_index;
        }
    };

    auto register_cycle = [&](Index cycle_id, bool orientation) {
        Index cycle_begin = contour_cycle_indices[cycle_id];
        Index cycle_end = contour_cycle_indices[cycle_id + 1];
        for (Index i = cycle_begin; i < cycle_end; i++) {
            SignedIndex si = contour_cycles[i];
            Index seg_id = index(si);
            assert(seg_id < num_contour_segments);
            bool seg_ori = orientation(si);

            if (seg_ori != orientation) {
                assert(positive_segment_map[seg_id] == invalid_signed_index);
                negative_segment_map[seg_id] = signed_index(cycle_id, orientation);
            } else {
                assert(negative_segment_map[seg_id] == invalid_signed_index);
                positive_segment_map[seg_id] = signed_index(cycle_id, orientation);
            }
        }
    };

    auto register_cycles = [&](Index cycle_index_begin, Index cycle_index_end, bool orientation) {
        for (Index ci = cycle_index_begin; ci < cycle_index_end; ci++) {
            assert(ci < num_cycles);
            cycle_orientations[ci] = orientation;
            register_cycle(ci, orientation);
        }
    };

    auto compute_components = [&](Index cycles_021_begin,
                                  Index cycles_021_end,
                                  Index cycles_123_begin,
                                  Index cycles_123_end,
                                  Index cycles_013_begin,
                                  Index cycles_013_end,
                                  Index cycles_032_begin,
                                  Index cycles_032_end,
                                  bool t021_ori,
                                  bool t123_ori,
                                  bool t013_ori,
                                  bool t032_ori) {
        clear_register_cycle(cycles_021_begin, cycles_021_end);
        clear_register_cycle(cycles_123_begin, cycles_123_end);
        clear_register_cycle(cycles_013_begin, cycles_013_end);
        clear_register_cycle(cycles_032_begin, cycles_032_end);

        register_cycles(cycles_021_begin, cycles_021_end, t021_ori);
        register_cycles(cycles_123_begin, cycles_123_end, t123_ori);
        register_cycles(cycles_013_begin, cycles_013_end, t013_ori);
        register_cycles(cycles_032_begin, cycles_032_end, t032_ori);

        extract_components(cycles_021_begin, cycles_021_end, t021_ori);
    };

    for (size_t ti = 0; ti < num_tets; ti++) {
        SignedIndex t021 = tetrahedra[ti * 4];
        SignedIndex t123 = tetrahedra[ti * 4 + 1];
        SignedIndex t013 = tetrahedra[ti * 4 + 2];
        SignedIndex t032 = tetrahedra[ti * 4 + 3];

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
    }

    return {
        std::move(polyhedra),
        std::move(polyhedron_indices),
        std::move(polyhedron_tet_indices)};
}

} // namespace mtetcol
