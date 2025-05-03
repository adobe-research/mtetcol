#include "utils.h"

#include <mtetcol/logger.h>
#include <mtetcol/disjoint_cycles.h>

#include <SmallVector.h>

#include <functional>
#include <stdexcept>
#include <array>

namespace mtetcol {

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
            idx_end - idx_begin
        );
        std::span<const Scalar> function_values_i(
            function_values.data() + idx_begin, 
            idx_end - idx_begin
        );

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
    std::vector<SignedIndex> contour_cycles;
    std::vector<Index> contour_cycle_indices;
    std::vector<Index> contour_cycle_triangle_indices;

    contour_cycle_indices.push_back(0);
    contour_cycle_triangle_indices.push_back(0);

    DisjointCycles cycle_engine(num_contour_vertices, contour_segments);

    auto register_edge = [&](Index e_id, bool ori) {
        Index seg_begin = contour_segment_indices[e_id];
        Index seg_end = contour_segment_indices[e_id + 1];
        assert(seg_begin < seg_end);
        assert((seg_end - seg_begin) % 2 == 0);

        for (Index i = seg_begin; i < seg_end; i+=2) {
            Index seg_id = i / 2;
            SignedIndex si = signed_index(seg_id, ori);
            cycle_engine.register_segment(si);
        }
    };

    assert(triangles.size() % 3 == 0);
    const size_t num_triangles = triangles.size() / 3;
    for (size_t ti = 0; ti < num_triangles; ti++) {
        logger().debug("Triangle {}", ti);
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

        register_edge(e01_id, e01_ori);
        register_edge(e12_id, e12_ori);
        register_edge(e20_id, e20_ori);

        cycle_engine.extract_cycles(contour_cycles, contour_cycle_indices);

        contour_cycle_triangle_indices.push_back(contour_cycle_indices.size() - 1);
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
    polyhedron_indices.push_back(0);
    polyhedron_tet_indices.push_back(0);

    // Map from segment index to the cycle containing the it and is conistently oriented.
    std::vector<SignedIndex> positive_segment_map(num_contour_segments, invalid_signed_index);
    // Map from segment index to the cycle containing the it and has opposite orientation.
    std::vector<SignedIndex> negative_segment_map(num_contour_segments, invalid_signed_index);
    // Tracking whether a cycle has been processed.
    std::vector<bool> involved(num_cycles, false);

    auto clear_register_cycle = [&](Index cycle_id) {
        assert(cycle_id < num_cycles);
        involved[cycle_id] = false;

        Index cycle_begin = contour_cycle_indices[cycle_id];
        Index cycle_end = contour_cycle_indices[cycle_id + 1];

        for (Index i = cycle_begin; i < cycle_end; i++) {
            SignedIndex si = contour_cycles[i];
            Index seg_id = index(si);
            positive_segment_map[seg_id] = invalid_signed_index;
            negative_segment_map[seg_id] = invalid_signed_index;
        }
    };

    auto clear_register_cycles = [&](Index cycle_index_begin, Index cycle_index_end) {
        for (Index ci = cycle_index_begin; ci < cycle_index_end; ci++) {
            clear_register_cycle(ci);
        }
    };

    auto register_cycle = [&](Index cycle_id, bool ori) {
        Index cycle_begin = contour_cycle_indices[cycle_id];
        Index cycle_end = contour_cycle_indices[cycle_id + 1];
        for (Index i = cycle_begin; i < cycle_end; i++) {
            SignedIndex si = contour_cycles[i];
            Index seg_id = index(si);
            assert(seg_id < num_contour_segments);
            bool seg_ori = orientation(si);

            if (seg_ori != ori) {
                assert(negative_segment_map[seg_id] == invalid_signed_index);
                negative_segment_map[seg_id] = signed_index(cycle_id, ori);
            } else {
                assert(positive_segment_map[seg_id] == invalid_signed_index);
                positive_segment_map[seg_id] = signed_index(cycle_id, ori);
            }
        }
    };

    auto register_cycles = [&](Index cycle_index_begin, Index cycle_index_end, bool ori) {
        for (Index ci = cycle_index_begin; ci < cycle_index_end; ci++) {
            assert(ci < num_cycles);
            register_cycle(ci, ori);
        }
    };

    std::function<void(Index, bool)> grow_polyhedron;
    grow_polyhedron = [&](Index cycle_id, bool ori) {
        assert(!involved[cycle_id]);
        involved[cycle_id] = true;
        polyhedra.push_back(signed_index(cycle_id, ori));

        Index cycle_begin = contour_cycle_indices[cycle_id];
        Index cycle_end = contour_cycle_indices[cycle_id + 1];
        for (Index ci = cycle_begin; ci < cycle_end; ci++) {
            SignedIndex si = contour_cycles[ci];
            Index seg_id = index(si);
            Index seg_ori = orientation(si);
            if (seg_ori == ori) {
                assert(positive_segment_map[seg_id] != invalid_signed_index);
                assert(index(positive_segment_map[seg_id]) == cycle_id);
                assert(positive_segment_map[seg_id] == signed_index(cycle_id, ori));
                SignedIndex adj_cycle = negative_segment_map[seg_id];
                assert(adj_cycle != invalid_signed_index);
                Index adj_cycle_id = index(adj_cycle);
                bool adj_cycle_ori = orientation(adj_cycle);

                if (!involved[adj_cycle_id]) {
                    grow_polyhedron(adj_cycle_id, adj_cycle_ori);
                }
            } else {
                assert(negative_segment_map[seg_id] != invalid_signed_index);
                assert(index(negative_segment_map[seg_id]) == cycle_id);
                assert(negative_segment_map[seg_id] == signed_index(cycle_id, ori));
                SignedIndex adj_cycle = positive_segment_map[seg_id];
                assert(adj_cycle != invalid_signed_index);
                Index adj_cycle_id = index(adj_cycle);
                bool adj_cycle_ori = orientation(adj_cycle);

                if (!involved[adj_cycle_id]) {
                    grow_polyhedron(adj_cycle_id, adj_cycle_ori);
                }
            }
        }
    };

    auto extract_components = [&](Index cycle_index_begin, Index cycle_index_end, bool ori) {
        for (Index ci = cycle_index_begin; ci < cycle_index_end; ci++) {
            if (!involved[ci]) {
                grow_polyhedron(ci, ori);
                polyhedron_indices.push_back(static_cast<Index>(polyhedra.size()));
            }
        }
    };

    auto compute_components =
        [&](SignedIndex t021, SignedIndex t123, SignedIndex t013, SignedIndex t032) {
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

            Index cycles_021_begin = contour_cycle_triangle_indices[t021_id];
            Index cycles_021_end = contour_cycle_triangle_indices[t021_id + 1];
            Index cycles_123_begin = contour_cycle_triangle_indices[t123_id];
            Index cycles_123_end = contour_cycle_triangle_indices[t123_id + 1];
            Index cycles_013_begin = contour_cycle_triangle_indices[t013_id];
            Index cycles_013_end = contour_cycle_triangle_indices[t013_id + 1];
            Index cycles_032_begin = contour_cycle_triangle_indices[t032_id];
            Index cycles_032_end = contour_cycle_triangle_indices[t032_id + 1];

            clear_register_cycles(cycles_021_begin, cycles_021_end);
            clear_register_cycles(cycles_123_begin, cycles_123_end);
            clear_register_cycles(cycles_013_begin, cycles_013_end);
            clear_register_cycles(cycles_032_begin, cycles_032_end);

            register_cycles(cycles_021_begin, cycles_021_end, t021_ori);
            register_cycles(cycles_123_begin, cycles_123_end, t123_ori);
            register_cycles(cycles_013_begin, cycles_013_end, t013_ori);
            register_cycles(cycles_032_begin, cycles_032_end, t032_ori);

            extract_components(cycles_021_begin, cycles_021_end, t021_ori);
            extract_components(cycles_123_begin, cycles_123_end, t123_ori);
            extract_components(cycles_013_begin, cycles_013_end, t013_ori);
            extract_components(cycles_032_begin, cycles_032_end, t032_ori);

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
