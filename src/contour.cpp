#include <mtetcol/contour.h>
#include <mtetcol/disjoint_cycles.h>

#include "hashmap.h"

#include <SmallVector.h>

namespace mtetcol {

namespace {

void print_polyhedron(const Contour<4>& contour, Index polyhedron_id)
{
    logger().info("Polyhedron {}", polyhedron_id);
    auto cycles = contour.get_polyhedron(polyhedron_id);
    size_t cycle_size = cycles.size();
    for (size_t i = 0; i < cycle_size; i++) {
        auto ci = cycles[i];
        Index cycle_id = index(ci);
        bool cycle_ori = orientation(ci);

        auto cycle = contour.get_cycle(cycle_id);
        std::string cycle_str;
        for (auto si : cycle) {
            cycle_str += std::to_string(value_of(si)) + " ";
        }
        logger().debug("Cycle {}: {}", i, cycle_str);
    }
}

std::vector<Index> triangulate(
    std::vector<Index>& segments,
    std::vector<SignedIndex>& cycles,
    std::vector<Index>& cycle_indices)
{
    size_t num_cycles = cycle_indices.size() - 1;

    OrientedEdgeMap diagonal_map;

    std::vector<SignedIndex> triangle_cycles;
    std::vector<Index> triangle_cycle_indices;
    std::vector<Index> cycle_to_triangle_map;

    triangle_cycles.reserve(cycles.size());
    triangle_cycle_indices.reserve(num_cycles + 1);
    triangle_cycle_indices.push_back(0);
    cycle_to_triangle_map.reserve(num_cycles + 1);
    cycle_to_triangle_map.push_back(0);

    auto get_segment = [&](SignedIndex si) -> std::array<Index, 2> {
        Index seg_id = index(si);
        bool seg_ori = orientation(si);

        assert(seg_id * 2 + 1 < segments.size());
        Index v0 = segments[seg_id * 2 + (seg_ori ? 0 : 1)];
        Index v1 = segments[seg_id * 2 + (seg_ori ? 1 : 0)];
        return {v0, v1};
    };

    auto add_diagnonal = [&](Index v0, Index v1) {
        Index num_segments = segments.size() / 2;
        auto [itr, inserted] = diagonal_map.try_emplace({v0, v1}, num_segments);
        Index diag_index = itr->second;
        if (inserted) {
            assert(diag_index == num_segments);
            segments.push_back(v0);
            segments.push_back(v1);
        }
        return diag_index;
    };

    for (size_t i = 0; i < num_cycles; i++) {
        auto start = cycle_indices[i];
        auto end = cycle_indices[i + 1];
        auto cycle = std::span<const SignedIndex>(cycles.data() + start, end - start);

        if (cycle.size() < 3) {
            // Cycle of size smaller than 3 is dropped.
        } else if (cycle.size() == 3) {
            // Check if the cycle is a triangle
            std::copy(cycle.begin(), cycle.end(), std::back_inserter(triangle_cycles));
            triangle_cycle_indices.push_back(triangle_cycles.size());
        } else {
            // Triangulate the cycle
            assert(cycle.size() > 3);
            size_t num_segments_in_cycle = end - start;

            // Pick the segment with minimum starting vertex index as the cycle start segment.
            // This ensures the same cycle from different polyhedra is triangulated consistently.
            size_t cycle_start_index = 0;
            for (size_t j = 0; j < num_segments_in_cycle; j++) {
                if (get_segment(cycle[j])[0] < get_segment(cycle[cycle_start_index])[0]) {
                    cycle_start_index = j;
                }
            }

            // Using the first vertex of cycle start segment as the vertex to generate triangle fan
            Index v0 = get_segment(cycle[cycle_start_index])[0];

            SignedIndex si_curr = invalid_signed_index;
            SignedIndex si_prev = invalid_signed_index;
            SignedIndex si_next = invalid_signed_index;

            for (size_t j = 1; j + 1 < num_segments_in_cycle; j++) {
                size_t idx = (cycle_start_index + j) % num_segments_in_cycle;
                auto si = cycle[idx];
                auto target_segment = get_segment(si);
                auto segment = get_segment(si);

                si_curr = si;

                if (j == 1) {
                    assert(v0 < segment[1]);
                    Index diag_index = add_diagnonal(v0, segment[1]);
                    si_prev = cycle[cycle_start_index];
                    si_next = signed_index(diag_index, false);
                } else if (j + 2 == num_segments_in_cycle) {
                    si_next = cycle
                        [(cycle_start_index + num_segments_in_cycle - 1) % num_segments_in_cycle];
                } else {
                    assert(v0 < segment[1]);
                    Index diag_index = add_diagnonal(v0, segment[1]);
                    si_next = signed_index(diag_index, false);
                }
                assert(si_prev != invalid_signed_index);
                assert(si_next != invalid_signed_index);
                triangle_cycles.push_back(si_prev);
                triangle_cycles.push_back(si_curr);
                triangle_cycles.push_back(si_next);
                triangle_cycle_indices.push_back(triangle_cycles.size());

                si_prev = -si_next;
            }
        }
        cycle_to_triangle_map.push_back(triangle_cycle_indices.size() - 1);
    }

    std::swap(cycles, triangle_cycles);
    std::swap(cycle_indices, triangle_cycle_indices);

    return cycle_to_triangle_map;
}

template <int dim>
std::vector<Index> compute_zero_crossing_vertices(
    const Contour<dim>& contour,
    std::span<const Scalar> function_values,
    Contour<dim>& result)
{
    static_assert(dim == 3 || dim == 4, "dim must be 3 or 4");
    size_t num_segments = contour.get_num_segments();
    std::vector<Index> zero_crossing_vertices(num_segments, invalid_index);

    auto interpolate_segment = [&](Index v0, Index v1) {
        Scalar t = function_values[v0] / (function_values[v0] - function_values[v1]);
        auto pos0 = contour.get_vertex(v0);
        auto pos1 = contour.get_vertex(v1);
        if constexpr (dim == 3) {
            result.add_vertex(
                {pos0[0] + t * (pos1[0] - pos0[0]),
                 pos0[1] + t * (pos1[1] - pos0[1]),
                 pos0[2] + t * (pos1[2] - pos0[2])});
        } else if constexpr (dim == 4) {
            result.add_vertex(
                {pos0[0] + t * (pos1[0] - pos0[0]),
                 pos0[1] + t * (pos1[1] - pos0[1]),
                 pos0[2] + t * (pos1[2] - pos0[2]),
                 pos0[3] + t * (pos1[3] - pos0[3])});
        }
        return static_cast<Index>(result.get_num_vertices() - 1);
    };

    for (size_t i = 0; i < num_segments; i++) {
        auto seg = contour.get_segment(i);
        Index v0 = seg[0];
        Index v1 = seg[1];
        if (function_values[v0] >= 0 && function_values[v1] < 0 ||
            function_values[v0] < 0 && function_values[v1] >= 0) {
            zero_crossing_vertices[i] = interpolate_segment(v0, v1);
        }
    }

    return zero_crossing_vertices;
}

template <int dim>
std::vector<Index> compute_zero_crossing_segments(
    const Contour<dim>& contour,
    std::span<const Scalar> function_values,
    const std::vector<Index>& zero_crossing_vertices,
    Contour<dim>& result)
{
    static_assert(dim == 3 || dim == 4, "dim must be 3 or 4");
    size_t num_cycles = contour.get_num_cycles();
    llvm_vecsmall::SmallVector<SignedIndex, 4> zero_crossing_segments;
    std::vector<Index> zero_crossing_segment_indices;
    zero_crossing_segment_indices.reserve(num_cycles + 1);
    zero_crossing_segment_indices.push_back(0);

    for (size_t i = 0; i < num_cycles; i++) {
        auto cycle = contour.get_cycle(i);
        size_t cycle_size = cycle.size();
        zero_crossing_segments.clear();

        for (SignedIndex sid : cycle) {
            Index seg_id = index(sid);
            if (zero_crossing_vertices[seg_id] != invalid_index) {
                zero_crossing_segments.push_back(sid);
                if (zero_crossing_segments.size() > 2) {
                    throw std::runtime_error("Cycle has more than 2 zero crossings, please call "
                                             "triangulate_cycles() first");
                }
            }
        }
        if (!zero_crossing_segments.empty()) {
            // Determine segment orientation
            assert(zero_crossing_segments.size() == 2);
            auto seg_id_0 = index(zero_crossing_segments[0]);
            bool seg_ori_0 = orientation(zero_crossing_segments[0]);
            auto seg_id_1 = index(zero_crossing_segments[1]);
            bool seg_ori_1 = orientation(zero_crossing_segments[1]);

            Index v0 = contour.get_segment(seg_id_0)[seg_ori_0 ? 0 : 1];
            if (function_values[v0] < 0) {
                result.add_segment(
                    zero_crossing_vertices[seg_id_1],
                    zero_crossing_vertices[seg_id_0]);
            } else {
                result.add_segment(
                    zero_crossing_vertices[seg_id_0],
                    zero_crossing_vertices[seg_id_1]);
            }
        }
        zero_crossing_segment_indices.push_back(result.get_num_segments());
    }
    return zero_crossing_segment_indices;
}

void map_cycle_regularity(
    const std::vector<Index>& cycle_to_triangle_map,
    std::vector<bool>& cycle_is_regular)
{
    std::vector<bool> updated_cycle_is_regular;
    updated_cycle_is_regular.reserve(cycle_to_triangle_map.back());

    size_t num_old_cycles = cycle_to_triangle_map.size() - 1;
    assert(num_old_cycles == cycle_is_regular.size());

    for (size_t i = 0; i < num_old_cycles; i++) {
        Index triangle_start = cycle_to_triangle_map[i];
        Index triangle_end = cycle_to_triangle_map[i + 1];
        updated_cycle_is_regular.insert(
            updated_cycle_is_regular.end(),
            triangle_end - triangle_start,
            cycle_is_regular[i]);
    }
    std::swap(cycle_is_regular, updated_cycle_is_regular);
}


} // namespace

template <>
void Contour<3>::triangulate_cycles()
{
    auto cycle_to_triangle_map = triangulate(m_segments, m_cycles, m_cycle_start_indices);

    // Update the cycle_is_regular vector
    map_cycle_regularity(cycle_to_triangle_map, m_cycle_is_regular);

#ifndef NDEBUG
    check_all_segments();
    check_all_cycles();
    check_all_polyhedra();
#endif
}

template <>
void Contour<4>::triangulate_cycles()
{
#ifndef NDEBUG
    check_all_segments();
    check_all_cycles();
    check_all_polyhedra();
#endif

    auto cycle_to_triangle_map = triangulate(m_segments, m_cycles, m_cycle_start_indices);

    // Update the cycle_is_regular vector
    map_cycle_regularity(cycle_to_triangle_map, m_cycle_is_regular);

#ifndef NDEBUG
    check_all_cycles();
#endif

    size_t num_polyhedra = get_num_polyhedra();
    std::vector<SignedIndex> updated_polyhedra;
    std::vector<Index> updated_polyhedron_start_indices;
    updated_polyhedra.reserve(m_polyhedra.size());
    updated_polyhedron_start_indices.reserve(m_polyhedron_start_indices.size());
    updated_polyhedron_start_indices.push_back(0);

    for (size_t i = 0; i < num_polyhedra; i++) {
        Index cycle_start = m_polyhedron_start_indices[i];
        Index cycle_end = m_polyhedron_start_indices[i + 1];
        for (size_t j = cycle_start; j < cycle_end; j++) {
            SignedIndex ci = m_polyhedra[j];
            Index cycle_id = index(ci);
            bool cycle_ori = orientation(ci);

            Index new_triangle_start = cycle_to_triangle_map[cycle_id];
            Index new_triangle_end = cycle_to_triangle_map[cycle_id + 1];

            for (Index k = new_triangle_start; k < new_triangle_end; k++) {
                updated_polyhedra.push_back(signed_index(k, cycle_ori));
            }
        }
        updated_polyhedron_start_indices.push_back(updated_polyhedra.size());
    }

    std::swap(m_polyhedra, updated_polyhedra);
    std::swap(m_polyhedron_start_indices, updated_polyhedron_start_indices);

    // Polyhedron regularity is unchnaged.

#ifndef NDEBUG
    check_all_segments();
    check_all_cycles();
    check_all_polyhedra();
#endif
}

template <>
Contour<3> Contour<3>::isocontour(std::span<Scalar> function_values) const
{
    assert(get_num_vertices() == function_values.size());
    Contour<3> result;

    std::vector<Index> zero_crossing_vertices =
        compute_zero_crossing_vertices(*this, function_values, result);

    compute_zero_crossing_segments(*this, function_values, zero_crossing_vertices, result);

#ifndef NDEBUG
    result.check_all_segments();
    result.check_all_cycles();
    result.check_all_polyhedra();
#endif

    return result;
}

template <>
Contour<4> Contour<4>::isocontour(std::span<Scalar> function_values) const
{
    assert(get_num_vertices() == function_values.size());
    Contour<4> result;

    std::vector<Index> zero_crossing_vertices =
        compute_zero_crossing_vertices(*this, function_values, result);

    std::vector<Index> zero_crossing_segment_indices =
        compute_zero_crossing_segments(*this, function_values, zero_crossing_vertices, result);

    size_t num_polyhedra = get_num_polyhedra();
    DisjointCycles disjoint_cycles(result.get_num_vertices(), result.m_segments);
    for (size_t i = 0; i < num_polyhedra; i++) {
        disjoint_cycles.clear();
        auto polyhedron = get_polyhedron(i);
        for (auto cid : polyhedron) {
            Index seg_start = zero_crossing_segment_indices[index(cid)];
            Index seg_end = zero_crossing_segment_indices[index(cid) + 1];
            bool cycle_ori = orientation(cid);

            for (Index seg_id = seg_start; seg_id < seg_end; seg_id++) {
                disjoint_cycles.register_segment(signed_index(seg_id, cycle_ori));
            }
        }
        size_t num_existing_cycles = result.m_cycle_start_indices.size();
        disjoint_cycles.extract_cycles(result.m_cycles, result.m_cycle_start_indices);

        for (size_t j = num_existing_cycles; j < result.m_cycle_start_indices.size(); j++) {
            result.m_cycle_is_regular.push_back(m_polyhedron_is_regular[i]);
        }
    }


#ifndef NDEBUG
    result.check_all_segments();
    result.check_all_cycles();
    result.check_all_polyhedra();
#endif

    return result;
}

} // namespace mtetcol
