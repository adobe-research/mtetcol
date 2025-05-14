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

    EdgeMap diagonal_map;

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

    /**
     * Add a new segment (v0, v1) to segments.
     * If ensure_unique is true, it will check if the segment has been added before.
     */
    auto add_diagonal = [&](Index v0, Index v1, bool ensure_unique) -> SignedIndex {
        Index num_segments = segments.size() / 2;
        Index diag_index = invalid_index;
        if (ensure_unique) {
            auto [itr, inserted] = diagonal_map.try_emplace({v0, v1}, num_segments);
            diag_index = itr->second;
            if (inserted) {
                assert(diag_index == num_segments);
                segments.push_back(v0);
                segments.push_back(v1);
            }
            bool diag_ori = segments[diag_index * 2] == v0;
            assert(!inserted || diag_ori);
            return signed_index(diag_index, diag_ori);
        } else {
            diag_index = num_segments;
            segments.push_back(v0);
            segments.push_back(v1);
            return signed_index(diag_index, true);
        }
    };

    for (size_t i = 0; i < num_cycles; i++) {
        auto start = cycle_indices[i];
        auto end = cycle_indices[i + 1];
        auto cycle = std::span<const SignedIndex>(cycles.data() + start, end - start);
        size_t num_segments_in_cycle = end - start;

        if (num_segments_in_cycle < 3) {
            // Cycle of size smaller than 3 is dropped.
        } else if (num_segments_in_cycle == 3) {
            // Check if the cycle is a triangle
            std::copy(cycle.begin(), cycle.end(), std::back_inserter(triangle_cycles));
            triangle_cycle_indices.push_back(triangle_cycles.size());
        } else if (num_segments_in_cycle == 4) {
            // Cycle is a quad, uniquely deterine a diagonal to insert
            std::array<Index, 4> cycle_vertices;
            for (size_t j = 0; j < 4; j++) {
                auto si = cycle[j];
                auto segment = get_segment(si);
                cycle_vertices[j] = segment[0];
            }

            if (std::min(cycle_vertices[0], cycle_vertices[2]) <
                std::min(cycle_vertices[1], cycle_vertices[3])) {
                // Insert diagonal 0-2
                SignedIndex diag_index = add_diagonal(cycle_vertices[0], cycle_vertices[2], true);

                // Add triangle (0, 1, 2)
                triangle_cycles.push_back(cycle[0]);
                triangle_cycles.push_back(cycle[1]);
                triangle_cycles.push_back(-diag_index);
                triangle_cycle_indices.push_back(triangle_cycles.size());

                // Add triangle (2, 3, 0)
                triangle_cycles.push_back(cycle[2]);
                triangle_cycles.push_back(cycle[3]);
                triangle_cycles.push_back(diag_index);
                triangle_cycle_indices.push_back(triangle_cycles.size());
            } else {
                // Insert diagonal 1-3
                SignedIndex diag_index = add_diagonal(cycle_vertices[1], cycle_vertices[3], true);

                // Add triangle (0, 1, 3)
                triangle_cycles.push_back(cycle[0]);
                triangle_cycles.push_back(diag_index);
                triangle_cycles.push_back(cycle[3]);
                triangle_cycle_indices.push_back(triangle_cycles.size());

                // Add triangle (1, 2, 3)
                triangle_cycles.push_back(cycle[1]);
                triangle_cycles.push_back(cycle[2]);
                triangle_cycles.push_back(-diag_index);
                triangle_cycle_indices.push_back(triangle_cycles.size());
            }
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
                    SignedIndex diag_index = add_diagonal(v0, segment[1], false);
                    si_prev = cycle[cycle_start_index];
                    si_next = -diag_index;
                } else if (j + 2 == num_segments_in_cycle) {
                    si_next = cycle
                        [(cycle_start_index + num_segments_in_cycle - 1) % num_segments_in_cycle];
                } else {
                    assert(v0 < segment[1]);
                    SignedIndex diag_index = add_diagonal(v0, segment[1], false);
                    si_next = -diag_index;
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

// compute a cubic root in [0,1] given values and gradients at two endpoints
Scalar get_cubic_root(Scalar val1, Scalar val2, Scalar g1, Scalar g2)
{
    if (val1 == 0) {
        return 0;
    }
    if (val2 == 0) {
        return 1;
    }
    assert(val1 * val2 < 0);

    // make sure val1 < 0 and val2 > 0
    if (val1 > 0) {
        val1 = -val1;
        val2 = -val2;
        g1 = -g1;
        g2 = -g2;
    }

    // compute the cubic function f(x) = a*x^3 + b*x^2 + c*x + d
    const Scalar a = g1 + g2 + 2 * (val1 - val2);
    const Scalar b = 3 * (val2 - val1) - 2 * g1 - g2;
    const Scalar c = g1;
    const Scalar d = val1;

    // initial guess: the linear root
    Scalar x = val1 / (val1 - val2);

    // root finding: combine Halley's method and bisect method
    // mostly a bisect method, but first find the next guess using Halley's method,
    // if the guess doesn't lie in the sign-changing interval, use the midpoint of that interval
    // terminate when the change in x is small
    Scalar xlo = 0;
    Scalar xhi = 1;
    constexpr Scalar x_tol = 1e-4;
    constexpr int max_iterations = 100; // Maximum number of iterations
    int iteration = 0; // Iteration counter
    while (true) {
        if (iteration++ >= max_iterations) {
            logger().warn("Maximum iterations reached in get_cubic_root");
            break;
        }
        Scalar f = d + x * (c + x * (b + x * a));
        if (f == 0) {
            break;
        }
        if (f < 0) {
            xlo = x;
        } else {
            xhi = x;
        }
        // f'(x) = 3*a*x^2 + 2*b*x + c
        // f''(x) = 6*a*x + 2*b
        const Scalar df = c + x * (2 * b + x * 3 * a);
        const Scalar ddf = 2 * b + x * 6 * a;
        const Scalar denominator = 2 * df * df - f * ddf;
        Scalar x_new;
        if (std::abs(denominator) < 1e-8) { // Safeguard against division by zero
            x_new = 0.5 * (xlo + xhi); // Use the midpoint of the interval
        } else {
            const Scalar dx = 2 * f * df / denominator;
            x_new = x - dx;
        }
        if (x_new <= xlo || x_new >= xhi) {
            x_new = 0.5 * (xlo + xhi);
        }
        if (std::abs(x_new - x) < x_tol) {
            x = x_new;
            break;
        }
        x = x_new;
    }

    return x;
}

template <int dim>
std::vector<Index> compute_zero_crossing_vertices(
    const Contour<dim>& contour,
    std::span<const Scalar> function_values,
    std::span<const Scalar> function_gradients,
    Contour<dim>& result)
{
    static_assert(dim == 3 || dim == 4, "dim must be 3 or 4");
    size_t num_segments = contour.get_num_segments();
    std::vector<Index> zero_crossing_vertices(num_segments, invalid_index);

    auto interpolate_segment = [&](Index v0, Index v1) {
        Scalar val0 = function_values[v0];
        Scalar val1 = function_values[v1];
        auto pos0 = contour.get_vertex(v0);
        auto pos1 = contour.get_vertex(v1);
        auto grad0 = function_gradients.subspan(v0 * dim, dim);
        auto grad1 = function_gradients.subspan(v1 * dim, dim);
        // directional derivatives
        Scalar g0, g1;
        // since the interval is scaled to (0,1), we don't normalize by dividing by the length of
        // p0p1
        if constexpr (dim == 3) {
            g0 = grad0[0] * (pos1[0] - pos0[0]) + grad0[1] * (pos1[1] - pos0[1]) +
                 grad0[2] * (pos1[2] - pos0[2]);
            g1 = grad1[0] * (pos1[0] - pos0[0]) + grad1[1] * (pos1[1] - pos0[1]) +
                 grad1[2] * (pos1[2] - pos0[2]);
        } else if constexpr (dim == 4) {
            g0 = grad0[0] * (pos1[0] - pos0[0]) + grad0[1] * (pos1[1] - pos0[1]) +
                 grad0[2] * (pos1[2] - pos0[2]) + grad0[3] * (pos1[3] - pos0[3]);
            g1 = grad1[0] * (pos1[0] - pos0[0]) + grad1[1] * (pos1[1] - pos0[1]) +
                 grad1[2] * (pos1[2] - pos0[2]) + grad1[3] * (pos1[3] - pos0[3]);
        }

        Scalar t = get_cubic_root(val0, val1, g0, g1);
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
Contour<3> Contour<3>::isocontour(
    std::span<Scalar> function_values,
    std::span<Scalar> function_gradients) const
{
    assert(get_num_vertices() == function_values.size());
    bool has_gradients = function_gradients.size() == get_num_vertices() * 3;
    Contour<3> result;

    std::vector<Index> zero_crossing_vertices;
    if (has_gradients) {
        zero_crossing_vertices =
            compute_zero_crossing_vertices(*this, function_values, function_gradients, result);
    } else {
        zero_crossing_vertices = compute_zero_crossing_vertices(*this, function_values, result);
    }
    compute_zero_crossing_segments(*this, function_values, zero_crossing_vertices, result);

#ifndef NDEBUG
    result.check_all_segments();
    result.check_all_cycles();
    result.check_all_polyhedra();
#endif

    return result;
}

template <>
Contour<4> Contour<4>::isocontour(
    std::span<Scalar> function_values,
    std::span<Scalar> function_gradients) const
{
    assert(get_num_vertices() == function_values.size());
    bool has_gradients = function_gradients.size() == get_num_vertices() * 4;
    Contour<4> result;

    std::vector<Index> zero_crossing_vertices;
    if (has_gradients) {
        zero_crossing_vertices =
            compute_zero_crossing_vertices(*this, function_values, function_gradients, result);
    } else {
        zero_crossing_vertices = compute_zero_crossing_vertices(*this, function_values, result);
    }

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
            assert(seg_end - seg_start < 2);
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
