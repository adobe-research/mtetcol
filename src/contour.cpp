#include <mtetcol/contour.h>
#include <mtetcol/disjoint_cycles.h>
#include <mtetcol/nondisjoint_cycles.h>

#include "hashmap.h"

#include <SmallVector.h>
#include <ankerl/unordered_dense.h>

namespace mtetcol {

template <int dim>
bool Contour<dim>::cycle_is_simple(Index cid) const
{
    auto cycle = this->get_cycle(cid);
    size_t cycle_size = cycle.size();

    ankerl::unordered_dense::map<Index, Index> vertex_count;

    for (size_t i = 0; i < cycle_size; i++) {
        auto seg = get_segment(index(cycle[i]));

        vertex_count[seg[0]]++;
        vertex_count[seg[1]]++;
    }

    for (auto [v, count] : vertex_count) {
        if (count != 2) {
            return false;
        }
    }
    return true;
}

llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<SignedIndex, 16>, 4> extract_simple_loops(
    std::span<const SignedIndex> cycle,
    std::function<std::array<Index, 2>(Index)> get_segment)
{
    using SmallVector = llvm_vecsmall::SmallVector<SignedIndex, 16>;
    using SubcycleList = llvm_vecsmall::SmallVector<SmallVector, 4>;

    if (cycle.empty()) {
        return {};
    }

    // Step 1: Identify and remove canceling segment pairs
    // Segments that connect the same two vertices in opposite directions cancel each other out
    // Build a map of directed edges: (v0, v1) -> list of segment indices
    ankerl::unordered_dense::map<std::pair<Index, Index>, SmallVector> edge_map;
    for (auto sid : cycle) {
        Index seg_idx = index(sid);
        bool seg_ori = orientation(sid);
        auto seg = get_segment(seg_idx);
        Index v0 = seg[seg_ori ? 0 : 1]; // start vertex
        Index v1 = seg[seg_ori ? 1 : 0]; // end vertex
        edge_map[{v0, v1}].push_back(sid);
    }

    // Find canceling pairs and mark segments to skip
    // Only segments connecting the same vertices in OPPOSITE directions cancel each other
    ankerl::unordered_dense::set<Index> canceled_segments;
    ankerl::unordered_dense::set<std::pair<Index, Index>> processed_edges;

    for (const auto& [edge, sids] : edge_map) {
        auto [v0, v1] = edge;

        // Skip if we already processed this edge pair
        if (processed_edges.count({v0, v1}) > 0 || processed_edges.count({v1, v0}) > 0) {
            continue;
        }
        processed_edges.insert({v0, v1});

        size_t forward_count = sids.size();
        size_t reverse_count = 0;
        SmallVector reverse_sids;

        // Check if reverse edge exists
        auto reverse_edge = std::make_pair(v1, v0);
        if (edge_map.count(reverse_edge) > 0) {
            reverse_sids = edge_map[reverse_edge];
            reverse_count = reverse_sids.size();
        }

        // Match and cancel pairs going in opposite directions
        size_t cancel_count = std::min(forward_count, reverse_count);
        for (size_t i = 0; i < cancel_count; i++) {
            canceled_segments.insert(index(sids[i]));
            canceled_segments.insert(index(reverse_sids[i]));
        }
    }

    // Step 2: Build segment array and use NonDisjointCycles for extraction
    // Collect non-canceled segments and their vertex pairs
    llvm_vecsmall::SmallVector<Index, 16> segment_vertices;
    SmallVector non_canceled_sids;
    ankerl::unordered_dense::map<Index, Index> old_to_new_index; // Map original index to new index

    Index new_idx = 0;
    for (auto sid : cycle) {
        Index seg_idx = index(sid);

        // Skip canceled segments
        if (canceled_segments.count(seg_idx) > 0) {
            continue;
        }

        auto seg = get_segment(seg_idx);
        segment_vertices.push_back(seg[0]);
        segment_vertices.push_back(seg[1]);

        old_to_new_index[seg_idx] = new_idx;
        non_canceled_sids.push_back(sid);
        new_idx++;
    }

    // If no non-canceled segments remain, return empty
    if (non_canceled_sids.empty()) {
        return {};
    }

    // Step 3: Use NonDisjointCycles to extract simple loops
    NonDisjointCycles cycle_extractor(segment_vertices);
    for (auto sid : non_canceled_sids) {
        Index old_idx = index(sid);
        Index mapped_idx = old_to_new_index[old_idx];
        bool ori = orientation(sid);
        cycle_extractor.register_segment(signed_index(mapped_idx, ori));
    }

    std::vector<SignedIndex> extracted_cycles;
    std::vector<Index> cycle_indices = {0}; // Initialize with 0
    cycle_extractor.extract_cycles(extracted_cycles, cycle_indices);

    // Step 4: Build reverse mapping for efficiency
    ankerl::unordered_dense::map<Index, Index> new_to_old_index;
    for (const auto& [old_idx, new_idx] : old_to_new_index) {
        new_to_old_index[new_idx] = old_idx;
    }

    // Convert the extracted cycles back to original indices
    SubcycleList subcycles;
    for (size_t i = 0; i < cycle_indices.size() - 1; i++) {
        Index start = cycle_indices[i];
        Index end = cycle_indices[i + 1];

        SmallVector subcycle;
        for (Index j = start; j < end; j++) {
            // Map back to original segment index
            Index mapped_idx = index(extracted_cycles[j]);
            bool ori = orientation(extracted_cycles[j]);
            Index original_idx = new_to_old_index[mapped_idx];

            subcycle.push_back(signed_index(original_idx, ori));
        }

        // Only include cycles with at least 3 segments
        if (subcycle.size() >= 3) {
            subcycles.push_back(subcycle);
        }
    }

    return subcycles;
}

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

template <int dim>
std::vector<Index> triangulate(
    std::vector<Scalar>& vertices,
    std::vector<Index>& segments,
    std::vector<SignedIndex>& cycles,
    std::vector<Index>& cycle_indices,
    bool optimal_triangulation)
{
    assert(vertices.size() % dim == 0);
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

            bool use_diag_02 = std::min(cycle_vertices[0], cycle_vertices[2]) <
                               std::min(cycle_vertices[1], cycle_vertices[3]);
            if (optimal_triangulation) {
                std::span<Scalar, dim> v0{vertices.data() + cycle_vertices[0] * dim, dim};
                std::span<Scalar, dim> v1{vertices.data() + cycle_vertices[1] * dim, dim};
                std::span<Scalar, dim> v2{vertices.data() + cycle_vertices[2] * dim, dim};
                std::span<Scalar, dim> v3{vertices.data() + cycle_vertices[3] * dim, dim};

                Scalar diag_02 = 0;
                Scalar diag_13 = 0;
                // Only spatial coordinates are used to compute the diagonal length.
                for (int d = 0; d + 1 < dim; d++) {
                    diag_02 += (v2[d] - v0[d]) * (v2[d] - v0[d]);
                    diag_13 += (v3[d] - v1[d]) * (v3[d] - v1[d]);
                }
                if (std::abs(diag_02 - diag_13) >= 1e-6) {
                    // Only pick a diagonal if the two diagonals are not equal.
                    // If almost equal, still pick the one based on vertex index to ensure
                    // uniqueness.
                    use_diag_02 = diag_02 < diag_13;
                }
            }

            if (use_diag_02) {
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
std::tuple<std::vector<Index>, std::vector<Index>, Index> compute_zero_crossing_vertices(
    const Contour<dim>& contour,
    std::span<const Scalar> function_values,
    std::span<const Scalar> function_gradients,
    bool use_snapping,
    Contour<dim>& result)
{
    static_assert(dim == 3 || dim == 4, "dim must be 3 or 4");
    size_t num_segments = contour.get_num_segments();
    std::vector<Index> zero_crossing_vertices(num_segments, invalid_index);

    bool has_gradients = function_gradients.size() == function_values.size() * dim;

    // For snapping
    constexpr Scalar snap_alpha = 0.1;
    std::vector<Index> snap_vertices(function_values.size(), invalid_index);
    std::vector<Scalar> snap_alphas(function_values.size(), 1);
    std::vector<bool> snap_segments(num_segments, false);
    std::vector<bool> snap_ends(num_segments, false);

    // First pass: compute all potential vertices
    for (size_t i = 0; i < num_segments; i++) {
        auto seg = contour.get_segment(i);
        Index v0 = seg[0];
        Index v1 = seg[1];

        if (function_values[v0] >= 0 && function_values[v1] < 0 ||
            function_values[v0] < 0 && function_values[v1] >= 0) {
            const Scalar val0 = function_values[v0];
            const Scalar val1 = function_values[v1];
            auto pos0 = contour.get_vertex(v0);
            auto pos1 = contour.get_vertex(v1);
            Scalar t;

            if (has_gradients) {
                auto grad0 = function_gradients.subspan(v0 * dim, dim);
                auto grad1 = function_gradients.subspan(v1 * dim, dim);
                // directional derivatives
                Scalar g0, g1;
                // since the interval is scaled to (0,1), we don't normalize by dividing by the
                // length of p0p1
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
                t = get_cubic_root(val0, val1, g0, g1);
            } else {
                t = val0 / (val0 - val1);
            }

            // Compute the position
            std::array<Scalar, dim> position;
            for (int d = 0; d < dim; d++) {
                position[d] = pos0[d] + t * (pos1[d] - pos0[d]);
            }
            result.add_vertex(position);

            // Record this as the current vertex for this segment
            zero_crossing_vertices[i] = static_cast<Index>(result.get_num_vertices() - 1);

            // Update snapping information if needed
            if (use_snapping) {
                snap_segments[i] = (t < snap_alpha || t > 1 - snap_alpha);
                snap_ends[i] = (t > 1 - snap_alpha);
                if (t < snap_alpha && t < snap_alphas[v0]) {
                    snap_vertices[v0] = zero_crossing_vertices[i];
                    snap_alphas[v0] = t;
                } else if (t > 1 - snap_alpha && (1 - t) < snap_alphas[v1]) {
                    snap_vertices[v1] = zero_crossing_vertices[i];
                    snap_alphas[v1] = 1 - t;
                }
            }
        }
    }

    // Second pass: compute snap map.
    ankerl::unordered_dense::map<Index, Index> snap_map;
    for (size_t i = 0; i < num_segments; i++) {
        if (zero_crossing_vertices[i] == invalid_index) {
            continue;
        }

        if (use_snapping && snap_segments[i]) {
            auto seg = contour.get_segment(i);
            Index v0 = seg[0];
            Index v1 = seg[1];
            if (snap_ends[i]) {
                assert(snap_vertices[v1] != invalid_index);
                snap_map[zero_crossing_vertices[i]] = snap_vertices[v1];
            } else {
                assert(snap_vertices[v0] != invalid_index);
                snap_map[zero_crossing_vertices[i]] = snap_vertices[v0];
            }
        }
    }

    // Remap vertex indices
    const Index num_vertices = result.get_num_vertices();
    std::vector<Index> vertex_map(num_vertices);
    std::iota(vertex_map.begin(), vertex_map.end(), 0);
    for (const auto& [key, value] : snap_map) {
        vertex_map[key] = value;
    }

    // Sort in-place and find unique (avoids copy)
    std::sort(vertex_map.begin(), vertex_map.end());
    auto last = std::unique(vertex_map.begin(), vertex_map.end());
    Index num_unique_vertices = static_cast<Index>(std::distance(vertex_map.begin(), last));

    // Build compact unique_map using hash map for sparse mapping
    ankerl::unordered_dense::map<Index, Index> unique_map;
    unique_map.reserve(num_unique_vertices);
    for (Index vid = 0; vid < num_unique_vertices; vid++) {
        unique_map[vertex_map[vid]] = vid;
    }

    // Rebuild vertex_map with final mapping
    std::iota(vertex_map.begin(), vertex_map.end(), 0);
    for (const auto& [key, value] : snap_map) {
        vertex_map[key] = value;
    }
    for (Index vid = 0; vid < num_vertices; vid++) {
        Index vnew = vertex_map[vid];
        assert(unique_map.count(vnew) > 0);
        vertex_map[vid] = unique_map[vnew];
    }

    return std::make_tuple(
        std::move(zero_crossing_vertices),
        std::move(vertex_map),
        num_unique_vertices);
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
    llvm_vecsmall::SmallVector<SignedIndex, 16> zero_crossing_segments;
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
            }
        }
        if (!zero_crossing_segments.empty()) {
            // Determine segment orientation
            SignedIndex first_segment = zero_crossing_segments[0];
            Index first_segment_id = index(first_segment);
            bool first_segment_ori = orientation(first_segment);
            Index first_vertex_id =
                contour.get_segment(first_segment_id)[first_segment_ori ? 0 : 1];
            bool first_vertex_sign = function_values[first_vertex_id] >= 0;
            size_t parity = first_vertex_sign ? 1 : 0;

            // Add diagonal segments so that the negative region stays connected.
            // The number of zero-crossing segments is assumed to be even because each crossing
            // corresponds to a transition between positive and negative regions, forming pairs.
            // Parity is used to determine the starting point for pairing segments, ensuring
            // consistent connectivity. The index manipulation below ensures that segments are
            // paired correctly based on the calculated parity.
            size_t num_zero_crossing_segments = zero_crossing_segments.size();
            assert(num_zero_crossing_segments % 2 == 0);
            for (size_t j = 0; j < num_zero_crossing_segments; j += 2) {
                size_t idx0 = (j + parity) % num_zero_crossing_segments;
                size_t idx1 = (j + 1 + parity) % num_zero_crossing_segments;

                auto seg_id_0 = index(zero_crossing_segments[idx0]);
                bool seg_ori_0 = orientation(zero_crossing_segments[idx0]);
                auto seg_id_1 = index(zero_crossing_segments[idx1]);
                bool seg_ori_1 = orientation(zero_crossing_segments[idx1]);

                // skip degenerate segment
                if (zero_crossing_vertices[seg_id_0] == zero_crossing_vertices[seg_id_1]) {
                    continue;
                }

                Index v0 = contour.get_segment(seg_id_0)[seg_ori_0 ? 1 : 0];
                assert(function_values[v0] >= 0);

                result.add_segment(
                    zero_crossing_vertices[seg_id_1],
                    zero_crossing_vertices[seg_id_0]);
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

template <int dim>
Contour<dim> remap_vertices(
    const Contour<dim>& input,
    std::span<const Index> vertex_map,
    Index num_snapped_vertices)
{
    assert(input.get_num_vertices() == vertex_map.size());
    Index num_vertices_before_snap = input.get_num_vertices();

    Contour<dim> result;

    // Map vertices
    result.resize_vertices(num_snapped_vertices);
    for (Index vid = 0; vid < num_vertices_before_snap; vid++) {
        Index vid_new = vertex_map[vid];
        assert(vid_new != invalid_index);
        assert(vid_new < num_snapped_vertices);
        auto old_pos = input.get_vertex(vid);
        auto new_pos = result.ref_vertex(vid_new);
        std::copy(old_pos.begin(), old_pos.end(), new_pos.begin());
    }

    // Map segments
    Index num_segments = input.get_num_segments();
    std::vector<Index> segment_map(num_segments, invalid_index);
    for (Index sid = 0; sid < num_segments; sid++) {
        auto seg = input.get_segment(sid);
        Index v0_new = vertex_map[seg[0]];
        Index v1_new = vertex_map[seg[1]];
        assert(v0_new != invalid_index && v1_new != invalid_index);
        assert(v0_new < num_snapped_vertices && v1_new < num_snapped_vertices);
        if (v0_new == v1_new) {
            // degenerate segment after snapping, skip
            continue;
        } else {
            result.add_segment(v0_new, v1_new);
            segment_map[sid] = result.get_num_segments() - 1;
        }
    }

    using SmallVector = llvm_vecsmall::SmallVector<SignedIndex, 16>;

    // Create a callback to get segment endpoints
    auto get_segment = [&](Index seg_idx) -> std::array<Index, 2> {
        auto seg = result.get_segment(seg_idx);
        return {seg[0], seg[1]};
    };

    // Map cycles
    // Each original cycle may map to multiple subcycles if it contains duplicate vertices
    Index num_cycles = input.get_num_cycles();
    using CycleIndexList = llvm_vecsmall::SmallVector<Index, 16>;
    std::vector<CycleIndexList> cycle_map(num_cycles);
    for (Index cid = 0; cid < num_cycles; cid++) {
        auto cycle = input.get_cycle(cid);
        SmallVector new_cycle;
        for (auto sid : cycle) {
            Index old_sid = index(sid);
            bool seg_ori = orientation(sid);
            Index new_sid = segment_map[old_sid];
            if (new_sid != invalid_index) {
                new_cycle.push_back(signed_index(new_sid, seg_ori));
            }
        }
        if (new_cycle.size() >= 3) {
            // Extract subcycles in case the cycle contains duplicate vertices
            auto subcycles = extract_simple_loops(
                std::span<const SignedIndex>(new_cycle.data(), new_cycle.size()),
                get_segment);
            // Add all subcycles and track their indices
            for (const auto& subcycle : subcycles) {
                result.add_cycle(
                    std::span<const SignedIndex>(subcycle.data(), subcycle.size()),
                    input.is_cycle_regular(cid));
                Index new_cid = result.get_num_cycles() - 1;
                cycle_map[cid].push_back(new_cid);
            }
        }
    }

    if constexpr (dim == 4) {
        // Map polyhedra
        Index num_polyhedra = input.get_num_polyhedra();
        for (Index pid = 0; pid < num_polyhedra; pid++) {
            auto polyhedron = input.get_polyhedron(pid);
            llvm_vecsmall::SmallVector<SignedIndex, 16> new_polyhedron;
            for (auto cid : polyhedron) {
                Index old_cid = index(cid);
                bool cycle_ori = orientation(cid);
                // Add all subcycles that the original cycle was split into
                for (Index new_cid : cycle_map[old_cid]) {
                    new_polyhedron.push_back(signed_index(new_cid, cycle_ori));
                }
            }
            if (new_polyhedron.size() >= 4) {
                result.add_polyhedron(
                    std::span<const SignedIndex>(new_polyhedron.data(), new_polyhedron.size()),
                    input.is_polyhedron_regular(pid));
            }
        }
    }

    return result;
}


} // namespace

template <>
void Contour<3>::triangulate_cycles(bool optimal_triangulation)
{
    auto cycle_to_triangle_map = triangulate<3>(
        m_vertices,
        m_segments,
        m_cycles,
        m_cycle_start_indices,
        optimal_triangulation);

    // Update the cycle_is_regular vector
    map_cycle_regularity(cycle_to_triangle_map, m_cycle_is_regular);

#ifndef NDEBUG
    check_all_segments();
    check_all_cycles();
    check_all_polyhedra();
#endif
}

template <>
void Contour<4>::triangulate_cycles(bool optimal_triangulation)
{
#ifndef NDEBUG
    check_all_segments();
    check_all_cycles();
    check_all_polyhedra();
#endif

    auto cycle_to_triangle_map = triangulate<4>(
        m_vertices,
        m_segments,
        m_cycles,
        m_cycle_start_indices,
        optimal_triangulation);

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
    std::span<Scalar> function_gradients,
    bool use_snapping) const
{
    assert(get_num_vertices() == function_values.size());
    Contour<3> result;

    auto [zero_crossing_vertices, vertex_map, num_snapped_vertices] =
        compute_zero_crossing_vertices(
            *this,
            function_values,
            function_gradients,
            use_snapping,
            result);
    compute_zero_crossing_segments(*this, function_values, zero_crossing_vertices, result);

    result = remap_vertices(result, vertex_map, num_snapped_vertices);

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
    std::span<Scalar> function_gradients,
    bool use_snapping) const
{
    assert(get_num_vertices() == function_values.size());
    Contour<4> result;

    auto [zero_crossing_vertices, vertex_map, num_snapped_vertices] =
        compute_zero_crossing_vertices(
            *this,
            function_values,
            function_gradients,
            use_snapping,
            result);

    std::vector<Index> zero_crossing_segment_indices =
        compute_zero_crossing_segments(*this, function_values, zero_crossing_vertices, result);

    size_t num_polyhedra = get_num_polyhedra();
    NonDisjointCycles cycles(result.m_segments);
    for (size_t i = 0; i < num_polyhedra; i++) {
        cycles.clear();
        auto polyhedron = get_polyhedron(i);
        for (auto cid : polyhedron) {
            Index seg_start = zero_crossing_segment_indices[index(cid)];
            Index seg_end = zero_crossing_segment_indices[index(cid) + 1];
            bool cycle_ori = orientation(cid);

            for (Index seg_id = seg_start; seg_id < seg_end; seg_id++) {
                cycles.register_segment(signed_index(seg_id, cycle_ori));
            }
        }
        size_t num_existing_cycles = result.m_cycle_start_indices.size();
        cycles.extract_cycles(result.m_cycles, result.m_cycle_start_indices);

        for (size_t j = num_existing_cycles; j < result.m_cycle_start_indices.size(); j++) {
            result.m_cycle_is_regular.push_back(m_polyhedron_is_regular[i]);
        }
    }

#ifndef NDEBUG
    result.check_all_segments();
    result.check_all_cycles();
    result.check_all_polyhedra();
#endif

    result = remap_vertices(result, vertex_map, num_snapped_vertices);

#ifndef NDEBUG
    result.check_all_segments();
    result.check_all_cycles();
    result.check_all_polyhedra();
#endif

    return result;
}

} // namespace mtetcol
