#include <ankerl/unordered_dense.h>
#include <mtetcol/simplicial_column.h>

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

namespace mtetcol {

namespace {
using Edge = std::array<Index, 2>; // [v0, v1]
using Triangle = std::array<Index, 3>; // [v0, v1, v2]

struct EdgeHash
{
    using is_transparent = void;
    using is_avalanching = void;

    [[nodiscard]] auto operator()(const Edge& edge) const noexcept -> uint64_t
    {
        ankerl::unordered_dense::hash<uint32_t> hash_fn;
        return hash_fn(edge[0]) ^ hash_fn(edge[1]);
    }
};

struct EdgeEqual
{
    using is_transparent = void;

    [[nodiscard]] bool operator()(const Edge& lhs, const Edge& rhs) const noexcept
    {
        return (lhs[0] == rhs[0] && lhs[1] == rhs[1]) || (lhs[0] == rhs[1] && lhs[1] == rhs[0]);
    }
};

struct TriangleHash
{
    using is_transparent = void;
    using is_avalanching = void;

    [[nodiscard]] auto operator()(const Triangle& tri) const noexcept -> uint64_t
    {
        ankerl::unordered_dense::hash<uint32_t> hash_fn;
        return hash_fn(tri[0]) ^ hash_fn(tri[1]) ^ hash_fn(tri[2]);
    }
};

struct TriangleEqual
{
    using is_transparent = void;

    [[nodiscard]] bool operator()(const Triangle& lhs, const Triangle& rhs) const noexcept
    {
        return (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2]) ||
               (lhs[0] == rhs[1] && lhs[1] == rhs[2] && lhs[2] == rhs[0]) ||
               (lhs[0] == rhs[2] && lhs[1] == rhs[0] && lhs[2] == rhs[1]) ||
               (lhs[0] == rhs[0] && lhs[1] == rhs[2] && lhs[2] == rhs[1]) ||
               (lhs[0] == rhs[1] && lhs[1] == rhs[0] && lhs[2] == rhs[2]) ||
               (lhs[0] == rhs[2] && lhs[1] == rhs[1] && lhs[2] == rhs[0]);
    }
};

/**
 * Mapping from [v0, v1] to edge index.
 */
using EdgeMap = ankerl::unordered_dense::map<Edge, Index, EdgeHash, EdgeEqual>;

/**
 * Mapping from [v0, v1, v2] to triangle index.
 */
using TriangleMap = ankerl::unordered_dense::map<Triangle, Index, TriangleHash, TriangleEqual>;

SignedIndex add_edge(const Edge& e, EdgeMap& edges)
{
    auto it = edges.find(e);
    if (it != edges.end()) {
        Index eid = it->second;
        return signed_index(eid, e[0] < e[1]);
    } else {
        Index eid = static_cast<Index>(edges.size());
        edges[e] = eid;
        return signed_index(eid, e[0] < e[1]);
    }
}

SignedIndex get_edge(const Edge& e, const EdgeMap& edges)
{
    auto it = edges.find(e);
    assert(it != edges.end());

    Index eid = it->second;
    return signed_index(eid, e[0] < e[1]);
}

SignedIndex add_triangle(const Triangle& t, TriangleMap& triangles)
{
    const int8_t sign = (t[0] < t[1] ? 1 : -1) + (t[1] < t[2] ? 1 : -1) + (t[2] < t[0] ? 1 : -1);
    assert(sign == 1 || sign == -1);

    auto it = triangles.find(t);
    if (it != triangles.end()) {
        Index tid = it->second;
        return signed_index(tid, sign == 1);
    } else {
        Index tid = static_cast<Index>(triangles.size());
        triangles[t] = tid;
        return signed_index(tid, sign == 1);
    }
}

template <int dim>
bool check_edges(const SimplicialColumn<dim>& columns)
{
    const size_t num_vertices = columns.get_num_spatial_vertices();
    const size_t num_edges = columns.get_num_spatial_edges();

    auto edges = columns.get_spatial_edges();
    for (size_t i = 0; i < num_edges; i++) {
        Index v0 = edges[i * 2];
        Index v1 = edges[i * 2 + 1];
        if (v0 >= v1) return false;
        if (v0 >= num_vertices) return false;
        if (v1 >= num_vertices) return false;
    }
    return true;
}

template <int dim>
bool check_triangles(const SimplicialColumn<dim>& columns)
{
    const size_t num_triangles = columns.get_num_spatial_triangles();
    auto edges = columns.get_spatial_edges();
    auto triangles = columns.get_spatial_triangles();

    auto check_adj_edges = [&](SignedIndex e0, SignedIndex e1) {
        Index e0_id = index(e0);
        Index e1_id = index(e1);
        bool e0_ori = orientation(e0);
        bool e1_ori = orientation(e1);

        // Local ids
        Index e0_li = e0_ori ? 1 : 0;
        Index e1_li = e1_ori ? 0 : 1;

        return edges[e0_id * 2 + e0_li] == edges[e1_id * 2 + e1_li];
    };

    for (size_t i = 0; i < num_triangles; i++) {
        SignedIndex e01 = triangles[i * 3];
        SignedIndex e12 = triangles[i * 3 + 1];
        SignedIndex e20 = triangles[i * 3 + 2];

        if (!check_adj_edges(e01, e12)) return false;
        if (!check_adj_edges(e12, e20)) return false;
        if (!check_adj_edges(e20, e01)) return false;
    }

    return true;
}

template <int dim>
bool check_tetrahedra(const SimplicialColumn<dim>& columns)
{
    const size_t num_tets = columns.get_num_spatial_tetrahedra();

    auto triangles = columns.get_spatial_triangles();
    auto tets = columns.get_spatial_tetrahedra();

    for (size_t i = 0; i < num_tets; i++) {
        SignedIndex t021 = tets[i * 4];
        SignedIndex t123 = tets[i * 4 + 1];
        SignedIndex t013 = tets[i * 4 + 2];
        SignedIndex t032 = tets[i * 4 + 3];

        Index t021_id = index(t021);
        Index t123_id = index(t123);
        Index t013_id = index(t013);
        Index t032_id = index(t032);

        int32_t t021_ori = orientation(t021) ? 1 : -1;
        int32_t t123_ori = orientation(t123) ? 1 : -1;
        int32_t t013_ori = orientation(t013) ? 1 : -1;
        int32_t t032_ori = orientation(t032) ? 1 : -1;

        int32_t sum = 0;

        for (size_t j = 0; j < 3; j++) {
            sum += value_of(triangles[t021_id * 3 + j]) * t021_ori;
            sum += value_of(triangles[t123_id * 3 + j]) * t123_ori;
            sum += value_of(triangles[t013_id * 3 + j]) * t013_ori;
            sum += value_of(triangles[t032_id * 3 + j]) * t032_ori;
        }

        assert(sum == 0);
        if (sum != 0) return false;
    }
    return true;
}

/**
 * @brief Extracts the zero-crossing times of a vertex.
 *
 * @note Function values that equals to `value` will be treated as have positive sign.
 *
 * @note When `cyclic` is false, we will add 0 as a zero-crossing time if the first function value
 * is non-negative, and add 1 as a zero-crossing time if the last function value is negative.
 *
 * @param[in] time_samples The time samples of the vertex.
 * @param[in] function_values The function values of the vertex.
 * @param[in] value The value to check for zero-crossing.
 * @param[in] cyclic Whether the time samples are cyclic.
 * @param[out] zero_crossing_times The output vector to store the zero-crossing times.
 */
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

template <int dim>
std::tuple<std::vector<Scalar>, std::vector<size_t>, std::vector<bool>> extract_contour_vertices(
    const std::span<Scalar> vertices,
    const std::vector<Scalar>& time_samples,
    const std::vector<Scalar>& function_values,
    const std::vector<Index>& vertex_start_indices,
    Scalar value,
    bool cyclic)
{
    assert(vertices.size() % (dim - 1) == 0);
    assert(vertices.size() / (dim - 1) + 1 == vertex_start_indices.size());
    size_t num_vertices = vertices.size() / (dim - 1);

    std::vector<Scalar> zero_crossing_times;
    std::vector<size_t> zero_crossing_indices;
    std::vector<bool> initial_signs;
    zero_crossing_times.reserve(time_samples.size());
    zero_crossing_indices.reserve(num_vertices + 1);
    zero_crossing_indices.push_back(0);
    initial_signs.reserve(num_vertices);
    for (size_t vi = 0; vi < num_vertices; vi++) {
        std::span<Scalar> position = vertices.subspan(vi * (dim - 1), dim - 1);
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
        std::cout << "edge: (" << v0 << ", " << v1 << ")" << std::endl;
    };

    auto print_segment = [&](Index eid, bool ori) {
        Index start = contour_segment_indices[eid];
        Index end = contour_segment_indices[eid + 1];

        Index num_segments_over_edge = (end - start) / 2;
        for (size_t si = 0; si < num_segments_over_edge; si++) {
            Index v0 = contour_segments[start + si * 2];
            Index v1 = contour_segments[start + si * 2 + 1];
            if (!ori) std::swap(v0, v1);

            std::cout << "segment: (" << v0 << ", " << v1 << ")" << std::endl;
        }
    };

    auto register_next_index = [&](Index eid, bool ori) {
        std::cout << "eid: " << eid << ", ori: " << ori << std::endl;
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
            std::cout << "Registering next index: (" << v0 << ", " << v1 << ") -> "
                      << next_index[v0] << std::endl;
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
        std::cout << "Triangle " << ti << ": " << std::endl;
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

} // namespace

template <>
void SimplicialColumn<4>::set_simplices(std::span<Index> simplices)
{
    assert(simplices.size() % 4 == 0);

    EdgeMap edges;
    TriangleMap triangles;

    size_t num_tets = simplices.size() / 4;
    m_tetrahedra.reserve(num_tets * 4);
    for (size_t i = 0; i < num_tets; i++) {
        auto tet = simplices.subspan(i * 4, 4);

        add_edge({tet[0], tet[1]}, edges);
        add_edge({tet[0], tet[2]}, edges);
        add_edge({tet[0], tet[3]}, edges);
        add_edge({tet[1], tet[2]}, edges);
        add_edge({tet[1], tet[3]}, edges);
        add_edge({tet[2], tet[3]}, edges);

        SignedIndex t021 = add_triangle({tet[0], tet[2], tet[1]}, triangles);
        SignedIndex t123 = add_triangle({tet[1], tet[2], tet[3]}, triangles);
        SignedIndex t013 = add_triangle({tet[0], tet[1], tet[3]}, triangles);
        SignedIndex t032 = add_triangle({tet[0], tet[3], tet[2]}, triangles);

        m_tetrahedra.push_back(t021);
        m_tetrahedra.push_back(t123);
        m_tetrahedra.push_back(t013);
        m_tetrahedra.push_back(t032);
    }

    m_edges.resize(edges.size() * 2, invalid_index);
    for (auto it = edges.begin(); it != edges.end(); ++it) {
        auto e = it->first;
        Index eid = it->second;
        m_edges[eid * 2] = e[0];
        m_edges[eid * 2 + 1] = e[1];
        if (e[0] >= e[1]) {
            std::swap(m_edges[eid * 2], m_edges[eid * 2 + 1]);
        }
    }

    m_triangles.resize(triangles.size() * 3, invalid_signed_index);
    for (auto it = triangles.begin(); it != triangles.end(); ++it) {
        auto t = it->first;
        Index tid = it->second;

        // Sort t to be ascending
        if (t[0] > t[1]) std::swap(t[0], t[1]);
        if (t[1] > t[2]) std::swap(t[1], t[2]);
        if (t[0] > t[1]) std::swap(t[0], t[1]);

        m_triangles[tid * 3] = get_edge({t[0], t[1]}, edges);
        m_triangles[tid * 3 + 1] = get_edge({t[1], t[2]}, edges);
        m_triangles[tid * 3 + 2] = get_edge({t[2], t[0]}, edges);
    }
    assert(check_edges(*this));
    assert(check_triangles(*this));
    assert(check_tetrahedra(*this));
}

template <>
void SimplicialColumn<3>::set_simplices(std::span<Index> simplices)
{
    assert(simplices.size() % 3 == 0);

    EdgeMap edges;

    size_t num_tris = simplices.size() / 3;
    m_triangles.reserve(num_tris * 3);
    for (size_t i = 0; i < num_tris; i++) {
        auto tri = simplices.subspan(i * 3, 3);

        SignedIndex e01 = add_edge({tri[0], tri[1]}, edges);
        SignedIndex e12 = add_edge({tri[1], tri[2]}, edges);
        SignedIndex e20 = add_edge({tri[2], tri[0]}, edges);

        m_triangles.push_back(e01);
        m_triangles.push_back(e12);
        m_triangles.push_back(e20);
    }

    m_edges.resize(edges.size() * 2, invalid_index);
    for (auto it = edges.begin(); it != edges.end(); ++it) {
        auto e = it->first;
        Index eid = it->second;
        m_edges[eid * 2] = e[0];
        m_edges[eid * 2 + 1] = e[1];
        if (e[0] >= e[1]) {
            std::swap(m_edges[eid * 2], m_edges[eid * 2 + 1]);
        }
    }

    assert(check_edges(*this));
    assert(check_triangles(*this));
}

template <>
Contour<4> SimplicialColumn<4>::extract_contour(Scalar value, bool cyclic) const
{
    auto [contour_times, contour_time_indices, initial_signs] = extract_contour_vertices<4>(
        m_vertices,
        m_time_samples,
        m_function_values,
        m_vertex_start_indices,
        value,
        cyclic);

    auto [contour_segments, contour_segment_indices] = extract_contour_segments(
        contour_times,
        contour_time_indices,
        initial_signs,
        m_edges,
        cyclic);

    auto [contour_cycles, contour_cycle_indices, contour_cycle_triangle_indices] =
        extract_contour_cycles(
            contour_times.size(),
            contour_segments,
            contour_segment_indices,
            m_edges,
            m_triangles);

    Contour<4> contour;

    size_t num_contour_vertices = contour_time_indices.size() - 1;
    for (size_t i = 0; i < num_contour_vertices; i++) {
        std::span<Scalar> position = m_vertices.subspan(i * 3, 3);
        std::span<Scalar> time_samples(
            contour_times.data() + contour_time_indices[i],
            contour_time_indices[i + 1] - contour_time_indices[i]);
        for (auto t : time_samples) {
            contour.add_vertex({position[0], position[1], position[2], t});
        }
    }

    size_t num_segments = contour_segments.size() / 2;
    for (size_t i = 0; i < num_segments; i++) {
        contour.add_segment(contour_segments[i * 2], contour_segments[i * 2 + 1]);
    }

    return contour;
}

template <>
Contour<3> SimplicialColumn<3>::extract_contour(Scalar value, bool cyclic) const
{
    // TODO
    return Contour<3>();
}

} // namespace mtetcol
