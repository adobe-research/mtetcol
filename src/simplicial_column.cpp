#include <ankerl/unordered_dense.h>
#include <mtetcol/simplicial_column.h>

#include <array>
#include <cassert>

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
    for (size_t i=0; i < num_edges; i++) {
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

    for (size_t i=0; i<num_triangles; i++) {
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

    for (size_t i=0; i<num_tets; i++) {
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

        for (size_t j=0; j<3; j++) {
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
Contour<4> SimplicialColumn<4>::extract_contour(Scalar value) const
{
    // TODO
    return Contour<4>();
}

template <>
Contour<3> SimplicialColumn<3>::extract_contour(Scalar value) const
{
    // TODO
    return Contour<3>();
}

} // namespace mtetcol
