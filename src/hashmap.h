#pragma once

#include <mtetcol/common.h>

#include <ankerl/unordered_dense.h>

#include <array>

namespace mtetcol {

using Edge = std::array<Index, 2>; // [v0, v1]
using Triangle = std::array<Index, 3>; // [v0, v1, v2]

struct EdgeHash
{
    using is_transparent = void;
    // using is_avalanching = void;

    [[nodiscard]] inline auto operator()(const Edge& edge) const noexcept -> uint64_t
    {
        ankerl::unordered_dense::hash<uint32_t> hash_fn;
        return hash_fn(edge[0]) ^ hash_fn(edge[1]);
    }
};

struct EdgeEqual
{
    using is_transparent = void;

    [[nodiscard]] inline bool operator()(const Edge& lhs, const Edge& rhs) const noexcept
    {
        return (lhs[0] == rhs[0] && lhs[1] == rhs[1]) || (lhs[0] == rhs[1] && lhs[1] == rhs[0]);
    }
};

struct TriangleHash
{
    using is_transparent = void;
    // using is_avalanching = void;

    [[nodiscard]] inline auto operator()(const Triangle& tri) const noexcept -> uint64_t
    {
        ankerl::unordered_dense::hash<uint32_t> hash_fn;
        return hash_fn(tri[0]) ^ hash_fn(tri[1]) ^ hash_fn(tri[2]);
    }
};

struct TriangleEqual
{
    using is_transparent = void;

    [[nodiscard]] inline bool operator()(const Triangle& lhs, const Triangle& rhs) const noexcept
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

inline SignedIndex add_edge(const Edge& e, EdgeMap& edges)
{
    auto [it, _] = edges.try_emplace(e, static_cast<Index>(edges.size()));
    return signed_index(it->second, e[0] < e[1]);
}

inline SignedIndex get_edge(const Edge& e, const EdgeMap& edges)
{
    auto it = edges.find(e);
    assert(it != edges.end());

    Index eid = it->second;
    return signed_index(eid, e[0] < e[1]);
}

inline SignedIndex add_triangle(const Triangle& t, TriangleMap& triangles)
{
    const int8_t sign = (t[0] < t[1] ? 1 : -1) + (t[1] < t[2] ? 1 : -1) + (t[2] < t[0] ? 1 : -1);
    assert(sign == 1 || sign == -1);

    auto [it, _] = triangles.try_emplace(t, static_cast<Index>(triangles.size()));
    return signed_index(it->second, sign == 1);
}

} // namespace mtetcol
