#include <mtetcol/simplicial_column.h>

#include "hashmap.h"
#include "logger.h"
#include "utils.h"

#include <spdlog/spdlog.h>

#include <array>
#include <cassert>
#include <span>

namespace mtetcol {

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
    auto [contour_times, contour_time_indices, initial_signs] = extract_contour_vertices(
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

    size_t num_cycles = contour_cycle_indices.size() - 1;
    for (size_t i = 0; i < num_cycles; i++) {
        contour.add_cycle(std::span<SignedIndex>(
            contour_cycles.data() + contour_cycle_indices[i],
            contour_cycle_indices[i + 1] - contour_cycle_indices[i]));
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
