#pragma once

#include <exception>
#include <span>
#include <vector>

namespace mtetcol {

using Scalar = double;
using Index = uint32_t;

class Contour
{
public:
    void add_vertex(Scalar x, Scalar y, Scalar z, Scalar t)
    {
        m_vertices.push_back(x);
        m_vertices.push_back(y);
        m_vertices.push_back(z);
        m_vertices.push_back(t);
    }
    std::span<const Scalar, 4> get_vertex(Index vid) const
    {
        return std::span<Scalar, 4>(m_vertices.data() + vid * 4, 4);
    }

    void add_segment(Index v0, Index v1)
    {
        m_segments.push_back(v0);
        m_segments.push_back(v1);
    }
    std::span<const Index, 2> get_segment(Index segid) const
    {
        return std::span<Index, 2>(m_segments.data() + segid * 2, 2);
    }

    void add_cycle(std::span<const Index> cycle, std::span<const bool> orientations)
    {
        if (cycle.size() != orientations.size()) {
            throw std::invalid_argument("Cycle and orientations must have the same size");
        }

        for (size_t i = 0; i < cycle.size(); ++i) {
            m_cycles.push_back(cycle[i]);
            m_cycle_segment_orientations.push_back(orientations[i]);
        }
        m_cycle_start_indices.push_back(m_cycles.size());
    }
    std::tuple<std::span<const Index>, std::span<const bool>> get_cycle(Index cid) const
    {
        if (cid >= m_cycle_start_indices.size() - 1) {
            throw std::out_of_range("Cycle ID out of range");
        }

        Index start = m_cycle_start_indices[cid];
        Index end = m_cycle_start_indices[cid + 1];

        return {
            std::span<Index>(m_cycles.data() + start, end - start),
            std::span<bool>(m_cycle_segment_orientations.data() + start, end - start)};
    }

    void add_polyhedron(std::span<const Index> polyhedron, std::span<const bool> orientations)
    {
        if (polyhedron.size() != orientations.size()) {
            throw std::invalid_argument("Polyhedron and orientations must have the same size");
        }

        for (size_t i = 0; i < polyhedron.size(); ++i) {
            m_polyhedra.push_back(polyhedron[i]);
            m_polyhedron_cycle_orientations.push_back(orientations[i]);
        }
        m_polyhedron_start_indices.push_back(m_polyhedra.size());
    }
    std::tuple<std::span<const Index>, std::span<const bool>> get_polyhedron(Index poly_id) const
    {
        if (poly_id >= m_polyhedron_start_indices.size() - 1) {
            throw std::out_of_range("Polyhedron ID out of range");
        }

        Index start = m_polyhedron_start_indices[poly_id];
        Index end = m_polyhedron_start_indices[poly_id + 1];

        return {
            std::span<Index>(m_polyhedra.data() + start, end - start),
            std::span<bool>(m_polyhedron_cycle_orientations.data() + start, end - start)};
    }


private:
    /**
     * @brief The vertices of the contour.
     *
     * Each vertex is represented by (x, y, z, t).
     */
    std::vector<Scalar> m_vertices;

    /**
     * @brief Contour edge segements.
     *
     * Each segment is represented as an oriented edge (v_0, v_1).
     */
    std::vector<Index> m_segments;

    /**
     * @brief Contour cycles.
     *
     * Each cycle is a chain of segments that form a closed loop.
     */
    std::vector<Index> m_cycles;
    std::vector<Index> m_cycle_start_indices = {0};
    std::vector<bool> m_cycle_segment_orientations;

    /**
     * @brief Contour polyhedra.
     *
     * Each polyhedron is a list of oriented cycles that forms the faces of the polyhedron.
     */
    std::vector<Index> m_polyhedra;
    std::vector<Index> m_polyhedron_start_indices = {0};
    std::vector<bool> m_polyhedron_cycle_orientations;
};

} // namespace mtetcol
