#pragma once

#include <algorithm>
#include <exception>
#include <iterator>
#include <span>
#include <vector>

#include <mtetcol/common.h>

namespace mtetcol {

/**
 * @brief A class representing a space-time contour composed of vertices, segments, cycles, and
 * polyhedra.
 *
 * @tparam dim The dimension of the ambient space. Default is 4 (3D space + 1D time).
 *
 * The Contour class provides a data structure for representing and manipulating 4D contours in
 * space-time. Each vertex has (x,y,t) or (x,y,z,t) coordinates, and these vertices are connected by
 * oriented segments. Segments can be grouped into cycles (closed loops), and cycles can be grouped
 * into polyhedra.
 */
template <int dim = 4>
class Contour
{
public:
    static_assert(dim == 3 || dim == 4, "Ambient dimension of contour must be 3 or 4");

    /**
     * @brief Adds a new vertex to the contour
     *
     * @param vertex A span of coordinates representing the vertex
     * @throws std::runtime_error if dim is not 4
     */
    void add_vertex(std::span<const Scalar, dim> vertex)
    {
        m_vertices.insert(m_vertices.end(), vertex.begin(), vertex.end());
    }

    /**
     * @brief Retrieves the coordinates of a vertex
     *
     * @param vid The index of the vertex to retrieve
     * @return std::span<const Scalar, 4> A span containing the (x,y,z,t) coordinates of the vertex
     * @throws std::out_of_range if vid is out of range
     */
    std::span<const Scalar, dim> get_vertex(Index vid) const
    {
        return std::span<Scalar, dim>(m_vertices.data() + vid * dim, dim);
    }

    /**
     * @brief Retrieves the number of vertices in the contour
     *
     * @return size_t The number of vertices
     */
    size_t get_num_vertices() const
    {
        return m_vertices.size() / dim;
    }

    /**
     * @brief Adds a new segment to the contour
     *
     * @param v0 Index of the first vertex
     * @param v1 Index of the second vertex
     */
    void add_segment(Index v0, Index v1)
    {
        m_segments.push_back(v0);
        m_segments.push_back(v1);
    }
    /**
     * @brief Retrieves the vertex indices of a segment
     *
     * @param segid The index of the segment to retrieve
     * @return std::span<const Index, 2> A span containing the two vertex indices
     * @throws std::out_of_range if segid is out of range
     */
    std::span<const Index, 2> get_segment(Index segid) const
    {
        return std::span<Index, 2>(m_segments.data() + segid * 2, 2);
    }

    /**
     * @brief Retrieves the number of segments in the contour
     *
     * @return size_t The number of segments
     */
    size_t get_num_segments() const
    {
        return m_segments.size() / 2;
    }

    /**
     * @brief Adds a new cycle to the contour
     *
     * A cycle is a closed loop of segments. Each segment in the cycle can be oriented
     * either in its original direction or reversed.
     *
     * @param cycle Span of segment signed indices that form the cycle
     */
    void add_cycle(std::span<const SignedIndex> cycle)
    {
        std::copy(cycle.begin(), cycle.end(), std::back_inserter(m_cycles));
        m_cycle_start_indices.push_back(m_cycles.size());
    }

    /**
     * @brief Retrieves a cycle from the contour
     *
     * @param cid The index of the cycle to retrieve
     *
     * @return    A span containing an ordered list of signed segment indices of the cycle
     *
     * @throws std::out_of_range if cid is out of range
     */
    std::span<const SignedIndex> get_cycle(Index cid) const
    {
        if (cid >= m_cycle_start_indices.size() - 1) {
            throw std::out_of_range("Cycle ID out of range");
        }

        Index start = m_cycle_start_indices[cid];
        Index end = m_cycle_start_indices[cid + 1];

        return std::span<SignedIndex>(m_cycles.data() + start, end - start);
    }

    /**
     * @brief Retrieves the number of cycles in the contour
     *
     * @return size_t The number of cycles
     */
    size_t get_num_cycles() const
    {
        assert(!m_cycle_start_indices.empty());
        return m_cycle_start_indices.size()  - 1;
    }

    /**
     * @brief Adds a new polyhedron to the contour
     *
     * A polyhedron is defined by a collection of oriented cycles that form its faces.
     *
     * @param polyhedron Span of cycle indices that form the polyhedron faces
     */
    void add_polyhedron(std::span<const SignedIndex> polyhedron)
    {
        std::copy(polyhedron.begin(), polyhedron.end(), std::back_inserter(m_polyhedra));
        m_polyhedron_start_indices.push_back(m_polyhedra.size());
    }

    /**
     * @brief Retrieves a polyhedron from the contour
     *
     * @param poly_id The index of the polyhedron to retrieve
     * @return A span containing a list of signed cycles of the polyhedron
     *
     * @throws std::out_of_range if poly_id is out of range
     */
    std::span<const SignedIndex> get_polyhedron(Index poly_id) const
    {
        if (poly_id >= m_polyhedron_start_indices.size() - 1) {
            throw std::out_of_range("Polyhedron ID out of range");
        }

        Index start = m_polyhedron_start_indices[poly_id];
        Index end = m_polyhedron_start_indices[poly_id + 1];

        return std::span<SignedIndex>(m_polyhedra.data() + start, end - start);
    }

    /**
     * @brief Retrieves the number of polyhedra in the contour
     *
     * @return size_t The number of polyhedra
     */
    size_t get_num_polyhedra() const
    {
        assert(!m_polyhedron_start_indices.empty());
        return m_polyhedron_start_indices.size() - 1;
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
    std::vector<SignedIndex> m_cycles;
    std::vector<Index> m_cycle_start_indices = {0};

    /**
     * @brief Contour polyhedra.
     *
     * Each polyhedron is a list of oriented cycles that forms the faces of the polyhedron.
     */
    std::vector<SignedIndex> m_polyhedra;
    std::vector<Index> m_polyhedron_start_indices = {0};
};

} // namespace mtetcol
