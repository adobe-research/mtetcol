#pragma once

#include <algorithm>
#include <cassert>
#include <exception>
#include <iterator>
#include <numeric>
#include <span>
#include <type_traits>
#include <vector>

#include <mtetcol/common.h>
#include <mtetcol/logger.h>

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
     */
    void add_vertex(std::span<const Scalar, dim> vertex)
    {
        m_vertices.insert(m_vertices.end(), vertex.begin(), vertex.end());
    }

    /**
     * @brief Adds a new vertex to the contour
     *
     * @param vertex A list of coordinates representing the vertex
     * @throws std::runtime_error if dim is not dim
     */
    void add_vertex(std::initializer_list<Scalar> vertex)
    {
        if (vertex.size() != dim) {
            throw std::runtime_error("Vertex size does not match the dimension of the contour");
        }
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
        return std::span<const Scalar, dim>(m_vertices.data() + vid * dim, dim);
    }

    /**
     * @brief Retrieves the number of vertices in the contour
     *
     * @return size_t The number of vertices
     */
    size_t get_num_vertices() const { return m_vertices.size() / dim; }

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
        return std::span<const Index, 2>(m_segments.data() + segid * 2, 2);
    }

    /**
     * @brief Retrieves the number of segments in the contour
     *
     * @return size_t The number of segments
     */
    size_t get_num_segments() const { return m_segments.size() / 2; }

    /**
     * @brief Adds a new cycle to the contour
     *
     * A cycle is a closed loop of segments. Each segment in the cycle can be oriented
     * either in its original direction or reversed.
     *
     * @param cycle Span of segment signed indices that form the cycle
     * @param is_regular True if the cycle is regular, false otherwise
     */
    void add_cycle(std::span<const SignedIndex> cycle, bool is_regular = true)
    {
        std::copy(cycle.begin(), cycle.end(), std::back_inserter(m_cycles));
        m_cycle_start_indices.push_back(m_cycles.size());
        assert(check_cycle(get_num_cycles() - 1));
        m_cycle_is_regular.push_back(is_regular);
    }

    /**
     * @brief Adds a new cycle to the contour
     *
     * A cycle is a closed loop of segments. Each segment in the cycle can be oriented
     * either in its original direction or reversed.
     *
     * @param cycle Span of segment signed indices that form the cycle
     * @param is_regular True if the cycle is regular, false otherwise
     */
    void add_cycle(std::initializer_list<SignedIndex> cycle, bool is_regular = true)
    {
        std::copy(cycle.begin(), cycle.end(), std::back_inserter(m_cycles));
        m_cycle_start_indices.push_back(m_cycles.size());
        assert(check_cycle(get_num_cycles() - 1));
        m_cycle_is_regular.push_back(is_regular);
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

        return std::span<const SignedIndex>(m_cycles.data() + start, end - start);
    }

    /**
     * @brief Check if a cycle is regular
     *
     * @param cid The index of the cycle to check
     *
     * @return bool True if the cycle is regular, false otherwise
     */
    bool is_cycle_regular(Index cid) const
    {
        if (cid >= m_cycle_is_regular.size()) {
            throw std::out_of_range("Cycle ID out of range");
        }
        return m_cycle_is_regular[cid];
    }

    /**
     * @brief Retrieves the regularity status of all cycles in the contour
     *
     * @return std::vector<bool> A vector indicating whether each cycle is regular
     */
    const std::vector<bool>& get_cycle_is_regular() const
    {
        return m_cycle_is_regular;
    }

    /**
     * @brief Retrieves the number of cycles in the contour
     *
     * @return size_t The number of cycles
     */
    size_t get_num_cycles() const
    {
        assert(!m_cycle_start_indices.empty());
        return m_cycle_start_indices.size() - 1;
    }

    /**
     * @brief Adds a new polyhedron to the contour
     *
     * A polyhedron is defined by a collection of oriented cycles that form its faces.
     *
     * @param polyhedron Span of cycle indices that form the polyhedron faces
     * @param is_regular True if the polyhedron is regular, false otherwise
     */
    void add_polyhedron(std::span<const SignedIndex> polyhedron, bool is_regular = true)
    {
        std::copy(polyhedron.begin(), polyhedron.end(), std::back_inserter(m_polyhedra));
        m_polyhedron_start_indices.push_back(m_polyhedra.size());
        assert(check_polyhedron(get_num_polyhedra() - 1));
        m_polyhedron_is_regular.push_back(is_regular);
    }

    /**
     * @brief Adds a new polyhedron to the contour
     *
     * A polyhedron is defined by a collection of oriented cycles that form its faces.
     *
     * @param polyhedron Span of cycle indices that form the polyhedron faces
     * @param is_regular True if the polyhedron is regular, false otherwise
     */
    void add_polyhedron(std::initializer_list<SignedIndex> polyhedron, bool is_regular = true)
    {
        std::copy(polyhedron.begin(), polyhedron.end(), std::back_inserter(m_polyhedra));
        m_polyhedron_start_indices.push_back(m_polyhedra.size());
        assert(check_polyhedron(get_num_polyhedra() - 1));
        m_polyhedron_is_regular.push_back(is_regular);
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

        return std::span<const SignedIndex>(m_polyhedra.data() + start, end - start);
    }

    /**
     * @brief Check if a polyhedron is regular
     *
     * @param poly_id The index of the polyhedron to check
     *
     * @return bool True if the polyhedron is regular, false otherwise
     */
    bool is_polyhedron_regular(Index poly_id) const
    {
        if (poly_id >= m_polyhedron_is_regular.size()) {
            throw std::out_of_range("Polyhedron ID out of range");
        }
        return m_polyhedron_is_regular[poly_id];
    }

    /**
     * @brief Retrieves the regularity status of all polyhedra in the contour
     *
     * @return std::vector<bool> A vector indicating whether each polyhedron is regular
     */
    const std::vector<bool>& get_polyhedron_is_regular() const
    {
        return m_polyhedron_is_regular;
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

    /**
     * @brief Update the contour so that all cycles are triangles.
     */
    void triangulate_cycles();

    /**
     * @brief Compute the isocontour of this contour.
     *
     * @param function_values A span of function values at the vertices of the contour.
     * @param function_gradients (optional) A span of function gradients at the vertices of the contour.
     *
     * @return A Contour object representing the isocontour.
     */
    Contour<dim> isocontour(
        std::span<Scalar> function_values,
        std::span<Scalar> function_gradients = {}) const;

private:
    void check_all_segments() const {
        const size_t num_vertices = get_num_vertices();
        for (auto vi : m_segments) {
            assert(vi >= 0 && vi < num_vertices);
        }
    }

    /**
     * @brief Check if the cycle is valid.
     *
     * A cycle is valid if each pair of consecutive segments share a vertex.
     * This ensures the cycle forms a closed loop.
     *
     * @param cid The index of the cycle to check.
     * @return True if the cycle is valid, false otherwise.
     */
    bool check_cycle(Index cid) const
    {
        auto cycle = get_cycle(cid);
        size_t cycle_size = cycle.size();

        for (size_t i = 0; i < cycle_size; i++) {
            size_t j = (i + 1) % cycle_size;

            auto seg_i = get_segment(index(cycle[i]));
            auto seg_j = get_segment(index(cycle[j]));

            auto ori_i = orientation(cycle[i]);
            auto ori_j = orientation(cycle[j]);

            Index li = ori_i ? 1 : 0;
            Index lj = ori_j ? 0 : 1;

            assert(seg_i[li] == seg_j[lj]);
            if (seg_i[li] != seg_j[lj]) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Check all cycles in the contour for validity.
     *
     * This method iterates through all cycles and calls check_cycle() on each one.
     * It is used for debugging purposes to ensure the contour's integrity.
     */
    void check_all_cycles() const
    {
        size_t num_cycles = get_num_cycles();
        assert(num_cycles == m_cycle_is_regular.size());
        for (size_t ci=0; ci < num_cycles; ci++) {
            check_cycle(ci);
        }
    }

    /**
     * @brief Check if the polyhedron is valid.
     *
     * A polyhedron is valid if the sum of signed segment indices of all cycles in the polyhedron is zero.
     * This ensures the polyhedron forms a closed surface.
     *
     * @param poly_id The index of the polyhedron to check.
     * @return True if the polyhedron is valid, false otherwise.
     */
    bool check_polyhedron(Index poly_id) const
    {
        logger().trace("### Check polyhedron {} ###", poly_id);
        int32_t sum = 0;
        auto poly = get_polyhedron(poly_id);
        for (auto ci : poly) {
            auto cycle = get_cycle(index(ci));
            bool cycle_ori = orientation(ci);
            int32_t sign = cycle_ori ? 1 : -1;
            logger().trace("Check cycle {} ({})", index(ci), sign);
            for (auto si : cycle) {
                Index seg_id = index(si);
                Index seg_ori = orientation(si);
                auto seg = get_segment(seg_id);
                Index v0 = seg[seg_ori == cycle_ori ? 0 : 1];
                Index v1 = seg[seg_ori == cycle_ori ? 1 : 0];

                logger().trace("Check segment {} ({} -> {})", value_of(si) * sign, v0, v1);
                sum += value_of(si) * sign;
            }
        }
        assert(sum == 0);
        return sum == 0;
    }

    /**
     * @brief Check all polyhedra in the contour for validity.
     *
     * This method iterates through all polyhedra and calls check_polyhedron() on each one.
     * It is used for debugging purposes to ensure the contour's integrity.
     */
    void check_all_polyhedra() const
    {
        size_t num_polyhedra = get_num_polyhedra();
        assert(num_polyhedra == m_polyhedron_is_regular.size());
        for (size_t pi=0; pi < num_polyhedra; pi++) {
            check_polyhedron(pi);
        }
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
    std::vector<bool> m_cycle_is_regular;

    /**
     * @brief Contour polyhedra.
     *
     * Each polyhedron is a list of oriented cycles that forms the faces of the polyhedron.
     */
    std::vector<SignedIndex> m_polyhedra;
    std::vector<Index> m_polyhedron_start_indices = {0};
    std::vector<bool> m_polyhedron_is_regular;
};

} // namespace mtetcol
