#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <exception>
#include <functional>
#include <iterator>
#include <numeric>
#include <span>
#include <type_traits>
#include <vector>

#include <mtetcol/common.h>
#include <mtetcol/logger.h>

#include <SmallVector.h>

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
        if (vid * dim + dim - 1 >= m_vertices.size()) {
            throw std::out_of_range("Vertex ID out of range");
        }
        return std::span<const Scalar, dim>(m_vertices.data() + vid * dim, dim);
    }

    /**
     * @brief Retrieves a mutable reference to the coordinates of a vertex
     *
     * @param vid The index of the vertex to retrieve
     *
     * @return std::span<Scalar, dim> A span containing the coordinates of the vertex (size depends
     * on dim)
     */
    std::span<Scalar, dim> ref_vertex(Index vid)
    {
        if (vid * dim + dim - 1 >= m_vertices.size()) {
            throw std::out_of_range("Vertex ID out of range");
        }
        return std::span<Scalar, dim>(m_vertices.data() + vid * dim, dim);
    }

    /**
     * @brief Resizes the vertex list to accommodate a specified number of vertices
     *
     * @param num_vertices The number of vertices to allocate space for
     */
    void resize_vertices(size_t num_vertices) { m_vertices.resize(num_vertices * dim); }

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
        if (segid * 2 + 1 >= m_segments.size()) {
            throw std::out_of_range("Segment ID out of range");
        }
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
        assert(cycle_is_valid(get_num_cycles() - 1));
        assert(cycle_is_simple(get_num_cycles() - 1));
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
        assert(cycle_is_valid(get_num_cycles() - 1));
        assert(cycle_is_simple(get_num_cycles() - 1));
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
    const std::vector<bool>& get_cycle_is_regular() const { return m_cycle_is_regular; }

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
    const std::vector<bool>& get_polyhedron_is_regular() const { return m_polyhedron_is_regular; }

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
     *
     * @param optimal_triangulation (optional) A boolean flag indicating whether to use optimal
     * triangulation. (Only apply to quad faces for now)
     */
    void triangulate_cycles(bool optimal_triangulation = true);

    /**
     * @brief Compute the isocontour of this contour.
     *
     * @param function_values A span of function values at the vertices of the contour.
     * @param function_gradients (optional) A span of function gradients at the vertices of the
     * contour.
     * @param use_snapping (optional) A boolean flag indicating whether to use snapping.
     *
     * @return A Contour object representing the isocontour.
     */
    Contour<dim> isocontour(
        std::span<Scalar> function_values,
        std::span<Scalar> function_gradients = {},
        bool use_snapping = false) const;

public:
    /**
     * @brief Check if the contour is manifold.
     *
     * @note This check is designed for 4D contours composed of polyhedra.
     *
     * A contour is manifold if no cycle is shared by more than two polyhedra.
     *
     * @return bool True if the contour is manifold, false otherwise.
     */
    bool is_manifold() const
    {
        size_t num_cycles = get_num_cycles();
        size_t num_polyhedra = get_num_polyhedra();

        if (num_polyhedra > 0) {
            std::vector<size_t> cycle_usage_count(num_cycles, 0);
            for (size_t pi = 0; pi < num_polyhedra; pi++) {
                auto poly = get_polyhedron(pi);
                for (auto ci : poly) {
                    Index cycle_id = index(ci);
                    cycle_usage_count[cycle_id]++;
                }
            }
            for (auto count : cycle_usage_count) {
                if (count > 2) {
                    return false;
                }
            }
        } else if (num_cycles > 0) {
            size_t num_segments = get_num_segments();
            std::vector<size_t> segment_usage_count(num_segments, 0);
            for (size_t ci = 0; ci < num_cycles; ci++) {
                auto cycle = get_cycle(ci);
                for (auto si : cycle) {
                    Index seg_id = index(si);
                    segment_usage_count[seg_id]++;
                }
            }
            for (auto count : segment_usage_count) {
                if (count > 2) {
                    return false;
                }
            }
        }

        return true;
    }


private:
    void check_all_segments() const
    {
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
    bool cycle_is_valid(Index cid) const
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
     * @brief Check if the cycle is simple.
     *
     * A cycle is simple if it does not self-interesect. I.e. there should not be any duplicate
     * vertices in the cycle.
     *
     * @param cid The index of the cycle to check.
     * @return True if the cycle is simple, false otherwise.
     */
    bool cycle_is_simple(Index cid) const;

    /**
     * @brief Check all cycles in the contour for validity.
     *
     * This method iterates through all cycles and ensures cycle is valid and simple.
     * It is used for debugging purposes to ensure the contour's integrity.
     */
    void check_all_cycles() const
    {
        size_t num_cycles = get_num_cycles();
        assert(num_cycles == m_cycle_is_regular.size());
        for (size_t ci = 0; ci < num_cycles; ci++) {
            assert(cycle_is_valid(ci));
            assert(cycle_is_simple(ci));
        }
    }

    /**
     * @brief Check if the polyhedron is valid.
     *
     * A polyhedron is valid if the sum of signed segment indices of all cycles in the polyhedron is
     * zero. This ensures the polyhedron forms a closed surface.
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
        for (size_t pi = 0; pi < num_polyhedra; pi++) {
            check_polyhedron(pi);
        }
        assert(is_manifold());
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

/**
 * @brief Extract simple loops from a cycle that may contain duplicate vertices.
 *
 * This function takes a cycle (represented as a sequence of signed segment indices) and
 * splits it into multiple simple loops if the cycle contains duplicate vertices (junction points).
 * Each junction vertex (appearing more than once in the cycle) acts as a splitting point.
 *
 * @param cycle A span of signed indices representing the cycle segments
 * @param get_segment A callback function that takes a segment index and returns the two vertex
 *                    indices of that segment as std::array<Index, 2>
 * @return A SmallVector of SmallVectors, where each inner SmallVector represents a simple loop
 *         (subcycle) extracted from the input cycle. If the input cycle has no duplicate vertices,
 *         returns a single-element vector containing the original cycle.
 *
 * @note The extraction algorithm prioritizes starting from junction vertices (valence > 1) to
 *       ensure proper splitting of cycles at duplicate vertices.
 */
llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<SignedIndex, 16>, 4> extract_simple_loops(
    std::span<const SignedIndex> cycle,
    std::function<std::array<Index, 2>(Index)> get_segment);

} // namespace mtetcol
