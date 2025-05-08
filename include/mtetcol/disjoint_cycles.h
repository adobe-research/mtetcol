#pragma once

#include <mtetcol/common.h>

#include <SmallVector.h>

#include <span>
#include <vector>

namespace mtetcol {

/**
 * @class DisjointCycles
 * @brief A class for managing and extracting disjoint cycles from a set of segments.
 *
 * This class provides functionality to:
 * - Register and unregister segments
 * - Extract cycles by chaining oriented segments
 * - Validate segment data
 *
 * The class maintains a graph-like structure where segments are edges between vertices,
 * and cycles are closed paths in this graph.
 */
class DisjointCycles
{
public:
    /**
     * @brief Constructs a DisjointCycles object with the given number of vertices and segments.
     *
     * @param num_vertices The number of vertices in the contour. This determines the size of
     *                     the internal vertex mapping array.
     * @param segments A Nx2 span of vertex indices representing the segments. Each segment
     *                 is defined by two vertex indices (start and end points).
     * @throws std::invalid_argument if the segments are not valid (e.g., invalid vertex indices)
     */
    DisjointCycles(size_t num_vertices, std::span<const Index> segments)
        : m_next_index(num_vertices, invalid_signed_index)
        , m_segments(segments)
    {
        check_segments();
    }

    /**
     * @brief Registers a segment to be included in the cycle extraction.
     *
     * @param si The signed index of the segment to register. A positive index indicates
     *           the segment's original orientation, while a negative index indicates
     *           the reverse orientation.
     */
    void register_segment(SignedIndex si) { m_active_segments.push_back(si); }

    /**
     * @brief Unregisters all segments from the cycle extraction.
     *
     * This clears all previously registered segments, allowing for a fresh start
     * in cycle extraction.
     */
    void clear();

    /**
     * @brief Extract cycles by chaining the registered oriented segments.
     *
     * This method processes all registered segments and identifies closed cycles.
     * Each cycle is represented as a sequence of signed segment indices.
     *
     * @param[out] cycles The vector to store the signed indices of segments that form the cycles.
     *                    Each cycle is stored as a contiguous sequence of segment indices.
     * @param[out] cycle_indices The vector to store the separating indices of the cycles in the
     *             cycles vector. Each element represents the start index of a new cycle in the
     *             cycles vector. The last element is the total number of segments in all cycles.
     */
    void extract_cycles(std::vector<SignedIndex>& cycles, std::vector<Index>& cycle_indices);

    /**
     * @brief Returns the number of registered segments.
     *
     * @return The number of currently registered segments.
     */
    size_t num_registered_segments() const { return m_active_segments.size(); }

private:
    /**
     * @brief Checks if the segments are valid.
     *
     * This method validates that:
     * - All vertex indices are within bounds
     * - Segments are properly formed
     *
     * @throws std::invalid_argument if the segments are not valid.
     */
    void check_segments() const;

private:
    std::vector<SignedIndex> m_next_index; ///< Maps each vertex to its next vertex in a cycle
    std::span<const Index> m_segments; ///< The input segments as vertex index pairs
    llvm_vecsmall::SmallVector<SignedIndex, 16>
        m_active_segments; ///< Currently active segments for cycle extraction
};

} // namespace mtetcol
