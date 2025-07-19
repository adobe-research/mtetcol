#pragma once

#include <mtetcol/common.h>

#include <SmallVector.h>

#include <span>
#include <stdexcept>

namespace mtetcol {

class NonDisjointCycles
{
public:
    /**
     * @brief Constructs a NonDisjointCycles object with the given segments.
     *
     * @param segments A span of vertex indices representing the segments. Each segment
     *                 is defined by a pair of vertex indices (start and end points).
     *
     * @throws std::invalid_argument if the segments are not valid (e.g., not an even number of
     * indices).
     */
    NonDisjointCycles(std::span<const Index> segments)
        : m_segments(segments)
    {
        if (m_segments.size() % 2 != 0) {
            throw std::invalid_argument("Segments must be a span of vertex index pairs");
        }
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
     * The method uses a depth-first search to find cycles in the graph of segments. Two cycles may
     * share vertices/segments.
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
    std::span<const Index> m_segments; ///< The input segments as vertex index pairs
    llvm_vecsmall::SmallVector<SignedIndex, 16>
        m_active_segments; ///< Currently active segments for cycle extraction
};

} // namespace mtetcol
