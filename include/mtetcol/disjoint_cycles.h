#pragma once

#include <mtetcol/common.h>

#include <SmallVector.h>

#include <span>
#include <vector>

namespace mtetcol {

class DisjointCycles
{
public:
    /**
     * @brief Constructs a DisjointCycles object with the given number of vertices and segments.
     *
     * @param num_vertices The number of vertices in the contour.
     * @param segments A Nx2 span of vertex indices representing the segments.
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
     * @param si The signed index of the segment to register.
     */
    void register_segment(SignedIndex si) { m_active_segments.push_back(si); }

    /**
     * @brief Unregisters all segments from the cycle extraction.
     */
    void clear();

    /**
     * @brief Extract cycles by chaining the registered oriented segments.
     *
     * @param[out] cycles The vector to store the signed index of segments that form the cycles.
     * @param[out] cycle_indices The vector to store the separating indices of the cycles in the
     *             cycles vector.
     */
    void extract_cycles(std::vector<SignedIndex>& cycles, std::vector<Index>& cycle_indices);

private:
    /**
     * @brief Checks if the segments are valid.
     *
     * @throws std::invalid_argument if the segments are not valid.
     */
    void check_segments() const;

private:
    std::vector<SignedIndex> m_next_index;
    std::span<const Index> m_segments;
    llvm_vecsmall::SmallVector<SignedIndex, 16> m_active_segments;
};

} // namespace mtetcol
