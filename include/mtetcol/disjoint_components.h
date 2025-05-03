#pragma once

#include <mtetcol/common.h>

#include <SmallVector.h>

#include <span>
#include <vector>

namespace mtetcol {

/**
 * @brief A class for managing and extracting disjoint components from cycles in a mesh.
 * 
 * This class maintains mappings between segments and cycles, allowing for the extraction
 * of polyhedral components from a set of cycles. It supports both positive and negative
 * segment mappings and can track active cycles for component extraction.
 */
class DisjointComponents
{
public:
    /**
     * @brief Constructs a DisjointComponents object with the given number of segments and cycle data.
     * 
     * @param num_segments The total number of segments in the mesh
     * @param cycles A span containing the cycle data
     * @param cycle_indices A span containing the indices into the cycles data
     */
    DisjointComponents(
        size_t num_segments,
        std::span<const SignedIndex> cycles,
        std::span<const Index> cycle_indices)
        : m_positive_segment_map(num_segments, invalid_signed_index)
        , m_negative_segment_map(num_segments, invalid_signed_index)
        , m_cycles(cycles)
        , m_cycle_indices(cycle_indices)
    {
        check_cycles();
    }

    /**
     * @brief Registers a cycle as active for component extraction.
     * 
     * @param ci The signed index of the cycle to register
     */
    void register_cycle(SignedIndex ci) { m_active_cycles.push_back(ci); }

    /**
     * @brief Clears all active cycles and resets the component state.
     */
    void clear();

    /**
     * @brief Extracts polyhedral components from the active cycles.
     * 
     * @param polyhedra Output vector to store the extracted polyhedra
     * @param polyhedron_indices Output vector to store the indices of the polyhedra
     */
    void extract_components(
        std::vector<SignedIndex>& polyhedra,
        std::vector<Index>& polyhedron_indices);

private:
    /**
     * @brief Validates the cycle data structure.
     * 
     * @throws std::runtime_error if the cycle data is invalid
     */
    void check_cycles() const;

private:
    std::vector<SignedIndex> m_positive_segment_map; ///< Maps segments to positive cycles
    std::vector<SignedIndex> m_negative_segment_map; ///< Maps segments to negative cycles
    std::span<const SignedIndex> m_cycles; ///< The complete set of cycles
    std::span<const Index> m_cycle_indices; ///< Indices into the cycles data
    llvm_vecsmall::SmallVector<SignedIndex, 16> m_active_cycles; ///< Currently active cycles for component extraction
};

} // namespace mtetcol
