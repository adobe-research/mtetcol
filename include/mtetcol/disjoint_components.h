#pragma once

#include <mtetcol/common.h>

#include <SmallVector.h>

#include <span>
#include <vector>

namespace mtetcol {

class DisjointComponents
{
public:
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

    void register_cycle(SignedIndex ci) { m_active_cycles.push_back(ci); }

    void clear();

    void extract_components(
        std::vector<SignedIndex>& polyhedra,
        std::vector<Index>& polyhedron_indices);

private:
    void check_cycles() const;

private:
    std::vector<SignedIndex> m_positive_segment_map;
    std::vector<SignedIndex> m_negative_segment_map;
    std::span<const SignedIndex> m_cycles;
    std::span<const Index> m_cycle_indices;
    llvm_vecsmall::SmallVector<SignedIndex, 16> m_active_cycles;
};

} // namespace mtetcol
