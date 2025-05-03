#include "disjoint_cycles.h"

#include <cassert>

namespace mtetcol {

bool DisjointCycles::check_segments() const
{
    size_t num_vertices = m_next_index.size();
    for (auto vi : m_segments) {
        assert(static_cast<size_t>(vi) < num_vertices);
        if (static_cast<size_t>(vi) >= num_vertices) {
            return false;
        }
    }
    return true;
}

void DisjointCycles::clear()
{
    assert(m_segments.size() % 2 == 0);
    const size_t num_segments = m_segments.size() / 2;
    for (auto si : m_active_segments) {
        Index seg_id = index(si);
        assert(seg_id < m_segments.size() / 2);
        Index v0 = m_segments[seg_id * 2];
        Index v1 = m_segments[seg_id * 2 + 1];
        m_next_index[v0] = invalid_signed_index;
        m_next_index[v1] = invalid_signed_index;
    }

    m_active_segments.clear();
}

void DisjointCycles::extract_cycles(std::vector<SignedIndex>& cycles, std::vector<Index>& cycle_indices)
{
    if (cycle_indices.empty() || cycle_indices.back() != cycles.size()) {
        throw std::invalid_argument("Cycle indices must be non-empty end with the size of cycles");
    }

    assert(m_segments.size() % 2 == 0);
    const size_t num_segments = m_segments.size() / 2;

    // Initialize vertex to next segment mapping
    for (auto si : m_active_segments) {
        Index seg_id = index(si);
        bool seg_ori = orientation(si);
        assert(seg_id < num_segments);
        Index v0 = m_segments[seg_id * 2];
        Index v1 = m_segments[seg_id * 2 + 1];
        if (seg_ori) {
            assert(m_next_index[v0] == invalid_signed_index);
            m_next_index[v0] = si;
        } else {
            assert(m_next_index[v1] == invalid_signed_index);
            m_next_index[v1] = si;
        }
    }

    // Grow the current cycle by one segment
    auto grow_cycle = [&](Index vid) -> Index {
        assert(m_next_index[vid] != invalid_signed_index);
        SignedIndex si = m_next_index[vid];
        m_next_index[vid] = invalid_signed_index;

        Index seg_id = index(si);
        bool seg_ori = orientation(si);
        assert(seg_id < num_segments);

        Index v0 = m_segments[seg_id * 2];
        Index v1 = m_segments[seg_id * 2 + 1];
        assert(seg_ori ? v0 == vid : v1 == vid);

        cycles.push_back(si);

        return seg_ori ? v1 : v0;
    };

    for (auto si : m_active_segments) {
        Index seg_id = index(si);
        assert(seg_id < num_segments);
        bool seg_ori = orientation(si);

        Index v0 = m_segments[seg_id * 2];
        Index v1 = m_segments[seg_id * 2 + 1];
        Index vid = seg_ori ? v0 : v1;

        if (m_next_index[vid] == invalid_signed_index) {
            continue;
        }

        while (m_next_index[vid] != invalid_signed_index) {
            vid = grow_cycle(vid);
        }

        cycle_indices.push_back(static_cast<Index>(cycles.size()));
    }
}

} // namespace mtetcol
