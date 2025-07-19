#include <mtetcol/logger.h>
#include <mtetcol/nondisjoint_cycles.h>

#include "hashmap.h"

#include <ankerl/unordered_dense.h>

#include <cassert>
#include <stdexcept>

namespace mtetcol {

void NonDisjointCycles::clear()
{
    m_active_segments.clear();
}

void NonDisjointCycles::extract_cycles(
    std::vector<SignedIndex>& cycles,
    std::vector<Index>& cycle_indices)
{
    // Validate input parameters
    if (cycle_indices.empty() || cycle_indices.back() != cycles.size()) {
        throw std::invalid_argument("Cycle indices must be non-empty and end with the size of cycles");
    }

    const size_t num_segments = m_segments.size() / 2;

    using IndexSet = ankerl::unordered_dense::set<Index>;
    using SignedIndexSet = ankerl::unordered_dense::set<SignedIndex, SignedIndexHash>;

    ankerl::unordered_dense::map<Index, SignedIndexSet> vertex_to_segments;
    vertex_to_segments.reserve(m_active_segments.size());

    // Initialize vertex to next segment mapping
    for (auto si : m_active_segments) {
        Index seg_id = index(si);
        bool seg_ori = orientation(si);
        assert(seg_id < num_segments);
        Index v0 = m_segments[seg_id * 2];
        Index v1 = m_segments[seg_id * 2 + 1];
        if (!seg_ori) {
            std::swap(v0, v1);
        }

        if (!vertex_to_segments.contains(v0)) {
            vertex_to_segments.insert({v0, {}});
        }
        vertex_to_segments[v0].insert(si);
    }

    std::vector<SignedIndex> segment_queue;
    IndexSet visited_vertices;

    auto extract_cycle = [&](Index vid) {
        size_t queue_size = segment_queue.size();
        assert(queue_size >= 2);
        bool cycle_found = false;
        for (size_t i=0; i<queue_size; i++) {
            SignedIndex si = segment_queue[i];
            Index seg_id = index(si);
            bool seg_ori = orientation(si);
            Index v0 = m_segments[seg_id * 2];
            Index v1 = m_segments[seg_id * 2 + 1];
            if (!seg_ori) {
                std::swap(v0, v1);
            }

            if (v0 == vid) {
                cycle_found = true;
            }
            if (cycle_found) {
                cycles.push_back(si);
            }
        }
    };

    std::function<void(Index)> dfs;
    dfs = [&](Index vid) {
        assert(!visited_vertices.contains(vid));
        visited_vertices.insert(vid);

        auto& segments = vertex_to_segments[vid];
        for (auto si : segments) {
            Index seg_id = index(si);
            bool seg_ori = orientation(si);
            Index v0 = m_segments[seg_id * 2];
            Index v1 = m_segments[seg_id * 2 + 1];
            Index next_vid = seg_ori ? v1 : v0;

            if (!visited_vertices.contains(next_vid)) {
                segment_queue.push_back(si);
                dfs(next_vid);
                segment_queue.pop_back();
            } else {
                // Cycle detected
                segment_queue.push_back(si);
                extract_cycle(next_vid);
                segment_queue.pop_back();
                cycle_indices.push_back(static_cast<Index>(cycles.size()));
            }
        }
    };

    for (auto si : m_active_segments) {
        Index seg_id = index(si);
        bool seg_ori = orientation(si);
        assert(seg_id < num_segments);
        Index v0 = m_segments[seg_id * 2];
        Index v1 = m_segments[seg_id * 2 + 1];
        if (!seg_ori) {
            std::swap(v0, v1);
        }

        if (visited_vertices.contains(v0)) {
            continue;
        }

        dfs(v0);
    }
}

} // namespace mtetcol
