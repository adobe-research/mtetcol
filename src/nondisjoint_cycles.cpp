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
        throw std::invalid_argument(
            "Cycle indices must be non-empty and end with the size of cycles");
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

    IndexSet used_segments; // Track segments that have been used in cycles

    // Extract cycles iteratively to ensure disjoint segments
    // Keep extracting until no more cycles can be found
    while (true) {
        std::vector<SignedIndex> segment_queue;
        IndexSet visited_vertices;
        bool cycle_extracted = false;

        auto extract_cycle = [&](Index vid) {
            size_t queue_size = segment_queue.size();
            assert(queue_size >= 2);
            bool reached_cycle_start = false;
            for (size_t i = 0; i < queue_size; i++) {
                SignedIndex si = segment_queue[i];
                Index seg_id = index(si);
                bool seg_ori = orientation(si);
                Index v0 = m_segments[seg_id * 2];
                Index v1 = m_segments[seg_id * 2 + 1];
                if (!seg_ori) {
                    std::swap(v0, v1);
                }

                if (v0 == vid) {
                    reached_cycle_start = true;
                }
                if (reached_cycle_start) {
                    cycles.push_back(si);
                    used_segments.insert(seg_id); // Mark segment as used
                }
            }
            cycle_extracted = true;
        };

        std::function<void(Index)> dfs;
        dfs = [&](Index vid) {
            if (visited_vertices.contains(vid)) {
                return;
            }
            visited_vertices.insert(vid);

            if (!vertex_to_segments.contains(vid)) {
                return;
            }

            auto& segments = vertex_to_segments[vid];
            for (auto si : segments) {
                Index seg_id = index(si);

                // Skip if this segment has already been used in another cycle
                if (used_segments.contains(seg_id)) {
                    continue;
                }

                bool seg_ori = orientation(si);
                Index v0 = m_segments[seg_id * 2];
                Index v1 = m_segments[seg_id * 2 + 1];
                Index next_vid = seg_ori ? v1 : v0;

                if (!visited_vertices.contains(next_vid)) {
                    segment_queue.push_back(si);
                    dfs(next_vid);
                    segment_queue.pop_back();
                } else if (!cycle_extracted) {
                    // Cycle detected - extract it and stop this DFS iteration
                    segment_queue.push_back(si);
                    extract_cycle(next_vid);
                    segment_queue.pop_back();
                    cycle_indices.push_back(static_cast<Index>(cycles.size()));
                    return; // Stop after finding one cycle
                }
            }
        };

        // Try to find one cycle starting from any unused segment
        bool started_dfs = false;
        for (auto si : m_active_segments) {
            Index seg_id = index(si);

            // Skip if this segment has already been used
            if (used_segments.contains(seg_id)) {
                continue;
            }

            bool seg_ori = orientation(si);
            Index v0 = m_segments[seg_id * 2];
            Index v1 = m_segments[seg_id * 2 + 1];
            if (!seg_ori) {
                std::swap(v0, v1);
            }

            // Start DFS from this vertex
            dfs(v0);
            started_dfs = true;

            // If we extracted a cycle, restart the search
            if (cycle_extracted) {
                break;
            }
        }

        // If no cycle was extracted or no DFS was started, we're done
        if (!cycle_extracted || !started_dfs) {
            break;
        }
    }
}

} // namespace mtetcol
