#include <mtetcol/disjoint_components.h>
#include <mtetcol/logger.h>

#include <cassert>
#include <functional>
#include <stdexcept>

namespace mtetcol {

void DisjointComponents::clear()
{
    for (auto cid : m_active_cycles) {
        Index cycle_id = index(cid);

        Index cycle_begin = m_cycle_indices[cycle_id];
        Index cycle_end = m_cycle_indices[cycle_id + 1];

        for (Index i = cycle_begin; i < cycle_end; i++) {
            SignedIndex si = m_cycles[i];
            Index seg_id = index(si);
            assert(seg_id < m_positive_segment_map.size());
            m_positive_segment_map[seg_id] = invalid_signed_index;
            m_negative_segment_map[seg_id] = invalid_signed_index;
        }
    }
    m_active_cycles.clear();
}

void DisjointComponents::extract_components(
    std::vector<SignedIndex>& polyhedra,
    std::vector<Index>& polyhedron_indices)
{
    for (auto cid : m_active_cycles) {
        Index cycle_id = index(cid);
        bool cycle_ori = orientation(cid);
        logger().trace("Cycle id: {}, orientation: {}", cycle_id, cycle_ori);

        Index cycle_begin = m_cycle_indices[cycle_id];
        Index cycle_end = m_cycle_indices[cycle_id + 1];
        if (cycle_end == cycle_begin + 2) {
            // 2 segment cycle (i.e. bubble), skipping
            SignedIndex si0 = m_cycles[cycle_begin];
            SignedIndex si1 = m_cycles[cycle_begin + 1];
            assert(si0 == -si1);
            continue;
        }
        for (Index i = cycle_begin; i < cycle_end; i++) {
            SignedIndex si = m_cycles[i];
            Index seg_id = index(si);
            bool seg_ori = orientation(si);
            assert(seg_id < m_positive_segment_map.size());

            if (cycle_ori == seg_ori) {
                assert(m_positive_segment_map[seg_id] == invalid_signed_index);
                m_positive_segment_map[seg_id] = signed_index(cycle_id, cycle_ori);
                logger().trace(
                    "Positive segment map: {} -> {}",
                    seg_id,
                    value_of(signed_index(cycle_id, cycle_ori)));
            } else {
                assert(m_negative_segment_map[seg_id] == invalid_signed_index);
                m_negative_segment_map[seg_id] = signed_index(cycle_id, cycle_ori);
                logger().trace(
                    "Negative segment map: {} -> {}",
                    seg_id,
                    value_of(signed_index(cycle_id, cycle_ori)));
            }
        }
    }

    size_t num_cycles = m_cycle_indices.size() - 1;
    std::vector<bool> involved(num_cycles, false);

    std::function<void(SignedIndex)> grow_component;
    grow_component = [&](SignedIndex cid) {
        Index cycle_id = index(cid);
        bool cycle_ori = orientation(cid);

        assert(!involved[cycle_id]);
        involved[cycle_id] = true;
        polyhedra.push_back(cid);

        Index cycle_begin = m_cycle_indices[cycle_id];
        Index cycle_end = m_cycle_indices[cycle_id + 1];
        if (cycle_end == cycle_begin + 2) {
            // 2 segment cycle (i.e. bubble), skipping
            return;
        }
        for (Index i = cycle_begin; i < cycle_end; i++) {
            SignedIndex si = m_cycles[i];
            Index seg_id = index(si);
            bool seg_ori = orientation(si);

            if (seg_ori == cycle_ori) {
                assert(m_positive_segment_map[seg_id] != invalid_signed_index);
                assert(index(m_positive_segment_map[seg_id]) == cycle_id);
                assert(m_positive_segment_map[seg_id] == signed_index(cycle_id, cycle_ori));

                SignedIndex adj_cycle = m_negative_segment_map[seg_id];
                if (adj_cycle == invalid_signed_index) {
                    throw std::runtime_error("Active cycle does not form closed component");
                }
                Index adj_cycle_id = index(adj_cycle);
                if (!involved[adj_cycle_id]) {
                    grow_component(adj_cycle);
                }
            } else {
                assert(m_negative_segment_map[seg_id] != invalid_signed_index);
                assert(index(m_negative_segment_map[seg_id]) == cycle_id);
                assert(m_negative_segment_map[seg_id] == signed_index(cycle_id, cycle_ori));

                SignedIndex adj_cycle = m_positive_segment_map[seg_id];
                if (adj_cycle == invalid_signed_index) {
                    throw std::runtime_error("Active cycle does not form closed component");
                }
                Index adj_cycle_id = index(adj_cycle);
                if (!involved[adj_cycle_id]) {
                    grow_component(adj_cycle);
                }
            }
        }
    };

    for (auto cid : m_active_cycles) {
        Index cycle_id = index(cid);
        bool cycle_ori = orientation(cid);
        if (involved[cycle_id]) {
            continue;
        }
        grow_component(cid);
        polyhedron_indices.push_back(static_cast<Index>(polyhedra.size()));
    }
}

void DisjointComponents::check_cycles() const
{
    if (m_cycle_indices.empty()) {
        throw std::invalid_argument("Cycle indices must not be empty");
    }
    size_t num_cycles = m_cycle_indices.size() - 1;
    size_t num_segments = m_positive_segment_map.size();

    for (auto cid : m_cycles) {
        Index id = index(cid);
        if (static_cast<size_t>(id) >= num_segments) {
            throw std::invalid_argument("Segment index out of bounds");
        }
    }
}

} // namespace mtetcol
