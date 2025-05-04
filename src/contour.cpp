#include <mtetcol/contour.h>

namespace mtetcol {

namespace {

std::vector<Index> triangulate(
    std::vector<Index>& segments,
    std::vector<SignedIndex>& cycles,
    std::vector<Index>& cycle_indices)
{
    size_t num_cycles = cycle_indices.size() - 1;

    std::vector<SignedIndex> triangle_cycles;
    std::vector<Index> triangle_cycle_indices;
    std::vector<Index> cycle_to_triangle_map;

    triangle_cycles.reserve(cycles.size());
    triangle_cycle_indices.reserve(num_cycles + 1);
    triangle_cycle_indices.push_back(0);
    cycle_to_triangle_map.reserve(num_cycles + 1);
    cycle_to_triangle_map.push_back(0);

    auto get_segment = [&](SignedIndex si) -> std::array<Index, 2> {
        Index seg_id = index(si);
        bool seg_ori = orientation(si);

        Index v0 = segments[seg_id * 2 + (seg_ori ? 0 : 1)];
        Index v1 = segments[seg_id * 2 + (seg_ori ? 1 : 0)];
        return {v0, v1};
    };

    for (size_t i = 0; i < num_cycles; i++) {
        auto start = cycle_indices[i];
        auto end = cycle_indices[i + 1];
        auto cycle = std::span<const SignedIndex>(cycles.data() + start, end - start);

        if (cycle.size() < 3) {
            // Cycle of size smaller than 3 is dropped.
        } else if (cycle.size() == 3) {
            // Check if the cycle is a triangle
            std::copy(cycle.begin(), cycle.end(), std::back_inserter(triangle_cycles));
            triangle_cycle_indices.push_back(triangle_cycles.size());
        } else {
            // Triangulate the cycle
            assert(cycle.size() > 3);
            size_t num_segments = end - start;
            auto first_segment = get_segment(cycle[0]);
            Index v0 = first_segment[0];

            SignedIndex si_curr = invalid_signed_index;
            SignedIndex si_prev = invalid_signed_index;
            SignedIndex si_next = invalid_signed_index;

            for (size_t j = 1; j + 1 < num_segments; j++) {
                auto si = cycle[j];
                auto target_segment = get_segment(si);
                auto segment = get_segment(cycle[j]);

                si_curr = si;

                if (j == 1) {
                    segments.push_back(v0);
                    segments.push_back(segment[1]);

                    si_prev = cycle[0];
                    si_next = signed_index(segments.size() / 2 - 1, false);

                } else if (j + 2 == num_segments) {
                    si_next = cycle.back();

                } else {
                    segments.push_back(v0);
                    segments.push_back(segment[1]);

                    si_next = signed_index(segments.size() / 2 - 1, false);
                }
                assert(si_prev != invalid_signed_index);
                assert(si_next != invalid_signed_index);
                triangle_cycles.push_back(si_prev);
                triangle_cycles.push_back(si_curr);
                triangle_cycles.push_back(si_next);
                triangle_cycle_indices.push_back(triangle_cycles.size());

                si_prev = -si_next;
            }
        }
        cycle_to_triangle_map.push_back(triangle_cycle_indices.size() - 1);
    }

    std::swap(cycles, triangle_cycles);
    std::swap(cycle_indices, triangle_cycle_indices);

    return cycle_to_triangle_map;
}

} // namespace

template <>
void Contour<3>::triangulate_cycles()
{
    [[maybe_unused]] auto cycle_to_triangle_map =
        triangulate(m_segments, m_cycles, m_cycle_start_indices);

    size_t num_cycles = get_num_cycles();
    for (size_t i = 0; i < num_cycles; i++) {
        check_cycle(i);
    }
}

template <>
void Contour<4>::triangulate_cycles()
{
    auto cycle_to_triangle_map = triangulate(m_segments, m_cycles, m_cycle_start_indices);

    std::vector<SignedIndex> updated_polyhedra;
    std::vector<Index> updated_polyhedron_start_indices;
    updated_polyhedra.reserve(m_polyhedra.size());
    updated_polyhedron_start_indices.reserve(m_polyhedron_start_indices.size());
    updated_polyhedron_start_indices.push_back(0);

    size_t num_polyhedra = get_num_polyhedra();

    for (size_t i = 0; i < num_polyhedra; i++) {
        Index cycle_start = m_polyhedron_start_indices[i];
        Index cycle_end = m_polyhedron_start_indices[i + 1];
        for (size_t j = cycle_start; j < cycle_end; j++) {
            SignedIndex ci = m_polyhedra[j];
            Index cycle_id = index(ci);
            bool cycle_ori = orientation(ci);

            Index new_triangle_start = cycle_to_triangle_map[cycle_id];
            Index new_triangle_end = cycle_to_triangle_map[cycle_id + 1];

            for (Index k = new_triangle_start; k < new_triangle_end; k++) {
                updated_polyhedra.push_back(signed_index(k, cycle_ori));
            }
        }
        updated_polyhedron_start_indices.push_back(updated_polyhedra.size());
    }

    std::swap(m_polyhedra, updated_polyhedra);
    std::swap(m_polyhedron_start_indices, updated_polyhedron_start_indices);

    num_polyhedra = get_num_polyhedra();
    for (size_t i = 0; i < num_polyhedra; i++) {
        check_polyhedron(i);
    }
}

template <>
Contour<3> Contour<3>::isocontour(std::span<Scalar> function_values) const
{
    // TODO
    return Contour<3>();
}

template <>
Contour<4> Contour<4>::isocontour(std::span<Scalar> function_values) const
{
    // TODO
    return Contour<4>();
}

} // namespace mtetcol
