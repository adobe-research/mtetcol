#pragma once

#include <span>
#include <vector>

#include <mtetcol/common.h>
#include <mtetcol/contour.h>

namespace mtetcol {

/**
 * @brief Simplicial column data structure for representing space-time domain.
 *
 * @tparam dim The dimension of the ambient space. Default is 4 (3D space + 1D time).
 *
 * `SimplicialColumn` is a data structure used to represent a given spacial-temporal scalar field.
 * It consists of a collection of simplices (triangles or tetrahedra) that discretize the space and
 * a sequence of ordered time/value samples associated with each vertex. It supports the extraction
 * of isocontours, where an isocontour is defined as a polyhedral mesh embedded in 4D space or a
 * polygonal mesh embedded in 3D space.
 */
template <int dim = 4>
class SimplicialColumn
{
public:
    static_assert(dim == 3 || dim == 4, "Ambient dimension of simplicial columns must be 3 or 4");

    void set_vertices(std::span<Scalar> vertices) { m_vertices = vertices; }

    void set_simplices(std::span<Index> simplices);

    /**
     * @brief Sets the condensed vertex time samples and function values.
     *
     * @param time_samples    A span of time samples associated with the vertices.
     * @param function_values A span of function values associated with the vertices.
     * @param vertex_start_indices A span of indices indicating the start of each vertex's time
     * samples.
     *
     * For example, to access the time samples and function values associated with vertex `i`:
     * ```
     *   Index vi_num_samples = vertex_start_indices[i + 1] - vertex_start_indices[i];
     *   span<Scalar> vi_times = time_samples.subspan(vertex_start_indices[i], vi_num_samples);
     *   span<Scalar> vi_values = function_values.subspan(vertex_start_indices[i], vi_num_samples);
     * ```
     */
    void set_time_samples(
        std::span<Scalar> time_samples,
        std::span<Scalar> function_values,
        std::span<Index> vertex_start_indices)
    {
        m_time_samples = time_samples;
        m_function_values = function_values;
        m_vertex_start_indices = vertex_start_indices;
    }

    /**
     * @brief Extract the isocountour from the simplicial column.
     *
     * @param value The isovalue for the contour extraction.
     *
     * @return A Contour object representing the extracted isocontour.
     */
    Contour<dim> extract_contour(Scalar value = 0) const;

private:
    std::span<Scalar> m_vertices;
    std::vector<Index> m_edges;
    std::vector<SignedIndex> m_triangles;
    std::vector<SignedIndex> m_tetrahedra;

    std::span<Scalar> m_time_samples;
    std::span<Scalar> m_function_values;
    std::span<Index> m_vertex_start_indices;
};

} // namespace mtetcol
