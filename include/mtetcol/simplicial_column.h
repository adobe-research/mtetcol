#pragma once

#include <cassert>
#include <functional>
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
 * `SimplicialColumn` is a data structure used to represent a given spatial-temporal scalar field.
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

    void set_vertices(std::span<Scalar> vertices)
    {
        assert(vertices.size() % (dim - 1) == 0);
        m_vertices = vertices;
    }

    void set_simplices(std::span<Index> simplices);

    /**
     * @brief Sets the condensed vertex time samples and function values.
     *
     * @note Time goes from 0 to 1, and all vertices should contain time samples at 0 and 1.
     * Time samples for each vertex should be stored in ascending order.
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
        m_time_samples = std::vector<Scalar>(time_samples.begin(), time_samples.end());
        m_function_values = std::vector<Scalar>(function_values.begin(), function_values.end());
        m_vertex_start_indices =
            std::vector<Index>(vertex_start_indices.begin(), vertex_start_indices.end());
    }

    /**
     * @brief Sets the time samples and function values using a function pointer.
     *
     * @note Time goes from 0 to 1, and all vertices should contain time samples at 0 and 1.
     * Time samples for each vertex should be stored in ascending order.
     *
     * @param time_samples    A function pointer to get the time samples associated with the
     *                        vertices. `time_samples(i)` returns a span of time samples for vertex
     *                        `i` in ascending order.
     * @param function_values A function pointer to get the function values associated with the
     *                        vertices. `function_values(i)` returns a span of function values for
     *                        vertex `i` in ascending order.
     */
    void set_time_samples(
        std::function<std::span<Scalar>(Index)> time_samples,
        std::function<std::span<Scalar>(Index)> function_values)
    {
        size_t num_vertices = get_num_spatial_vertices();

        m_time_samples.clear();
        m_function_values.clear();
        m_vertex_start_indices.clear();

        m_time_samples.reserve(num_vertices * 20);
        m_function_values.reserve(num_vertices * 20);
        m_vertex_start_indices.reserve(num_vertices + 1);

        m_vertex_start_indices = {0};
        for (Index vi = 0; vi < static_cast<Index>(num_vertices); ++vi) {
            auto times = time_samples(vi);
            auto values = function_values(vi);

            assert(times.size() == values.size());
            assert(times.size() > 0);

            m_time_samples.insert(std::back_inserter(m_time_samples), times.begin(), times.end());
            m_function_values.insert(
                std::back_inserter(m_function_values),
                values.begin(),
                values.end());
            m_vertex_start_indices.push_back(m_time_samples.size());
        }
    }

    /**
     * @brief Extract the isocountour from the simplicial column.
     *
     * @param value The isovalue for the contour extraction.
     *
     * @return A Contour object representing the extracted isocontour.
     */
    Contour<dim> extract_contour(Scalar value = 0, bool cyclic = false) const;

public:
    /**
     * @brief Get the number of spatial vertices.
     *
     * @return The number of vertices.
     */
    [[nodiscard]] size_t get_num_spatial_vertices() const { return m_vertices.size() / (dim - 1); }

    /**
     * @brief Get spatial vertex coordinates.
     *
     * @return A span of spatial vertex buffer.
     */
    [[nodiscard]] std::span<const Scalar> get_spatial_vertices() const
    {
        return std::span<const Scalar>(m_vertices.data(), m_vertices.size());
    }

    /**
     * @brief Get the number of spatial edges.
     *
     * @return The number of spatial edges.
     */
    [[nodiscard]] size_t get_num_spatial_edges() const { return m_edges.size() / 2; }

    /**
     * @brief Get the spatial edge buffer.
     *
     * @return The spatial edge buffer.
     */
    [[nodiscard]] std::span<const Index> get_spatial_edges() const
    {
        return std::span<const Index>(m_edges.data(), m_edges.size());
    }

    /**
     * @brief Get the number of spatial triangles.
     *
     * @return The number of spatial triangles.
     */
    [[nodiscard]] size_t get_num_spatial_triangles() const { return m_triangles.size() / 3; }

    /**
     * @brief Get the spatial triangle buffer.
     *
     * @return The spatial triangle buffer.
     */
    [[nodiscard]] std::span<const SignedIndex> get_spatial_triangles() const
    {
        return std::span<const SignedIndex>(m_triangles.data(), m_triangles.size());
    }

    /**
     * @brief Get the number of spatial tetrahedra.
     *
     * @return The number of spatial tetrahedra.
     */
    [[nodiscard]] size_t get_num_spatial_tetrahedra() const { return m_tetrahedra.size() / 4; }

    /**
     * @brief Get the spatial tetrahedra buffer.
     *
     * @return The spatial tetrahedra buffer.
     */
    [[nodiscard]] std::span<const SignedIndex> get_spatial_tetrahedra() const
    {
        return std::span<const SignedIndex>(m_tetrahedra.data(), m_tetrahedra.size());
    }

private:
    std::span<Scalar> m_vertices;
    std::vector<Index> m_edges;
    std::vector<SignedIndex> m_triangles;
    std::vector<SignedIndex> m_tetrahedra;

    std::vector<Scalar> m_time_samples;
    std::vector<Scalar> m_function_values;
    std::vector<Index> m_vertex_start_indices;
};

} // namespace mtetcol
