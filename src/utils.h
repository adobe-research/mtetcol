#pragma once

#include <mtetcol/common.h>
#include <mtetcol/logger.h>
#include <mtetcol/simplicial_column.h>

#include <span>
#include <tuple>
#include <vector>

namespace mtetcol {

/**
 * @brief Checks if the spatial edges of a simplicial column are valid.
 *
 * The spatial edges are valid if their end vertices are in ascending order and
 * the vertex indices are within the range of the number of spatial vertices.
 *
 * @param columns The simplicial column to check.
 * @return True if the edges are valid, false otherwise.
 */
template <int dim>
bool check_edges(const SimplicialColumn<dim>& columns)
{
    const size_t num_vertices = columns.get_num_spatial_vertices();
    const size_t num_edges = columns.get_num_spatial_edges();

    auto edges = columns.get_spatial_edges();
    for (size_t i = 0; i < num_edges; i++) {
        Index v0 = edges[i * 2];
        Index v1 = edges[i * 2 + 1];
        if (v0 >= v1) return false;
        if (v0 >= num_vertices) return false;
        if (v1 >= num_vertices) return false;
    }
    return true;
}

/**
 * @brief Checks if the spatial triangles of a simplicial column are valid.
 *
 * The spatial triangles are valid if the oriented edges of each triangle forms a closed loop.
 *
 * @param columns The simplicial column to check.
 * @return True if the triangles are valid, false otherwise.
 */
template <int dim>
bool check_triangles(const SimplicialColumn<dim>& columns)
{
    const size_t num_triangles = columns.get_num_spatial_triangles();
    auto edges = columns.get_spatial_edges();
    auto triangles = columns.get_spatial_triangles();

    auto check_adj_edges = [&](SignedIndex e0, SignedIndex e1) {
        Index e0_id = index(e0);
        Index e1_id = index(e1);
        bool e0_ori = orientation(e0);
        bool e1_ori = orientation(e1);

        // Local ids
        Index e0_li = e0_ori ? 1 : 0;
        Index e1_li = e1_ori ? 0 : 1;

        return edges[e0_id * 2 + e0_li] == edges[e1_id * 2 + e1_li];
    };

    for (size_t i = 0; i < num_triangles; i++) {
        SignedIndex e01 = triangles[i * 3];
        SignedIndex e12 = triangles[i * 3 + 1];
        SignedIndex e20 = triangles[i * 3 + 2];

        if (!check_adj_edges(e01, e12)) return false;
        if (!check_adj_edges(e12, e20)) return false;
        if (!check_adj_edges(e20, e01)) return false;
    }

    return true;
}

/**
 * @brief Checks if the spatial tetrahedra of a simplicial column are valid.
 *
 * The spatial tetrahedra are valid if all edges within the tet are consistently oriented.
 *
 * @param columns The simplicial column to check.
 * @return True if the tetrahedra are valid, false otherwise.
 */
template <int dim>
bool check_tetrahedra(const SimplicialColumn<dim>& columns)
{
    const size_t num_tets = columns.get_num_spatial_tetrahedra();

    auto triangles = columns.get_spatial_triangles();
    auto tets = columns.get_spatial_tetrahedra();

    for (size_t i = 0; i < num_tets; i++) {
        SignedIndex t021 = tets[i * 4];
        SignedIndex t123 = tets[i * 4 + 1];
        SignedIndex t013 = tets[i * 4 + 2];
        SignedIndex t032 = tets[i * 4 + 3];

        Index t021_id = index(t021);
        Index t123_id = index(t123);
        Index t013_id = index(t013);
        Index t032_id = index(t032);

        int32_t t021_ori = orientation(t021) ? 1 : -1;
        int32_t t123_ori = orientation(t123) ? 1 : -1;
        int32_t t013_ori = orientation(t013) ? 1 : -1;
        int32_t t032_ori = orientation(t032) ? 1 : -1;

        int32_t sum = 0;

        for (size_t j = 0; j < 3; j++) {
            sum += value_of(triangles[t021_id * 3 + j]) * t021_ori;
            sum += value_of(triangles[t123_id * 3 + j]) * t123_ori;
            sum += value_of(triangles[t013_id * 3 + j]) * t013_ori;
            sum += value_of(triangles[t032_id * 3 + j]) * t032_ori;
        }

        assert(sum == 0);
        if (sum != 0) return false;
    }
    return true;
}

/**
 * @brief Extracts zero-crossing times from a vertex's time series data.
 *
 * This function finds all points in time where a function's value crosses a given threshold.
 * It supports both cyclic and non-cyclic time series data.
 *
 * @param time_samples A span of time points, must start at 0 and end at 1
 * @param function_values A span of function values corresponding to the time points
 * @param value The threshold value to detect crossings
 * @param cyclic Whether the time series is cyclic (wraps around)
 * @param zero_crossing_times Output vector to store the crossing times
 * @throws std::invalid_argument if input validation fails
 */
void extract_vertex_zero_crossing(
    std::span<const Scalar> time_samples,
    std::span<const Scalar> function_values,
    Scalar value,
    bool cyclic,
    std::vector<Scalar>& zero_crossing_times);

/**
 * @brief Extracts contour vertices from a set of time series data.
 *
 * This function processes multiple vertices' time series data to find contour vertices
 * where the function values cross a given threshold. It handles both cyclic and non-cyclic data.
 *
 * @param time_samples Vector of time points for all vertices
 * @param function_values Vector of function values for all vertices
 * @param vertex_start_indices Indices marking the start of each vertex's data
 * @param value The threshold value to detect crossings
 * @param cyclic Whether the time series are cyclic
 * @return Tuple containing:
 *         - Vector of zero-crossing times
 *         - Vector of indices separating each vertex's crossings
 *         - Vector of initial signs for each vertex
 * @throws std::invalid_argument if input validation fails
 */
std::tuple<std::vector<Scalar>, std::vector<size_t>, std::vector<bool>> extract_contour_vertices(
    const std::vector<Scalar>& time_samples,
    const std::vector<Scalar>& function_values,
    const std::vector<Index>& vertex_start_indices,
    Scalar value,
    bool cyclic);

/**
 * @brief Extracts the contour segments from the contour vertices.
 *
 * @param contour_times The zero-crossing times of the vertices.
 * @param contour_time_indices Indices separating the zero-crossing times of each spatial vertex.
 * @param initial_signs The initial function signs for each spatial vertex.
 * @param edges The spatial edges of the simplicial column.
 * @param cyclic Whether the time samples are cyclic.
 *
 * @return A tuple containing:
 * - A vector of contour segments. (Each segment is represented by two consecutive indices.)
 * - A vector of oriented contour segment indices over each edge relative to the edge orinetation.
 * - A vector of indices separating the contour segments over edges of each spatial edge.
 * - A vector of booleans indicating if each spatial edge is __simple__.
 * - A vector of [-1, 0, 1] indicating the cyclic shift direction (for cyclic domain and simple edge
 * only)
 */
std::tuple<
    std::vector<Index>,
    std::vector<SignedIndex>,
    std::vector<Index>,
    std::vector<bool>,
    std::vector<int8_t>>
extract_contour_segments(
    const std::vector<Scalar>& contour_times,
    const std::vector<size_t>& contour_time_indices,
    const std::vector<bool>& initial_signs,
    const std::vector<Index>& edges,
    bool cyclic);

/**
 * @brief Extracts the contour cycles from the contour segments.
 *
 * @param num_contour_vertices The number of contour vertices.
 * @param contour_segments The contour segments.
 * @param contour_segment_over_edges The oriented contour segment indices over each edge.
 * @param contour_segment_over_edges_indices Indices separating the contour segments over edges of
 * each spatial edge.
 * @param edge_is_simple A vector indicating if each spatial edge is simple.
 * @param edge_shift A vector indicating the cyclic shift direction of the time sampels for the
 * first vertex of the edge.
 * @param edges The spatial edges of the simplicial column.
 * @param triangles The spatial triangles of the simplicial column.
 *
 * @return A tuple containing:
 * - A vector of signed indices of contour segments representing the contour cycles.
 * - A vector of indices separating the individual cycles.
 * - A vector of indices separating the contour cycles of each spatial triangle.
 * - A vector of booleans indicating if each spatial triangle is simple.
 */
auto extract_contour_cycles(
    const size_t num_contour_vertices,
    const std::vector<Index>& contour_segments,
    const std::vector<SignedIndex>& contour_segment_over_edges,
    const std::vector<Index>& contour_segment_over_edges_indices,
    const std::vector<bool>& edge_is_simple,
    const std::vector<int8_t>& edge_shift,
    const std::vector<Index>& edges,
    const std::vector<SignedIndex>& triangles) -> std::
    tuple<std::vector<SignedIndex>, std::vector<Index>, std::vector<Index>, std::vector<bool>>;

/**
 * @brief Extracts the contour polyhedra from the contour cycles.
 *
 * @param num_contour_segments The number of contour segments.
 * @param contour_cycles The contour cycles.
 * @param contour_cycle_indices Indices separating the individual cycles.
 * @param contour_cycle_triangle_indices Indices separating the contour cycles of each spatial
 * triangle.
 * @param triangle_is_simple A vector indicating if each spatial triangle is simple.
 * @param tetrahedra The spatial tetrahedra of the simplicial column.
 *
 * @return A tuple containing:
 * - A vector of signed indices of contour cycles representing the contour polyhedra.
 * - A vector of indices separating the individual polyhedra.
 * - A vector of indices separating the contour polyhedra of each spatial tetrahedron.
 */
std::tuple<std::vector<SignedIndex>, std::vector<Index>, std::vector<Index>>
extract_contour_polyhedra(
    size_t num_contour_segments,
    const std::vector<SignedIndex>& contour_cycles,
    const std::vector<Index>& contour_cycle_indices,
    const std::vector<Index>& contour_cycle_triangle_indices,
    const std::vector<bool>& triangle_is_simple,
    const std::vector<SignedIndex>& tetrahedra);

} // namespace mtetcol
