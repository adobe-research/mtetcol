#pragma once

#include <mtetcol/contour.h>
#include <span>
#include <string_view>

namespace mtetcol {

/**
 * @brief Saves a 3D contour to a file.
 *
 * @param filename The path to the output file where the contour will be saved
 * @param contour The 3D contour to be saved
 * @param function_values Optional span of scalar function values associated with the contour
 */
void save_contour(
    std::string_view filename,
    const Contour<3>& contour,
    std::span<Scalar> function_values = {});

/**
 * @brief Saves a 4D contour to a file, optionally with function values.
 *
 * @param filename The path to the output file where the contour will be saved
 * @param contour The 4D contour to be saved
 * @param function_values Optional span of scalar function values associated with the contour
 */
void save_contour(
    std::string_view filename,
    const Contour<4>& contour,
    std::span<Scalar> function_values = {});


/**
 * @brief Saves a polyhedron from a contour to a file (mostly for debugging purposes).
 *
 * @param filename The path to the output file where the polyhedron will be saved
 * @param contour The 4D contour containing the polyhedron
 * @param polyhedron_id The ID of the polyhedron to be saved
 */
void save_polyhedron(std::string_view filename, const Contour<4>& contour, Index polyhedron_id);

} // namespace mtetcol
