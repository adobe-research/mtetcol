#pragma once

#include <mtetcol/contour.h>

#include <string_view>

namespace mtetcol {

void save_contour(std::string_view filename, const Contour<3>& contour);
void save_contour(
    std::string_view filename,
    const Contour<4>& contour,
    std::span<Scalar> function_values = {});

} // namespace mtetcol
