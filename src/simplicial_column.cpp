#include <mtetcol/simplicial_column.h>

namespace mtetcol {

template <>
void SimplicialColumn<4>::set_simplices(std::span<Index> simplices)
{
    // TODO
}

template <>
void SimplicialColumn<3>::set_simplices(std::span<Index> simplices)
{
    // TODO
}

template <>
Contour<4> SimplicialColumn<4>::extract_contour(Scalar value) const
{
    // TODO
    return Contour<4>();
}

template <>
Contour<3> SimplicialColumn<3>::extract_contour(Scalar value) const
{
    // TODO
    return Contour<3>();
}

} // namespace mtetcol
