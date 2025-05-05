#pragma once

#include <mtetcol/common.h>

#include <array>

namespace mtetcol {

template <int dim>
class ImplicitFunction
{
public:
    virtual Scalar value(std::array<Scalar, dim> pos) const = 0;
    virtual std::array<Scalar, dim> gradient(std::array<Scalar, dim> pos) const = 0;
};

class ImplicitSphere : public ImplicitFunction<3>
{
public:
    ImplicitSphere(Scalar radius, std::array<Scalar, 3> center)
        : m_radius(radius)
        , m_center(center)
    {}

    Scalar value(std::array<Scalar, 3> pos) const override
    {
        return std::sqrt(
                   (pos[0] - m_center[0]) * (pos[0] - m_center[0]) +
                   (pos[1] - m_center[1]) * (pos[1] - m_center[1]) +
                   (pos[2] - m_center[2]) * (pos[2] - m_center[2])) -
               m_radius;
    }

    std::array<Scalar, 3> gradient(std::array<Scalar, 3> pos) const override
    {
        Scalar r = std::sqrt(
            (pos[0] - m_center[0]) * (pos[0] - m_center[0]) +
            (pos[1] - m_center[1]) * (pos[1] - m_center[1]) +
            (pos[2] - m_center[2]) * (pos[2] - m_center[2]));
        if (r == 0) return {0, 0, 0};

        return {(pos[0] - m_center[0]) / r, (pos[1] - m_center[1]) / r, (pos[2] - m_center[2]) / r};
    }

private:
    Scalar m_radius;
    std::array<Scalar, 3> m_center;
};

} // namespace mtetcol
