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

class ImplicitTorus : public ImplicitFunction<3>
{
public:
    ImplicitTorus(Scalar R, Scalar r, std::array<Scalar, 3> center)
        : m_R(R)
        , m_r(r)
        , m_center(center)
    {}

    Scalar value(std::array<Scalar, 3> pos) const override
    {
        Scalar x = pos[0] - m_center[0];
        Scalar y = pos[1] - m_center[1];
        Scalar z = pos[2] - m_center[2];

        return std::sqrt(z * z + (std::sqrt(x * x + y * y) - m_R) * (std::sqrt(x * x + y * y) - m_R)) - m_r;
    }

    std::array<Scalar, 3> gradient(std::array<Scalar, 3> pos) const override
    {
        Scalar x = pos[0] - m_center[0];
        Scalar y = pos[1] - m_center[1];
        Scalar z = pos[2] - m_center[2];

        Scalar r = std::sqrt(x * x + y * y);
        Scalar r2 = r * r;
        Scalar r3 = r2 * r;

        return {x / r3, y / r3, z / (r2 + (r - m_R) * (r - m_R))};
    }

private:
    Scalar m_R;
    Scalar m_r;
    std::array<Scalar, 3> m_center;
}

} // namespace mtetcol
