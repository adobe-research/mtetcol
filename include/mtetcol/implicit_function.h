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

class ImplicitCircle : public ImplicitFunction<2>
{
public:
    ImplicitCircle(Scalar radius, std::array<Scalar, 2> center)
        : m_radius(radius)
        , m_center(center)
    {}

    Scalar value(std::array<Scalar, 2> pos) const override
    {
        return std::sqrt(
                   (pos[0] - m_center[0]) * (pos[0] - m_center[0]) +
                   (pos[1] - m_center[1]) * (pos[1] - m_center[1])) -
               m_radius;
    }

    std::array<Scalar, 2> gradient(std::array<Scalar, 2> pos) const override
    {
        Scalar r = std::sqrt(
            (pos[0] - m_center[0]) * (pos[0] - m_center[0]) +
            (pos[1] - m_center[1]) * (pos[1] - m_center[1]));
        if (r == 0) return {0, 0};

        return {(pos[0] - m_center[0]) / r, (pos[1] - m_center[1]) / r};
    }

private:
    Scalar m_radius;
    std::array<Scalar, 2> m_center;
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
        Scalar len_xy = std::sqrt(x * x + y * y);

        return std::sqrt(z * z + (len_xy - m_R) * (len_xy - m_R)) - m_r;
    }

    std::array<Scalar, 3> gradient(std::array<Scalar, 3> pos) const override
    {
        Scalar x = pos[0] - m_center[0];
        Scalar y = pos[1] - m_center[1];
        Scalar z = pos[2] - m_center[2];

        Scalar len_xy = std::sqrt(x * x + y * y);

        // Avoid division by zero (if point is at z-axis)
        if (len_xy < 1e-6f) {
            return {0, 0, static_cast<Scalar>(z >= 0 ? 1 : -1)};
        }

        Scalar a = len_xy - m_R;
        Scalar q_len = std::sqrt(a * a + z * z);

        // Again avoid division by zero if exactly on torus surface
        if (q_len < 1e-6f) {
            return {0, 0, 0}; // Undefined
        }

        Scalar dx = (a / q_len) * (x / len_xy);
        Scalar dy = (a / q_len) * (y / len_xy);
        Scalar dz = z / q_len;

        return {dx, dy, dz};
    }

private:
    Scalar m_R;
    Scalar m_r;
    std::array<Scalar, 3> m_center;
};

} // namespace mtetcol
