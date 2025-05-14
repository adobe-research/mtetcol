#pragma once

#include <mtetcol/common.h>

#include <array>

namespace mtetcol {

/**
 * @brief Base class for implicit functions in N-dimensional space.
 *
 * An implicit function defines a surface as the zero level set of a scalar function.
 * The function returns positive values outside the surface, negative values inside,
 * and zero on the surface.
 *
 * @tparam dim The dimension of the space (2 for 2D, 3 for 3D)
 */
template <int dim>
class ImplicitFunction
{
public:
    /**
     * @brief Evaluates the implicit function at a given position.
     *
     * @param pos The position to evaluate at
     * @return Scalar The signed distance to the surface (positive outside, negative inside)
     */
    virtual Scalar value(std::array<Scalar, dim> pos) const = 0;

    /**
     * @brief Computes the gradient of the implicit function at a given position.
     *
     * @param pos The position to evaluate at
     * @return std::array<Scalar, dim> The normalized gradient vector
     */
    virtual std::array<Scalar, dim> gradient(std::array<Scalar, dim> pos) const = 0;
};

/**
 * @brief Implicit function representing a 2D circle.
 *
 * The circle is defined by its radius and center point.
 */
class ImplicitCircle : public ImplicitFunction<2>
{
public:
    /**
     * @brief Constructs a new implicit circle.
     *
     * @param radius The radius of the circle
     * @param center The center point of the circle
     */
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
    Scalar m_radius; ///< The radius of the circle
    std::array<Scalar, 2> m_center; ///< The center point of the circle
};

/**
 * @brief Implicit function representing a 3D sphere.
 *
 * The sphere is defined by its radius and center point.
 */
class ImplicitSphere : public ImplicitFunction<3>
{
public:
    /**
     * @brief Constructs a new implicit sphere.
     *
     * @param radius The radius of the sphere
     * @param center The center point of the sphere
     */
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
    Scalar m_radius; ///< The radius of the sphere
    std::array<Scalar, 3> m_center; ///< The center point of the sphere
};

/**
 * @brief Implicit function representing a 3D torus.
 *
 * The torus is defined by two radii (R and r) and a center point.
 * R is the distance from the center of the tube to the center of the torus,
 * and r is the radius of the tube.
 */
class ImplicitTorus : public ImplicitFunction<3>
{
public:
    /**
     * @brief Constructs a new implicit torus.
     *
     * @param R The major radius (distance from center of tube to center of torus)
     * @param r The minor radius (radius of the tube)
     * @param center The center point of the torus
     */
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
    Scalar m_R; ///< The major radius of the torus
    Scalar m_r; ///< The minor radius of the torus
    std::array<Scalar, 3> m_center; ///< The center point of the torus
};

} // namespace mtetcol
