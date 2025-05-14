#pragma once

#include <mtetcol/common.h>

#include <array>
#include <span>

namespace mtetcol {

/**
 * @brief Base class for geometric transformations in n-dimensional space.
 *
 * This abstract class defines the interface for geometric transformations that can be applied
 * to points in n-dimensional space. It provides methods for both position transformation
 * and velocity calculation.
 *
 * @tparam dim The dimensionality of the space (2D or 3D)
 */
template <int dim>
class Transform
{
public:
    /**
     * @brief Transforms a point in space according to the transformation rules.
     *
     * @param pos The input position to transform
     * @param t The time parameter for time-dependent transformations
     * @return std::array<Scalar, dim> The transformed position
     */
    virtual std::array<Scalar, dim> transform(std::array<Scalar, dim> pos, Scalar t) const = 0;

    /**
     * @brief Calculates the velocity of a point under the transformation.
     *
     * @param pos The position at which to calculate the velocity
     * @param t The time parameter for time-dependent transformations
     * @return std::array<Scalar, dim> The velocity vector
     */
    virtual std::array<Scalar, dim> velocity(std::array<Scalar, dim> pos, Scalar t) const = 0;
};

/**
 * @brief A translation transformation that moves points along a constant vector.
 *
 * This class implements a translation transformation where points are moved along
 * a constant vector scaled by time.
 *
 * @tparam dim The dimensionality of the space (2D or 3D)
 */
template <int dim>
class Translation : public Transform<dim>
{
public:
    /**
     * @brief Constructs a translation transformation.
     *
     * @param translation The translation vector that defines the direction and magnitude
     *                    of the translation
     */
    Translation(std::array<Scalar, dim> translation)
        : m_translation(translation)
    {}

    std::array<Scalar, dim> transform(std::array<Scalar, dim> pos, Scalar t) const override
    {
        for (int i = 0; i < dim; ++i) {
            pos[i] += m_translation[i] * t;
        }
        return pos;
    }

    std::array<Scalar, dim> velocity(std::array<Scalar, dim> pos, Scalar t) const override
    {
        return m_translation;
    }

private:
    std::array<Scalar, dim> m_translation;
};

/**
 * @brief A rotation transformation around an axis in 3D or around a point in 2D.
 *
 * This class implements rotation transformations. In 3D, it rotates points around
 * a specified axis through a center point. In 2D, it rotates points around a center point.
 *
 * @tparam dim The dimensionality of the space (2D or 3D)
 */
template <int dim>
class Rotation : public Transform<dim>
{
public:
    /**
     * @brief Constructs a rotation transformation.
     *
     * @param center The center point of rotation
     * @param axis The rotation axis (only used in 3D)
     * @param angle The total angle of rotation in degrees (default: 360)
     */
    Rotation(std::array<Scalar, dim> center, std::array<Scalar, dim> axis, Scalar angle = 360)
        : m_axis(axis)
        , m_center(center)
        , m_angle(angle)
    {}

    std::array<Scalar, dim> transform(std::array<Scalar, dim> pos, Scalar t) const override
    {
        if constexpr (dim == 3) {
            // Convert angle to radians
            Scalar angle = t * m_angle * M_PI / 180.0;

            // Normalize the axis
            Scalar axis_length = 0;
            for (int i = 0; i < dim; ++i) {
                axis_length += m_axis[i] * m_axis[i];
            }
            axis_length = std::sqrt(axis_length);

            std::array<Scalar, dim> normalized_axis;
            for (int i = 0; i < dim; ++i) {
                normalized_axis[i] = m_axis[i] / axis_length;
            }

            // Rodrigues' rotation formula
            Scalar cos_angle = std::cos(angle);
            Scalar sin_angle = std::sin(angle);

            // Translate point to origin
            for (int i = 0; i < dim; ++i) {
                pos[i] -= m_center[i];
            }

            // Apply rotation
            std::array<Scalar, dim> result;
            result[0] =
                pos[0] * (cos_angle + normalized_axis[0] * normalized_axis[0] * (1 - cos_angle)) +
                pos[1] * (normalized_axis[0] * normalized_axis[1] * (1 - cos_angle) -
                          normalized_axis[2] * sin_angle) +
                pos[2] * (normalized_axis[0] * normalized_axis[2] * (1 - cos_angle) +
                          normalized_axis[1] * sin_angle);

            result[1] =
                pos[0] * (normalized_axis[1] * normalized_axis[0] * (1 - cos_angle) +
                          normalized_axis[2] * sin_angle) +
                pos[1] * (cos_angle + normalized_axis[1] * normalized_axis[1] * (1 - cos_angle)) +
                pos[2] * (normalized_axis[1] * normalized_axis[2] * (1 - cos_angle) -
                          normalized_axis[0] * sin_angle);

            result[2] =
                pos[0] * (normalized_axis[2] * normalized_axis[0] * (1 - cos_angle) -
                          normalized_axis[1] * sin_angle) +
                pos[1] * (normalized_axis[2] * normalized_axis[1] * (1 - cos_angle) +
                          normalized_axis[0] * sin_angle) +
                pos[2] * (cos_angle + normalized_axis[2] * normalized_axis[2] * (1 - cos_angle));

            // Translate back from origin
            for (int i = 0; i < dim; ++i) {
                result[i] += m_center[i];
            }

            return result;
        } else {
            static_assert(dim == 2, "Rotation is only implemented for 2D and 3d");

            // Convert angle to radians
            Scalar angle = t * m_angle * M_PI / 180.0;

            pos[0] -= m_center[0];
            pos[1] -= m_center[1];

            std::array<Scalar, dim> result;
            result[0] = pos[0] * std::cos(angle) - pos[1] * std::sin(angle) + m_center[0];
            result[1] = pos[0] * std::sin(angle) + pos[1] * std::cos(angle) + m_center[1];

            return result;
        }
    }

    std::array<Scalar, dim> velocity(std::array<Scalar, dim> pos, Scalar t) const override
    {
        if constexpr (dim == 3) {
            // Normalize the axis
            Scalar axis_length = 0;
            for (int i = 0; i < dim; ++i) {
                axis_length += m_axis[i] * m_axis[i];
            }
            axis_length = std::sqrt(axis_length);

            std::array<Scalar, dim> normalized_axis;
            for (int i = 0; i < dim; ++i) {
                normalized_axis[i] = m_axis[i] / axis_length;
            }

            pos = transform(pos, t);
            // Translate point to origin
            for (int i = 0; i < dim; ++i) {
                pos[i] -= m_center[i];
            }

            // Cross product of axis and position gives the velocity direction
            std::array<Scalar, dim> velocity;
            velocity[0] = (normalized_axis[1] * pos[2] - normalized_axis[2] * pos[1]) * m_angle *
                          M_PI / 180.0;
            velocity[1] = (normalized_axis[2] * pos[0] - normalized_axis[0] * pos[2]) * m_angle *
                          M_PI / 180.0;
            velocity[2] = (normalized_axis[0] * pos[1] - normalized_axis[1] * pos[0]) * m_angle *
                          M_PI / 180.0;

            return velocity;
        } else {
            static_assert(dim == 2, "Rotation is only implemented for 2D and 3d");

            pos = transform(pos, t);
            pos[0] -= m_center[0];
            pos[1] -= m_center[1];

            return {
                -pos[1] * m_angle * M_PI / 180.0,
                pos[0] * m_angle * M_PI / 180.0,
            };
        }
    }

    /**
     * @brief Calculates velocity using finite difference approximation.
     *
     * This is a helper method that computes velocity using central difference
     * approximation. It's used for verification purposes.
     *
     * @param pos The position at which to calculate the velocity
     * @param t The time parameter
     * @return std::array<Scalar, dim> The approximated velocity vector
     */
    std::array<Scalar, dim> finite_difference(std::array<Scalar, dim> pos, Scalar t) const
    {
        constexpr Scalar delta = 1e-6;
        std::array<Scalar, dim> pos_plus_delta = transform(pos, t + delta);
        std::array<Scalar, dim> pos_minus_delta = transform(pos, t - delta);
        std::array<Scalar, dim> velocity;

        for (int i = 0; i < dim; ++i) {
            velocity[i] = (pos_plus_delta[i] - pos_minus_delta[i]) / (2 * delta);
        }

        return velocity;
    }

private:
    std::array<Scalar, dim> m_center; ///< Center point of rotation
    std::array<Scalar, dim> m_axis; ///< Rotation axis (3D only)
    Scalar m_angle; ///< Total rotation angle in degrees
};

/**
 * @brief Composes two transformations by applying them in sequence.
 *
 * This class combines two transformations by applying them one after another.
 * The first transformation is applied to the input position, and then the second
 * transformation is applied to the result.
 *
 * @tparam dim The dimensionality of the space (2D or 3D)
 */
template <int dim>
class Compose : public Transform<dim>
{
public:
    /**
     * @brief Constructs a composition of two transformations.
     *
     * @param transform1 The first transformation to apply
     * @param transform2 The second transformation to apply
     */
    Compose(Transform<dim>& transform1, Transform<dim>& transform2)
        : m_transform1(transform1)
        , m_transform2(transform2)
    {}

    std::array<Scalar, dim> transform(std::array<Scalar, dim> pos, Scalar t) const override
    {
        auto intermediate = m_transform1.transform(pos, t);
        return m_transform2.transform(intermediate, t);
    }

    std::array<Scalar, dim> velocity(std::array<Scalar, dim> pos, Scalar t) const override
    {
        // Using finite difference to calculate the velocity
        constexpr Scalar delta = 1e-6;
        auto value_prev = transform(pos, t - delta);
        auto value_next = transform(pos, t + delta);
        std::array<Scalar, dim> velocity;
        for (int i = 0; i < dim; ++i) {
            velocity[i] = (value_next[i] - value_prev[i]) / (2 * delta);
        }
        return velocity;
    }

private:
    Transform<dim>& m_transform1; ///< First transformation
    Transform<dim>& m_transform2; ///< Second transformation
};

} // namespace mtetcol
