#pragma once

#include <mtetcol/common.h>

#include <array>
#include <nonstd/indirect_value.hpp>
#include <span>

namespace mtetcol {

template <int dim>
class Transform
{
public:
    virtual std::array<Scalar, dim> transform(std::array<Scalar, dim> pos, Scalar t) const = 0;
    virtual std::array<Scalar, dim> velocity(std::array<Scalar, dim> pos, Scalar t) const = 0;
};

template <int dim>
class Translation : public Transform<dim>
{
public:
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

template <int dim>
class Rotation : public Transform<dim>
{
public:
    Rotation(std::array<Scalar, dim> center, std::array<Scalar, dim> axis)
        : m_axis(axis)
        , m_center(center)
    {}

    std::array<Scalar, dim> transform(std::array<Scalar, dim> pos, Scalar t) const override
    {
        static_assert(dim == 3, "Rotation is only implemented for 3D");

        // Convert angle to radians
        Scalar angle = t * 2 * M_PI; // Full rotation in 1 second

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
    }

    std::array<Scalar, dim> velocity(std::array<Scalar, dim> pos, Scalar t) const override
    {
        static_assert(dim == 3, "Rotation is only implemented for 3D");

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
        velocity[0] = (normalized_axis[1] * pos[2] - normalized_axis[2] * pos[1]) * 2 * M_PI;
        velocity[1] = (normalized_axis[2] * pos[0] - normalized_axis[0] * pos[2]) * 2 * M_PI;
        velocity[2] = (normalized_axis[0] * pos[1] - normalized_axis[1] * pos[0]) * 2 * M_PI;

        return velocity;
    }

    std::array<Scalar, dim> finite_difference(
        std::array<Scalar, dim> pos, Scalar t) const
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
    std::array<Scalar, dim> m_center;
    std::array<Scalar, dim> m_axis;
};


} // namespace mtetcol
