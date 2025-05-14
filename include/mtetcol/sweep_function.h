#pragma once

#include <mtetcol/common.h>

#include "implicit_function.h"
#include "transform.h"

#include <array>
#include <stdexcept>

namespace mtetcol {

template <int dim>
class SpaceTimeFunction
{
public:
    virtual ~SpaceTimeFunction() = default;
    virtual Scalar value(std::array<Scalar, dim> pos, Scalar t) const = 0;
    virtual Scalar time_derivative(std::array<Scalar, dim> pos, Scalar t) const = 0;
    // Gradient w.r.t. (x, t)
    virtual std::array<Scalar, dim + 1> gradient(std::array<Scalar, dim> pos, Scalar t) const = 0;
};

template <int dim>
class ExplicitForm : public SpaceTimeFunction<dim>
{
public:
    ExplicitForm(
        std::function<Scalar(std::array<Scalar, dim>, Scalar)> func,
        std::function<Scalar(std::array<Scalar, dim>, Scalar)> time_derivative = nullptr,
        std::function<std::array<Scalar, dim + 1>(std::array<Scalar, dim>, Scalar)> gradient = nullptr)
        : m_function(func)
        , m_time_derivative(time_derivative)
        , m_gradient(gradient)
    {
        assert(m_function != nullptr);
    }

    virtual Scalar value(std::array<Scalar, dim> pos, Scalar t) const override
    {
        return m_function(pos, t);
    }

    virtual Scalar time_derivative(std::array<Scalar, dim> pos, Scalar t) const override
    {
        if (m_time_derivative == nullptr) {
            // Finite difference
            auto delta_t = 1e-6;
            auto value1 = m_function(pos, t);
            auto value2 = m_function(pos, t + delta_t);
            return (value2 - value1) / delta_t;
        } else {
            return m_time_derivative(pos, t);
        }
    }

    virtual std::array<Scalar, dim + 1> gradient(std::array<Scalar, dim> pos, Scalar t) const override
    {
        if (m_gradient == nullptr) {
            // Finite difference
            auto delta = 1e-6;
            std::array<Scalar, dim + 1> gradient;
            for (int i = 0; i < dim; ++i) {
                auto pos_delta = pos;
                pos_delta[i] += delta;
                auto value1 = m_function(pos, t);
                auto value2 = m_function(pos_delta, t);
                gradient[i] = (value2 - value1) / delta;
            }
            gradient[dim] = time_derivative(pos, t);
            return gradient;
        } else {
            return m_gradient(pos, t);
        }
    }

private:
    std::function<Scalar(std::array<Scalar, dim>, Scalar)> m_function;
    std::function<Scalar(std::array<Scalar, dim>, Scalar)> m_time_derivative;
    std::function<std::array<Scalar, dim + 1>(std::array<Scalar, dim>, Scalar)> m_gradient;
};

template <int dim>
class SweepFunction : public SpaceTimeFunction<dim>
{
public:
    SweepFunction(ImplicitFunction<dim>& implicit_function, Transform<dim>& transform)
        : m_implicit_function(&implicit_function)
        , m_transform(&transform)
    {}

    virtual Scalar value(std::array<Scalar, dim> pos, Scalar t) const override
    {
        assert(m_implicit_function != nullptr);
        assert(m_transform != nullptr);
        auto transformed_pos = m_transform->transform(pos, t);
        return m_implicit_function->value(transformed_pos);
    }

    virtual Scalar time_derivative(std::array<Scalar, dim> pos, Scalar t) const override
    {
        assert(m_implicit_function != nullptr);
        assert(m_transform != nullptr);
        auto transformed_pos = m_transform->transform(pos, t);
        auto velocity = m_transform->velocity(pos, t);
        auto spacial_grad = m_implicit_function->gradient(transformed_pos);
        if constexpr (dim == 2) {
            return spacial_grad[0] * velocity[0] + spacial_grad[1] * velocity[1];
        } else if constexpr (dim == 3) {
            return spacial_grad[0] * velocity[0] + spacial_grad[1] * velocity[1] +
                   spacial_grad[2] * velocity[2];
        } else {
            throw std::runtime_error("Unsupported dimension");
        }
    }

    virtual std::array<Scalar, dim + 1> gradient(std::array<Scalar, dim> pos, Scalar t)
        const override
    {
        assert(m_implicit_function != nullptr);
        assert(m_transform != nullptr);

        const auto transformed_pos = m_transform->transform(pos, t);
        const auto g_f = m_implicit_function->gradient(transformed_pos);
        const auto J = m_transform->position_Jacobian(pos, t);

        /* spatial part  ∇_x F = Jᵀ ∇f */
        std::array<Scalar, dim + 1> grad{};
        for (int i = 0; i < dim; ++i) {
            Scalar sum = 0;
            for (int k = 0; k < dim; ++k) sum += J[k][i] * g_f[k];
            grad[i] = sum;
        }

        /* time component */
        grad[dim] = time_derivative(pos, t);

        return grad;
    }

private:
    ImplicitFunction<dim>* m_implicit_function = nullptr;
    Transform<dim>* m_transform = nullptr;
};

} // namespace mtetcol
