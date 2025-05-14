#pragma once

#include <mtetcol/common.h>

#include "implicit_function.h"
#include "transform.h"

#include <array>
#include <stdexcept>

namespace mtetcol {

/**
 * @brief Abstract base class for space-time functions
 *
 * This class defines the interface for functions that depend on both space and time.
 * It provides methods to evaluate the function value, its time derivative, and gradient.
 *
 * @tparam dim The spatial dimension of the function
 */
template <int dim>
class SpaceTimeFunction
{
public:
    virtual ~SpaceTimeFunction() = default;

    /**
     * @brief Evaluate the function at a given position and time
     *
     * @param pos The spatial position as an array of coordinates
     * @param t The time value
     * @return Scalar The function value at the given position and time
     */
    virtual Scalar value(std::array<Scalar, dim> pos, Scalar t) const = 0;

    /**
     * @brief Compute the time derivative of the function
     *
     * @param pos The spatial position as an array of coordinates
     * @param t The time value
     * @return Scalar The time derivative at the given position and time
     */
    virtual Scalar time_derivative(std::array<Scalar, dim> pos, Scalar t) const = 0;

    /**
     * @brief Compute the gradient of the function with respect to both space and time
     *
     * The gradient is returned as an array of size dim+1, where the first dim elements
     * represent the spatial gradient and the last element represents the time derivative.
     *
     * @param pos The spatial position as an array of coordinates
     * @param t The time value
     * @return std::array<Scalar, dim + 1> The gradient vector
     */
    virtual std::array<Scalar, dim + 1> gradient(std::array<Scalar, dim> pos, Scalar t) const = 0;
};

/**
 * @brief Concrete implementation of SpaceTimeFunction using explicit function definitions
 *
 * This class allows creating space-time functions from explicit function definitions
 * for the value, time derivative, and gradient. If the time derivative or gradient
 * are not provided, they are computed using finite differences.
 *
 * @tparam dim The spatial dimension of the function
 */
template <int dim>
class ExplicitForm : public SpaceTimeFunction<dim>
{
public:
    /**
     * @brief Construct a new ExplicitForm object
     *
     * @param func The function defining the value
     * @param time_derivative Optional function defining the time derivative
     * @param gradient Optional function defining the gradient
     */
    ExplicitForm(
        std::function<Scalar(std::array<Scalar, dim>, Scalar)> func,
        std::function<Scalar(std::array<Scalar, dim>, Scalar)> time_derivative = nullptr,
        std::function<std::array<Scalar, dim + 1>(std::array<Scalar, dim>, Scalar)> gradient =
            nullptr)
        : m_function(func)
        , m_time_derivative(time_derivative)
        , m_gradient(gradient)
    {
        assert(m_function != nullptr);
    }

    /**
     * @brief Evaluate the function at a given position and time
     *
     * Evaluates the function using the provided function object. This is a direct
     * evaluation of the function without any numerical approximation.
     *
     * @param pos The spatial position as an array of coordinates
     * @param t The time value
     * @return Scalar The function value at the given position and time
     */
    virtual Scalar value(std::array<Scalar, dim> pos, Scalar t) const override
    {
        return m_function(pos, t);
    }

    /**
     * @brief Compute the time derivative of the function
     *
     * If a time derivative function was provided during construction, it is used directly.
     * Otherwise, the time derivative is approximated using a forward finite difference
     * with a small time step (1e-6).
     *
     * @param pos The spatial position as an array of coordinates
     * @param t The time value
     * @return Scalar The time derivative at the given position and time
     */
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

    /**
     * @brief Compute the gradient of the function
     *
     * If a gradient function was provided during construction, it is used directly.
     * Otherwise, the gradient is approximated using forward finite differences with
     * a small step size (1e-6) for each spatial dimension. The time component of
     * the gradient is computed using the time_derivative method.
     *
     * @param pos The spatial position as an array of coordinates
     * @param t The time value
     * @return std::array<Scalar, dim + 1> The gradient vector, where the first dim
     *         elements represent the spatial gradient and the last element represents
     *         the time derivative
     */
    virtual std::array<Scalar, dim + 1> gradient(std::array<Scalar, dim> pos, Scalar t)
        const override
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
    std::function<Scalar(std::array<Scalar, dim>, Scalar)>
        m_function; ///< The function defining the value
    std::function<Scalar(std::array<Scalar, dim>, Scalar)>
        m_time_derivative; ///< Optional function defining the time derivative
    std::function<std::array<Scalar, dim + 1>(std::array<Scalar, dim>, Scalar)>
        m_gradient; ///< Optional function defining the gradient
};

/**
 * @brief Space-time function created by sweeping an implicit function through space
 *
 * This class represents a space-time function created by applying a transformation
 * to an implicit function. The transformation can be time-dependent, allowing the
 * implicit function to move through space.
 *
 * The swept function F(x,t) is defined as F(x,t) = f(T(x,t)), where f is the implicit
 * function and T is the transformation. The time derivative and gradient are computed
 * using the chain rule, taking into account both the spatial gradient of the implicit
 * function and the properties of the transformation.
 *
 * @tparam dim The spatial dimension of the function (2 or 3)
 */
template <int dim>
class SweepFunction : public SpaceTimeFunction<dim>
{
public:
    /**
     * @brief Construct a new SweepFunction object
     *
     * Creates a new swept function by combining an implicit function with a transformation.
     * The implicit function will be evaluated at positions transformed by the given
     * transformation.
     *
     * @param implicit_function The implicit function to be swept through space
     * @param transform The transformation to apply to the implicit function
     */
    SweepFunction(ImplicitFunction<dim>& implicit_function, Transform<dim>& transform)
        : m_implicit_function(&implicit_function)
        , m_transform(&transform)
    {}

    /**
     * @brief Evaluate the swept function at a given position and time
     *
     * Computes F(x,t) = f(T(x,t)), where f is the implicit function and T is the
     * transformation. The function first transforms the input position using the
     * transformation, then evaluates the implicit function at the transformed position.
     *
     * @param pos The spatial position as an array of coordinates
     * @param t The time value
     * @return Scalar The function value at the given position and time
     */
    virtual Scalar value(std::array<Scalar, dim> pos, Scalar t) const override
    {
        assert(m_implicit_function != nullptr);
        assert(m_transform != nullptr);
        auto transformed_pos = m_transform->transform(pos, t);
        return m_implicit_function->value(transformed_pos);
    }

    /**
     * @brief Compute the time derivative of the swept function
     *
     * The time derivative is computed using the chain rule: ∂F/∂t = ∇f · ∂T/∂t,
     * where ∇f is the spatial gradient of the implicit function and ∂T/∂t is the
     * velocity of the transformation. This represents how the function value changes
     * as the implicit function moves through space.
     *
     * @param pos The spatial position as an array of coordinates
     * @param t The time value
     * @return Scalar The time derivative at the given position and time
     */
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

    /**
     * @brief Compute the gradient of the swept function
     *
     * The gradient is computed using the chain rule, taking into account both
     * the spatial gradient of the implicit function and the Jacobian of the
     * transformation. The spatial part of the gradient is computed as ∇_x F = J^T ∇f,
     * where J is the position Jacobian of the transformation and ∇f is the gradient
     * of the implicit function. The time component is computed using the time_derivative
     * method.
     *
     * @param pos The spatial position as an array of coordinates
     * @param t The time value
     * @return std::array<Scalar, dim + 1> The gradient vector, where the first dim
     *         elements represent the spatial gradient and the last element represents
     *         the time derivative
     */
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
    ImplicitFunction<dim>* m_implicit_function = nullptr; ///< The implicit function being swept
    Transform<dim>* m_transform = nullptr; ///< The transformation applied to the implicit function
};

} // namespace mtetcol
