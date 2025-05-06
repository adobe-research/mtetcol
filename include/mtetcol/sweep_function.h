#pragma once

#include <mtetcol/common.h>

#include "implicit_function.h"
#include "transform.h"

#include <array>
#include <stdexcept>

namespace mtetcol {

template <int dim>
class SweepFunction
{
public:
    SweepFunction(ImplicitFunction<dim>& implicit_function, Transform<dim>& transform)
        : m_implicit_function(&implicit_function)
        , m_transform(&transform)
    {}

    Scalar value(std::array<Scalar, dim> pos, Scalar t) const
    {
        assert(m_implicit_function != nullptr);
        assert(m_transform != nullptr);
        auto transformed_pos = m_transform->transform(pos, t);
        return m_implicit_function->value(transformed_pos);
    }

    Scalar time_derivative(std::array<Scalar, dim> pos, Scalar t) const
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

private:
    ImplicitFunction<dim>* m_implicit_function = nullptr;
    Transform<dim>* m_transform = nullptr;
};

} // namespace mtetcol
