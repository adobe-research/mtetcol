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


} // namespace mtetcol
