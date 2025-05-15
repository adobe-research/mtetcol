#pragma once

#include <mtetcol/common.h>

#include <array>
#include <cmath>
#include <numbers>

namespace mtetcol {

inline std::pair<Scalar, std::array<Scalar, 4>> raw_flipping_donut(std::array<Scalar, 4> inputs)
{
    Scalar value;
    std::array<Scalar, 4> gradient;
    // Unpack the inputs
    Scalar xx = inputs[0];
    Scalar yy = inputs[1];
    Scalar zz = inputs[2];
    Scalar tt = inputs[3];

    // Constants
    constexpr Scalar pi = std::numbers::pi;

    // Precomputed terms
    Scalar cos_pi_tt = std::cos(pi * tt);
    Scalar sin_pi_tt = std::sin(pi * tt);
    Scalar term_zz = -0.51 - 0.01 * tt + zz;
    Scalar term_tt = -0.5 - 0.01 * (1 - tt);
    Scalar term_tt2 = -0.25 - 0.01 * (1 - tt) - 0.51 * tt;

    Scalar sqrt_inner = std::sqrt(
        0.0 + term_zz * term_zz +
        std::pow(
            -0.01 + term_tt * cos_pi_tt + yy * cos_pi_tt + term_tt2 * sin_pi_tt + xx * sin_pi_tt,
            2));
    Scalar sqrt_outer = -0.2 + sqrt_inner;

    Scalar term1 =
        -0.01 + term_tt2 * cos_pi_tt + xx * cos_pi_tt - term_tt * sin_pi_tt - yy * sin_pi_tt;

    // Compute scalar value
    value = -0.0025 + std::pow(term1, 2) + std::pow(sqrt_outer, 2);

    // Compute gradient
    // Gradient w.r.t. xx
    gradient[0] =
        2 * cos_pi_tt * term1 +
        (2 * sin_pi_tt *
         (-0.01 + term_tt * cos_pi_tt + yy * cos_pi_tt + term_tt2 * sin_pi_tt + xx * sin_pi_tt) *
         sqrt_outer) /
            sqrt_inner;

    // Gradient w.r.t. yy
    gradient[1] =
        -2 * sin_pi_tt * term1 +
        (2 * cos_pi_tt *
         (-0.01 + term_tt * cos_pi_tt + yy * cos_pi_tt + term_tt2 * sin_pi_tt + xx * sin_pi_tt) *
         sqrt_outer) /
            sqrt_inner;

    // Gradient w.r.t. zz
    gradient[2] = (2 * term_zz * sqrt_outer) / sqrt_inner;

    // Gradient w.r.t. tt
    gradient[3] = 2 *
                      (-0.5 * cos_pi_tt - pi * term_tt * cos_pi_tt - pi * yy * cos_pi_tt -
                       0.01 * sin_pi_tt - pi * term_tt2 * sin_pi_tt - pi * xx * sin_pi_tt) *
                      term1 +
                  ((-0.02 * term_zz +
                    2 *
                        (-0.01 + term_tt * cos_pi_tt + yy * cos_pi_tt + term_tt2 * sin_pi_tt +
                         xx * sin_pi_tt) *
                        (0.01 * cos_pi_tt + pi * term_tt2 * cos_pi_tt + pi * xx * cos_pi_tt -
                         0.5 * sin_pi_tt - pi * term_tt * sin_pi_tt - pi * yy * sin_pi_tt)) *
                   sqrt_outer) /
                      sqrt_inner;
    return {value, gradient};
}

} // namespace mtetcol
