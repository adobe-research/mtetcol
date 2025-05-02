#pragma once

#include <cassert>
#include <cmath>
#include <limits>

#include <strong_type/strong_type.hpp>


namespace mtetcol {

using Scalar = double;
using Index = uint32_t;

/**
 * Signed index to represent both orientation and index.
 */
using SignedIndex = strong::type<int32_t, struct SignedIndexTag, strong::equality>;

constexpr Index invalid_index = std::numeric_limits<Index>::max();
constexpr SignedIndex invalid_signed_index = SignedIndex(0);

/**
 * @brief Get the value of a signed index.
 *
 * @param signed_index The signed index to get the value from.
 *
 * @return The value of the signed index.
 */
[[nodiscard]] inline SignedIndex signed_index(Index index, bool orientation)
{
    if (orientation) {
        return SignedIndex(static_cast<int32_t>(index + 1));
    } else {
        return SignedIndex(-static_cast<int32_t>(index + 1));
    }
}

/**
 * @brief Get the unsigned index value of a signed index.
 *
 * @param signed_index The signed index to get the value from.
 *
 * @return The value of the unsigned index.
 */
[[nodiscard]] inline Index index(SignedIndex signed_index)
{
    const int32_t value = value_of(signed_index);
    assert(value != 0);
    if (value >= 0) {
        return static_cast<Index>(value - 1);
    } else {
        return static_cast<Index>(-value - 1);
    }
}

/**
 * @brief Get the orientation of a signed index.
 *
 * @param signed_index The signed index to get the orientation from.
 *
 * @return True if the signed index is positive, false otherwise.
 */
[[nodiscard]] inline bool orientation(SignedIndex signed_index)
{
    const int32_t value = value_of(signed_index);
    assert(value != 0);
    return value > 0;
}

[[nodiscard]] inline SignedIndex operator-(SignedIndex signed_index)
{
    return SignedIndex(-value_of(signed_index));
}

} // namespace mtetcol

