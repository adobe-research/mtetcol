#pragma once

#include <nonstd/indirect_value.hpp>
#include <span>

namespace mtetcol {

using Scalar = double;
using Index = uint32_t;
class MTetColImpl;


struct Contour
{
    /**
     * @brief The vertices of the contour.
     *
     * Each vertex is represented by (x, y, z, t).
     */
    std::vector<Scalar> vertices;

    /**
     * @brief Contour edge segements.
     *
     * Each segment is represented as an oriented edge (v_0, v_1).
     */
    std::vector<Index> segments;

    /**
     * @brief Contour cycles.
     *
     * Each cycle is a chain of segments that form a closed loop.
     */
    std::vector<Index> cycles;
    std::vector<Index> cycle_indices;
};

class MTetCol
{
public:
    MTetCol();
    ~MTetCol();

    MTetCol(const MTetCol&) noexcept;
    MTetCol& operator=(const MTetCol&) noexcept;

    MTetCol(MTetCol&&) noexcept;
    MTetCol& operator=(MTetCol&&) noexcept;

public:
    void add_vertices(std::span<const Scalar> vertices);
    void add_tets(std::span<const Index> tets);
    void add_zero_crossing(size_t vid, Scalar t);

    void set_cyclic(bool cyclic);
    void get_cyclic(bool& cyclic) const;

private:
    nonstd::indirect_value<MTetColImpl> m_impl;
};

} // namespace mtetcol
