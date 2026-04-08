#pragma once

// EMBER Integer Exact CSG - Core Types
// Uses geometry<26, 55> from Nehring-Wirxel et al. 2021:
//   26-bit input positions, 55-bit normals
//   All intermediate results stay within 256 bits

#include <integer-plane-geometry/geometry.hh>
#include <integer-plane-geometry/plane.hh>
#include <integer-plane-geometry/point.hh>
#include <integer-plane-geometry/line.hh>
#include <integer-plane-geometry/intersect.hh>
#include <integer-plane-geometry/classify.hh>
#include <integer-plane-geometry/any_point.hh>
#include <integer-plane-geometry/are_parallel.hh>
#include <integer-plane-geometry/integer_math.hh>

#include <cstdint>
#include <array>

namespace ember
{
// Primary geometry configuration: 26-bit positions, 55-bit normals
// This gives ~15nm accuracy for 1m^3 scenes (almost 8 decimal digits)
using geometry_t = ipg::geometry<26, 55>;

using plane_t = ipg::plane<geometry_t>;
using point4_t = ipg::point4<geometry_t>;
using line_t = ipg::line<geometry_t>;
using pos_t = geometry_t::pos_t;
using pos_scalar_t = geometry_t::pos_scalar_t;
using normal_scalar_t = geometry_t::normal_scalar_t;
using plane_d_t = geometry_t::plane_d_t;
using det_abc_t = geometry_t::determinant_abc_t;
using det_xxd_t = geometry_t::determinant_xxd_t;

// Integer AABB for subdivision
struct IAABB
{
    pos_t min;
    pos_t max;

    pos_scalar_t extent(int axis) const
    {
        return (&max.x)[axis] - (&min.x)[axis];
    }

    int longest_axis() const
    {
        auto ex = extent(0);
        auto ey = extent(1);
        auto ez = extent(2);
        if (ex >= ey && ex >= ez) return 0;
        if (ey >= ez) return 1;
        return 2;
    }

    pos_scalar_t midpoint(int axis) const
    {
        // Integer midpoint (rounds toward zero)
        return ((&min.x)[axis] + (&max.x)[axis]) / 2;
    }

    // Create an axis-aligned splitting plane at position 'val' along 'axis'
    // Plane normal points in positive axis direction: n[axis]=1, d=-val
    // classify(p, plane) > 0 means p is on the positive side (p[axis] > val)
    plane_t splitting_plane(int axis, pos_scalar_t val) const
    {
        plane_t p;
        p.a = (axis == 0) ? normal_scalar_t(1) : normal_scalar_t(0);
        p.b = (axis == 1) ? normal_scalar_t(1) : normal_scalar_t(0);
        p.c = (axis == 2) ? normal_scalar_t(1) : normal_scalar_t(0);
        // d = -val (so that ax+by+cz+d = x[axis] - val = 0 at the split)
        p.d = plane_d_t(-val);
        return p;
    }

    // Split this AABB into two halves along the given axis at the given value
    IAABB left_half(int axis, pos_scalar_t val) const
    {
        IAABB r = *this;
        (&r.max.x)[axis] = val;
        return r;
    }

    IAABB right_half(int axis, pos_scalar_t val) const
    {
        IAABB r = *this;
        (&r.min.x)[axis] = val;
        return r;
    }

    bool contains(pos_t const& p) const
    {
        return p.x >= min.x && p.x <= max.x &&
               p.y >= min.y && p.y <= max.y &&
               p.z >= min.z && p.z <= max.z;
    }
};

// Maximum input coordinate magnitude: 2^26 - 2 = 67108862
// Leaves room for clipping planes at AABB boundaries
static constexpr int32_t MAX_COORD = (1 << 25) - 1;

} // namespace ember
