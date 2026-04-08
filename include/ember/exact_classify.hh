#pragma once

// EMBER - Exact point4-vs-plane classify that correctly handles sign extension.
//
// The upstream ipg::classify has a bug: mul<253>(det_abc_t, plane_d_t) calls
// imul<4>(i192, i128) which doesn't sign-extend the 3-word i192 to 4 words
// before multiplying. For negative w values, the high bits of the product
// are wrong, causing classify to return incorrect signs.
//
// This file provides a correct implementation by promoting all operands to
// i256 before multiplication.

#include <ember/types.hh>

namespace ember
{

// Correctly classify a point4 against a plane, handling sign extension.
// Returns -1, 0, or +1.
inline tg::i8 exact_classify(point4_t const& pt, plane_t const& s)
{
    // Promote ALL operands to i256 to avoid sign-extension bugs in imul.
    // The computation: d = x*a + y*b + z*c + w*d_plane
    // where x,y,z are already i256, w is i192 (needs promotion), and
    // a,b,c are i64, d_plane is i128.

    using wide_t = det_xxd_t; // i256

    // Promote w to i256 (sign-extends correctly via constructor)
    wide_t w_wide = wide_t(pt.w);

    // Promote plane coefficients to i256
    wide_t a_wide = wide_t(int64_t(s.a));
    wide_t b_wide = wide_t(int64_t(s.b));
    wide_t c_wide = wide_t(int64_t(s.c));
    wide_t d_wide = wide_t(s.d); // plane_d_t (i128) → i256

    // Compute dot product entirely in i256.
    // Each product is i256 * i256 → need at most 512 bits for exact product,
    // but the VALUES are bounded: x (196 bits) * a (55 bits) = 251 bits,
    // w (168 bits) * d (83 bits) = 251 bits. Sum fits in 253 bits < 256 bits.
    // So i256 arithmetic is sufficient for the VALUES even though it truncates
    // the TYPE-based product width.

    // Use operator* directly (i256 * i256 → i256, which truncates to 256 bits
    // but is correct since the values fit in 253 bits).
    wide_t result = pt.x * a_wide + pt.y * b_wide + pt.z * c_wide + w_wide * d_wide;

    return tg::sign(result) * tg::sign(pt.w);
}

// Overload for pos_t (i32) - delegates to the original (no sign-extension issue)
inline tg::i8 exact_classify(pos_t const& pt, plane_t const& p)
{
    return ipg::classify<geometry_t>(pt, p);
}

} // namespace ember
