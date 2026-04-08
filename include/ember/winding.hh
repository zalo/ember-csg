#pragma once

// EMBER Integer Exact CSG - Winding Number Types
// Winding Number Vectors (WNV) and Transition Vectors (WNTV)
// Boolean operation indicator functions

#include <vector>
#include <cstdint>
#include <cstdio>
#include <functional>

namespace ember
{

// Winding Number Vector: one integer per input mesh
// w_i counts how often mesh i was "entered" at a point
using WNV = std::vector<int>;

// Winding Number Transition Vector: change in WNV when crossing a polygon
// For a polygon belonging to mesh j with outward-facing normal:
//   delta_w = (0, ..., 0, 1, 0, ..., 0) with 1 at position j
using WNTV = std::vector<int>;

// Front and back winding numbers for a classified polygon
struct WindingPair
{
    WNV w_front; // WNV on the front side (normal direction)
    WNV w_back;  // WNV on the back side (opposite normal)
};

// Boolean operation types
enum class BooleanOp
{
    Union,              // A ∪ B
    Intersection,       // A ∩ B
    Difference,         // A - B
    SymmetricDifference // A ⊕ B
};

// Indicator function: given a WNV, returns true if the point is "inside" the result
using IndicatorFn = std::function<bool(WNV const&)>;

inline IndicatorFn make_indicator(BooleanOp op, int num_meshes)
{
    switch (op)
    {
    case BooleanOp::Union:
        return [](WNV const& w) {
            for (auto wi : w)
                if (wi != 0) return true;
            return false;
        };
    case BooleanOp::Intersection:
        return [](WNV const& w) {
            for (auto wi : w)
                if (wi == 0) return false;
            return true;
        };
    case BooleanOp::Difference:
        // A - B: inside A and outside B (and outside C, D, ...)
        return [](WNV const& w) {
            if (w[0] == 0) return false;
            for (size_t i = 1; i < w.size(); i++)
                if (w[i] != 0) return false;
            return true;
        };
    case BooleanOp::SymmetricDifference:
        return [](WNV const& w) {
            int count = 0;
            for (auto wi : w)
                if (wi != 0) count++;
            return (count % 2) == 1;
        };
    }
    return [](WNV const&) { return false; };
}

// Classify a polygon given its front and back WNVs:
//   Returns: +1 = emit as-is (out→in transition in result)
//            -1 = emit inverted (in→out transition, flip winding)
//             0 = skip (both sides same classification)
inline int classify_polygon_output(WNV const& w_front, WNV const& w_back, IndicatorFn const& indicator)
{
    bool front_in = indicator(w_front);
    bool back_in = indicator(w_back);

    if (!front_in && back_in)
        return +1; // out→in: emit as-is
    if (front_in && !back_in)
        return -1; // in→out: emit inverted
    return 0;      // in→in or out→out: skip
}

// Early termination: check if any reachable WNV from a reference point
// could produce a useful polygon classification.
// Returns true if the subproblem can be discarded.
inline bool can_early_terminate(WNV const& ref_wnv,
                                std::vector<WNTV> const& available_wntvs,
                                IndicatorFn const& indicator)
{
    // BFS/DFS over reachable WNVs
    // Start from ref_wnv, apply each available WNTV (add or subtract)
    // If all reachable (w_front, w_back) pairs result in skip → can terminate

    // Simple implementation: check a bounded neighborhood
    // For typical cases (2-mesh booleans), this is very fast

    std::vector<WNV> frontier = {ref_wnv};
    std::vector<WNV> visited = {ref_wnv};

    // Search deep enough to find all reachable WNV states.
    // For n meshes with bounded winding numbers, the reachable set is finite.
    // Use a generous limit.
    constexpr int MAX_DEPTH = 8;

    for (int depth = 0; depth < MAX_DEPTH && !frontier.empty(); depth++)
    {
        std::vector<WNV> next_frontier;
        for (auto const& w : frontier)
        {
            for (auto const& dw : available_wntvs)
            {
                // Apply +dw (crossing front-to-back) and -dw (crossing back-to-front)
                for (int sign = -1; sign <= 1; sign += 2)
                {
                    WNV w_next = w;
                    for (size_t i = 0; i < w.size() && i < dw.size(); i++)
                        w_next[i] += sign * dw[i];

                    // Check if this transition produces useful output
                    if (classify_polygon_output(w, w_next, indicator) != 0)
                        return false;
                    if (classify_polygon_output(w_next, w, indicator) != 0)
                        return false;

                    // Add to frontier if not visited
                    bool found = false;
                    for (auto const& v : visited)
                        if (v == w_next) { found = true; break; }
                    if (!found)
                    {
                        visited.push_back(w_next);
                        next_frontier.push_back(w_next);
                    }
                }
            }
        }
        frontier = std::move(next_frontier);
    }

    return true; // All reachable transitions produce skip → can discard
}

// Propagate WNV along a segment from point x to point y
// T contains polygons intersected by the segment
// For each polygon t in T: w_y += sign(dot(n_t, x-y)) * delta_w_t
inline WNV propagate_wnv(WNV const& w_x,
                          int sign_direction, // +1 or -1 based on dot(normal, segment_dir)
                          WNTV const& delta_w)
{
    WNV result = w_x;
    for (size_t i = 0; i < result.size() && i < delta_w.size(); i++)
        result[i] += sign_direction * delta_w[i];
    return result;
}

} // namespace ember
