#pragma once

// EMBER Integer Exact CSG - Main Public API
//
// Usage:
//   1. Create InputMesh objects with positions and triangles
//   2. Call ember::boolean_operation() with the meshes and desired operation
//   3. Receive BooleanResult with classified output polygons
//
// The implementation follows the EMBER paper:
//   "Exact Mesh Booleans via Efficient and Robust Local Arrangements"
//   Nehring-Wirxel et al.
//
// Key guarantees:
//   - All internal computation is exact using integer homogeneous coordinates
//   - Output polygons are convex with well-defined winding numbers
//   - Handles coplanar faces, exact edge hits, self-intersections
//   - Inaccuracies only at import/export (floating-point ↔ integer conversion)

#include <ember/clip.hh>
#include <ember/intersect_polygons.hh>
#include <ember/local_bsp.hh>
#include <ember/mesh.hh>
#include <ember/resolve_tjunctions.hh>
#include <ember/polygon.hh>
#include <ember/segment_trace.hh>
#include <ember/subdivision.hh>
#include <ember/types.hh>
#include <ember/exact_classify.hh>
#include <ember/winding.hh>

#include <vector>

namespace ember
{

// Configuration for the EMBER algorithm
struct EmberConfig
{
    int leaf_threshold = LEAF_THRESHOLD; // Polygon count threshold for leaf processing
    int max_depth = MAX_DEPTH;          // Maximum subdivision depth
    bool use_early_termination = true;  // Enable WNV reachability early-out
    bool use_cog_splitting = true;      // Use center-of-gravity splitting heuristic
    bool assume_nsi = false;            // Assume no self-intersections for all meshes
    bool assume_nnc = false;            // Assume no nested components for all meshes
};

// Perform a boolean operation on two or more input meshes
inline BooleanResult boolean_operation(
    std::vector<InputMesh> const& meshes,
    BooleanOp op,
    EmberConfig const& config = {})
{
    BooleanResult result;

    if (meshes.empty())
        return result;

    // Step 1: Prepare input - convert to integer polygon soup
    PolygonSoup soup = prepare_input(meshes);

    if (soup.polygons.empty())
        return result;

    // Apply config overrides
    if (config.assume_nsi)
        for (auto& p : soup.polygons) p.no_self_intersections = true;
    if (config.assume_nnc)
        for (auto& p : soup.polygons) p.no_nested_components = true;

    // Step 2: Create indicator function for the boolean operation
    IndicatorFn indicator = make_indicator(op, soup.num_meshes);

    // Step 3: Set up initial reference point
    // The reference point should be outside all input meshes.
    // Use a corner of the AABB (which is padded beyond all mesh vertices).
    pos_t ref_point = soup.bounds.min;

    // Initial WNV: all zeros (outside all meshes)
    WNV ref_wnv(soup.num_meshes, 0);

    // Step 4: Create initial subdivision task
    SubdivisionTask initial_task;
    initial_task.polygons = std::move(soup.polygons);
    initial_task.bounds = soup.bounds;
    initial_task.root_bounds = soup.bounds;
    initial_task.ref_point = ref_point;
    initial_task.ref_wnv = ref_wnv;
    initial_task.depth = 0;

    // Step 5: Run subdivision and classification
    std::vector<ClassifiedPolygon> classified;
    subdivide(std::move(initial_task), indicator, classified);

    // Step 6: Assemble output
    result.output.num_meshes = soup.num_meshes;
    result.output.bounds = soup.bounds;
    result.output.transform = soup.transform;
    result.transform = soup.transform;

    // Build AABB boundary planes for validity checking.
    // Valid BSP leaf vertices must lie within the root AABB because all input
    // vertices and all clipping planes are within bounds.
    // Vertices outside indicate degenerate BSP leaves from nearly-parallel planes.
    plane_t aabb_planes[6];
    for (int a = 0; a < 3; a++)
    {
        // x[a] >= min[a]  →  x[a] - min[a] >= 0  →  classify > 0 means OK
        aabb_planes[a*2].a = (a == 0) ? normal_scalar_t(1) : normal_scalar_t(0);
        aabb_planes[a*2].b = (a == 1) ? normal_scalar_t(1) : normal_scalar_t(0);
        aabb_planes[a*2].c = (a == 2) ? normal_scalar_t(1) : normal_scalar_t(0);
        aabb_planes[a*2].d = plane_d_t(-int64_t((&soup.bounds.min.x)[a]));
        // x[a] <= max[a]  →  -x[a] + max[a] >= 0
        aabb_planes[a*2+1].a = (a == 0) ? normal_scalar_t(-1) : normal_scalar_t(0);
        aabb_planes[a*2+1].b = (a == 1) ? normal_scalar_t(-1) : normal_scalar_t(0);
        aabb_planes[a*2+1].c = (a == 2) ? normal_scalar_t(-1) : normal_scalar_t(0);
        aabb_planes[a*2+1].d = plane_d_t(int64_t((&soup.bounds.max.x)[a]));
    }

    int filtered = 0;
    for (auto& cp : classified)
    {
        if (cp.classification == -1)
            cp.polygon = cp.polygon.inverted();

        // Exact validity check: all vertices must be within root AABB
        bool valid = true;
        for (int vi = 0; vi < cp.polygon.vertex_count() && valid; vi++)
        {
            auto pt = cp.polygon.vertex(vi);
            if (!pt.is_valid()) { valid = false; break; }

            (void)vi; // vertex check is done by classify below

            for (int p = 0; p < 6; p++)
            {
                if (exact_classify(pt, aabb_planes[p]) < 0)
                { valid = false; break; }
            }
        }
        if (!valid) { filtered++; continue; }

        result.output.polygons.push_back(std::move(cp.polygon));
        result.classifications.push_back(cp.classification);
    }
    if (filtered > 0)
        std::fprintf(stderr, "[ember] filtered %d/%zu degenerate polygons\n",
                     filtered, classified.size());

    return result;
}

// Convenience wrappers for common operations
inline BooleanResult boolean_union(InputMesh const& a, InputMesh const& b, EmberConfig const& config = {})
{
    return boolean_operation({a, b}, BooleanOp::Union, config);
}

inline BooleanResult boolean_intersection(InputMesh const& a, InputMesh const& b, EmberConfig const& config = {})
{
    return boolean_operation({a, b}, BooleanOp::Intersection, config);
}

inline BooleanResult boolean_difference(InputMesh const& a, InputMesh const& b, EmberConfig const& config = {})
{
    return boolean_operation({a, b}, BooleanOp::Difference, config);
}

// Extract output polygons as floating-point vertices (in original coordinate space)
inline std::vector<OutputPolygon> extract_output(BooleanResult const& result)
{
    std::vector<OutputPolygon> out;
    out.reserve(result.output.polygons.size());
    auto const& xf = result.transform;

    for (auto const& poly : result.output.polygons)
    {
        OutputPolygon op;
        op.source_mesh = poly.mesh_index;
        op.source_polygon = poly.polygon_index;

        int n = poly.vertex_count();
        op.vertices.resize(n);
        for (int i = 0; i < n; i++)
        {
            auto dp = poly.vertex_dpos(i);
            op.vertices[i] = xf.to_float(dp.x, dp.y, dp.z);
        }

        out.push_back(std::move(op));
    }

    return out;
}

inline TriangleSoup triangulate_output(BooleanResult const& result)
{
    TriangleSoup soup;
    auto const& xf = result.transform;

    for (auto const& poly : result.output.polygons)
    {
        int n = poly.vertex_count();
        if (n < 3) continue;

        // Fan triangulation from vertex 0
        int base = static_cast<int>(soup.vertices.size());

        // Add vertices (converted back to original coordinate space)
        for (int i = 0; i < n; i++)
        {
            auto dp = poly.vertex_dpos(i);
            soup.vertices.push_back(xf.to_float(dp.x, dp.y, dp.z));
        }

        // Add fan triangles
        for (int i = 1; i < n - 1; i++)
        {
            soup.triangles.push_back({base, base + i, base + i + 1});
        }
    }

    return soup;
}

// Triangulate and resolve T-junctions for manifold output
inline TriangleSoup triangulate_and_resolve(BooleanResult const& result)
{
    return resolve_tjunctions(triangulate_output(result));
}

// Write output as OBJ format string
inline std::string to_obj(TriangleSoup const& soup)
{
    std::string obj;
    obj.reserve(soup.vertices.size() * 40 + soup.triangles.size() * 20);

    obj += "# EMBER Integer Exact CSG output\n";

    for (auto const& v : soup.vertices)
    {
        obj += "v ";
        obj += std::to_string(v.x) + " ";
        obj += std::to_string(v.y) + " ";
        obj += std::to_string(v.z) + "\n";
    }

    for (auto const& tri : soup.triangles)
    {
        obj += "f ";
        obj += std::to_string(tri[0] + 1) + " ";
        obj += std::to_string(tri[1] + 1) + " ";
        obj += std::to_string(tri[2] + 1) + "\n";
    }

    return obj;
}

} // namespace ember
