// Trace every polygon through the EMBER pipeline and identify exactly why
// each one is emitted or dropped. For a single-mesh "identity" boolean
// (union with empty mesh), every input polygon should be emitted.

#include <ember/ember.hh>
#include <cstdio>
#include <cmath>

int main()
{
    auto mesh_a = ember::load_obj("sphere_a.obj", true, true);
    if (mesh_a.triangles.empty()) { std::printf("sphere_a.obj not found\n"); return 1; }
    std::printf("Loaded sphere: %zu tris\n", mesh_a.triangles.size());

    // Test 1: Single mesh union with itself should output ALL triangles
    // (union of A with nothing = A)
    ember::EmberConfig config;
    config.assume_nsi = true;
    config.assume_nnc = true;

    // Create a tiny degenerate mesh for "B" so we have 2 meshes
    ember::InputMesh mesh_b;
    mesh_b.nsi = true;
    mesh_b.nnc = true;
    mesh_b.positions = {{100, 100, 100}, {100.001, 100, 100}, {100, 100.001, 100}};
    mesh_b.triangles = {{0, 1, 2}};

    auto result = ember::boolean_union(mesh_a, mesh_b, config);
    int mesh_a_out = 0, mesh_b_out = 0;
    for (auto const& p : result.output.polygons)
    {
        if (p.mesh_index == 0) mesh_a_out++;
        else mesh_b_out++;
    }
    std::printf("\nTest 1: Union(sphere, tiny_far_triangle)\n");
    std::printf("  Input: %zu from A, 1 from B\n", mesh_a.triangles.size());
    std::printf("  Output: %d from A, %d from B\n", mesh_a_out, mesh_b_out);
    std::printf("  Missing from A: %zu\n", mesh_a.triangles.size() - mesh_a_out);

    // Test 2: Trace classification for EVERY polygon of A
    // Manually run the pipeline stages to see where polygons are lost
    std::printf("\n=== Test 2: Manual pipeline trace ===\n");

    auto soup = ember::prepare_input({mesh_a, mesh_b});
    auto indicator = ember::make_indicator(ember::BooleanOp::Union, 2);

    // Count polygons per mesh
    int n_a = 0, n_b = 0;
    for (auto const& p : soup.polygons)
    {
        if (p.mesh_index == 0) n_a++;
        else n_b++;
    }
    std::printf("  Soup: %d from A, %d from B\n", n_a, n_b);

    // Manually classify each polygon (no subdivision, no BSP - just direct WNV)
    int emit_a = 0, skip_a = 0, emit_b = 0, skip_b = 0;
    int wrong_wnv = 0;

    for (size_t i = 0; i < soup.polygons.size(); i++)
    {
        auto const& poly = soup.polygons[i];

        // Build temp poly for find_probe_point
        ember::ConvexPolygon tmp;
        tmp.support = poly.support;
        tmp.edges = poly.edges;

        ember::pos_t probe;
        int probe_side;
        if (!ember::find_probe_point(tmp, probe, probe_side))
        {
            if (i < 5) std::printf("  poly %zu: probe FAILED\n", i);
            skip_a++;
            continue;
        }

        // Classify using the robust method
        ember::WNV w_front = ember::classify_leaf_polygon(
            poly.support, poly.edges, soup.bounds.min, ember::WNV(2, 0),
            soup.polygons, soup.bounds, poly.delta_w);
        ember::WNV w_back = ember::propagate_wnv(w_front, 1, poly.delta_w);
        int cls = ember::classify_polygon_output(w_front, w_back, indicator);

        // For union with a far-away tiny B: all A polygons should be emitted
        // because they're outside B, so w_front should be (0,0) or (1,0) but NOT (x,1)
        bool expected_emit = (poly.mesh_index == 0); // all A should emit

        if (poly.mesh_index == 0)
        {
            if (cls != 0) emit_a++;
            else
            {
                skip_a++;
                // Why was it skipped?
                if (skip_a <= 10)
                {
                    std::printf("  DROPPED A poly %zu: w_front=(%d,%d) w_back=(%d,%d) cls=%d side=%d\n",
                        i, w_front[0], w_front[1], w_back[0], w_back[1], cls, probe_side);

                    // Sanity: is this probe actually inside A?
                    // The probe should be right next to the polygon surface.
                    // For side>0 (front), if normal is outward, probe is outside A.
                    // For union: outside A + outside B → w_front=(0,0) → w_back=(1,0) → out→in → emit
                    // If w_front=(1,0) → w_back=(2,0) → in→in → skip. But WNV=2 is impossible!

                    // Check: what does classify at the sphere center give?
                    double cx = (0.0 - soup.transform.center_x) * soup.transform.scale;
                    double cy = (0.0 - soup.transform.center_y) * soup.transform.scale;
                    double cz = (0.0 - soup.transform.center_z) * soup.transform.scale;
                    auto wnv_center = ember::point_in_meshes_robust(cx, cy, cz, soup.polygons, 2);
                    if (skip_a <= 3)
                        std::printf("    (center of A: wnv=(%d,%d))\n", wnv_center[0], wnv_center[1]);

                    if (w_front[0] != 0 && w_front[1] == 0 && probe_side > 0)
                    {
                        // w_front=(1,0) with probe on front: means probe is INSIDE A
                        // but if normal is outward, probe should be OUTSIDE A
                        // This means the normal is INWARD!
                        wrong_wnv++;
                    }
                }
            }
        }
        else
        {
            if (cls != 0) emit_b++;
            else skip_b++;
        }
    }

    std::printf("\n  Classification results:\n");
    std::printf("    A: %d emitted, %d skipped (expected %d emitted)\n", emit_a, skip_a, n_a);
    std::printf("    B: %d emitted, %d skipped\n", emit_b, skip_b);
    std::printf("    Wrong WNV (probe inside own mesh): %d\n", wrong_wnv);

    // Test 3: Check normal directions
    std::printf("\n=== Test 3: Normal direction check ===\n");
    int inward = 0, outward = 0;
    for (size_t i = 0; i < soup.polygons.size() && i < 512; i++)
    {
        auto const& poly = soup.polygons[i];
        if (poly.mesh_index != 0) continue;

        // The support normal
        auto sn = poly.support.normal();
        // The polygon center in float space
        auto center = poly.center_of_gravity();
        auto fc = soup.transform.to_float(center.x, center.y, center.z);

        // For sphere at (0,0,0), outward normal has positive dot with position
        double dot = double(sn.x) * (fc.x) + double(sn.y) * (fc.y) + double(sn.z) * (fc.z);
        if (dot > 0) outward++;
        else inward++;
    }
    std::printf("  Sphere A normals: %d outward, %d inward\n", outward, inward);
    if (inward > outward)
        std::printf("  WARNING: majority of normals point INWARD!\n");

    // Test 4: Two-sphere union - trace dropped polygons
    std::printf("\n=== Test 4: Two-sphere union trace ===\n");
    auto mesh_b2 = ember::load_obj("sphere_b.obj", true, true);
    if (mesh_b2.triangles.empty()) { std::printf("sphere_b.obj not found\n"); return 1; }

    auto soup2 = ember::prepare_input({mesh_a, mesh_b2});
    auto ind_union = ember::make_indicator(ember::BooleanOp::Union, 2);
    auto ind_isect = ember::make_indicator(ember::BooleanOp::Intersection, 2);

    std::printf("  Soup: %zu total polygons\n", soup2.polygons.size());

    // Classify each INPUT polygon (before BSP) directly
    int cls_counts[2][3] = {}; // [mesh][emit/skip/fail]
    int wnv_dist[2][2][2] = {}; // [mesh][w0][w1] distribution

    for (size_t i = 0; i < soup2.polygons.size(); i++)
    {
        auto const& poly = soup2.polygons[i];
        int m = poly.mesh_index;

        ember::WNV w_front = ember::classify_leaf_polygon(
            poly.support, poly.edges, soup2.bounds.min, ember::WNV(2, 0),
            soup2.polygons, soup2.bounds, poly.delta_w);
        ember::WNV w_back = ember::propagate_wnv(w_front, 1, poly.delta_w);
        int cls_u = ember::classify_polygon_output(w_front, w_back, ind_union);

        if (cls_u != 0) cls_counts[m][0]++;
        else cls_counts[m][1]++;

        // Track WNV distribution (clamp to [0,1] for display)
        int w0 = std::clamp(w_front[0], 0, 1);
        int w1 = std::clamp(w_front[1], 0, 1);
        wnv_dist[m][w0][w1]++;
    }

    std::printf("  Union classification (no BSP, direct per-polygon):\n");
    std::printf("    Mesh A: %d emit, %d skip\n", cls_counts[0][0], cls_counts[0][1]);
    std::printf("    Mesh B: %d emit, %d skip\n", cls_counts[1][0], cls_counts[1][1]);

    std::printf("  w_front distribution:\n");
    for (int m = 0; m < 2; m++)
    {
        std::printf("    Mesh %d: (0,0)=%d (0,1)=%d (1,0)=%d (1,1)=%d\n",
            m, wnv_dist[m][0][0], wnv_dist[m][0][1], wnv_dist[m][1][0], wnv_dist[m][1][1]);
    }

    // For union: emit if w_front=(0,0)→w_back=(1,0) or (0,1) → out→in
    // Skip if w_front=(0,1)→w_back=(1,1) → in→in (inside other mesh)
    // or w_front=(1,0)→w_back=(2,0) → ???
    std::printf("\n  Expected for union:\n");
    std::printf("    Mesh A: w_front=(0,0)→emit, w_front=(0,1)→skip (inside B)\n");
    std::printf("    Mesh B: w_front=(0,0)→emit, w_front=(1,0)→skip (inside A)\n");

    // Test 5: Count BSP leaves per polygon to find over-splitting
    std::printf("\n=== Test 5: BSP leaf count per polygon ===\n");
    {
        int single_leaf = 0, multi_leaf = 0, max_leaves = 0;
        int total_leaves = 0;
        int total_segs = 0;

        for (size_t i = 0; i < soup2.polygons.size(); i++)
        {
            auto const& poly_i = soup2.polygons[i];
            ember::LocalBSP bsp;
            bsp.init(poly_i);
            int segs = 0;

            for (size_t j = 0; j < soup2.polygons.size(); j++)
            {
                if (i == j) continue;
                if (poly_i.no_self_intersections &&
                    soup2.polygons[j].mesh_index == poly_i.mesh_index)
                    continue;

                auto isect = ember::intersect_polygons(poly_i, soup2.polygons[j], (int)j);
                if (isect.type == ember::PairwiseIntersection::Type::Segment)
                {
                    bsp.add_segment(isect.segment.v0, isect.segment.v1, isect.segment.split_plane);
                    segs++;
                }
                else if (isect.type == ember::PairwiseIntersection::Type::Overlap)
                {
                    bsp.add_overlap(soup2.polygons[j], soup2.polygons[j].polygon_index);
                }
            }

            std::vector<ember::BSPLeaf*> leaves;
            bsp.collect_leaves(leaves);
            int nl = (int)leaves.size();
            total_leaves += nl;
            total_segs += segs;

            if (nl <= 1) single_leaf++;
            else multi_leaf++;
            if (nl > max_leaves) max_leaves = nl;

            if (nl > 5 && multi_leaf <= 3)
                std::printf("  poly %zu (mesh %d): %d leaves, %d segments\n",
                    i, poly_i.mesh_index, nl, segs);
        }

        std::printf("  Single-leaf polys: %d (no intersection)\n", single_leaf);
        std::printf("  Multi-leaf polys: %d (BSP-split)\n", multi_leaf);
        std::printf("  Max leaves per poly: %d\n", max_leaves);
        std::printf("  Total leaves: %d, total segments: %d\n", total_leaves, total_segs);
    }

    // Now also run the ACTUAL boolean and count output
    ember::EmberConfig cfg2;
    cfg2.assume_nsi = true;
    cfg2.assume_nnc = true;
    cfg2.leaf_threshold = 100000; // disable subdivision
    auto result2 = ember::boolean_union(mesh_a, mesh_b2, cfg2);
    int out_a = 0, out_b = 0;
    for (auto const& p : result2.output.polygons)
    {
        if (p.mesh_index == 0) out_a++;
        else out_b++;
    }
    std::printf("\n  Actual boolean output: %d from A, %d from B, %zu total\n",
        out_a, out_b, result2.output.polygons.size());
    std::printf("  Expected: ~400 from each (total ~800), with BSP splits at intersection\n");

    return 0;
}
