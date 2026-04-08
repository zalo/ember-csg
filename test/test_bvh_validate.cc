// Validate that BVH-accelerated queries produce IDENTICAL results to brute-force.
// Tests: trace_segment, classify_leaf_polygon
// Any difference is a bug.

#include <ember/ember.hh>
#include <cstdio>

int main()
{
    auto mesh_a = ember::load_obj("sphere_a.obj", true, true);
    auto mesh_b = ember::load_obj("sphere_b.obj", true, true);
    if (mesh_a.triangles.empty()) { std::printf("sphere_a.obj not found\n"); return 1; }

    auto soup = ember::prepare_input({mesh_a, mesh_b});
    auto const& polys = soup.polygons;
    int n = (int)polys.size();
    std::printf("Polygons: %d\n", n);

    // Build BVH
    ember::BVH bvh;
    bvh.build(polys);
    std::printf("BVH: %zu nodes\n", bvh.nodes.size());

    // Test 1: trace_segment BVH vs brute-force
    std::printf("\n=== Test 1: trace_segment BVH vs brute-force ===\n");
    int seg_tests = 0, seg_mismatches = 0;

    // Test from AABB corner to various polygon centers
    ember::pos_t ref = soup.bounds.min;
    ember::WNV ref_wnv(2, 0);

    for (int i = 0; i < std::min(n, 100); i++)
    {
        ember::pos_t interior;
        ember::ConvexPolygon tmp;
        tmp.support = polys[i].support;
        tmp.edges = polys[i].edges;
        if (!ember::find_interior_point(tmp, interior)) continue;

        // Test X segment
        ember::pos_t p1(interior.x, ref.y, ref.z);
        auto bf_result = ember::trace_segment(ref, p1, ref_wnv, polys, polys[i].support, polys[i].mesh_index);
        auto bvh_result = ember::trace_segment(ref, p1, ref_wnv, polys, polys[i].support, polys[i].mesh_index, &bvh);

        seg_tests++;
        if (bf_result != bvh_result)
        {
            seg_mismatches++;
            if (seg_mismatches <= 3)
                std::printf("  MISMATCH seg X poly %d: bf=(%d,%d) bvh=(%d,%d)\n",
                    i, bf_result[0], bf_result[1], bvh_result[0], bvh_result[1]);
        }

        // Test Y segment
        ember::pos_t p2(interior.x, interior.y, ref.z);
        auto bf2 = ember::trace_segment(p1, p2, bf_result, polys, polys[i].support, polys[i].mesh_index);
        auto bvh2 = ember::trace_segment(p1, p2, bvh_result, polys, polys[i].support, polys[i].mesh_index, &bvh);
        seg_tests++;
        if (bf2 != bvh2)
        {
            seg_mismatches++;
            if (seg_mismatches <= 3)
                std::printf("  MISMATCH seg Y poly %d: bf=(%d,%d) bvh=(%d,%d)\n",
                    i, bf2[0], bf2[1], bvh2[0], bvh2[1]);
        }
    }
    std::printf("  trace_segment: %d tests, %d mismatches\n", seg_tests, seg_mismatches);

    // Test 2: classify_leaf_polygon BVH vs brute-force
    std::printf("\n=== Test 2: classify_leaf_polygon BVH vs brute-force ===\n");
    int cls_tests = 0, cls_mismatches = 0;

    for (int i = 0; i < std::min(n, 200); i++)
    {
        auto const& poly = polys[i];
        auto bf_wf = ember::classify_leaf_polygon(
            poly.support, poly.edges, soup.bounds.min, ref_wnv,
            polys, soup.bounds, poly.delta_w, poly.mesh_index);
        auto bvh_wf = ember::classify_leaf_polygon(
            poly.support, poly.edges, soup.bounds.min, ref_wnv,
            polys, soup.bounds, poly.delta_w, poly.mesh_index, &bvh);

        cls_tests++;
        if (bf_wf != bvh_wf)
        {
            cls_mismatches++;
            if (cls_mismatches <= 5)
                std::printf("  MISMATCH poly %d (mesh %d): bf=(%d,%d) bvh=(%d,%d)\n",
                    i, poly.mesh_index, bf_wf[0], bf_wf[1], bvh_wf[0], bvh_wf[1]);
        }
    }
    std::printf("  classify_leaf_polygon: %d tests, %d mismatches\n", cls_tests, cls_mismatches);

    // Test 1b: Validate dual-BVH finds ALL overlapping cross-mesh pairs
    std::printf("\n=== Test 1b: dual-BVH pair completeness ===\n");
    {
        // Build per-mesh BVHs (same as process_leaf does)
        std::vector<int> mesh0_ids, mesh1_ids;
        for (int i = 0; i < n; i++)
        {
            if (polys[i].mesh_index == 0) mesh0_ids.push_back(i);
            else mesh1_ids.push_back(i);
        }
        std::vector<ember::ConvexPolygon> mesh0_polys, mesh1_polys;
        for (int i : mesh0_ids) mesh0_polys.push_back(polys[i]);
        for (int i : mesh1_ids) mesh1_polys.push_back(polys[i]);

        ember::BVH bvh0, bvh1;
        bvh0.build(mesh0_polys);
        bvh1.build(mesh1_polys);

        // BVH pairs
        std::set<std::pair<int,int>> bvh_pairs;
        bvh0.intersect_pairs(bvh1, [&](int li, int lj) {
            bvh_pairs.insert({mesh0_ids[li], mesh1_ids[lj]});
        });

        // Brute-force pairs (check AABB overlap)
        std::set<std::pair<int,int>> bf_pairs;
        for (int i : mesh0_ids)
        {
            for (int j : mesh1_ids)
            {
                if (bvh.prim_aabbs[i].intersects(bvh.prim_aabbs[j]))
                    bf_pairs.insert({i, j});
            }
        }

        // Count actual geometric intersections in each set
        int bvh_geom = 0, bf_geom = 0;
        for (auto [i, j] : bvh_pairs)
        {
            auto isect = ember::intersect_polygons(polys[i], polys[j], j);
            if (isect.type == ember::PairwiseIntersection::Type::Segment) bvh_geom++;
        }
        for (auto [i, j] : bf_pairs)
        {
            auto isect = ember::intersect_polygons(polys[i], polys[j], j);
            if (isect.type == ember::PairwiseIntersection::Type::Segment) bf_geom++;
        }

        int missing = 0;
        for (auto& p : bf_pairs)
            if (bvh_pairs.find(p) == bvh_pairs.end()) missing++;

        std::printf("  BVH pairs: %zu (geom intersections: %d)\n", bvh_pairs.size(), bvh_geom);
        std::printf("  Brute-force pairs: %zu (geom intersections: %d)\n", bf_pairs.size(), bf_geom);
        std::printf("  Missing from BVH: %d\n", missing);
    }

    // Test 2b: Check WNV distribution for intersection
    std::printf("\n=== Test 2b: WNV distribution ===\n");
    auto ind_isect = ember::make_indicator(ember::BooleanOp::Intersection, 2);
    int cls_dist[3] = {}; // emit+1, skip, emit-1
    int wnv_dist[2][2] = {};
    for (int i = 0; i < std::min(n, 256); i++)
    {
        auto const& poly = polys[i];
        auto wf = ember::classify_leaf_polygon(
            poly.support, poly.edges, soup.bounds.min, ref_wnv,
            polys, soup.bounds, poly.delta_w, poly.mesh_index, &bvh);
        auto wb = ember::propagate_wnv(wf, 1, poly.delta_w);
        int cls = ember::classify_polygon_output(wf, wb, ind_isect);
        cls_dist[cls + 1]++;
        int w0 = std::clamp(wf[0], 0, 1), w1 = std::clamp(wf[1], 0, 1);
        wnv_dist[w0][w1]++;
        if (i < 5)
            std::printf("  poly %d (mesh %d): wf=(%d,%d) wb=(%d,%d) cls=%d\n",
                i, poly.mesh_index, wf[0], wf[1], wb[0], wb[1], cls);
    }
    std::printf("  WNV: (0,0)=%d (0,1)=%d (1,0)=%d (1,1)=%d\n",
        wnv_dist[0][0], wnv_dist[0][1], wnv_dist[1][0], wnv_dist[1][1]);
    std::printf("  Classification: emit+1=%d skip=%d emit-1=%d\n", cls_dist[2], cls_dist[1], cls_dist[0]);

    // Test 3: Full boolean result BVH vs brute-force
    std::printf("\n=== Test 3: Full boolean union result ===\n");
    ember::EmberConfig cfg;
    cfg.assume_nsi = true;
    cfg.assume_nnc = true;

    auto result = ember::boolean_union(mesh_a, mesh_b, cfg);
    auto tsoup = ember::triangulate_and_resolve(result);
    double vol = 0;
    for (auto& t : tsoup.triangles)
    {
        auto& v0 = tsoup.vertices[t[0]]; auto& v1 = tsoup.vertices[t[1]]; auto& v2 = tsoup.vertices[t[2]];
        vol += (v0.x*(v1.y*v2.z-v1.z*v2.y)+v0.y*(v1.z*v2.x-v1.x*v2.z)+v0.z*(v1.x*v2.y-v1.y*v2.x))/6;
    }
    std::printf("  Union: %zu polys, %zu tris, vol=%.4f\n",
        result.output.polygons.size(), tsoup.triangles.size(), std::abs(vol));

    auto result_i = ember::boolean_intersection(mesh_a, mesh_b, cfg);
    auto tsoup_i = ember::triangulate_and_resolve(result_i);
    double vol_i = 0;
    for (auto& t : tsoup_i.triangles)
    {
        auto& v0 = tsoup_i.vertices[t[0]]; auto& v1 = tsoup_i.vertices[t[1]]; auto& v2 = tsoup_i.vertices[t[2]];
        vol_i += (v0.x*(v1.y*v2.z-v1.z*v2.y)+v0.y*(v1.z*v2.x-v1.x*v2.z)+v0.z*(v1.x*v2.y-v1.y*v2.x))/6;
    }
    std::printf("  Intersection: %zu polys, %zu tris, vol=%.4f\n",
        result_i.output.polygons.size(), tsoup_i.triangles.size(), std::abs(vol_i));

    int total_errors = seg_mismatches + cls_mismatches;
    std::printf("\nTotal errors: %d\n", total_errors);
    return total_errors > 0 ? 1 : 0;
}
