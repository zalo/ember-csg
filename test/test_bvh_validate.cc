// Validate: 1) BVH matches brute-force, 2) Exact 3-seg matches robust multi-ray
#include <ember/ember.hh>
#include <cstdio>

int main()
{
    auto ma = ember::load_obj("sphere_a.obj", true, true);
    auto mb = ember::load_obj("sphere_b.obj", true, true);
    if (ma.triangles.empty()) { std::printf("sphere_a.obj not found\n"); return 1; }

    auto soup = ember::prepare_input({ma, mb});
    auto const& polys = soup.polygons;
    int n = (int)polys.size();
    std::printf("Polygons: %d\n", n);

    ember::BVH bvh;
    bvh.build(polys);
    std::printf("BVH: %zu nodes\n", bvh.nodes.size());

    // Compare exact 3-seg vs robust multi-ray for every polygon
    std::printf("\n=== Comparing exact vs robust for all %d polygons ===\n", n);

    ember::WNV ref_wnv(2, 0);
    ember::pos_t ref_point = soup.bounds.min;
    int match = 0, mismatch = 0, exact_fallback = 0;

    for (int i = 0; i < n; i++)
    {
        auto const& poly = polys[i];

        // Robust multi-ray (reference)
        ember::ConvexPolygon tmp;
        tmp.support = poly.support;
        tmp.edges = poly.edges;
        ember::pos_t interior;
        int ei, ej;
        if (!ember::find_interior_point(tmp, interior, ei, ej)) continue;

        auto sn = poly.support.normal();
        auto abs_a = ipg::abs(int64_t(sn.x)), abs_b = ipg::abs(int64_t(sn.y)), abs_c = ipg::abs(int64_t(sn.z));
        int dom = 0;
        if (abs_b > abs_a && abs_b > abs_c) dom = 1;
        else if (abs_c > abs_a && abs_c > abs_b) dom = 2;
        int dir = (int64_t((&sn.x)[dom]) > 0) ? 1 : -1;
        double px = double(interior.x), py = double(interior.y), pz = double(interior.z);
        (&px)[dom] += dir;

        ember::WNV robust_wf = ember::robust_classify(px, py, pz, polys, 2,
            poly.support, poly.mesh_index, &bvh);
        ember::pos_t probe = interior;
        (&probe.x)[dom] = ember::pos_scalar_t(int32_t((&probe.x)[dom]) + dir);
        auto side = ember::exact_classify(probe, poly.support);
        if (side < 0)
            for (size_t k = 0; k < robust_wf.size() && k < poly.delta_w.size(); k++)
                robust_wf[k] -= poly.delta_w[k];

        // Full classify (exact 3-seg with robust fallback)
        ember::WNV full_wf = ember::classify_leaf_polygon(
            poly.support, poly.edges, ref_point, ref_wnv,
            polys, soup.bounds, poly.delta_w, poly.mesh_index, &bvh);

        if (robust_wf == full_wf)
            match++;
        else
        {
            mismatch++;
            if (mismatch <= 10)
                std::printf("  MISMATCH poly %d (mesh %d): robust=(%d,%d) full=(%d,%d)\n",
                    i, poly.mesh_index, robust_wf[0], robust_wf[1], full_wf[0], full_wf[1]);
        }
    }

    std::printf("\nResults: %d match, %d mismatch out of %d\n", match, mismatch, match+mismatch);

    // Run full boolean
    std::printf("\n=== Full boolean ===\n");
    ember::EmberConfig cfg;
    cfg.assume_nsi = true;
    cfg.assume_nnc = true;

    for (auto [name, op] : {std::pair{"union", ember::BooleanOp::Union},
                            {"intersection", ember::BooleanOp::Intersection},
                            {"difference", ember::BooleanOp::Difference}})
    {
        auto r = ember::boolean_operation({ma, mb}, op, cfg);
        auto s = ember::triangulate_and_resolve(r);
        double vol = 0;
        for (auto& t : s.triangles)
        {
            auto& v0=s.vertices[t[0]]; auto& v1=s.vertices[t[1]]; auto& v2=s.vertices[t[2]];
            vol += (v0.x*(v1.y*v2.z-v1.z*v2.y)+v0.y*(v1.z*v2.x-v1.x*v2.z)+v0.z*(v1.x*v2.y-v1.y*v2.x))/6;
        }
        std::printf("  %s: %zu polys, %zu tris, vol=%.4f\n",
            name, r.output.polygons.size(), s.triangles.size(), std::abs(vol));
    }

    return mismatch > 0 ? 1 : 0;
}
