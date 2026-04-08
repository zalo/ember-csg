// Validate BVH: compare BVH-accelerated vs brute-force classification
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
    std::printf("BVH built: %d prims\n", bvh.prim_count);

    // Compare BVH vs brute-force classify for every polygon
    ember::WNV ref_wnv(2, 0);
    int match = 0, mismatch = 0;
    for (int i = 0; i < std::min(n, 200); i++)
    {
        auto const& poly = polys[i];
        auto bf = ember::classify_leaf_polygon(
            poly.support, poly.edges, soup.bounds.min, ref_wnv,
            polys, soup.bounds, poly.delta_w, poly.mesh_index, nullptr);
        auto bvh_r = ember::classify_leaf_polygon(
            poly.support, poly.edges, soup.bounds.min, ref_wnv,
            polys, soup.bounds, poly.delta_w, poly.mesh_index, &bvh);
        if (bf == bvh_r) match++;
        else {
            mismatch++;
            if (mismatch <= 5)
                std::printf("  MISMATCH poly %d: bf=(%d,%d) bvh=(%d,%d)\n",
                    i, bf[0], bf[1], bvh_r[0], bvh_r[1]);
        }
    }
    std::printf("classify: %d match, %d mismatch\n", match, mismatch);

    // Full boolean
    ember::EmberConfig cfg; cfg.assume_nsi=true; cfg.assume_nnc=true;
    for (auto [name, op] : {std::pair{"union", ember::BooleanOp::Union},
        {"intersection", ember::BooleanOp::Intersection},
        {"difference", ember::BooleanOp::Difference}}) {
        auto r = ember::boolean_operation({ma, mb}, op, cfg);
        auto s = ember::triangulate_and_resolve(r);
        double vol=0;
        for(auto&t:s.triangles){auto&v0=s.vertices[t[0]];auto&v1=s.vertices[t[1]];auto&v2=s.vertices[t[2]];
            vol+=(v0.x*(v1.y*v2.z-v1.z*v2.y)+v0.y*(v1.z*v2.x-v1.x*v2.z)+v0.z*(v1.x*v2.y-v1.y*v2.x))/6;}
        std::printf("  %s: %zu polys %zu tris vol=%.4f\n", name, r.output.polygons.size(), s.triangles.size(), std::abs(vol));
    }
    return mismatch > 0 ? 1 : 0;
}
