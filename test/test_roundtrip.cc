// Test: roundtrip a single mesh through EMBER's internal representation.
// Load OBJ → prepare_input (scale to integers, build plane-based polygons) → export OBJ
// No CSG, no BSP, no subdivision. Just test the representation fidelity.

#include <ember/ember.hh>
#include <cstdio>
#include <cmath>

int main()
{
    // Load the manifold3d sphere
    auto mesh = ember::load_obj("sphere_a.obj", true, true);
    std::printf("Loaded: %zu verts, %zu tris\n", mesh.positions.size(), mesh.triangles.size());

    if (mesh.triangles.empty())
    {
        std::printf("ERROR: sphere_a.obj not found. Run: python3 ../../compare/export_spheres.py\n");
        return 1;
    }

    // Convert to polygon soup (integer plane-based representation)
    auto soup = ember::prepare_input({mesh});
    std::printf("Polygon soup: %zu polygons\n", soup.polygons.size());
    std::printf("Transform: scale=%.4f center=(%.6f, %.6f, %.6f)\n",
        soup.transform.scale, soup.transform.center_x,
        soup.transform.center_y, soup.transform.center_z);
    std::printf("AABB: (%d,%d,%d) to (%d,%d,%d)\n",
        int(soup.bounds.min.x), int(soup.bounds.min.y), int(soup.bounds.min.z),
        int(soup.bounds.max.x), int(soup.bounds.max.y), int(soup.bounds.max.z));

    // Validate each polygon: check vertices are sane
    int bad_verts = 0, bad_polys = 0, total_verts = 0;
    double max_err = 0;

    for (size_t pi = 0; pi < soup.polygons.size(); pi++)
    {
        auto const& poly = soup.polygons[pi];
        int nv = poly.vertex_count();

        // Get the original triangle vertices
        auto const& tri = mesh.triangles[pi];
        tg::dpos3 orig[3] = {
            mesh.positions[tri.v0],
            mesh.positions[tri.v1],
            mesh.positions[tri.v2]
        };

        bool poly_bad = false;
        for (int vi = 0; vi < nv; vi++)
        {
            auto dp = poly.vertex_dpos(vi);
            // Convert back to original space
            auto fp = soup.transform.to_float(dp.x, dp.y, dp.z);

            // Find closest original vertex
            double min_dist = 1e30;
            for (int oi = 0; oi < 3; oi++)
            {
                double dx = fp.x - orig[oi].x;
                double dy = fp.y - orig[oi].y;
                double dz = fp.z - orig[oi].z;
                double d = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (d < min_dist) min_dist = d;
            }

            if (min_dist > max_err) max_err = min_dist;
            total_verts++;

            if (min_dist > 0.01)
            {
                bad_verts++;
                poly_bad = true;
                if (bad_verts <= 5)
                {
                    std::printf("  BAD vertex: poly %zu v%d: dpos=(%.2f,%.2f,%.2f) "
                        "float=(%.6f,%.6f,%.6f) dist=%.6f\n",
                        pi, vi, dp.x, dp.y, dp.z, fp.x, fp.y, fp.z, min_dist);
                    std::printf("    original: (%.6f,%.6f,%.6f) (%.6f,%.6f,%.6f) (%.6f,%.6f,%.6f)\n",
                        orig[0].x, orig[0].y, orig[0].z,
                        orig[1].x, orig[1].y, orig[1].z,
                        orig[2].x, orig[2].y, orig[2].z);

                    // Also print the raw point4
                    auto pt = poly.vertex(vi);
                    std::printf("    point4 w_valid=%d\n", pt.is_valid() ? 1 : 0);
                }
            }
        }
        if (poly_bad) bad_polys++;
    }

    // Direct test: check exact_classify on a known point
    if (soup.polygons.size() > 4)
    {
        auto& poly4 = soup.polygons[4];
        auto pt4 = poly4.vertex(0);
        std::printf("\nDirect test poly 4 vertex 0:\n");
        std::printf("  point4 valid=%d\n", pt4.is_valid() ? 1 : 0);

        // Check what exact_classify gives for axis test planes at known positions
        // The original vertex should be at the integer-scaled position of the triangle vertex
        auto const& tri4 = mesh.triangles[4];
        auto orig_v = mesh.positions[tri4.v1]; // vertex(0) = p1 for make_triangle convention
        double scale = soup.transform.scale;
        double cx = soup.transform.center_x;
        int32_t expected_x = int32_t(std::round((orig_v.x - cx) * scale));
        std::printf("  expected integer x: %d\n", expected_x);

        // Test classify at expected_x and expected_x-1
        ember::plane_t test_plane{};
        test_plane.a = ember::normal_scalar_t(1);

        test_plane.d = ember::plane_d_t(-int64_t(expected_x));
        auto c_at = ember::exact_classify(pt4, test_plane);
        test_plane.d = ember::plane_d_t(-int64_t(expected_x - 1));
        auto c_below = ember::exact_classify(pt4, test_plane);
        test_plane.d = ember::plane_d_t(-int64_t(expected_x + 1));
        auto c_above = ember::exact_classify(pt4, test_plane);

        std::printf("  classify at x=%d: %d (expected 0)\n", expected_x, c_at);
        std::printf("  classify at x=%d: %d (expected +1)\n", expected_x - 1, c_below);
        std::printf("  classify at x=%d: %d (expected -1)\n", expected_x + 1, c_above);

        auto dpos = poly4.vertex_dpos(0);
        std::printf("  vertex_dpos: (%.2f, %.2f, %.2f)\n", dpos.x, dpos.y, dpos.z);
        std::printf("  expected dpos: (%.2f, %.2f, %.2f)\n",
            (orig_v.x - cx) * scale, (orig_v.y - soup.transform.center_y) * scale,
            (orig_v.z - soup.transform.center_z) * scale);
    }

    std::printf("\nRoundtrip results:\n");
    std::printf("  Total vertices: %d\n", total_verts);
    std::printf("  Bad vertices (error > 0.01): %d (%.1f%%)\n",
        bad_verts, 100.0 * bad_verts / std::max(total_verts, 1));
    std::printf("  Bad polygons: %d / %zu (%.1f%%)\n",
        bad_polys, soup.polygons.size(), 100.0 * bad_polys / std::max(soup.polygons.size(), size_t(1)));
    std::printf("  Max vertex error: %.8f\n", max_err);

    // Write roundtrip output
    ember::BooleanResult fake_result;
    fake_result.output = soup;
    fake_result.transform = soup.transform;
    auto tri_soup = ember::triangulate_output(fake_result);
    auto obj = ember::to_obj(tri_soup);

    if (auto* f = std::fopen("roundtrip_sphere.obj", "w"))
    {
        std::fputs(obj.c_str(), f);
        std::fclose(f);
        std::printf("\nWritten roundtrip_sphere.obj (%zu verts, %zu tris)\n",
            tri_soup.vertices.size(), tri_soup.triangles.size());
    }

    // Compute volume of roundtrip output
    double vol = 0;
    for (auto const& t : tri_soup.triangles)
    {
        auto& v0 = tri_soup.vertices[t[0]];
        auto& v1 = tri_soup.vertices[t[1]];
        auto& v2 = tri_soup.vertices[t[2]];
        double cx = v1.y*v2.z - v1.z*v2.y;
        double cy = v1.z*v2.x - v1.x*v2.z;
        double cz = v1.x*v2.y - v1.y*v2.x;
        vol += (v0.x*cx + v0.y*cy + v0.z*cz) / 6.0;
    }
    std::printf("  Roundtrip volume: %.6f (expected ~4.095 for unit sphere)\n", std::abs(vol));

    return (bad_verts > 0) ? 1 : 0;
}
