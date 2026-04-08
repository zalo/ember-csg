// Validate WNV classification by comparing the raycast-based method against
// a brute-force reference: cast 7 rays in random directions and take the
// majority vote for each mesh's inside/outside status.

#include <ember/ember.hh>
#include <cstdio>
#include <cmath>
#include <cstdlib>

// Brute-force: Moller-Trumbore ray-triangle test in double precision
static double ray_tri(double ox, double oy, double oz,
                      double dx, double dy, double dz,
                      double v0x, double v0y, double v0z,
                      double v1x, double v1y, double v1z,
                      double v2x, double v2y, double v2z)
{
    double e1x = v1x-v0x, e1y = v1y-v0y, e1z = v1z-v0z;
    double e2x = v2x-v0x, e2y = v2y-v0y, e2z = v2z-v0z;
    double hx = dy*e2z - dz*e2y, hy = dz*e2x - dx*e2z, hz = dx*e2y - dy*e2x;
    double a = e1x*hx + e1y*hy + e1z*hz;
    if (std::abs(a) < 1e-30) return -1;
    double f = 1.0/a;
    double sx = ox-v0x, sy = oy-v0y, sz = oz-v0z;
    double u = f * (sx*hx + sy*hy + sz*hz);
    if (u < -1e-8 || u > 1.0+1e-8) return -1;
    double qx = sy*e1z - sz*e1y, qy = sz*e1x - sx*e1z, qz = sx*e1y - sy*e1x;
    double v = f * (dx*qx + dy*qy + dz*qz);
    if (v < -1e-8 || u+v > 1.0+1e-8) return -1;
    return f * (e2x*qx + e2y*qy + e2z*qz);
}

// Brute-force point-in-mesh: cast a ray and count parity of crossings
// with triangles of a specific mesh_index. Uses original float positions.
static int brute_force_inside(double px, double py, double pz,
                               double dx, double dy, double dz,
                               ember::InputMesh const& mesh)
{
    int crossings = 0;
    for (auto const& tri : mesh.triangles)
    {
        auto& v0 = mesh.positions[tri.v0];
        auto& v1 = mesh.positions[tri.v1];
        auto& v2 = mesh.positions[tri.v2];
        double t = ray_tri(px, py, pz, dx, dy, dz,
                           v0.x, v0.y, v0.z, v1.x, v1.y, v1.z, v2.x, v2.y, v2.z);
        if (t > 1e-6) crossings++;
    }
    return crossings & 1; // parity: 1 = inside, 0 = outside
}

// Robust brute-force: cast 7 rays in different directions, majority vote
static int robust_inside(double px, double py, double pz, ember::InputMesh const& mesh)
{
    static const double dirs[][3] = {
        { 0.8507, 0.5257, 0.0353},
        {-0.3891, 0.9213, 0.0123},
        { 0.1234,-0.4567, 0.8809},
        {-0.7071,-0.7071, 0.0100},
        { 0.5774, 0.5774, 0.5774},
        {-0.2345, 0.1234,-0.9642},
        { 0.9511, 0.0312,-0.3073},
    };
    int votes = 0;
    for (auto& d : dirs)
        votes += brute_force_inside(px, py, pz, d[0], d[1], d[2], mesh);
    return (votes > 3) ? 1 : 0; // majority
}

int main()
{
    auto mesh_a = ember::load_obj("sphere_a.obj", true, true);
    auto mesh_b = ember::load_obj("sphere_b.obj", true, true);
    if (mesh_a.triangles.empty() || mesh_b.triangles.empty())
    {
        std::printf("ERROR: sphere OBJ not found\n");
        return 1;
    }

    // Prepare the polygon soup (same as boolean_operation does)
    std::vector<ember::InputMesh> meshes = {mesh_a, mesh_b};
    auto soup = ember::prepare_input(meshes);
    auto const& xf = soup.transform;

    std::printf("Soup: %zu polygons, scale=%.1f center=(%.3f,%.3f,%.3f)\n",
        soup.polygons.size(), xf.scale, xf.center_x, xf.center_y, xf.center_z);

    // Direct sanity check: test a known-inside point in both spaces
    {
        // Center of sphere A is at (0,0,0) in float space
        // In integer space: (0 - 0.5) * scale = -11184809
        double int_center_x = (0.0 - xf.center_x) * xf.scale;
        double int_center_y = (0.0 - xf.center_y) * xf.scale;
        double int_center_z = (0.0 - xf.center_z) * xf.scale;
        std::printf("Sphere A center in integer space: (%.0f, %.0f, %.0f)\n",
            int_center_x, int_center_y, int_center_z);

        // Test: is center of sphere A inside sphere A?
        int bf_at_center = robust_inside(0, 0, 0, mesh_a);
        std::printf("Brute-force at sphere A center (float): inside_A=%d\n", bf_at_center);

        // EMBER ray cast at center in integer space
        auto ember_at_center = ember::point_in_meshes_robust(
            int_center_x, int_center_y, int_center_z, soup.polygons, 2);
        std::printf("EMBER at sphere A center (integer): (%d, %d)\n",
            ember_at_center[0], ember_at_center[1]);

        // Also test: cast rays in float space against polygon vertex_dpos positions
        // (this checks if the polygon geometry is correct in float space)
        int crossings_float = 0;
        double dir[3] = {0.85065, 0.52573, 0.03532};
        for (auto const& poly : soup.polygons)
        {
            if (poly.mesh_index != 0) continue;
            auto v0 = poly.vertex_dpos(0);
            auto v1 = poly.vertex_dpos(1);
            auto v2 = poly.vertex_dpos(2);
            // Convert from integer space to float space for comparison
            auto fv0 = xf.to_float(v0.x, v0.y, v0.z);
            auto fv1 = xf.to_float(v1.x, v1.y, v1.z);
            auto fv2 = xf.to_float(v2.x, v2.y, v2.z);
            double t = ray_tri(0, 0, 0, dir[0], dir[1], dir[2],
                fv0.x, fv0.y, fv0.z, fv1.x, fv1.y, fv1.z, fv2.x, fv2.y, fv2.z);
            if (t > 1e-6) crossings_float++;
        }
        std::printf("Float-space ray crossings at (0,0,0) with EMBER polys: %d (parity=%d)\n",
            crossings_float, crossings_float & 1);

        // And in integer space directly:
        int crossings_int = 0;
        for (auto const& poly : soup.polygons)
        {
            if (poly.mesh_index != 0) continue;
            auto v0 = poly.vertex_dpos(0);
            auto v1 = poly.vertex_dpos(1);
            auto v2 = poly.vertex_dpos(2);
            double t = ray_tri(int_center_x, int_center_y, int_center_z,
                dir[0], dir[1], dir[2],
                v0.x, v0.y, v0.z, v1.x, v1.y, v1.z, v2.x, v2.y, v2.z);
            if (t > 1e-6) crossings_int++;
        }
        std::printf("Integer-space ray crossings at center with EMBER polys: %d (parity=%d)\n",
            crossings_int, crossings_int & 1);
    }

    // Direct test at poly 0's probe position
    {
        auto const& poly0 = soup.polygons[0];
        ember::pos_t probe;
        int probe_side;
        ember::find_probe_point(poly0, probe, probe_side);
        double ppx = double(probe.x), ppy = double(probe.y), ppz = double(probe.z);
        std::printf("\nPoly 0 probe integer: (%.0f, %.0f, %.0f) side=%d\n", ppx, ppy, ppz, probe_side);

        // Direct 7-ray test at this point in integer space
        auto wnv_int = ember::point_in_meshes_robust(ppx, ppy, ppz, soup.polygons, 2);
        std::printf("EMBER at poly 0 probe (int space): (%d, %d)\n", wnv_int[0], wnv_int[1]);

        // Same point in float space, test against original mesh
        auto fp = xf.to_float(ppx, ppy, ppz);
        int bf_at_probe = robust_inside(fp.x, fp.y, fp.z, mesh_a);
        std::printf("Brute-force at poly 0 probe (float space): inside_A=%d\n", bf_at_probe);
        std::printf("Probe float: (%.6f, %.6f, %.6f)\n", fp.x, fp.y, fp.z);

        // Single ray in integer space: detailed crossing count
        double dir[3] = {0.85065, 0.52573, 0.03532};
        int cross_m0 = 0;
        for (auto const& p : soup.polygons)
        {
            if (p.mesh_index != 0) continue;
            int nv = p.vertex_count();
            auto v0 = p.vertex_dpos(0);
            for (int k = 1; k < nv - 1; k++)
            {
                auto v1 = p.vertex_dpos(k);
                auto v2 = p.vertex_dpos(k+1);
                double t = ray_tri(ppx, ppy, ppz, dir[0], dir[1], dir[2],
                    v0.x, v0.y, v0.z, v1.x, v1.y, v1.z, v2.x, v2.y, v2.z);
                if (t > 1e-6) cross_m0++;
            }
        }
        std::printf("Single ray crossings mesh 0 at probe (int space): %d (parity=%d)\n",
            cross_m0, cross_m0 & 1);
    }

    // Test the WNV at various probe points and compare with brute-force
    // Sample probe points: polygon centers (offset along normal)
    int total = 0, mismatches = 0;
    int match_00 = 0, match_01 = 0, match_10 = 0, match_11 = 0;

    int num_to_test = std::min(int(soup.polygons.size()), 200);
    for (int pi = 0; pi < num_to_test; pi++)
    {
        auto const& poly = soup.polygons[pi];

        // Find probe point (same as classify_leaf_polygon)
        ember::pos_t probe;
        int probe_side;
        if (!ember::find_probe_point(poly, probe, probe_side))
            continue;

        // Convert probe to float space for brute-force
        auto fp = xf.to_float(double(probe.x), double(probe.y), double(probe.z));

        // Brute-force WNV at probe
        int bf_w0 = robust_inside(fp.x, fp.y, fp.z, mesh_a);
        int bf_w1 = robust_inside(fp.x, fp.y, fp.z, mesh_b);

        // EMBER raycast WNV at probe
        ember::WNV ember_wnv = ember::classify_probe_by_raycast(
            probe, 0, true, soup.bounds, soup.polygons, 2);

        int em_w0 = ember_wnv[0], em_w1 = ember_wnv[1];

        total++;
        bool ok = (bf_w0 == em_w0 && bf_w1 == em_w1);
        if (!ok) mismatches++;

        if (bf_w0 == 0 && bf_w1 == 0) match_00++;
        if (bf_w0 == 0 && bf_w1 == 1) match_01++;
        if (bf_w0 == 1 && bf_w1 == 0) match_10++;
        if (bf_w0 == 1 && bf_w1 == 1) match_11++;

        if (!ok && mismatches <= 10)
        {
            std::printf("  MISMATCH poly %d: brute=(%d,%d) ember=(%d,%d) probe=(%+.4f,%+.4f,%+.4f) side=%d\n",
                pi, bf_w0, bf_w1, em_w0, em_w1, fp.x, fp.y, fp.z, probe_side);
        }
    }

    std::printf("\nWNV validation: %d tested, %d mismatches (%.1f%%)\n",
        total, mismatches, 100.0 * mismatches / std::max(total, 1));
    std::printf("Brute-force WNV distribution: (0,0)=%d (0,1)=%d (1,0)=%d (1,1)=%d\n",
        match_00, match_01, match_10, match_11);

    // Also test: what does the FULL classify_leaf_polygon return vs brute-force?
    std::printf("\n=== Full classification comparison (union) ===\n");
    auto indicator = ember::make_indicator(ember::BooleanOp::Union, 2);
    int cls_match = 0, cls_total = 0, cls_mismatch = 0;
    for (int pi = 0; pi < num_to_test; pi++)
    {
        auto const& poly = soup.polygons[pi];

        // EMBER's full classification
        ember::WNV w_front = ember::classify_leaf_polygon(
            poly.support, poly.edges, soup.bounds.min, ember::WNV(2, 0),
            soup.polygons, soup.bounds, poly.delta_w);
        ember::WNV w_back = ember::propagate_wnv(w_front, 1, poly.delta_w);
        int ember_cls = ember::classify_polygon_output(w_front, w_back, indicator);

        // Brute-force classification at probe
        ember::pos_t probe;
        int probe_side;
        if (!ember::find_probe_point(poly, probe, probe_side))
            continue;
        auto fp = xf.to_float(double(probe.x), double(probe.y), double(probe.z));
        int bf_w0 = robust_inside(fp.x, fp.y, fp.z, mesh_a);
        int bf_w1 = robust_inside(fp.x, fp.y, fp.z, mesh_b);

        // Determine expected WNV at front side
        ember::WNV bf_wnv_at_probe = {bf_w0, bf_w1};
        ember::WNV bf_front, bf_back;
        if (probe_side > 0) {
            bf_front = bf_wnv_at_probe;
        } else {
            bf_front = bf_wnv_at_probe;
            for (size_t k = 0; k < bf_front.size() && k < poly.delta_w.size(); k++)
                bf_front[k] -= poly.delta_w[k];
        }
        bf_back = ember::propagate_wnv(bf_front, 1, poly.delta_w);
        int bf_cls = ember::classify_polygon_output(bf_front, bf_back, indicator);

        cls_total++;
        if (ember_cls == bf_cls) cls_match++;
        else {
            cls_mismatch++;
            if (cls_mismatch <= 5)
                std::printf("  CLS MISMATCH poly %d (mesh %d): ember=%d bf=%d "
                    "ember_wf=(%d,%d) bf_wf=(%d,%d) side=%d\n",
                    pi, poly.mesh_index, ember_cls, bf_cls,
                    w_front[0], w_front[1], bf_front[0], bf_front[1], probe_side);
        }
    }

    std::printf("\nClassification: %d/%d match (%.1f%%), %d mismatch\n",
        cls_match, cls_total, 100.0*cls_match/std::max(cls_total,1), cls_mismatch);

    return (mismatches > 0) ? 1 : 0;
}
