// Diagnose non-manifold edges in EMBER sphere output.
// Find boundary edges (shared by only 1 triangle) and trace back to
// the source polygons to understand why they don't match up.

#include <ember/ember.hh>
#include <ember/resolve_tjunctions.hh>
#include <cstdio>
#include <cmath>
#include <map>
#include <set>
#include <algorithm>

struct Edge { int a, b; };
bool operator<(Edge const& l, Edge const& r) {
    return std::tie(l.a, l.b) < std::tie(r.a, r.b);
}
Edge make_edge(int a, int b) { return a < b ? Edge{a, b} : Edge{b, a}; }

int main()
{
    // Single sphere - roundtrip should be manifold
    auto mesh_a = ember::load_obj("sphere_a.obj", true, true);
    if (mesh_a.triangles.empty()) { std::printf("sphere_a.obj not found\n"); return 1; }

    ember::EmberConfig cfg;
    cfg.assume_nsi = true;
    cfg.assume_nnc = true;

    // Test: two non-overlapping spheres. Union should be both spheres intact.
    // sphere_a at (0,0,0), sphere_b at (1,1,1) with radius 1 → they DO overlap.
    // So for a non-overlapping test, use the actual sphere_b.
    // Actually, let's test the ACTUAL two-sphere union since that's what matters.
    auto mesh_b = ember::load_obj("sphere_b.obj", true, true);
    if (mesh_b.triangles.empty()) { std::printf("sphere_b.obj not found\n"); return 1; }

    auto result = ember::boolean_union(mesh_a, mesh_b, cfg);
    auto raw_soup = ember::triangulate_output(result);
    auto soup = ember::resolve_tjunctions(raw_soup);

    std::printf("Before T-junction resolve: %zu verts, %zu tris\n",
        raw_soup.vertices.size(), raw_soup.triangles.size());

    std::printf("Single sphere passthrough: %zu verts, %zu tris\n",
        soup.vertices.size(), soup.triangles.size());

    // Merge vertices within tolerance
    double tol = 1e-4;
    struct V3 { double x, y, z; };
    auto vkey = [&](V3 const& v) -> std::tuple<int,int,int> {
        return {int(std::round(v.x / tol)), int(std::round(v.y / tol)), int(std::round(v.z / tol))};
    };

    std::map<std::tuple<int,int,int>, int> vert_map;
    std::vector<V3> merged_verts;
    auto get_vert = [&](ember::OutputVertex const& ov) -> int {
        V3 v{ov.x, ov.y, ov.z};
        auto k = vkey(v);
        auto it = vert_map.find(k);
        if (it != vert_map.end()) return it->second;
        int id = (int)merged_verts.size();
        vert_map[k] = id;
        merged_verts.push_back(v);
        return id;
    };

    std::vector<std::array<int,3>> merged_tris;
    for (auto const& t : soup.triangles)
    {
        int a = get_vert(soup.vertices[t[0]]);
        int b = get_vert(soup.vertices[t[1]]);
        int c = get_vert(soup.vertices[t[2]]);
        if (a != b && b != c && a != c)
            merged_tris.push_back({a, b, c});
    }

    std::printf("After merge: %zu verts, %zu tris\n", merged_verts.size(), merged_tris.size());

    // Count edge usage
    std::map<Edge, std::vector<int>> edge_faces;
    for (int fi = 0; fi < (int)merged_tris.size(); fi++)
    {
        auto& t = merged_tris[fi];
        edge_faces[make_edge(t[0], t[1])].push_back(fi);
        edge_faces[make_edge(t[1], t[2])].push_back(fi);
        edge_faces[make_edge(t[2], t[0])].push_back(fi);
    }

    int boundary = 0, manifold = 0, nonmanifold = 0;
    for (auto& [e, faces] : edge_faces)
    {
        if (faces.size() == 1) boundary++;
        else if (faces.size() == 2) manifold++;
        else nonmanifold++;
    }

    std::printf("Edges: %zu total, %d boundary, %d manifold, %d non-manifold\n",
        edge_faces.size(), boundary, manifold, nonmanifold);

    // Write the resolved mesh
    {
        auto obj = ember::to_obj(soup);
        if (auto* f = std::fopen("sphere_union_resolved.obj", "w"))
        { std::fputs(obj.c_str(), f); std::fclose(f); }
        std::printf("Written sphere_union_resolved.obj\n");
    }

    if (boundary == 0 && nonmanifold == 0)
    {
        std::printf("MANIFOLD!\n");
        return 0;
    }

    // Show boundary edges and their incident faces
    std::printf("\nBoundary edges (first 10):\n");
    int shown = 0;
    for (auto& [e, faces] : edge_faces)
    {
        if (faces.size() != 1) continue;
        if (shown++ >= 10) break;

        auto& va = merged_verts[e.a];
        auto& vb = merged_verts[e.b];
        double len = std::sqrt(
            (va.x-vb.x)*(va.x-vb.x) +
            (va.y-vb.y)*(va.y-vb.y) +
            (va.z-vb.z)*(va.z-vb.z));

        std::printf("  edge v%d-v%d: (%.6f,%.6f,%.6f)-(%.6f,%.6f,%.6f) len=%.6f face=%d\n",
            e.a, e.b, va.x, va.y, va.z, vb.x, vb.y, vb.z, len, faces[0]);

        // Find nearby edges that SHOULD match but don't
        for (auto& [e2, f2] : edge_faces)
        {
            if (e2.a == e.a && e2.b == e.b) continue;
            if (f2.size() != 1) continue;

            // Check if any vertex of e2 is close to any vertex of e
            auto& v2a = merged_verts[e2.a];
            auto& v2b = merged_verts[e2.b];

            auto dist = [](V3 const& a, V3 const& b) {
                return std::sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
            };

            double d_aa = dist(va, v2a), d_ab = dist(va, v2b);
            double d_ba = dist(vb, v2a), d_bb = dist(vb, v2b);
            double min_d = std::min({d_aa, d_ab, d_ba, d_bb});

            if (min_d < 0.01 && min_d > tol)
            {
                std::printf("    NEAR-MISS with v%d-v%d: dist=%.8f (gap between BSP vertices)\n",
                    e2.a, e2.b, min_d);
                break;
            }
        }
    }

    // Check: how many source polygons have > 3 edges (BSP-split)?
    int original = 0, split = 0;
    for (auto const& p : result.output.polygons)
    {
        if (p.vertex_count() == 3) original++;
        else split++;
    }
    std::printf("\nOutput polygons: %d original triangles, %d BSP-split (>3 edges)\n", original, split);

    return 1;
}
