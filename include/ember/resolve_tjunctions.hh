#pragma once

// Post-process EMBER polygon soup to produce a manifold triangle mesh.
//
// Strategy:
//   1. Merge exact-duplicate vertices (tight tolerance)
//   2. Find boundary edges (1 incident face)
//   3. Stitch boundary edges by merging nearby vertices from different faces
//   4. Resolve remaining T-junctions by splitting edges at intervening vertices
//   5. Fix winding consistency via signed volume

#include <ember/mesh.hh>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <vector>

namespace ember
{

inline TriangleSoup resolve_tjunctions(TriangleSoup const& input)
{
    auto const& verts = input.vertices;
    auto const& tris = input.triangles;
    if (tris.empty()) return input;

    // Step 1: Merge exact-duplicate vertices (tight tolerance for rounding)
    int nv = (int)verts.size();
    double tight_tol = 1e-5;
    double inv_tol = 1.0 / tight_tol;

    struct VKey { int x, y, z;
        bool operator<(VKey const& o) const { return std::tie(x,y,z) < std::tie(o.x,o.y,o.z); }};
    std::map<VKey, int> vmap;
    std::vector<OutputVertex> mv;
    std::vector<int> remap(nv);

    for (int i = 0; i < nv; i++)
    {
        VKey k{(int)std::round(verts[i].x * inv_tol),
               (int)std::round(verts[i].y * inv_tol),
               (int)std::round(verts[i].z * inv_tol)};
        auto it = vmap.find(k);
        if (it != vmap.end()) remap[i] = it->second;
        else { remap[i] = (int)mv.size(); vmap[k] = remap[i]; mv.push_back(verts[i]); }
    }

    std::vector<std::array<int,3>> mt;
    for (auto const& t : tris)
    {
        int a = remap[t[0]], b = remap[t[1]], c = remap[t[2]];
        if (a != b && b != c && a != c) mt.push_back({a, b, c});
    }

    // Step 2: Match boundary edge PAIRS and merge their vertices.
    // Adjacent BSP-split polygons create near-coincident but not identical
    // vertices at the intersection curve. We match pairs of boundary edges
    // by proximity and merge corresponding endpoints.
    double stitch_tol = 1e-4; // Must be > merge_tol to catch T-junction near-misses

    auto vdist = [&](int a, int b) {
        double dx = mv[a].x-mv[b].x, dy = mv[a].y-mv[b].y, dz = mv[a].z-mv[b].z;
        return std::sqrt(dx*dx+dy*dy+dz*dz);
    };
    auto vmid = [&](int a) -> std::tuple<double,double,double> {
        return {mv[a].x, mv[a].y, mv[a].z};
    };

    for (int pass = 0; pass < 5; pass++)
    {
        struct Edge { int a, b;
            bool operator<(Edge const& o) const { return std::tie(a,b) < std::tie(o.a,o.b); }};
        auto mk = [](int a, int b) -> Edge { return a<b ? Edge{a,b} : Edge{b,a}; };

        std::map<Edge, int> edge_count;
        for (auto& t : mt)
        {
            edge_count[mk(t[0],t[1])]++;
            edge_count[mk(t[1],t[2])]++;
            edge_count[mk(t[2],t[0])]++;
        }

        // Collect boundary edges (exactly 1 face)
        struct BEdge { int a, b, face; };
        std::vector<BEdge> bedges;
        for (auto& [e, c] : edge_count)
            if (c == 1) bedges.push_back({e.a, e.b, -1});
        if (bedges.empty()) break;

        // Match boundary edge pairs: for each boundary edge, find the nearest
        // boundary edge (from a DIFFERENT face) and merge vertices pairwise
        std::vector<int> merge_to(mv.size());
        for (int i = 0; i < (int)mv.size(); i++) merge_to[i] = i;
        bool any_merged = false;

        std::vector<bool> matched(bedges.size(), false);
        for (int i = 0; i < (int)bedges.size(); i++)
        {
            if (matched[i]) continue;
            auto& ei = bedges[i];
            double best_d = stitch_tol * 2;
            int best_j = -1;
            bool flip = false;

            for (int j = i+1; j < (int)bedges.size(); j++)
            {
                if (matched[j]) continue;
                auto& ej = bedges[j];

                // Try both orientations
                double d_same = vdist(ei.a, ej.a) + vdist(ei.b, ej.b);
                double d_flip = vdist(ei.a, ej.b) + vdist(ei.b, ej.a);
                double d = std::min(d_same, d_flip);
                if (d < best_d) { best_d = d; best_j = j; flip = (d_flip < d_same); }
            }

            if (best_j >= 0)
            {
                matched[i] = matched[best_j] = true;
                auto& ej = bedges[best_j];
                if (flip)
                {
                    if (ei.a != ej.b) { merge_to[ej.b] = ei.a; any_merged = true; }
                    if (ei.b != ej.a) { merge_to[ej.a] = ei.b; any_merged = true; }
                }
                else
                {
                    if (ei.a != ej.a) { merge_to[ej.a] = ei.a; any_merged = true; }
                    if (ei.b != ej.b) { merge_to[ej.b] = ei.b; any_merged = true; }
                }
            }
        }
        if (!any_merged) break;

        // Resolve transitive merges
        for (int i = 0; i < (int)mv.size(); i++)
        {
            int r = i; while (merge_to[r] != r) r = merge_to[r]; merge_to[i] = r;
        }
        for (auto& t : mt) { t[0] = merge_to[t[0]]; t[1] = merge_to[t[1]]; t[2] = merge_to[t[2]]; }

        // Remove degenerate
        std::vector<std::array<int,3>> clean;
        for (auto& t : mt)
            if (t[0] != t[1] && t[1] != t[2] && t[0] != t[2]) clean.push_back(t);
        mt = std::move(clean);
    }

    // Step 3: Remove duplicate and near-duplicate faces
    {
        auto face_key = [](std::array<int,3> t) -> std::tuple<int,int,int> {
            if (t[0] > t[1]) std::swap(t[0], t[1]);
            if (t[1] > t[2]) std::swap(t[1], t[2]);
            if (t[0] > t[1]) std::swap(t[0], t[1]);
            return {t[0], t[1], t[2]};
        };
        std::set<std::tuple<int,int,int>> seen;
        std::vector<std::array<int,3>> unique;
        for (auto& t : mt)
        {
            auto k = face_key(t);
            if (seen.insert(k).second) unique.push_back(t);
        }
        mt = std::move(unique);
    }

    // Step 4: T-junction resolution on boundary edges
    auto dist_to_seg = [&](int vi, int ea, int eb, double& t_out) -> double {
        double dx = mv[eb].x-mv[ea].x, dy = mv[eb].y-mv[ea].y, dz = mv[eb].z-mv[ea].z;
        double len2 = dx*dx + dy*dy + dz*dz;
        if (len2 < 1e-30) { t_out = 0; return 1e30; }
        double px = mv[vi].x-mv[ea].x, py = mv[vi].y-mv[ea].y, pz = mv[vi].z-mv[ea].z;
        t_out = (px*dx + py*dy + pz*dz) / len2;
        double t = std::clamp(t_out, 0.0, 1.0);
        double qx = mv[ea].x+t*dx-mv[vi].x, qy = mv[ea].y+t*dy-mv[vi].y, qz = mv[ea].z+t*dz-mv[vi].z;
        return std::sqrt(qx*qx + qy*qy + qz*qz);
    };

    struct SEdge { int a, b;
        bool operator<(SEdge const& o) const { return std::tie(a,b) < std::tie(o.a,o.b); }};
    auto mke = [](int a, int b) -> SEdge { return a<b ? SEdge{a,b} : SEdge{b,a}; };

    for (int iter = 0; iter < 8; iter++)
    {
        std::map<SEdge, std::vector<int>> ef;
        for (int fi = 0; fi < (int)mt.size(); fi++)
        {
            auto& t = mt[fi];
            ef[mke(t[0],t[1])].push_back(fi);
            ef[mke(t[1],t[2])].push_back(fi);
            ef[mke(t[2],t[0])].push_back(fi);
        }

        bool any_split = false;
        std::set<int> to_remove;
        std::vector<std::array<int,3>> new_tris;

        for (auto& [edge, faces] : ef)
        {
            if (faces.size() != 1) continue;
            std::vector<std::pair<double,int>> on_edge;
            for (int vi = 0; vi < (int)mv.size(); vi++)
            {
                if (vi == edge.a || vi == edge.b) continue;
                double t;
                double d = dist_to_seg(vi, edge.a, edge.b, t);
                if (d < stitch_tol && t > 0.001 && t < 0.999)
                    on_edge.push_back({t, vi});
            }
            if (on_edge.empty()) continue;
            std::sort(on_edge.begin(), on_edge.end());

            for (int fi : faces)
            {
                if (to_remove.count(fi)) continue;
                auto& tri = mt[fi];
                for (int ei = 0; ei < 3; ei++)
                {
                    int ea = tri[ei], eb = tri[(ei+1)%3], ec = tri[(ei+2)%3];
                    if (mke(ea,eb).a != edge.a || mke(ea,eb).b != edge.b) continue;
                    to_remove.insert(fi);
                    any_split = true;
                    std::vector<int> chain = {ea};
                    for (auto& [t, vi] : on_edge) chain.push_back(vi);
                    chain.push_back(eb);
                    for (int k = 0; k+1 < (int)chain.size(); k++)
                        if (chain[k] != ec && chain[k+1] != ec && chain[k] != chain[k+1])
                            new_tris.push_back({chain[k], chain[k+1], ec});
                    break;
                }
            }
        }
        if (!any_split) break;
        std::vector<std::array<int,3>> kept;
        for (int i = 0; i < (int)mt.size(); i++)
            if (!to_remove.count(i)) kept.push_back(mt[i]);
        for (auto& t : new_tris) kept.push_back(t);
        mt = std::move(kept);
    }

    // Step 5: Fix winding
    double vol = 0;
    for (auto& t : mt)
    {
        auto& v0 = mv[t[0]]; auto& v1 = mv[t[1]]; auto& v2 = mv[t[2]];
        vol += (v0.x*(v1.y*v2.z-v1.z*v2.y)+v0.y*(v1.z*v2.x-v1.x*v2.z)+v0.z*(v1.x*v2.y-v1.y*v2.x))/6;
    }
    if (vol < 0) for (auto& t : mt) std::swap(t[0], t[1]);

    TriangleSoup out;
    out.vertices.resize(mv.size());
    for (int i = 0; i < (int)mv.size(); i++) out.vertices[i] = mv[i];
    out.triangles = std::move(mt);
    return out;
}

} // namespace ember
