#pragma once

// EMBER Integer Exact CSG - Face-Local BSP Tree
// Per the paper (Section 4.3):
//   For each polygon t, build a local BSP by adding intersection segments.
//   Each BSP leaf is a convex sub-polygon of t.
//   After all segments are added, no BSP leaf has an interior intersection
//   with any other polygon.

#include <ember/intersect_polygons.hh>
#include <ember/polygon.hh>
#include <ember/types.hh>
#include <ember/exact_classify.hh>

#include <vector>

namespace ember
{

// BSP leaf: a convex sub-polygon of the host polygon
struct BSPLeaf
{
    std::vector<plane_t> edges;
    bool enabled = true;
};

// BSP node: either a leaf or an inner node with a splitting plane.
// Children are raw pointers into the owning LocalBSP's node pool.
struct BSPNode
{
    bool is_leaf = true;
    BSPLeaf leaf;

    plane_t split_plane;
    BSPNode* negative = nullptr;
    BSPNode* positive = nullptr;
};

// Face-local BSP tree for a single polygon.
// All nodes are pool-allocated in `nodes` — no per-node heap allocation.
struct LocalBSP
{
    plane_t support;
    int host_polygon_index;

    // Node pool: all BSPNodes live here. Pointers are stable because
    // we reserve enough capacity upfront and never exceed it.
    std::vector<BSPNode> nodes;
    BSPNode* root = nullptr;

    BSPNode* alloc_leaf(std::vector<plane_t> const& edges)
    {
        nodes.emplace_back();
        auto* n = &nodes.back();
        n->is_leaf = true;
        n->leaf.edges = edges;
        return n;
    }
    BSPNode* alloc_leaf(std::vector<plane_t>&& edges)
    {
        nodes.emplace_back();
        auto* n = &nodes.back();
        n->is_leaf = true;
        n->leaf.edges = std::move(edges);
        return n;
    }

    void init(ConvexPolygon const& poly)
    {
        support = poly.support;
        host_polygon_index = poly.polygon_index;
        // Reserve generous capacity to prevent reallocation (invalidates pointers)
        nodes.reserve(64);
        root = alloc_leaf(poly.edges);
    }

    void add_segment(point4_t const& v0, point4_t const& v1, plane_t const& s)
    {
        if (root)
            add_segment_recursive(root, v0, v1, s);
    }

    void add_overlap(ConvexPolygon const& other, int other_poly_idx)
    {
        for (auto const& edge : other.edges)
        {
            if (root)
                add_plane_split_recursive(root, edge);
        }
        if (host_polygon_index > other_poly_idx)
            mark_overlapping_leaves(root, other);
    }

    void collect_leaves(std::vector<BSPLeaf*>& out)
    {
        if (root)
            collect_leaves_recursive(root, out);
    }

private:
    void add_segment_recursive(BSPNode* node, point4_t const& v0, point4_t const& v1, plane_t const& s)
    {
        if (node->is_leaf)
        {
            // Paper Section 4.3: split the leaf by the segment's plane.
            // split_leaf checks if the plane actually divides the leaf.
            split_leaf(node, s);
            return;
        }

        auto c0 = exact_classify(v0, node->split_plane);
        auto c1 = exact_classify(v1, node->split_plane);

        if (c0 == 0 && c1 == 0) return; // Segment on splitting plane

        if (c0 <= 0 && c1 <= 0)
        {
            add_segment_recursive(node->negative, v0, v1, s);
        }
        else if (c0 >= 0 && c1 >= 0)
        {
            add_segment_recursive(node->positive, v0, v1, s);
        }
        else
        {
            auto int_line = ipg::intersect<geometry_t>(support, s);
            auto v_mid = ipg::intersect(int_line, node->split_plane);

            if (c0 < 0)
            {
                add_segment_recursive(node->negative, v0, v_mid, s);
                add_segment_recursive(node->positive, v_mid, v1, s);
            }
            else
            {
                add_segment_recursive(node->positive, v0, v_mid, s);
                add_segment_recursive(node->negative, v_mid, v1, s);
            }
        }
    }

    void split_leaf(BSPNode* node, plane_t const& s)
    {
        auto& old_edges = node->leaf.edges;
        bool was_enabled = node->leaf.enabled;

        int n = static_cast<int>(old_edges.size());
        if (n < 3) return;

        // Classify vertices against the split plane
        // Use a small fixed-size buffer to avoid allocation for typical polygons
        tg::i8 cls_buf[16];
        std::vector<tg::i8> cls_vec;
        tg::i8* cls;
        if (n <= 16) { cls = cls_buf; }
        else { cls_vec.resize(n); cls = cls_vec.data(); }

        bool has_pos = false, has_neg = false;
        for (int i = 0; i < n; i++)
        {
            auto v = ipg::intersect(support, old_edges[i], old_edges[(i + 1) % n]);
            cls[i] = exact_classify(v, s);
            if (cls[i] > 0) has_pos = true;
            if (cls[i] < 0) has_neg = true;
        }

        if (!has_pos || !has_neg) return; // Plane doesn't divide this leaf

        // Build edge lists (same algorithm as clip_polygon)
        plane_t s_inv = s.inverted();
        std::vector<plane_t> neg_edges, pos_edges;
        neg_edges.reserve(n + 2);
        pos_edges.reserve(n + 2);

        for (int i = 0; i < n; i++)
        {
            int j = (i + 1) % n;
            plane_t seg_edge = old_edges[j];

            if (cls[i] <= 0 && cls[j] <= 0)
            {
                neg_edges.push_back(seg_edge);
            }
            else if (cls[i] <= 0 && cls[j] > 0)
            {
                neg_edges.push_back(seg_edge);
                neg_edges.push_back(s);
                pos_edges.push_back(seg_edge);
            }
            else if (cls[i] > 0 && cls[j] <= 0)
            {
                pos_edges.push_back(seg_edge);
                pos_edges.push_back(s_inv);
                neg_edges.push_back(seg_edge);
            }
            else
            {
                pos_edges.push_back(seg_edge);
            }
        }

        // Convert this leaf to an inner node
        node->is_leaf = false;
        node->split_plane = s;
        node->leaf.edges.clear(); // Free memory

        // Ensure pool capacity before allocating (prevents reallocation)
        if (nodes.size() + 2 > nodes.capacity())
            nodes.reserve(nodes.capacity() * 2);

        node->negative = alloc_leaf(std::move(neg_edges));
        node->negative->leaf.enabled = was_enabled;
        node->positive = alloc_leaf(std::move(pos_edges));
        node->positive->leaf.enabled = was_enabled;
    }

    void add_plane_split_recursive(BSPNode* node, plane_t const& s)
    {
        if (node->is_leaf)
        {
            split_leaf(node, s);
            return;
        }
        add_plane_split_recursive(node->negative, s);
        add_plane_split_recursive(node->positive, s);
    }

    void mark_overlapping_leaves(BSPNode* node, ConvexPolygon const& other)
    {
        if (node->is_leaf)
        {
            if (!node->leaf.enabled) return;
            int n = static_cast<int>(node->leaf.edges.size());
            if (n < 3) return;

            auto test_pt = ipg::intersect(support, node->leaf.edges[0], node->leaf.edges[1]);
            bool inside = true;
            for (auto const& e : other.edges)
            {
                if (exact_classify(test_pt, e) >= 0) { inside = false; break; }
            }
            if (inside) node->leaf.enabled = false;
            return;
        }
        mark_overlapping_leaves(node->negative, other);
        mark_overlapping_leaves(node->positive, other);
    }

    void collect_leaves_recursive(BSPNode* node, std::vector<BSPLeaf*>& out)
    {
        if (node->is_leaf)
        {
            if (node->leaf.enabled) out.push_back(&node->leaf);
            return;
        }
        collect_leaves_recursive(node->negative, out);
        collect_leaves_recursive(node->positive, out);
    }
};

} // namespace ember
