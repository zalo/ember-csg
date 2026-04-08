#pragma once

// EMBER Integer Exact CSG - Face-Local BSP Tree
// Per the paper (Section 4.3):
//   For each polygon t, build a local BSP by adding intersection segments.
//   Each BSP leaf is a convex sub-polygon of t.
//   After all segments are added, no BSP leaf has an interior intersection
//   with any other polygon.
//
// BSP "add segment" operation:
//   Recursively split BSP nodes at the segment's splitting plane.
//   At leaf nodes, the splitting plane creates a new inner node.
//
// Overlap handling (Section 4.3, Fig. 8):
//   When two coplanar polygons overlap, add all edges of the other polygon.
//   Mark overlapping leaves as "disabled" for the higher-indexed polygon.

#include <ember/intersect_polygons.hh>
#include <ember/polygon.hh>
#include <ember/types.hh>
#include <ember/exact_classify.hh>

#include <memory>
#include <vector>

namespace ember
{

// Forward declarations
struct BSPNode;

// BSP leaf: a convex sub-polygon of the host polygon
struct BSPLeaf
{
    std::vector<plane_t> edges; // Edge planes of this sub-polygon
    bool enabled = true;        // False if disabled due to overlap resolution
    int overlap_priority = -1;  // If overlapping, the polygon index that "owns" this region
};

// BSP node: either a leaf or an inner node with a splitting plane
struct BSPNode
{
    bool is_leaf = true;

    // Leaf data
    BSPLeaf leaf;

    // Inner node data
    plane_t split_plane;                    // The plane that splits this node
    std::unique_ptr<BSPNode> negative;      // Child on negative side (classify <= 0)
    std::unique_ptr<BSPNode> positive;      // Child on positive side (classify > 0)

    // Create a leaf node from a polygon's edge planes
    static std::unique_ptr<BSPNode> make_leaf(std::vector<plane_t> const& edges)
    {
        auto node = std::make_unique<BSPNode>();
        node->is_leaf = true;
        node->leaf.edges = edges;
        return node;
    }
};

// Face-local BSP tree for a single polygon
struct LocalBSP
{
    plane_t support;        // Supporting plane of the host polygon
    int host_polygon_index; // Index of the host polygon (for overlap resolution)
    std::unique_ptr<BSPNode> root;

    // Initialize BSP from a convex polygon
    void init(ConvexPolygon const& poly)
    {
        support = poly.support;
        host_polygon_index = poly.polygon_index;
        root = BSPNode::make_leaf(poly.edges);
    }

    // Add an intersection segment to the BSP
    // The segment is defined by two endpoints v0, v1 and a splitting plane s
    // (The segment lies on the line intersect(support, s))
    void add_segment(point4_t const& v0, point4_t const& v1, plane_t const& s)
    {
        if (root)
            add_segment_recursive(root.get(), v0, v1, s);
    }

    // Add all edges of an overlapping polygon (for C4 overlap handling)
    void add_overlap(ConvexPolygon const& other, int other_poly_idx)
    {
        // Add each edge of the other polygon as a splitting plane
        for (auto const& edge : other.edges)
        {
            // For each edge, we need to construct a "segment" that covers the overlap.
            // Since the edge plane extends infinitely on the supporting plane,
            // we add it as a full-plane split rather than a bounded segment.
            if (root)
                add_plane_split_recursive(root.get(), edge);
        }

        // Mark overlap regions: leaves that are inside the other polygon
        // should be disabled if this polygon has higher index
        if (host_polygon_index > other_poly_idx)
        {
            mark_overlapping_leaves(root.get(), other);
        }
    }

    // Enumerate all enabled BSP leaves
    void collect_leaves(std::vector<BSPLeaf*>& out)
    {
        if (root)
            collect_leaves_recursive(root.get(), out);
    }

    // Get the edge planes for a leaf, including inherited edge planes from ancestors
    // This reconstructs the full convex polygon for each leaf
    void collect_leaf_polygons(std::vector<std::vector<plane_t>>& out)
    {
        std::vector<plane_t> accumulated_edges;
        if (root)
            collect_leaf_polygons_recursive(root.get(), accumulated_edges, out);
    }

private:
    // Recursive add_segment: split BSP nodes at the segment's plane
    void add_segment_recursive(BSPNode* node, point4_t const& v0, point4_t const& v1, plane_t const& s)
    {
        if (node->is_leaf)
        {
            // At leaf: the segment's splitting plane divides this leaf into two
            // Check if the segment actually intersects this leaf
            // by classifying v0 and v1 against the leaf's edge planes

            bool v0_inside = true, v1_inside = true;
            for (auto const& e : node->leaf.edges)
            {
                if (exact_classify(v0, e) > 0) v0_inside = false;
                if (exact_classify(v1, e) > 0) v1_inside = false;
            }

            // Also check if the splitting plane actually cuts this leaf
            // by classifying leaf vertices against the splitting plane
            bool has_pos = false, has_neg = false;
            int n = static_cast<int>(node->leaf.edges.size());
            for (int i = 0; i < n; i++)
            {
                auto v = ipg::intersect(support, node->leaf.edges[i], node->leaf.edges[(i + 1) % n]);
                auto c = exact_classify(v, s);
                if (c > 0) has_pos = true;
                if (c < 0) has_neg = true;
            }

            if (!has_pos || !has_neg)
            {
                // Splitting plane doesn't cut this leaf, or segment is entirely outside
                return;
            }

            // Split the leaf at plane s
            split_leaf(node, s);
            return;
        }

        // Inner node: classify v0 and v1 against the splitting plane
        auto c0 = exact_classify(v0, node->split_plane);
        auto c1 = exact_classify(v1, node->split_plane);

        if (c0 == 0 && c1 == 0)
        {
            // Segment lies on the splitting plane → stop (no further subdivision needed)
            return;
        }

        if (c0 <= 0 && c1 <= 0)
        {
            // Both on negative side
            add_segment_recursive(node->negative.get(), v0, v1, s);
        }
        else if (c0 >= 0 && c1 >= 0)
        {
            // Both on positive side
            add_segment_recursive(node->positive.get(), v0, v1, s);
        }
        else
        {
            // Straddling: split the segment at the splitting plane
            // New vertex: intersect of the segment's line with the splitting plane
            auto int_line = ipg::intersect<geometry_t>(support, s);
            auto v_mid = ipg::intersect(int_line, node->split_plane);

            if (c0 < 0)
            {
                // v0 on negative, v1 on positive
                add_segment_recursive(node->negative.get(), v0, v_mid, s);
                add_segment_recursive(node->positive.get(), v_mid, v1, s);
            }
            else
            {
                // v0 on positive, v1 on negative
                add_segment_recursive(node->positive.get(), v0, v_mid, s);
                add_segment_recursive(node->negative.get(), v_mid, v1, s);
            }
        }
    }

    // Split a leaf node at a plane, creating an inner node
    void split_leaf(BSPNode* node, plane_t const& s)
    {
        auto old_edges = std::move(node->leaf.edges);
        bool was_enabled = node->leaf.enabled;

        node->is_leaf = false;
        node->split_plane = s;

        // Create two new leaf nodes by clipping the old edges against s
        std::vector<plane_t> neg_edges, pos_edges;
        plane_t s_inv = s.inverted();

        int n = static_cast<int>(old_edges.size());

        // Classify vertices of the old leaf
        std::vector<tg::i8> classifications(n);
        for (int i = 0; i < n; i++)
        {
            auto v = ipg::intersect(support, old_edges[i], old_edges[(i + 1) % n]);
            classifications[i] = exact_classify(v, s);
        }

        bool has_pos = false, has_neg = false;
        for (auto c : classifications)
        {
            if (c > 0) has_pos = true;
            if (c < 0) has_neg = true;
        }

        if (!has_pos || !has_neg)
        {
            // Degenerate split: plane doesn't actually divide the leaf
            // Restore as leaf
            node->is_leaf = true;
            node->leaf.edges = std::move(old_edges);
            node->leaf.enabled = was_enabled;
            return;
        }

        // Build edge lists using the same algorithm as clip_polygon
        for (int i = 0; i < n; i++)
        {
            int j = (i + 1) % n;
            plane_t seg_edge = old_edges[j];

            if (classifications[i] <= 0 && classifications[j] <= 0)
            {
                neg_edges.push_back(seg_edge);
            }
            else if (classifications[i] <= 0 && classifications[j] > 0)
            {
                neg_edges.push_back(seg_edge);
                neg_edges.push_back(s);
                pos_edges.push_back(seg_edge);
            }
            else if (classifications[i] > 0 && classifications[j] <= 0)
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

        node->negative = BSPNode::make_leaf(neg_edges);
        node->negative->leaf.enabled = was_enabled;
        node->positive = BSPNode::make_leaf(pos_edges);
        node->positive->leaf.enabled = was_enabled;
    }

    // Add a plane split to the entire BSP (for overlap handling)
    void add_plane_split_recursive(BSPNode* node, plane_t const& s)
    {
        if (node->is_leaf)
        {
            // Check if the plane actually splits this leaf
            int n = static_cast<int>(node->leaf.edges.size());
            bool has_pos = false, has_neg = false;
            for (int i = 0; i < n; i++)
            {
                auto v = ipg::intersect(support, node->leaf.edges[i], node->leaf.edges[(i + 1) % n]);
                auto c = exact_classify(v, s);
                if (c > 0) has_pos = true;
                if (c < 0) has_neg = true;
            }

            if (has_pos && has_neg)
                split_leaf(node, s);
            return;
        }

        // Recurse into both children
        add_plane_split_recursive(node->negative.get(), s);
        add_plane_split_recursive(node->positive.get(), s);
    }

    // Mark BSP leaves that overlap with the given polygon as disabled
    void mark_overlapping_leaves(BSPNode* node, ConvexPolygon const& other)
    {
        if (node->is_leaf)
        {
            if (!node->leaf.enabled) return;

            // Check if this leaf is inside the other polygon
            // Use the centroid of the leaf for the test
            int n = static_cast<int>(node->leaf.edges.size());
            if (n < 3) return;

            // Check if ANY vertex of the leaf is STRICTLY inside the other polygon.
            // Strict test (< 0) prevents falsely disabling adjacent non-overlapping
            // coplanar triangles that share an edge (shared-edge vertices classify
            // as 0 which should NOT be treated as overlapping).
            auto test_pt = ipg::intersect(support, node->leaf.edges[0], node->leaf.edges[1]);
            bool inside = true;
            for (auto const& e : other.edges)
            {
                if (exact_classify(test_pt, e) >= 0)
                {
                    inside = false;
                    break;
                }
            }

            if (inside)
                node->leaf.enabled = false;
            return;
        }

        mark_overlapping_leaves(node->negative.get(), other);
        mark_overlapping_leaves(node->positive.get(), other);
    }

    // Collect all enabled leaves
    void collect_leaves_recursive(BSPNode* node, std::vector<BSPLeaf*>& out)
    {
        if (node->is_leaf)
        {
            if (node->leaf.enabled)
                out.push_back(&node->leaf);
            return;
        }
        collect_leaves_recursive(node->negative.get(), out);
        collect_leaves_recursive(node->positive.get(), out);
    }

    // Collect full edge plane lists for each enabled leaf
    void collect_leaf_polygons_recursive(BSPNode* node,
                                         std::vector<plane_t>& accumulated,
                                         std::vector<std::vector<plane_t>>& out)
    {
        if (node->is_leaf)
        {
            if (node->leaf.enabled)
                out.push_back(node->leaf.edges);
            return;
        }

        // Negative child: add split_plane constraint
        accumulated.push_back(node->split_plane);
        collect_leaf_polygons_recursive(node->negative.get(), accumulated, out);
        accumulated.pop_back();

        // Positive child: add inverted split_plane constraint
        accumulated.push_back(node->split_plane.inverted());
        collect_leaf_polygons_recursive(node->positive.get(), accumulated, out);
        accumulated.pop_back();
    }
};

} // namespace ember
