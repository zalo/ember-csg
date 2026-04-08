// EMBER Integer Exact CSG - Test: Boolean operations on two cubes
//
// Creates two overlapping unit cubes and tests union, intersection, and difference.
// This validates the complete EMBER pipeline:
//   input preparation → subdivision → local BSP → classification → output

#include <ember/ember.hh>

#include <cstdio>
#include <cmath>

// Helper: create a cube InputMesh centered at (cx, cy, cz) with half-extent s
ember::InputMesh make_cube(double cx, double cy, double cz, double s)
{
    ember::InputMesh mesh;
    mesh.nsi = true;
    mesh.nnc = true;

    // 8 vertices of the cube
    mesh.positions = {
        {cx - s, cy - s, cz - s}, // 0: ---
        {cx + s, cy - s, cz - s}, // 1: +--
        {cx + s, cy + s, cz - s}, // 2: ++-
        {cx - s, cy + s, cz - s}, // 3: -+-
        {cx - s, cy - s, cz + s}, // 4: --+
        {cx + s, cy - s, cz + s}, // 5: +-+
        {cx + s, cy + s, cz + s}, // 6: +++
        {cx - s, cy + s, cz + s}, // 7: -++
    };

    // 12 triangles (2 per face, outward-facing CCW)
    mesh.triangles = {
        // Front face (z = cz+s): 4, 5, 6, 7
        {4, 5, 6},
        {4, 6, 7},
        // Back face (z = cz-s): 0, 3, 2, 1
        {0, 3, 2},
        {0, 2, 1},
        // Right face (x = cx+s): 1, 2, 6, 5
        {1, 2, 6},
        {1, 6, 5},
        // Left face (x = cx-s): 0, 4, 7, 3
        {0, 4, 7},
        {0, 7, 3},
        // Top face (y = cy+s): 3, 7, 6, 2
        {3, 7, 6},
        {3, 6, 2},
        // Bottom face (y = cy-s): 0, 1, 5, 4
        {0, 1, 5},
        {0, 5, 4},
    };

    return mesh;
}

void print_result(const char* name, ember::BooleanResult const& result)
{
    auto polys = ember::extract_output(result);
    auto soup = ember::triangulate_output(result);

    std::printf("\n=== %s ===\n", name);
    std::printf("Output polygons: %zu\n", polys.size());
    std::printf("Triangulated: %zu vertices, %zu triangles\n",
                soup.vertices.size(), soup.triangles.size());

    // Print first few polygon vertex counts
    if (!polys.empty())
    {
        std::printf("Polygon vertex counts:");
        for (size_t i = 0; i < std::min(polys.size(), size_t(10)); i++)
            std::printf(" %zu", polys[i].vertices.size());
        if (polys.size() > 10)
            std::printf(" ...");
        std::printf("\n");
    }

    // Write OBJ to file
    auto obj = ember::to_obj(soup);
    char filename[256];
    std::snprintf(filename, sizeof(filename), "%s.obj", name);
    if (auto* f = std::fopen(filename, "w"))
    {
        std::fputs(obj.c_str(), f);
        std::fclose(f);
        std::printf("Written to %s\n", filename);
    }
}

void validate_polygon(ember::ConvexPolygon const& poly, const char* label)
{
    int n = poly.vertex_count();
    // Check that the centroid is on the supporting plane and inside all edges
    tg::dpos3 center = {0, 0, 0};
    for (int i = 0; i < n; i++)
    {
        auto v = poly.vertex_dpos(i);
        center.x += v.x;
        center.y += v.y;
        center.z += v.z;
    }
    center.x /= n;
    center.y /= n;
    center.z /= n;

    // Check each vertex is on the support plane and on exactly two edge planes
    for (int i = 0; i < n; i++)
    {
        auto v = poly.vertex(i);
        auto cs = ipg::classify<ember::geometry_t>(v, poly.support);
        if (cs != 0)
            std::printf("  [%s] WARNING: vertex %d not on support plane (classify=%d)\n", label, i, cs);

        // Each vertex should satisfy classify <= 0 for all edges
        for (int j = 0; j < n; j++)
        {
            auto ce = ipg::classify<ember::geometry_t>(v, poly.edges[j]);
            if (ce > 0)
                std::printf("  [%s] WARNING: vertex %d outside edge %d (classify=%d)\n", label, i, j, ce);
        }
    }
}

int main()
{
    std::printf("EMBER Integer Exact CSG - Cube Boolean Test\n");
    std::printf("============================================\n");

    // Create two overlapping cubes
    // Cube A centered at origin, half-extent 1.0
    // Cube B centered at (0.5, 0.5, 0.5), half-extent 1.0
    auto cube_a = make_cube(0.0, 0.0, 0.0, 1.0);
    auto cube_b = make_cube(0.5, 0.5, 0.5, 1.0);

    std::printf("\nCube A: center (0,0,0), half-extent 1.0 → 12 triangles\n");
    std::printf("Cube B: center (0.5,0.5,0.5), half-extent 1.0 → 12 triangles\n");

    // Validate: convert to polygon soup and check edge orientations
    {
        auto soup = ember::prepare_input({cube_a});
        std::printf("\nValidating cube A polygons (%zu polygons):\n", soup.polygons.size());
        for (size_t i = 0; i < std::min(soup.polygons.size(), size_t(3)); i++)
        {
            char label[32];
            std::snprintf(label, sizeof(label), "poly_%zu", i);
            validate_polygon(soup.polygons[i], label);
        }
        std::printf("Validation complete.\n");
    }

    ember::EmberConfig config;
    config.assume_nsi = true;
    config.assume_nnc = true;

    // Test Union
    std::printf("\nComputing Union (A ∪ B)...\n");
    auto result_union = ember::boolean_union(cube_a, cube_b, config);
    print_result("union", result_union);

    // Test Intersection
    std::printf("\nComputing Intersection (A ∩ B)...\n");
    auto result_isect = ember::boolean_intersection(cube_a, cube_b, config);
    print_result("intersection", result_isect);

    // Test Difference
    std::printf("\nComputing Difference (A - B)...\n");
    auto result_diff = ember::boolean_difference(cube_a, cube_b, config);
    print_result("difference", result_diff);

    // ==========================================
    // Test with sphere meshes for more complex validation
    // ==========================================
    std::printf("\n\n===== SPHERE TESTS =====\n");

    // Load manifold3d-exported spheres (run compare/export_spheres.py first)
    auto sphere_a = ember::load_obj("sphere_a.obj", true, true);
    auto sphere_b = ember::load_obj("sphere_b.obj", true, true);

    if (sphere_a.triangles.empty() || sphere_b.triangles.empty())
    {
        std::printf("ERROR: Could not load sphere_a.obj / sphere_b.obj\n");
        std::printf("Run: python3 ../compare/export_spheres.py\n");
        return 1;
    }

    std::printf("Sphere A: %zu verts, %zu triangles (loaded from manifold3d export)\n",
                sphere_a.positions.size(), sphere_a.triangles.size());
    std::printf("Sphere B: %zu verts, %zu triangles (loaded from manifold3d export)\n",
                sphere_b.positions.size(), sphere_b.triangles.size());

    // Debug: check the transform for sphere inputs
    {
        auto test_soup = ember::prepare_input({sphere_a, sphere_b});
        std::printf("Sphere transform: scale=%.4f center=(%.4f, %.4f, %.4f)\n",
            test_soup.transform.scale, test_soup.transform.center_x,
            test_soup.transform.center_y, test_soup.transform.center_z);
    }

    // Sphere Union
    std::printf("\nComputing Sphere Union...\n");
    auto sphere_union = ember::boolean_union(sphere_a, sphere_b, config);
    print_result("sphere_union", sphere_union);

    // Sphere Intersection
    std::printf("\nComputing Sphere Intersection...\n");
    auto sphere_isect = ember::boolean_intersection(sphere_a, sphere_b, config);
    print_result("sphere_intersection", sphere_isect);

    // Sphere Difference
    std::printf("\nComputing Sphere Difference (A - B)...\n");
    auto sphere_diff = ember::boolean_difference(sphere_a, sphere_b, config);
    print_result("sphere_difference", sphere_diff);

    std::printf("\nAll tests completed.\n");
    return 0;
}
