#pragma once

// EMBER Integer Exact CSG - Input/Output Mesh (Polygon Soup)
// Converts triangle meshes to polygon soup with WNTVs
// Handles coordinate scaling to 26-bit integer range

#include <ember/polygon.hh>
#include <ember/types.hh>
#include <ember/winding.hh>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace ember
{

// Input triangle: three vertex indices into a position array
struct Triangle
{
    int v0, v1, v2;
};

// Input mesh: positions + triangles
struct InputMesh
{
    std::vector<tg::dpos3> positions; // Floating-point positions
    std::vector<Triangle> triangles;

    bool nsi = false; // No Self-Intersections flag
    bool nnc = false; // No Nested Components flag
};

// Output vertex in original floating-point space
struct OutputVertex
{
    double x, y, z;
};

// Coordinate transform: maps between integer and original floating-point space
struct CoordinateTransform
{
    double center_x = 0, center_y = 0, center_z = 0;
    double scale = 1.0;

    OutputVertex to_float(double ix, double iy, double iz) const
    {
        return {ix / scale + center_x, iy / scale + center_y, iz / scale + center_z};
    }
};

// Polygon soup: the working representation for EMBER
struct PolygonSoup
{
    std::vector<ConvexPolygon> polygons;
    IAABB bounds;
    int num_meshes = 0;
    CoordinateTransform transform;

    // Compute AABB from all polygon vertices
    void compute_bounds()
    {
        if (polygons.empty()) return;

        auto first = polygons[0].vertex_dpos(0);
        auto min_x = first.x, min_y = first.y, min_z = first.z;
        auto max_x = first.x, max_y = first.y, max_z = first.z;

        for (auto const& poly : polygons)
        {
            for (int i = 0; i < poly.vertex_count(); i++)
            {
                auto p = poly.vertex_dpos(i);
                min_x = std::min(min_x, p.x);
                min_y = std::min(min_y, p.y);
                min_z = std::min(min_z, p.z);
                max_x = std::max(max_x, p.x);
                max_y = std::max(max_y, p.y);
                max_z = std::max(max_z, p.z);
            }
        }

        // Pad by 1 to ensure no vertex is exactly on the boundary
        bounds.min = pos_t(pos_scalar_t(int32_t(std::floor(min_x)) - 1),
                           pos_scalar_t(int32_t(std::floor(min_y)) - 1),
                           pos_scalar_t(int32_t(std::floor(min_z)) - 1));
        bounds.max = pos_t(pos_scalar_t(int32_t(std::ceil(max_x)) + 1),
                           pos_scalar_t(int32_t(std::ceil(max_y)) + 1),
                           pos_scalar_t(int32_t(std::ceil(max_z)) + 1));
    }
};


// Result of a boolean operation
struct BooleanResult
{
    PolygonSoup output;
    CoordinateTransform transform; // To convert output back to original coords

    // For each output polygon, the classification:
    //   +1 = emitted as-is (out→in)
    //   -1 = emitted inverted (in→out)
    std::vector<int> classifications;
};

// Scale and convert input meshes to integer polygon soup
// All meshes are combined into a single PolygonSoup with mesh indices
inline PolygonSoup prepare_input(std::vector<InputMesh> const& meshes)
{
    PolygonSoup soup;
    soup.num_meshes = static_cast<int>(meshes.size());

    // Find global AABB in floating-point
    double gmin_x = 1e30, gmin_y = 1e30, gmin_z = 1e30;
    double gmax_x = -1e30, gmax_y = -1e30, gmax_z = -1e30;

    for (auto const& mesh : meshes)
    {
        for (auto const& p : mesh.positions)
        {
            gmin_x = std::min(gmin_x, p.x);
            gmin_y = std::min(gmin_y, p.y);
            gmin_z = std::min(gmin_z, p.z);
            gmax_x = std::max(gmax_x, p.x);
            gmax_y = std::max(gmax_y, p.y);
            gmax_z = std::max(gmax_z, p.z);
        }
    }

    // Compute scale factor to fit into 26-bit integer range
    // Leave some headroom for clipping planes
    double extent = std::max({gmax_x - gmin_x, gmax_y - gmin_y, gmax_z - gmin_z, 1e-10});
    double center_x = (gmin_x + gmax_x) * 0.5;
    double center_y = (gmin_y + gmax_y) * 0.5;
    double center_z = (gmin_z + gmax_z) * 0.5;

    // Scale to fit in [-MAX_COORD+2, MAX_COORD-2] with margin
    double scale = (MAX_COORD - 4) * 2.0 / extent;

    soup.transform.center_x = center_x;
    soup.transform.center_y = center_y;
    soup.transform.center_z = center_z;
    soup.transform.scale = scale;

    int poly_index = 0;
    for (int mesh_idx = 0; mesh_idx < static_cast<int>(meshes.size()); mesh_idx++)
    {
        auto const& mesh = meshes[mesh_idx];

        // Convert positions to integer
        std::vector<pos_t> int_positions(mesh.positions.size());
        for (size_t i = 0; i < mesh.positions.size(); i++)
        {
            auto const& p = mesh.positions[i];
            int32_t ix = static_cast<int32_t>(std::round((p.x - center_x) * scale));
            int32_t iy = static_cast<int32_t>(std::round((p.y - center_y) * scale));
            int32_t iz = static_cast<int32_t>(std::round((p.z - center_z) * scale));

            // Clamp to valid range
            ix = std::clamp(ix, -MAX_COORD + 2, MAX_COORD - 2);
            iy = std::clamp(iy, -MAX_COORD + 2, MAX_COORD - 2);
            iz = std::clamp(iz, -MAX_COORD + 2, MAX_COORD - 2);

            int_positions[i] = pos_t(pos_scalar_t(ix), pos_scalar_t(iy), pos_scalar_t(iz));
        }

        // Convert triangles to polygons
        for (auto const& tri : mesh.triangles)
        {
            auto const& p0 = int_positions[tri.v0];
            auto const& p1 = int_positions[tri.v1];
            auto const& p2 = int_positions[tri.v2];

            // Skip degenerate triangles
            auto poly = make_triangle(p0, p1, p2, mesh_idx, poly_index);
            if (!poly.support.is_valid())
                continue;

            // Set WNTV: standard case for closed mesh
            // delta_w has 1 at mesh_index position
            poly.delta_w.resize(meshes.size(), 0);
            poly.delta_w[mesh_idx] = 1;

            poly.no_self_intersections = mesh.nsi;
            poly.no_nested_components = mesh.nnc;

            soup.polygons.push_back(std::move(poly));
            poly_index++;
        }
    }

    soup.compute_bounds();
    return soup;
}

struct OutputPolygon
{
    std::vector<OutputVertex> vertices;
    int source_mesh;    // Which input mesh this came from
    int source_polygon; // Original polygon index
};

// Load an OBJ file into an InputMesh
inline InputMesh load_obj(std::string const& path, bool nsi = false, bool nnc = false)
{
    InputMesh mesh;
    mesh.nsi = nsi;
    mesh.nnc = nnc;

    std::ifstream file(path);
    if (!file.is_open())
        return mesh;

    std::string line;
    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string token;
        iss >> token;

        if (token == "v")
        {
            double x, y, z;
            iss >> x >> y >> z;
            mesh.positions.push_back({x, y, z});
        }
        else if (token == "f")
        {
            std::vector<int> indices;
            std::string vert;
            while (iss >> vert)
            {
                // Handle "v", "v/vt", "v/vt/vn", "v//vn" formats
                int idx = std::stoi(vert.substr(0, vert.find('/')));
                indices.push_back(idx - 1); // OBJ is 1-indexed
            }
            // Fan-triangulate polygons with > 3 vertices
            for (size_t i = 1; i + 1 < indices.size(); i++)
            {
                mesh.triangles.push_back({indices[0], indices[i], indices[i + 1]});
            }
        }
    }

    return mesh;
}

} // namespace ember
