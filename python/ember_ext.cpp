// EMBER CSG — nanobind Python bindings
//
// Exposes the EMBER boolean CSG API to Python so benchmarks can call
// EMBER in-process (no subprocess / OBJ round-trip overhead).

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/ndarray.h>

#include <ember/ember.hh>

#include <chrono>

namespace nb = nanobind;
using namespace nb::literals;

// ---------------------------------------------------------------------------
// helpers — build InputMesh from numpy arrays
// ---------------------------------------------------------------------------

static ember::InputMesh mesh_from_arrays(
    nb::ndarray<double, nb::shape<-1, 3>> verts,
    nb::ndarray<int32_t, nb::shape<-1, 3>> tris,
    bool nsi, bool nnc)
{
    ember::InputMesh mesh;
    mesh.nsi = nsi;
    mesh.nnc = nnc;

    size_t nv = verts.shape(0);
    mesh.positions.resize(nv);
    for (size_t i = 0; i < nv; i++)
        mesh.positions[i] = {verts(i, 0), verts(i, 1), verts(i, 2)};

    size_t nt = tris.shape(0);
    mesh.triangles.resize(nt);
    for (size_t i = 0; i < nt; i++)
        mesh.triangles[i] = {tris(i, 0), tris(i, 1), tris(i, 2)};

    return mesh;
}

// Return (vertices Nx3, triangles Mx3) as vectors of arrays
using VertList = std::vector<std::array<double, 3>>;
using TriList  = std::vector<std::array<int32_t, 3>>;

static std::tuple<VertList, TriList> soup_to_arrays(ember::TriangleSoup const& soup)
{
    VertList verts(soup.vertices.size());
    for (size_t i = 0; i < soup.vertices.size(); i++)
        verts[i] = {soup.vertices[i].x, soup.vertices[i].y, soup.vertices[i].z};

    TriList tris(soup.triangles.size());
    for (size_t i = 0; i < soup.triangles.size(); i++)
        tris[i] = {soup.triangles[i][0], soup.triangles[i][1], soup.triangles[i][2]};

    return {std::move(verts), std::move(tris)};
}

// ---------------------------------------------------------------------------
// core boolean that returns (verts, tris, time_ms)
// ---------------------------------------------------------------------------

static std::tuple<VertList, TriList, double> run_boolean(
    nb::ndarray<double, nb::shape<-1, 3>> verts_a,
    nb::ndarray<int32_t, nb::shape<-1, 3>> tris_a,
    nb::ndarray<double, nb::shape<-1, 3>> verts_b,
    nb::ndarray<int32_t, nb::shape<-1, 3>> tris_b,
    ember::BooleanOp op,
    bool nsi, bool nnc)
{
    auto a = mesh_from_arrays(verts_a, tris_a, nsi, nnc);
    auto b = mesh_from_arrays(verts_b, tris_b, nsi, nnc);

    ember::EmberConfig cfg;
    cfg.assume_nsi = nsi;
    cfg.assume_nnc = nnc;

    auto t0 = std::chrono::high_resolution_clock::now();
    auto result = ember::boolean_operation({a, b}, op, cfg);
    auto soup   = ember::triangulate_and_resolve(result);
    auto t1 = std::chrono::high_resolution_clock::now();

    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    auto [v, t] = soup_to_arrays(soup);
    return {std::move(v), std::move(t), ms};
}

// ---------------------------------------------------------------------------
// profiled boolean — returns per-stage timings
// ---------------------------------------------------------------------------

struct ProfileResult {
    VertList verts;
    TriList  tris;
    double ms_input;       // mesh_from_arrays + prepare_input
    double ms_subdivide;   // subdivide (the core algorithm)
    double ms_triangulate; // triangulate_and_resolve
    double ms_output;      // soup_to_arrays conversion
    double ms_total;
};

static ProfileResult run_boolean_profiled(
    nb::ndarray<double, nb::shape<-1, 3>> verts_a,
    nb::ndarray<int32_t, nb::shape<-1, 3>> tris_a,
    nb::ndarray<double, nb::shape<-1, 3>> verts_b,
    nb::ndarray<int32_t, nb::shape<-1, 3>> tris_b,
    ember::BooleanOp op,
    bool nsi, bool nnc)
{
    using Clock = std::chrono::high_resolution_clock;
    auto wall0 = Clock::now();

    // Stage 1: Input preparation
    auto t0 = Clock::now();
    auto a = mesh_from_arrays(verts_a, tris_a, nsi, nnc);
    auto b = mesh_from_arrays(verts_b, tris_b, nsi, nnc);

    ember::EmberConfig cfg;
    cfg.assume_nsi = nsi;
    cfg.assume_nnc = nnc;

    ember::PolygonSoup soup = ember::prepare_input({a, b});
    if (cfg.assume_nsi)
        for (auto& p : soup.polygons) p.no_self_intersections = true;
    if (cfg.assume_nnc)
        for (auto& p : soup.polygons) p.no_nested_components = true;

    ember::IndicatorFn indicator = ember::make_indicator(op, soup.num_meshes);
    auto t1 = Clock::now();
    double ms_input = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // Stage 2: Subdivision (core algorithm)
    t0 = Clock::now();
    ember::pos_t ref_point = soup.bounds.min;
    ember::WNV ref_wnv(soup.num_meshes, 0);

    ember::SubdivisionTask initial_task;
    initial_task.polygons = std::move(soup.polygons);
    initial_task.bounds = soup.bounds;
    initial_task.root_bounds = soup.bounds;
    initial_task.ref_point = ref_point;
    initial_task.ref_wnv = ref_wnv;
    initial_task.depth = 0;

    std::vector<ember::ClassifiedPolygon> classified;
    ember::subdivide(std::move(initial_task), indicator, classified);
    t1 = Clock::now();
    double ms_subdivide = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // Assemble BooleanResult
    ember::BooleanResult result;
    result.output.num_meshes = soup.num_meshes;
    result.output.bounds = soup.bounds;
    result.output.transform = soup.transform;
    result.transform = soup.transform;
    for (auto& cp : classified) {
        if (cp.classification == -1) cp.polygon = cp.polygon.inverted();
        result.output.polygons.push_back(std::move(cp.polygon));
        result.classifications.push_back(cp.classification);
    }

    // Stage 3: Triangulate + T-junction resolve
    t0 = Clock::now();
    auto tri_soup = ember::triangulate_and_resolve(result);
    t1 = Clock::now();
    double ms_triangulate = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // Stage 4: Output conversion
    t0 = Clock::now();
    auto [v, t] = soup_to_arrays(tri_soup);
    t1 = Clock::now();
    double ms_output = std::chrono::duration<double, std::milli>(t1 - t0).count();

    auto wall1 = Clock::now();
    double ms_total = std::chrono::duration<double, std::milli>(wall1 - wall0).count();

    return {std::move(v), std::move(t), ms_input, ms_subdivide, ms_triangulate, ms_output, ms_total};
}

// ---------------------------------------------------------------------------
// NB_MODULE
// ---------------------------------------------------------------------------

NB_MODULE(ember_ext, m) {
    m.doc() = "EMBER Integer Exact CSG — Python bindings via nanobind";

    // Expose BooleanOp enum
    nb::enum_<ember::BooleanOp>(m, "BooleanOp")
        .value("Union",               ember::BooleanOp::Union)
        .value("Intersection",        ember::BooleanOp::Intersection)
        .value("Difference",          ember::BooleanOp::Difference)
        .value("SymmetricDifference", ember::BooleanOp::SymmetricDifference);

    // General boolean: (verts_a, tris_a, verts_b, tris_b, op, nsi, nnc) -> (verts, tris, ms)
    m.def("boolean", [](
        nb::ndarray<double, nb::shape<-1, 3>>  va,
        nb::ndarray<int32_t, nb::shape<-1, 3>> ta,
        nb::ndarray<double, nb::shape<-1, 3>>  vb,
        nb::ndarray<int32_t, nb::shape<-1, 3>> tb,
        ember::BooleanOp op,
        bool nsi, bool nnc) {
        return run_boolean(va, ta, vb, tb, op, nsi, nnc);
    }, "verts_a"_a, "tris_a"_a, "verts_b"_a, "tris_b"_a,
       "op"_a, "nsi"_a = true, "nnc"_a = true,
       "Perform a CSG boolean and return (vertices, triangles, elapsed_ms).");

    // Convenience wrappers
    m.def("boolean_union", [](
        nb::ndarray<double, nb::shape<-1, 3>>  va,
        nb::ndarray<int32_t, nb::shape<-1, 3>> ta,
        nb::ndarray<double, nb::shape<-1, 3>>  vb,
        nb::ndarray<int32_t, nb::shape<-1, 3>> tb,
        bool nsi, bool nnc) {
        return run_boolean(va, ta, vb, tb, ember::BooleanOp::Union, nsi, nnc);
    }, "verts_a"_a, "tris_a"_a, "verts_b"_a, "tris_b"_a,
       "nsi"_a = true, "nnc"_a = true,
       "CSG union. Returns (vertices, triangles, elapsed_ms).");

    m.def("boolean_intersection", [](
        nb::ndarray<double, nb::shape<-1, 3>>  va,
        nb::ndarray<int32_t, nb::shape<-1, 3>> ta,
        nb::ndarray<double, nb::shape<-1, 3>>  vb,
        nb::ndarray<int32_t, nb::shape<-1, 3>> tb,
        bool nsi, bool nnc) {
        return run_boolean(va, ta, vb, tb, ember::BooleanOp::Intersection, nsi, nnc);
    }, "verts_a"_a, "tris_a"_a, "verts_b"_a, "tris_b"_a,
       "nsi"_a = true, "nnc"_a = true,
       "CSG intersection. Returns (vertices, triangles, elapsed_ms).");

    m.def("boolean_difference", [](
        nb::ndarray<double, nb::shape<-1, 3>>  va,
        nb::ndarray<int32_t, nb::shape<-1, 3>> ta,
        nb::ndarray<double, nb::shape<-1, 3>>  vb,
        nb::ndarray<int32_t, nb::shape<-1, 3>> tb,
        bool nsi, bool nnc) {
        return run_boolean(va, ta, vb, tb, ember::BooleanOp::Difference, nsi, nnc);
    }, "verts_a"_a, "tris_a"_a, "verts_b"_a, "tris_b"_a,
       "nsi"_a = true, "nnc"_a = true,
       "CSG difference. Returns (vertices, triangles, elapsed_ms).");

    // Load OBJ file directly
    m.def("load_obj", [](const std::string& path, bool nsi, bool nnc)
        -> std::tuple<VertList, TriList>
    {
        auto mesh = ember::load_obj(path, nsi, nnc);
        VertList verts(mesh.positions.size());
        for (size_t i = 0; i < mesh.positions.size(); i++)
            verts[i] = {mesh.positions[i].x, mesh.positions[i].y, mesh.positions[i].z};

        TriList tris(mesh.triangles.size());
        for (size_t i = 0; i < mesh.triangles.size(); i++)
            tris[i] = {mesh.triangles[i].v0, mesh.triangles[i].v1, mesh.triangles[i].v2};

        return {std::move(verts), std::move(tris)};
    }, "path"_a, "nsi"_a = true, "nnc"_a = true,
       "Load an OBJ file. Returns (vertices, triangles).");

    // Compute volume from triangle soup (divergence theorem)
    m.def("volume", [](
        nb::ndarray<double, nb::shape<-1, 3>>  verts,
        nb::ndarray<int32_t, nb::shape<-1, 3>> tris) -> double
    {
        double vol = 0;
        for (size_t i = 0; i < (size_t)tris.shape(0); i++) {
            int i0 = tris(i, 0), i1 = tris(i, 1), i2 = tris(i, 2);
            double v0x = verts(i0, 0), v0y = verts(i0, 1), v0z = verts(i0, 2);
            double v1x = verts(i1, 0), v1y = verts(i1, 1), v1z = verts(i1, 2);
            double v2x = verts(i2, 0), v2y = verts(i2, 1), v2z = verts(i2, 2);
            vol += (v0x * (v1y * v2z - v1z * v2y)
                  + v0y * (v1z * v2x - v1x * v2z)
                  + v0z * (v1x * v2y - v1y * v2x)) / 6.0;
        }
        return std::abs(vol);
    }, "verts"_a, "tris"_a,
       "Compute signed volume of a closed triangle mesh.");

    // Profiled boolean: returns dict with per-stage timings
    m.def("boolean_profiled", [](
        nb::ndarray<double, nb::shape<-1, 3>>  va,
        nb::ndarray<int32_t, nb::shape<-1, 3>> ta,
        nb::ndarray<double, nb::shape<-1, 3>>  vb,
        nb::ndarray<int32_t, nb::shape<-1, 3>> tb,
        ember::BooleanOp op,
        bool nsi, bool nnc) {
        auto r = run_boolean_profiled(va, ta, vb, tb, op, nsi, nnc);
        return std::make_tuple(
            std::move(r.verts), std::move(r.tris),
            r.ms_input, r.ms_subdivide, r.ms_triangulate, r.ms_output, r.ms_total);
    }, "verts_a"_a, "tris_a"_a, "verts_b"_a, "tris_b"_a,
       "op"_a, "nsi"_a = true, "nnc"_a = true,
       "Profiled boolean. Returns (verts, tris, ms_input, ms_subdivide, ms_triangulate, ms_output, ms_total).");

    m.attr("__version__") = "0.1.0";
}
