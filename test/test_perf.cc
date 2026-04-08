// Performance test for EMBER - simple benchmarks, under 15 seconds total
#include <ember/ember.hh>
#include <cstdio>
#include <chrono>

using Clock = std::chrono::high_resolution_clock;

ember::InputMesh make_cube(double cx, double cy, double cz, double s)
{
    ember::InputMesh m;
    m.nsi = true; m.nnc = true;
    m.positions = {
        {cx-s,cy-s,cz-s},{cx+s,cy-s,cz-s},{cx+s,cy+s,cz-s},{cx-s,cy+s,cz-s},
        {cx-s,cy-s,cz+s},{cx+s,cy-s,cz+s},{cx+s,cy+s,cz+s},{cx-s,cy+s,cz+s}};
    m.triangles = {{4,5,6},{4,6,7},{0,3,2},{0,2,1},{1,2,6},{1,6,5},
                   {0,4,7},{0,7,3},{3,7,6},{3,6,2},{0,1,5},{0,5,4}};
    return m;
}

int main()
{
    ember::EmberConfig cfg;
    cfg.assume_nsi = true;
    cfg.assume_nnc = true;

    std::printf("EMBER Performance Test\n");
    std::printf("======================\n\n");

    // Cube-Cube (should be instant)
    {
        auto a = make_cube(0, 0, 0, 1);
        auto b = make_cube(0.5, 0.5, 0.5, 1);
        auto t0 = Clock::now();
        auto r = ember::boolean_union(a, b, cfg);
        auto t1 = Clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1-t0).count();
        auto s = ember::triangulate_output(r);
        std::printf("  Cube-Cube union:         %7.1f ms  %4zu tris\n", ms, s.triangles.size());
    }

    // Sphere-Sphere (low res)
    auto sa = ember::load_obj("sphere_a.obj", true, true);
    auto sb = ember::load_obj("sphere_b.obj", true, true);
    if (!sa.triangles.empty())
    {
        auto t0 = Clock::now();
        auto r = ember::boolean_union(sa, sb, cfg);
        auto t1 = Clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1-t0).count();
        auto s = ember::triangulate_output(r);
        std::printf("  Sphere-Sphere union:     %7.1f ms  %4zu tris (%zu+%zu input)\n",
            ms, s.triangles.size(), sa.triangles.size(), sb.triangles.size());

        t0 = Clock::now();
        r = ember::boolean_intersection(sa, sb, cfg);
        t1 = Clock::now();
        ms = std::chrono::duration<double, std::milli>(t1-t0).count();
        s = ember::triangulate_output(r);
        std::printf("  Sphere-Sphere intersect: %7.1f ms  %4zu tris\n", ms, s.triangles.size());

        t0 = Clock::now();
        r = ember::boolean_difference(sa, sb, cfg);
        t1 = Clock::now();
        ms = std::chrono::duration<double, std::milli>(t1-t0).count();
        s = ember::triangulate_output(r);
        std::printf("  Sphere-Sphere diff:      %7.1f ms  %4zu tris\n", ms, s.triangles.size());
    }

    // Rounded cube (cube ^ sphere)
    auto cc = ember::load_obj("complex_cube.obj", true, true);
    auto cs = ember::load_obj("complex_sphere_big.obj", true, true);
    if (!cc.triangles.empty() && !cs.triangles.empty())
    {
        auto t0 = Clock::now();
        auto r = ember::boolean_intersection(cc, cs, cfg);
        auto t1 = Clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1-t0).count();
        auto s = ember::triangulate_output(r);
        std::printf("  Rounded cube (cube^sph): %7.1f ms  %4zu tris\n", ms, s.triangles.size());
        auto obj = ember::to_obj(s);
        if (auto* f = std::fopen("ember_rounded_cube.obj", "w")) { std::fputs(obj.c_str(), f); std::fclose(f); }
    }

    std::printf("\nDone.\n");
    return 0;
}
