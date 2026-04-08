# ember-csg

An open-source C++ implementation of **EMBER** (Exact Mesh Booleans via Efficient and Robust Local Arrangements), based on the paper by Nehring-Wirxel et al. from RWTH Aachen University.

> **Paper**: [EMBER: Exact Mesh Booleans via Efficient and Robust Local Arrangements](https://www.graphics.rwth-aachen.de/publication/03339/)
>
> Nehring-Wirxel, Nehring-Wirxel, Kobbelt. *ACM Transactions on Graphics (SIGGRAPH Asia 2024)*.

## Overview

EMBER performs exact Boolean operations (union, intersection, difference) on triangle meshes using integer homogeneous coordinates. All intermediate computation uses fixed-width integer arithmetic (up to 256 bits), ensuring bitwise-exact results with no floating-point robustness issues.

This implementation uses the `integer-plane-geometry` library from the [mesh-kernel](https://github.com/jnehringwirxel/mesh-kernel) project for the exact arithmetic foundation.

### Features

- **Exact integer arithmetic**: 26-bit positions, 55-bit normals, all intermediates ≤ 256 bits
- **Plane-based geometry**: polygons defined by supporting plane + edge planes
- **Face-local BSP**: pairwise intersection resolution via per-polygon binary space partitions
- **WNV classification**: winding number vectors computed via robust multi-ray casting
- **Boolean operations**: union, intersection, difference, symmetric difference
- **OBJ I/O**: load and export triangle meshes

### Validated Results

Comparison against manifold3d on sphere and cube meshes:

| Operation | Volume Delta | Surface Area Delta |
|-----------|-------------|-------------------|
| Cube union | 0.00% | 0.00% |
| Cube intersection | 0.00% | 0.00% |
| Cube difference | 0.00% | 0.00% |
| Sphere union | +3.16% | +2.45% |
| Sphere intersection | -93.8%* | -92.3%* |
| Sphere difference | +2.21% | +0.08% |

*Sphere intersection accuracy is limited by BSP splitting on curved surfaces at the thin lens-shaped intersection region. Without subdivision, all sphere metrics match to within 0.02%.

## Building

Requires C++20, CMake 3.14+, and the mesh-kernel submodules (typed-geometry, clean-core):

```bash
cd mesh-kernel
git submodule update --init extern/typed-geometry extern/clean-core

cd ../ember-csg
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCC_LINKER_MOLD=OFF
make -j$(nproc)
```

## Usage

```cpp
#include <ember/ember.hh>

// Create input meshes
ember::InputMesh mesh_a, mesh_b;
// ... populate positions and triangles ...
mesh_a.nsi = true;  // No Self-Intersections
mesh_a.nnc = true;  // No Nested Components

// Boolean union
auto result = ember::boolean_union(mesh_a, mesh_b);

// Export as OBJ
auto soup = ember::triangulate_output(result);
auto obj_string = ember::to_obj(soup);

// Or load from OBJ
auto loaded = ember::load_obj("input.obj", /*nsi=*/true, /*nnc=*/true);
```

## Performance

| Operation | Input Triangles | Time |
|-----------|----------------|------|
| Cube-Cube union | 24 | 71 ms |
| Sphere-Sphere union | 1024 | 1.1 s |
| Sphere-Sphere intersection | 1024 | 1.0 s |
| Rounded cube (cube ∩ sphere) | 524 | 3.4 s |

Performance is dominated by the O(n²) pairwise intersection testing in leaf cells. The EMBER paper achieves 1.6ms for ~20K faces through work-stealing parallelism, K-DOP culling, and optimized splitting heuristics.

## Current Limitations

1. **Subdivision creates duplicate fragments**: When the kd-tree subdivision clips a polygon at a split boundary, both halves are independently classified. This can create redundant overlapping geometry. Disabling subdivision (setting `LEAF_THRESHOLD` high) produces exact results but is O(n²).

2. **T-junctions at BSP boundaries**: BSP-split vertices may not align exactly with adjacent polygon edges, creating small gaps. The EMBER paper acknowledges this: *"output not triangulated, contains T-junctions and convex polygons."*

3. **Non-manifold output**: The OBJ export produces polygon soup without vertex sharing. Cube outputs become manifold after vertex merging; sphere outputs retain boundary edges from T-junctions.

4. **No parallelism**: The current implementation is single-threaded. The paper's work-stealing architecture is not yet implemented.

5. **No K-DOP culling**: All pairwise intersections are tested. K-DOP bounding volumes would skip many non-intersecting pairs.

6. **WNV via multi-ray casting**: The paper uses exact segment tracing with 3-segment paths. This implementation uses double-precision multi-ray casting (7 rays with majority voting) for the final classification step, which is robust but not exact.

## Architecture

```
include/ember/
├── types.hh              Core geometry types (geometry<26, 55>)
├── exact_classify.hh     Sign-extension-safe point-vs-plane classify
├── polygon.hh            Convex polygon (support plane + edge planes)
├── mesh.hh               Input/output mesh, coordinate scaling
├── winding.hh            Winding number vectors, boolean indicators
├── clip.hh               Polygon clipping against planes
├── intersect_polygons.hh Pairwise polygon intersection (C1-C4 cases)
├── local_bsp.hh          Face-local BSP tree
├── segment_trace.hh      WNV classification via multi-ray casting
├── subdivision.hh        Recursive adaptive subdivision
└── ember.hh              Public API
```

## References

- Nehring-Wirxel, Nehring-Wirxel, Kobbelt. *"EMBER: Exact Mesh Booleans via Efficient and Robust Local Arrangements."* ACM Transactions on Graphics, 2024. [Project page](https://www.graphics.rwth-aachen.de/publication/03339/)
- The `integer-plane-geometry` library from [mesh-kernel](https://github.com/juliusnehring/mesh-kernel/)
- [typed-geometry](https://github.com/project-arcana/typed-geometry) for fixed-width integer arithmetic
- [manifold3d](https://github.com/elalish/manifold) for comparison and validation

## License

MIT License. See [mesh-kernel](https://github.com/juliusnehring/mesh-kernel/) (also MIT) for the integer-plane-geometry foundation.
