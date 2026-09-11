# libgeo
> Spatial database library using Morton code indexing for efficient multi-dimensional coordinate queries.

A small library for spatial/geographic databases. Store and query data indexed by multi-dimensional coordinates (typically 3D) with efficient range queries powered by Morton codes (Z-order space-filling curves).

## Key Features

- **Morton Code Indexing**: Efficient spatial-to-linear mapping preserving locality
- **Multi-Value Cells**: Several values can share one grid cell (append/read-all/replace/clear)
- **Fast Range Queries**: Query rectangular regions with O(log n + k) complexity
- **Z-Interval Skip**: Box walks jump whole empty aligned cubes instead of decoding every key
- **Kernel-Form Fill**: `Point3_2.fill_bbox()` streams a box straight into a sealed recall candidate set
- **File Persistence**: Optional disk storage via libqmap
- **Flexible Dimensions**: Optimized for 3D, designed for N-dimensional support
- **Simple C API**: Minimal, easy-to-use interface
- **Sparse Data Efficient**: Only stores occupied coordinates

## Quick Start

```c
#include <ttypt/pointcfg.h>

int main() {
    // Initialize (required)
    geo_init();
    
    // Create database (256-entry initial table, auto-grows)
    uint32_t db = geo_open(NULL, NULL, 0xFF);
    
    // Store value at 3D coordinate
    int16_t pos[3] = {10, 20, 30};
    Point3_2.put(db, pos, 42);
    
    // Retrieve value
    uint32_t value = Point3_2.get(db, pos);
    if (value != GEO_MISS) {
        printf("Value: %u\n", value);
    }
    
    // Query rectangular region
    int16_t start[3] = {0, 0, 0};
    uint16_t lengths[3] = {50, 50, 50};
    uint32_t iter = Point3_2.iter(db, start, lengths);
    
    int16_t point[3];
    uint32_t ref;
    while (Point3_2.next(point, &ref, iter)) {
        printf("Point (%d,%d,%d) = %u\n",
               point[0], point[1], point[2], ref);
    }
    
    return 0;
}
```

**Compile:**
```sh
cc -o myapp myapp.c -lgeo -lqmap -lqsys -lxxhash
```

## Installation

Check out [these instructions](https://github.com/tty-pt/ci/blob/main/docs/install.md#install-ttypt-packages) and use "libgeo" as the package name.

**Dependencies:**
- libqmap >= 0.6.0
- libqsys
- libxxhash

## Coordinate System

Two lane widths, picked per config object (`Point<D>_<B>`, B = bytes/lane):

- **2-byte lanes** (`Point1_2`..`Point4_2`): `int16_t`, -32768 to 32767
  per dimension, 1..4 dimensions (optimized for 3D)
- **4-byte lanes** (`Point2_4`): `int32_t`, -2147483648 to 2147483647
  per dimension, 2 dimensions, dense codec filling all 64 key bits
- **Use Cases**: Game worlds, voxel engines, spatial simulations, particle systems

Coordinates are internally converted to Morton codes (Z-order) for efficient storage and spatial queries.

## Use Cases

- **Voxel/Block Games**: Minecraft-like worlds, efficient chunk management
- **Spatial Databases**: Store and query location-based data
- **Particle Systems**: Manage thousands of particles with spatial queries
- **Level-of-Detail**: Spatial queries for LOD selection
- **Collision Detection**: Broad-phase spatial partitioning
- **Any Application**: Needing fast N-dimensional range queries

## Examples

The `examples/` directory contains complete, documented programs:

- **basic.c**: Fundamental operations (put, get, del)
- **iteration.c**: Spatial queries and region iteration
- **persistence.c**: File-backed databases
- **3d_world.c**: Realistic voxel world example

Build and run examples:
```sh
cd examples
make
./basic
./iteration
./persistence
./3d_world
```

See `examples/README.md` for detailed descriptions.

## Documentation

Use the man pages for complete API documentation:
```sh
man geo_open    # Open/create databases
man geo_iter    # Spatial iteration
man geo_get     # Retrieve values
man geo_put     # Store values
man morton_set  # Morton code encoding
man point_add   # Point utilities
```

**Generate man pages:**
```sh
make docs
```

Man pages are generated from Doxygen comments in header files:
- `include/ttypt/geo.h` - Main API
- `include/ttypt/pointcfg.h` - Config objects (recommended surface)
- `include/ttypt/morton.h` - Morton code utilities
- `include/ttypt/point.h` - Point arithmetic

## Performance Notes

- **Capacity**: The mask sizes the initial table (capacity = mask + 1);
  the map **auto-grows by doubling** when it fills — no termination, no resizing API
- **Mask Choice**: use 2^n - 1 near your expected entry count (e.g., 0xFF = 256)
- **Morton Codes**: Provide good spatial locality for cache-friendly access
- **Range Queries**: O(log n) seek + one decode per examined entry + O(k) visits;
  the Z-interval skip jumps whole empty aligned cubes, so dense intervals
  decode far fewer entries than they span
- **Sparse Data**: Only allocated coordinates consume memory
- **Memory**: Inherits qmap overhead (~32 bytes/entry + key/value sizes)

## Thread Safety

⚠️ **NOT thread-safe** - Inherits libqmap's global state limitations.

Use external synchronization (mutexes) if accessing from multiple threads.

## File Persistence

File-backed databases (when filename is provided to `geo_open()`):
- **Automatic Loading**: Data loads from disk on open
- **Automatic Saving**: Data saves to disk at process exit
- **Manual Save**: Call `qmap_save()` for mid-execution persistence
- **Multiple Databases**: Store multiple logical databases in one file

Example:
```c
// Open persistent database
uint32_t db = geo_open("world.db", "main", 0xFFFF);

// ... modify data ...

// Optional: save before exit
qmap_save();

// Automatic: saves on process exit anyway
```

## API Overview

The recommended surface is one config object per point type
(`include/ttypt/pointcfg.h`) — a struct of static methods, one per
operation, so a language server autocompletes the whole set off a single
symbol. `<D>` = dimensions, `<B>` = bytes per lane (2 = int16, 4 = int32):

```c
#include <ttypt/pointcfg.h>

int16_t p[3] = {10, 20, 30};
Point3_2.put(db, p, 42);             // store 42 at the cell
uint32_t v = Point3_2.get(db, p);    // read it back

int16_t s[3] = {0, 0, 0};
uint16_t l[3] = {4, 4, 4};
uint32_t cur = Point3_2.iter(db, s, l);
while (Point3_2.next(p, &v, cur))    // advance the box iterator
    ;
```

| Object | Config | Lanes | Box lengths | Advance |
|--------|--------|-------|-------------|---------|
| `Point1_2` | 1D, 2-byte lanes | `int16_t[1]` | `uint16_t[1]` | `.next` = `geo_next` |
| `Point2_2` | 2D, 2-byte lanes | `int16_t[2]` | `uint16_t[2]` | `.next` = `geo_next` |
| `Point3_2` | 3D, 2-byte lanes | `int16_t[3]` | `uint16_t[3]` | `.next` = `geo_next` |
| `Point4_2` | 4D, 2-byte lanes | `int16_t[4]` | `uint16_t[4]` | `.next` = `geo_next` |
| `Point2_4` | 2D, 4-byte lanes (dense full-key codec) | `int32_t[2]` | `int32_t[2]` | `.next` = `geo_next32` |

Each object has the same ~20 members: `morton_set`, `morton_get`
(codec), `add`, `sub`, `min`, `max`, `copy`, `vol`, `set` (broadcast),
`debug`, `idx` (point utils), `put`, `get`, `replace`, `del`, `del_all`,
`cell_count`, `get_multi` (+ shared `geo_cell_next()` to drain it), `iter`
(+ fused `.next` to advance), `fill_bbox` (sealed recall candidate set).

Member names mirror the flat functions one-to-one (`put` = `geo_put_N`,
`idx` = `point_idx_N`, …) with a single deliberate exception: the DB
replace op is `.replace`, because `.set` is already the point broadcast
(`point_set_N`) — the flat name for it is `geo_set_N`. The full mapping
lives in `include/ttypt/pointcfg.h`.

Shared setup/globals (not per-config):
`geo_init()`, `geo_open()`, `geo_last_scan_count()`, `geo_cell_next()`,
`morton_set_bulk()` / `morton_set_bulk4()` (SIMD batch encode).

One config per database: the int16 configs (`Point1_2..Point4_2`) and
`Point2_4` share the uint64 key space with different layouts, and the file
carries no config tag — never reopen a database with another config's
functions.

Low-level flat functions (ABI + tight-loop fast path, same behavior):
`morton_set_1..4()` / `morton_get_1..4()`, `point_*_1..4()`,
`geo_put/set/get/del/del_all/cell_count/get_multi/iter/fill_..._1..4()`,
the `geo_*_2_32()` 32-bit family, and the `geo_ops[1..4]` runtime-dim table.
Use these directly in hot per-cell loops: the config objects add one
indirect call per member (measured ~6.5x slower on a bare codec
round-trip; see docs/PERF.md), while scatter DB ops are qmap-dominated
and unaffected in practice.

## Kernel Form (recall composition)

`Point3_2.fill_bbox()` is the space-axis adapter for recall-style
composition: it streams every value in a bounding box into a
`rec_set_t` candidate set and seals it (sorted + deduplicated), without
ever materializing the box volume:

```c
#include <ttypt/pointcfg.h>
#include <ttypt/rec.h>

rec_set_t *cands = rec_set_new();
int16_t s[3] = {0, 0, 0};
uint16_t l[3] = {16, 16, 16};
if (Point3_2.fill_bbox(db, s, l, cands) == 0) {
    size_t n = rec_set_count(cands);      /* distinct refs, sorted */
    const rec_ref_t *refs = rec_set_at(cands);
    /* ... intersect/union with other axes, rank, fetch ... */
}
rec_set_free(cands);
```

Boxes larger than `GEO_FILL_MAX_VOL` (1M cells) are rejected with `-1`.
Requires libqmap >= 0.8.0 (multi-value chains + `rec.h`).

This follows the recall-kernel adapter contract
(`docs/RECALL-KERNEL.md` in libqmap): one `int rec_axis_fill_*(params,
rec_set_t *out)` that streams matches into the set and seals it, plain
`int` return (0 ok / -1 error), additive — the standard
`geo_iter`/`geo_next` cursor and `geo_search` raw path remain. The ref is
libgeo's stored `uint32` cell value widened to `rec_ref_t`; the kernel
never interprets it. Because both the raw iterator and the fill share
`geo_box_visit`, the Z-interval skip speeds the fill too, and fills never
build the box-volume array `geo_iter` collects.

## Building from Source

```sh
git clone <repository-url>
cd geo
make
sudo make install
```

Or install to custom prefix:
```sh
make PREFIX=$HOME/.local install
```

## Testing

The test suite includes 89+ tests covering:
- Unit tests (65 tests)
- Integration tests (13 tests)
- Property-based tests (12K+ assertions)
- Stress tests (11 tests)
- Performance benchmarks (4)

Run all tests:
```sh
cd tests
make clean test
```

Run benchmarks:
```sh
make bench
```

Run fuzz tests (requires libFuzzer):
```sh
make fuzz-standalone  # Standalone mode
```

See `tests/TESTING.md` for detailed documentation.

## Project Structure

```
geo/
├── include/ttypt/    # Public headers
│   ├── geo.h        # Main API
│   ├── morton.h     # Morton code utilities
│   └── point.h      # Point arithmetic
├── src/             # Implementation
│   └── libgeo.c
├── examples/        # Example programs
├── man/             # Generated man pages (via make docs)
├── CHANGELOG.md     # Version history
└── README.md        # This file
```

## Version

Current version: **0.4.1**

See `CHANGELOG.md` for version history and changes.

## License

See LICENSE file for details.

## Contributing

Issues and pull requests welcome. Please maintain code style and add tests for new features.

## See Also

- **libqmap**: Underlying hash table and persistence layer
- **Morton Codes**: https://en.wikipedia.org/wiki/Z-order_curve
- **Spatial Indexing**: Multi-dimensional range query paper referenced in source

## Contact

For bugs, questions, or suggestions, please open an issue on the project repository.
