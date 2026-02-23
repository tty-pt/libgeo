# libgeo
> Spatial database library using Morton code indexing for efficient multi-dimensional coordinate queries.

A small library for spatial/geographic databases. Store and query data indexed by multi-dimensional coordinates (typically 3D) with efficient range queries powered by Morton codes (Z-order space-filling curves).

## Key Features

- **Morton Code Indexing**: Efficient spatial-to-linear mapping preserving locality
- **Fast Range Queries**: Query rectangular regions with O(log n + k) complexity
- **File Persistence**: Optional disk storage via libqmap
- **Flexible Dimensions**: Optimized for 3D, designed for N-dimensional support
- **Simple C API**: Minimal, easy-to-use interface
- **Sparse Data Efficient**: Only stores occupied coordinates

## Quick Start

```c
#include <ttypt/geo.h>

int main() {
    // Initialize (required)
    geo_init();
    
    // Create database (256 entry capacity)
    uint32_t db = geo_open(NULL, NULL, 0xFF);
    
    // Store value at 3D coordinate
    int16_t pos[3] = {10, 20, 30};
    geo_put(db, pos, 42, 3);
    
    // Retrieve value
    uint32_t value = geo_get(db, pos, 3);
    if (value != GEO_MISS) {
        printf("Value: %u\n", value);
    }
    
    // Query rectangular region
    int16_t start[3] = {0, 0, 0};
    uint16_t lengths[3] = {50, 50, 50};
    uint32_t iter = geo_iter(db, start, lengths, 3);
    
    int16_t point[3];
    uint32_t ref;
    while (geo_next(point, &ref, iter)) {
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

- **Type**: `int16_t` (signed 16-bit integers)
- **Range**: -32768 to 32767 per dimension
- **Dimensions**: Variable (2D, 3D, etc.), optimized for 3D
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
- `include/ttypt/morton.h` - Morton code utilities
- `include/ttypt/point.h` - Point arithmetic

## Performance Notes

- **Capacity**: Fixed at creation (mask parameter), no dynamic resizing
- **Mask Calculation**: capacity = mask + 1, use 2^n - 1 (e.g., 0xFF = 256)
- **Morton Codes**: Provide good spatial locality for cache-friendly access
- **Range Queries**: O(log n + k) where k is the number of results
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

| Function | Purpose |
|----------|---------|
| `geo_init()` | Initialize libgeo (call first) |
| `geo_open()` | Create/open spatial database |
| `geo_put()` | Store value at coordinate |
| `geo_get()` | Retrieve value from coordinate |
| `geo_del()` | Delete entry at coordinate |
| `geo_iter()` | Create region iterator |
| `geo_next()` | Advance iterator, get next point |
| `morton_set()` | Encode coordinate to Morton code |
| `morton_get()` | Decode Morton code to coordinate |
| `point_*()` | Vector/point utility functions |

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

Current version: **0.4.0**

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
