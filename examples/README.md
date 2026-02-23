# Libgeo Examples

This directory contains example programs demonstrating various features of libgeo.

## Building Examples

From this directory:
```sh
make
```

Or build individual examples:
```sh
make basic
make iteration
make persistence
make 3d_world
```

## Example Programs

### 1. basic.c - Fundamental Operations

**What it demonstrates:**
- Initializing libgeo with `geo_init()`
- Creating in-memory databases with `geo_open()`
- Storing values with `geo_put()`
- Retrieving values with `geo_get()`
- Deleting entries with `geo_del()`
- Handling missing entries (`GEO_MISS`)
- Using both 2D and 3D coordinates

**Run:**
```sh
./basic
```

**Expected output:** Shows basic CRUD operations on spatial coordinates.

---

### 2. iteration.c - Spatial Queries

**What it demonstrates:**
- Populating a spatial database with multiple points
- Creating rectangular region queries with `geo_iter()`
- Iterating through results with `geo_next()`
- Morton code (Z-order) traversal effects
- Querying empty regions
- Handling sparse data efficiently
- 2D and 3D iteration examples

**Run:**
```sh
./iteration
```

**Expected output:** Various spatial queries showing iteration patterns and results.

---

### 3. persistence.c - File-Backed Databases

**What it demonstrates:**
- Creating file-backed databases for persistence
- Automatic data loading on reopen
- Explicit saving with `qmap_save()`
- Automatic saving at process exit
- Storing multiple logical databases in one file
- Verifying data persistence across runs

**Run:**
```sh
./persistence
```

**Expected output:** Creates `example_world.db` and `example_multi.db` files,
demonstrates data persisting across database close/reopen operations.

**Generated files:**
- `example_world.db` - Single database example
- `example_multi.db` - Multiple databases in one file

**Cleanup:**
```sh
rm example_*.db
```

---

### 4. 3d_world.c - Realistic Voxel World

**What it demonstrates:**
- Chunk-based terrain generation
- Different block types (grass, dirt, stone, water, etc.)
- View frustum / render distance queries
- Block placement and destruction
- World statistics gathering
- Practical use case for game development

**Run:**
```sh
./3d_world
```

**Expected output:** Generates a 3×3 chunk voxel world, demonstrates player views,
block manipulation, and shows performance with realistic data.

**Generated files:**
- `world.db` - Persistent voxel world database

**Cleanup:**
```sh
rm world.db
```

---

## Prerequisites

These examples require:
- libgeo (built from parent directory)
- libqmap >= 0.6.0
- libqsys
- libxxhash

Ensure you've built and installed the main libgeo library before compiling examples:

```sh
cd ..
make
sudo make install  # Or install to PREFIX if needed
cd examples
make
```

## Learning Path

Recommended order for learning:

1. **basic.c** - Start here to learn fundamental operations
2. **iteration.c** - Learn spatial queries and region traversal
3. **persistence.c** - Understand file-backed databases
4. **3d_world.c** - See realistic application patterns

## Tips

- All examples include detailed comments explaining each step
- Run with valgrind to verify no memory leaks:
  ```sh
  valgrind --leak-check=full ./basic
  ```
- Modify the examples to experiment with different parameters
- Check the man pages for detailed API documentation:
  ```sh
  man geo_open
  man geo_iter
  ```

## Common Issues

**Undefined reference errors:**
- Make sure libgeo is built: `cd .. && make`
- Check library paths in `Makefile` match your installation

**Segmentation faults:**
- Did you call `geo_init()` before other geo functions?
- Check that coordinate arrays have sufficient space for dimensions

**File permission errors:**
- Ensure write permissions in the current directory for persistence examples

## Further Reading

- Main README: `../README.md`
- API Documentation: `man geo_open`, `man geo_iter`, etc.
- Header files: `../include/ttypt/geo.h`, `morton.h`, `point.h`

## License

These examples are provided with the same license as libgeo.
Feel free to use them as starting points for your own projects.
