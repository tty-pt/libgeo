#ifndef GEO_H
#define GEO_H

/**
 * @file geo.h
 * @brief Public API for libgeo — spatial indexing with Morton codes.
 *
 * Libgeo provides efficient spatial database operations using Morton codes
 * (Z-order space-filling curves) for multi-dimensional coordinate indexing.
 * Built on top of libqmap for persistence and hash table operations.
 *
 * Coordinates are signed 16-bit integers (int16_t) ranging from -32768 to 32767,
 * suitable for game worlds, voxel engines, and spatial simulations.
 *
 * @note Depends on libqmap >= 0.6.0 and libqsys.
 */

#include <stdint.h>
#include <ttypt/qmap.h>

#include "morton.h"

/** @defgroup geo_core Geo core API
 *  @brief Spatial map operations using Morton code indexing.
 *
 *  Libgeo provides a spatial database that maps multi-dimensional coordinates
 *  to arbitrary uint32_t values. Internally, coordinates are converted to
 *  Morton codes (Z-order) for efficient spatial queries and storage.
 *
 *  Coordinate System:
 *  - Type: int16_t (signed 16-bit integers)
 *  - Range: -32768 to 32767 per dimension
 *  - Currently optimized for 3D (dim=3), designed for future N-dimensional support
 *
 *  @note Thread Safety: Libgeo inherits libqmap's thread-safety properties.
 *        It uses global state and is NOT thread-safe. Use external
 *        synchronization if accessing from multiple threads.
 *
 *  @note Capacity: Map capacity is fixed at creation via the mask parameter.
 *        Capacity = mask + 1. Attempting to exceed capacity terminates the
 *        process. There is no dynamic resizing.
 *
 *  @note Memory: Malloc failures terminate the process immediately via CBUG().
 *        There is no graceful error handling for out-of-memory conditions.
 *
 *  @warning File Persistence: File-backed databases are automatically saved
 *           at process exit via libqmap. Explicit qmap_save() calls are only
 *           needed for mid-execution checkpointing.
 *
 *  @see geo_morton
 *  @see geo_point
 *  @{
 */

/**
 * @brief Sentinel value returned by geo_get when no entry exists.
 *
 * This value (130056652770671ULL) is chosen to be outside the valid
 * Morton code space, ensuring it never conflicts with actual stored values.
 * Different from QM_MISS (UINT32_MAX) because Morton codes are uint64_t.
 *
 * Usage:
 * @code
 * uint32_t value = geo_get(db, point, 3);
 * if (value == GEO_MISS) {
 *     // No entry at this coordinate
 * }
 * @endcode
 */
#define GEO_MISS (130056652770671ULL)

/**
 * @brief Initialize the geo subsystem.
 *
 * Registers custom types with libqmap (uint64_t for Morton codes, uint32_t
 * for values) and sets up the Morton code comparator for sorted iteration.
 * Also initializes the internal IDM (ID Manager) for iterator handles.
 *
 * @warning Must be called before any other geo functions. Calling other
 *          functions without initialization results in undefined behavior.
 *
 * @note This function can be called multiple times safely (idempotent if
 *       qmap types are already registered).
 *
 * Example:
 * @code
 * int main() {
 *     geo_init();  // Always call first
 *     uint32_t db = geo_open(NULL, NULL, 0xFF);
 *     // ... use database
 *     return 0;
 * }
 * @endcode
 *
 * @see geo_open
 */
void geo_init(void);

/**
 * @brief Open or create a spatial database.
 *
 * Creates a spatial map backed by libqmap with QM_SORTED for efficient
 * range queries. Internally uses Morton codes (uint64_t) as keys and
 * uint32_t for values. If a filename is provided, data is automatically
 * loaded from disk if the file exists.
 *
 * @param[in] filename File path for persistence, or NULL for in-memory only.
 *                     Example: "world.db"
 * @param[in] database Logical database name within the file, or NULL to
 *                     skip file association. Multiple databases can share
 *                     one file. Example: "main", "players", "chunks"
 * @param[in] mask     Hash table mask (capacity = mask + 1). Must be 2^n - 1
 *                     for optimal performance. Examples: 0xFF (256 entries),
 *                     0xFFF (4096 entries), 0xFFFF (65536 entries)
 *
 * @return Database handle for use with other geo functions.
 *
 * @note File Persistence: File-backed databases automatically load existing
 *       data when opened. Data is automatically saved at process exit via
 *       libqmap's destructor. Call qmap_save() explicitly only if you need
 *       mid-execution checkpointing.
 *
 * @note Capacity: The mask parameter determines maximum capacity. This cannot
 *       be changed after creation. Choose a mask that accommodates your
 *       expected data size. Exceeding capacity terminates the process.
 *
 * @note The database is created with QM_SORTED flag, enabling ordered
 *       iteration by Morton code for efficient spatial range queries.
 *
 * @warning The filename and database parameters are not copied. If providing
 *          string literals, they must remain valid for the database lifetime.
 *
 * Example (in-memory):
 * @code
 * geo_init();
 * uint32_t db = geo_open(NULL, NULL, 0xFF);  // 256-entry in-memory DB
 * @endcode
 *
 * Example (persistent):
 * @code
 * geo_init();
 * uint32_t db = geo_open("world.db", "main", 0xFFF);  // 4096-entry persistent DB
 * @endcode
 *
 * @see geo_init
 * @see qmap_open
 * @see qmap_save
 * @see qmap_close
 */
uint32_t geo_open(char *filename, char *database, uint32_t mask);

/**
 * @brief Create an iterator for all points within a rectangular bounding box.
 *
 * Queries the spatial database for all entries within the axis-aligned
 * bounding box defined by start point 's' and lengths 'l'. The iterator
 * pre-allocates and collects all matching results using Morton code range
 * queries internally.
 *
 * @param[in] pdb_hd Database handle from geo_open().
 * @param[in] s      Start point (minimum corner of bounding box).
 *                   Array of int16_t with at least 'dim' elements.
 * @param[in] l      Lengths of bounding box per dimension (unsigned).
 *                   Array of uint16_t with at least 'dim' elements.
 *                   The box extends from s[i] to s[i]+l[i] (exclusive end).
 * @param[in] dim    Number of dimensions. Currently optimized for 3D,
 *                   designed for future multi-dimensional support.
 *
 * @return Iterator handle for use with geo_next(). The handle is
 *         automatically freed when geo_next() returns 0.
 *
 * @note Results are collected in Morton code order (Z-order space-filling
 *       curve), NOT spatial order. Nearby points may not be adjacent in
 *       iteration sequence.
 *
 * @note Memory: Allocates an array sized to the bounding box volume
 *       (product of all lengths). Large volumes may consume significant
 *       memory. The volume is calculated as l[0] * l[1] * ... * l[dim-1].
 *
 * @warning The volume calculation can overflow for very large boxes.
 *          Reasonable box sizes are recommended (< 1 million points).
 *
 * Example (2D region):
 * @code
 * int16_t start[2] = {0, 0};
 * uint16_t lengths[2] = {100, 100};  // 100x100 region
 * uint32_t iter = geo_iter(db, start, lengths, 2);
 * @endcode
 *
 * Example (3D region):
 * @code
 * int16_t start[3] = {-10, -10, -10};
 * uint16_t lengths[3] = {20, 20, 20};  // 20x20x20 cube (8000 points)
 * uint32_t iter = geo_iter(db, start, lengths, 3);
 * @endcode
 *
 * @see geo_next
 * @see morton_set
 */
uint32_t geo_iter(uint32_t pdb_hd, int16_t *s,
		uint16_t *l, uint8_t dim);

/**
 * @brief Advance iterator and retrieve the next point/value pair.
 *
 * Retrieves the next entry from the iterator created by geo_iter().
 * Skips empty cells in the bounding box. When all entries have been
 * returned (or the box was empty), returns 0 and automatically frees
 * the iterator's internal memory.
 *
 * @param[out] p   Output point array. Must have space for at least 'dim'
 *                 int16_t elements (where dim was passed to geo_iter()).
 *                 Filled with the coordinates of the next point.
 * @param[out] ref Output pointer for the stored value (uint32_t).
 *                 Filled with the value stored at this point.
 * @param[in]  cur Iterator handle from geo_iter().
 *
 * @return 1 if a point/value pair was retrieved (output written to p and ref).
 *         0 when iteration is complete (no more entries). The iterator is
 *         automatically freed on return 0.
 *
 * @note The iterator is stateful and maintains position. Each call advances
 *       to the next entry.
 *
 * @note Empty cells in the bounding box are automatically skipped.
 *       Only coordinates with stored values are returned.
 *
 * @warning After this function returns 0, the iterator handle becomes invalid.
 *          Do not call geo_next() again with the same handle.
 *
 * @warning The iterator's internal memory is freed when this returns 0.
 *          There is no separate "close" or "free" function.
 *
 * Typical usage pattern:
 * @code
 * int16_t start[3] = {0, 0, 0};
 * uint16_t lengths[3] = {10, 10, 10};
 * uint32_t iter = geo_iter(db, start, lengths, 3);
 * 
 * int16_t point[3];
 * uint32_t value;
 * while (geo_next(point, &value, iter)) {
 *     printf("Point (%d,%d,%d) = %u\n",
 *            point[0], point[1], point[2], value);
 * }
 * // Iterator is automatically freed after loop
 * @endcode
 *
 * @see geo_iter
 */
int geo_next(int16_t *p, uint32_t *ref, uint32_t cur);

/**
 * @brief Delete the entry at a spatial coordinate.
 *
 * Removes the value stored at the given point from the database.
 * Internally converts the coordinate to a Morton code and calls qmap_del().
 * If no entry exists at the coordinate, this is a no-op (safe to call).
 *
 * @param[in] pdb_hd Database handle from geo_open().
 * @param[in] p      Point coordinate. Array of int16_t with at least
 *                   'dim' elements. Coordinates are signed 16-bit integers.
 * @param[in] dim    Number of dimensions.
 *
 * @note This operation invalidates any pointers obtained from geo_get()
 *       or geo_next() that refer to this coordinate.
 *
 * @note Safe to call on non-existent coordinates (no error, no effect).
 *
 * Example:
 * @code
 * int16_t pos[3] = {10, 20, 30};
 * geo_del(db, pos, 3);  // Remove entry at (10,20,30)
 * @endcode
 *
 * @see geo_get
 * @see geo_put
 * @see qmap_del
 */
static inline void
geo_del(uint32_t pdb_hd, int16_t *p, uint8_t dim)
{
	uint64_t at = morton_set(p, dim);
	qmap_del(pdb_hd, &at);
}

/**
 * @brief Retrieve the value stored at a spatial coordinate.
 *
 * Looks up the value at the given point in the database. Internally
 * converts the coordinate to a Morton code and calls qmap_get().
 *
 * @param[in] pdb_hd Database handle from geo_open().
 * @param[in] p      Point coordinate. Array of int16_t with at least
 *                   'dim' elements. Coordinates are signed 16-bit integers
 *                   ranging from -32768 to 32767.
 * @param[in] dim    Number of dimensions.
 *
 * @return The stored uint32_t value at this coordinate, or GEO_MISS
 *         (130056652770671ULL) if no entry exists at this point.
 *
 * @note GEO_MISS is different from QM_MISS because Morton codes use
 *       uint64_t space. The value 130056652770671ULL is guaranteed not
 *       to conflict with valid uint32_t stored values.
 *
 * @note The returned value is a copy, not a pointer. Unlike qmap_get()
 *       which returns a pointer, geo_get() returns the actual uint32_t value.
 *
 * Example:
 * @code
 * int16_t pos[3] = {10, 20, 30};
 * uint32_t value = geo_get(db, pos, 3);
 * if (value == GEO_MISS) {
 *     printf("No entry at (10,20,30)\n");
 * } else {
 *     printf("Value at (10,20,30) = %u\n", value);
 * }
 * @endcode
 *
 * @see geo_put
 * @see geo_del
 * @see GEO_MISS
 * @see qmap_get
 */
static inline uint32_t
geo_get(uint32_t pdb_hd, int16_t *p, uint8_t dim)
{
	uint64_t at = morton_set(p, dim);
	const void *ref = qmap_get(pdb_hd, &at);

	if (ref)
		return * (uint32_t *) ref;

	return QM_MISS;
}

/**
 * @brief Store a value at a spatial coordinate.
 *
 * Inserts or updates the value at the given point in the database.
 * Internally converts the coordinate to a Morton code and calls qmap_put().
 * If an entry already exists at this coordinate, it is replaced.
 *
 * @param[in] pdb_hd Database handle from geo_open().
 * @param[in] p      Point coordinate. Array of int16_t with at least
 *                   'dim' elements. Coordinates are signed 16-bit integers
 *                   ranging from -32768 to 32767.
 * @param[in] thing  Value to store (uint32_t). Can be any 32-bit value
 *                   including 0. Avoid using QM_MISS (0xFFFFFFFF) as it
 *                   may cause confusion, though it's technically valid.
 * @param[in] dim    Number of dimensions.
 *
 * @note If the database is file-backed (filename provided to geo_open()),
 *       changes are automatically saved at process exit. Call qmap_save()
 *       explicitly for mid-execution persistence.
 *
 * @note Replacing an existing entry may invalidate pointers if the internal
 *       qmap allocation changes (though v0.6.0+ has allocation reuse).
 *
 * Example (store single value):
 * @code
 * int16_t pos[3] = {10, 20, 30};
 * geo_put(db, pos, 42, 3);  // Store value 42 at (10,20,30)
 * @endcode
 *
 * Example (update existing value):
 * @code
 * int16_t pos[3] = {10, 20, 30};
 * uint32_t old = geo_get(db, pos, 3);
 * if (old != GEO_MISS) {
 *     geo_put(db, pos, old + 1, 3);  // Increment
 * }
 * @endcode
 *
 * Example (populate a grid):
 * @code
 * for (int16_t x = 0; x < 10; x++) {
 *     for (int16_t y = 0; y < 10; y++) {
 *         int16_t pos[2] = {x, y};
 *         geo_put(db, pos, x * 10 + y, 2);
 *     }
 * }
 * @endcode
 *
 * @see geo_get
 * @see geo_del
 * @see qmap_put
 * @see qmap_save
 */
static inline void
geo_put(uint32_t pdb_hd, int16_t *p,
		uint32_t thing, uint8_t dim)
{
	uint64_t code = morton_set(p, dim);
	qmap_put(pdb_hd, &code, &thing);
}

/** @} */

#endif
