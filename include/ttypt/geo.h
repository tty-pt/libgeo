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
 * @note Depends on libqmap >= 0.8.0 (per-key multi-value chains,
 *       qmap_get_multi, rec.h) and libqsys.
 */

#include <stdint.h>
#include <ttypt/qmap.h>
#include <ttypt/rec.h>

/* Optimization tunable — a 0/1 flag. GEO_SIMD_MORTON gates the batch
 * encoders (morton_set_bulk / morton_set_bulk4); it defaults to 1.
 * Override with -DGEO_SIMD_MORTON=0/1 (note: -U does NOT work — this
 * block re-defines an undefined macro to its default). The former
 * per-path tunables (INLINE/HOIST/CLZ/UNROLL) were promoted to
 * unconditional after measurement showed the generic fallbacks never
 * win; see docs/PERF.md. */
#ifndef GEO_SIMD_MORTON
#define GEO_SIMD_MORTON 1
#endif

#include "morton.h"

/** @defgroup geo_core Geo core API
 *  @brief Spatial map operations using Morton code indexing.
 *
 *  Libgeo provides a spatial database that maps multi-dimensional coordinates
 *  to arbitrary uint32_t values. Internally, coordinates are converted to
 *  Morton codes (Z-order) for efficient spatial queries and storage.
 *
 *  Value Semantics (multi-value cells):
 *  - A grid cell may hold MULTIPLE values (QM_MULTIVALUE map). geo_put()
 *    appends; geo_get() returns the first value; geo_get_multi() iterates
 *    all values in insertion order; geo_del() removes the first value;
 *    geo_del_all() removes every value; geo_set() replaces all values
 *    with a single new one.
 *
 *  Coordinate System:
 *  - Type: int16_t (signed 16-bit integers)
 *  - Range: -32768 to 32767 per dimension
 *  - Dimensions 1..4 supported; optimized for 3D (dim=3) with 4D
 *    (dim=4) support. 1D/2D/3D codes are bit-identical to v0.5.0;
 *    4D codes densely fill all 64 key bits.
 *
 *  @note Keyspace: all dimensions share the uint64 Morton key space.
 *        Use one dimension per database.
 *
 *  @note Thread Safety: Libgeo inherits libqmap's thread-safety properties.
 *        It uses global state and is NOT thread-safe. Use external
 *        synchronization if accessing from multiple threads.
 *
 *  @note Capacity: The mask parameter sizes the initial hash table
 *        (capacity = mask + 1). The map AUTO-GROWS by doubling when it
 *        fills (inherited from libqmap; libgeo never passes QM_NOGROW),
 *        so exceeding the initial capacity is safe. Choose a mask near
 *        your expected entry count to avoid early regrows.
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
 * Sentinel value returned by geo_get when no entry exists.
 */
#define GEO_MISS UINT32_MAX

/**
 * Maximum bounding-box volume (in cells) accepted by rec_axis_fill_bbox().
 * Boxes whose volume exceeds this are rejected with -1: they are almost
 * certainly a caller bug, and the sealed set would be huge. Raw
 * geo_iter() carries no such cap — it only allocates for entries found.
 */
#define GEO_FILL_MAX_VOL 1048576u

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
 * @param[in] mask     Hash table mask (initial capacity = mask + 1).
 *                     Must be 2^n - 1 for optimal performance. Examples:
 *                     0xFF (256 entries), 0xFFF (4096 entries),
 *                     0xFFFF (65536 entries). The map auto-grows past this;
 *                     pick a mask near your expected entry count.
 *
 * @return Database handle for use with other geo functions.
 *
 * @note File Persistence: File-backed databases automatically load existing
 *       data when opened. Data is automatically saved at process exit via
 *       libqmap's destructor. Call qmap_save() explicitly only if you need
 *       mid-execution checkpointing.
 *
 * @note Capacity: The mask sizes the initial table; the map auto-grows
 *       by doubling when it fills. Nothing terminates on growth.
 *
 * @note The database is created with QM_SORTED | QM_MULTIVALUE: ordered
 *       iteration by Morton code for efficient spatial range queries,
 *       with multi-value cells (see the Value Semantics note above).
 *
 * @warning The filename and database parameters are not copied. If providing
 *          string literals, they must remain valid for the database lifetime.
 *
 * Example (in-memory):
 * @code
 * geo_init();
 * uint32_t db = geo_open(NULL, NULL, 0xFF);  // 256-entry initial table
 * @endcode
 *
 * Example (persistent):
 * @code
 * geo_init();
 * uint32_t db = geo_open("world.db", "main", 0xFFF);  // 4096-entry initial table
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
 *                   The box covers s[i]..s[i]+l[i] inclusive per dimension.
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
 * @note Every stored (point, value) pair is returned exactly once, including
 *       multiple values stored at the same cell. Results are NOT deduplicated.
 *
 * @note Memory: Allocates a flat list that grows with the number of matching
 *       entries (not the bounding box volume), so sparse queries over large
 *       boxes stay cheap.
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
 * @brief Delete the first value stored at a spatial coordinate.
 *
 * Removes the first value stored at the given point from the database.
 * Internally converts the coordinate to a Morton code and calls qmap_del().
 * If no entry exists at the coordinate, this is a no-op (safe to call).
 * When several values share the cell, only the earliest-inserted one is
 * removed; use geo_del_all() to clear the cell.
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
 * @brief Delete every value stored at a spatial coordinate.
 *
 * Removes all values stored at the given point from the database.
 * Internally converts the coordinate to a Morton code and calls
 * qmap_del_all(). If no entry exists at the coordinate, this is a no-op
 * (safe to call, returns 0).
 *
 * @param[in] pdb_hd Database handle from geo_open().
 * @param[in] p      Point coordinate. Array of int16_t with at least
 *                   'dim' elements.
 * @param[in] dim    Number of dimensions.
 *
 * @return Number of values removed (0 if the cell was empty).
 *
 * Example:
 * @code
 * int16_t pos[3] = {10, 20, 30};
 * uint32_t n = geo_del_all(db, pos, 3);  // Clear the cell
 * @endcode
 *
 * @see geo_del
 * @see geo_set
 * @see qmap_del_all
 */
static inline uint32_t
geo_del_all(uint32_t pdb_hd, int16_t *p, uint8_t dim)
{
	uint64_t at = morton_set(p, dim);
	uint32_t n = qmap_count(pdb_hd, &at);

	if (n)
		qmap_del_all(pdb_hd, &at);

	return n;
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
 * @return The first stored uint32_t value at this coordinate, or GEO_MISS
 *         (UINT32_MAX) if no entry exists at this point. When several
 *         values share the cell, the earliest-inserted one is returned;
 *         use geo_get_multi() to retrieve them all.
 *
 * @note GEO_MISS equals UINT32_MAX (0xFFFFFFFF), the same as QM_MISS.
 *       This is the standard sentinel value for missing entries.
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

	return GEO_MISS;
}

/**
 * @brief Count the values stored at a spatial coordinate.
 *
 * @param[in] pdb_hd Database handle from geo_open().
 * @param[in] p      Point coordinate. Array of int16_t with at least
 *                   'dim' elements.
 * @param[in] dim    Number of dimensions.
 *
 * @return Number of values stored at this coordinate (0 if empty).
 *
 * @see geo_get
 * @see geo_get_multi
 * @see qmap_count
 */
static inline uint32_t
geo_cell_count(uint32_t pdb_hd, int16_t *p, uint8_t dim)
{
	uint64_t at = morton_set(p, dim);

	return qmap_count(pdb_hd, &at);
}

/**
 * @brief Diagnostic: index entries examined by the last box walk.
 *
 * Returns the number of index entries examined by the most recent
 * geo_iter()/rec_axis_fill_bbox() box walk in this process. Tests and
 * benchmarks use it to prove the Z-interval skip engages (entries
 * examined well below the morton-interval width on sparse boxes).
 *
 * @return Entry-examination count of the last box walk (0 if none ran yet).
 */
uint32_t geo_last_scan_count(void);

/**
 * @brief Create an iterator for all values stored at one coordinate.
 *
 * Opens a chain-aware cursor over every value stored at the given point,
 * in insertion order. Use geo_cell_next() to retrieve the values.
 *
 * @param[in] pdb_hd Database handle from geo_open().
 * @param[in] p      Point coordinate. Array of int16_t with at least
 *                   'dim' elements.
 * @param[in] dim    Number of dimensions.
 *
 * @return Iterator handle for use with geo_cell_next(), or QM_MISS when
 *         the cell holds no values.
 *
 * Example:
 * @code
 * int16_t pos[3] = {10, 20, 30};
 * uint32_t cur = geo_get_multi(db, pos, 3);
 * if (cur != QM_MISS) {
 *     uint32_t value;
 *     while (geo_cell_next(&value, cur))
 *         printf("value %u\n", value);
 * }
 * @endcode
 *
 * @see geo_cell_next
 * @see geo_cell_count
 */
uint32_t geo_get_multi(uint32_t pdb_hd, int16_t *p, uint8_t dim);

/**
 * @brief Advance a cell iterator and retrieve the next value.
 *
 * Retrieves the next value from the iterator created by geo_get_multi().
 * When all values have been returned, returns 0 and automatically frees
 * the iterator.
 *
 * @param[out] ref Output pointer for the stored value (uint32_t).
 * @param[in]  cur Iterator handle from geo_get_multi().
 *
 * @return 1 if a value was retrieved. 0 when iteration is complete; the
 *         iterator is automatically freed on return 0.
 *
 * @warning After this function returns 0, the iterator handle becomes invalid.
 *
 * @see geo_get_multi
 */
int geo_cell_next(uint32_t *ref, uint32_t cur);

/**
 * @brief Store a value at a spatial coordinate (append).
 *
 * Appends the value at the given point in the database. Internally
 * converts the coordinate to a Morton code and calls qmap_put().
 * If entries already exist at this coordinate, the new value is ADDED
 * alongside them (multi-value cell) — nothing is replaced. Use geo_set()
 * for replace semantics, geo_get_multi() to read all values back.
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
 *     geo_set(db, pos, old + 1, 3);  // Increment (replace)
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

/**
 * @brief Store a value at a spatial coordinate (replace).
 *
 * Replaces every value at the given point with a single new value:
 * geo_del_all() followed by geo_put(). This preserves the historical
 * single-value overwrite convenience on top of multi-value cells.
 *
 * @param[in] pdb_hd Database handle from geo_open().
 * @param[in] p      Point coordinate. Array of int16_t with at least
 *                   'dim' elements.
 * @param[in] thing  Value to store (uint32_t).
 * @param[in] dim    Number of dimensions.
 *
 * Example:
 * @code
 * int16_t pos[3] = {10, 20, 30};
 * geo_put(db, pos, 42, 3);
 * geo_set(db, pos, 99, 3);  // Cell now holds exactly {99}
 * @endcode
 *
 * @see geo_put
 * @see geo_del_all
 * @see geo_get
 */
static inline void
geo_set(uint32_t pdb_hd, int16_t *p,
		uint32_t thing, uint8_t dim)
{
	uint64_t code = morton_set(p, dim);

	qmap_del_all(pdb_hd, &code);
	qmap_put(pdb_hd, &code, &thing);
}

/**
 * @brief Fill a recall candidate set with the values in a bounding box.
 *
 * Kernel-form space-axis adapter: streams every value stored in the
 * axis-aligned box [s, s+l] into a recall candidate set, then seals it
 * (sort + dedup). Multiple values sharing one cell each enter the set;
 * sealing collapses exact duplicates. The walk never materializes the
 * box volume — sparse queries over large boxes stay cheap.
 *
 * @param[in] pdb_hd Database handle from geo_open().
 * @param[in] s      Start point (minimum corner). Array of int16_t.
 * @param[in] l      Lengths per dimension. Array of uint16_t.
 * @param[in] dim    Number of dimensions (1..4).
 * @param[out] out   Caller-owned candidate set; appended to, then sealed.
 *                   May already hold refs (result is the union, sealed).
 *
 * @return 0 on success (out sealed). -1 when out is NULL, dim is not
 *         1..4, or the box volume exceeds GEO_FILL_MAX_VOL (note: 4D
 *         volumes are 4-way products, so keep each side small).
 *
 * Example:
 * @code
 * rec_set_t *cands = rec_set_new();
 * int16_t s[3] = {0, 0, 0};
 * uint16_t l[3] = {16, 16, 16};
 * if (rec_axis_fill_bbox(db, s, l, 3, cands) == 0) {
 *     // rec_set_count(cands) distinct refs, sorted
 * }
 * rec_set_free(cands);
 * @endcode
 *
 * @see GEO_FILL_MAX_VOL
 * @see geo_iter
 */
int rec_axis_fill_bbox(uint32_t pdb_hd, int16_t *s,
		uint16_t *l, uint8_t dim, rec_set_t *out);

#if GEO_SIMD_MORTON
/**
 * @brief Batch-encode multiple 3D points to Morton codes using SIMD.
 *
 * Encodes up to 4 points (AVX2) or 2 points (NEON) in parallel,
 * falling back to scalar for remaining points. The output array must
 * have space for at least @p n elements.
 *
 * @param[out] out     Output Morton codes. Array of uint64_t with 'n' elements.
 * @param[in]  points  Input points. Array of int16_t[3] with 'n' entries.
 * @param[in]  n       Number of points to encode (0..UINT32_MAX).
 *
 * @return Number of points encoded (always == n).
 */
uint32_t morton_set_bulk(uint64_t *out, int16_t points[][3], uint32_t n);

/**
 * @brief Batch-encode multiple 4D points to Morton codes using SIMD.
 *
 * 4D analogue of morton_set_bulk(): encodes 4 points at a time under
 * AVX2 (dense stride-4 packing), scalar tail + fallback otherwise.
 * Output codes match morton_set() with dim=4 exactly.
 *
 * @param[out] out     Output Morton codes. Array of uint64_t with 'n' elements.
 * @param[in]  points  Input points. Array of int16_t[4] with 'n' entries.
 * @param[in]  n       Number of points to encode (0..UINT32_MAX).
 *
 * @return Number of points encoded (always == n).
 */
uint32_t morton_set_bulk4(uint64_t *out, int16_t points[][4], uint32_t n);
#endif

/** @} */

#endif
