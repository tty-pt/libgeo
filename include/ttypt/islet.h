#ifndef ISLET_H
#define ISLET_H

/**
 * @file islet.h
 * @brief Public API for libislet — spatial indexing with Morton codes.
 *
 * Islet provides efficient spatial database operations using Morton codes
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

/* Optimization tunable — a 0/1 flag. ISLET_SIMD_MORTON gates the batch
 * encoders (morton_set_bulk / morton_set_bulk4); it defaults to 1.
 * Override with -DISLET_SIMD_MORTON=0/1 (note: -U does NOT work — this
 * block re-defines an undefined macro to its default). The former
 * per-path tunables (INLINE/HOIST/CLZ/UNROLL) were promoted to
 * unconditional after measurement showed the generic fallbacks never
 * win; see docs/PERF.md. */
#ifndef ISLET_SIMD_MORTON
#define ISLET_SIMD_MORTON 1
#endif

#include "morton.h"

/** @defgroup islet_core Islet core API
 *  @brief Spatial map operations using Morton code indexing.
 *
 *  Islet provides a spatial database that maps multi-dimensional coordinates
 *  to arbitrary uint32_t values. Internally, coordinates are converted to
 *  Morton codes (Z-order) for efficient spatial queries and storage.
 *
 *  Value Semantics (multi-value cells):
 *  - A grid cell may hold MULTIPLE values (QM_MULTIVALUE map). islet_put_N()
 *    appends; islet_get_N() returns the first value; islet_get_multi_N()
 *    iterates all values in insertion order; islet_del_N() removes the
 *    first value; islet_del_all_N() removes every value; islet_set_N()
 *    replaces all values with a single new one (N = 1..4).
 *
 *  Coordinate System:
 *  - Type: int16_t (signed 16-bit integers)
 *  - Range: -32768 to 32767 per dimension
 *  - Dimensions 1..4 supported; optimized for 3D (dim=3) with 4D
 *    (dim=4) support. 1D/2D/3D codes are bit-identical to v0.5.0;
 *    4D codes densely fill all 64 key bits.
 *
 *  The 2D x 32-bit config (int32_t lanes, islet_*_2_32 / morton_*_2_32)
 *  is a second dense family member: 2 x 32 = 64 key bits, no reserved
 *  bits, lane range -2147483648..2147483647.
 *
 *  @note Keyspace: all dimensions and configs share the uint64 Morton
 *        key space. Use one dimension AND one config per database, and
 *        never reopen a database with a different config's functions:
 *        the file carries no config tag, so a cross-config reopen
 *        silently decodes garbage.
 *
 *  @note islet_ops[] covers the int16-lane configs (1..4) only; the
 *        2D x 32-bit config is reached through its suffixed monomorphs
 *        (islet_put_2_32(), islet_iter_2_32(), ...) directly.
 *
 *  @note Preferred surface: the config objects in pointcfg.h
 *        (Point1_2..Point4_2, Point2_4) group everything above per
 *        config as static-method structs. This header's flat functions
 *        remain the ABI and the tight-loop fast path.
 *
 *  @note Thread Safety: Islet inherits libqmap's thread-safety properties.
 *        It uses global state and is NOT thread-safe. Use external
 *        synchronization if accessing from multiple threads.
 *
 *  @note Capacity: The mask parameter sizes the initial hash table
 *        (capacity = mask + 1). The map AUTO-GROWS by doubling when it
 *        fills (inherited from libqmap; libislet never passes QM_NOGROW),
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
 *  @see islet_morton
 *  @see islet_point
 *  @{
 */

/**
 * Sentinel value returned by islet_get when no entry exists.
 */
#define ISLET_MISS UINT32_MAX

/**
 * Maximum bounding-box volume (in cells) accepted by the rec_axis_fill_bbox
 * family (including rec_axis_fill_bbox_2_32).
 * Boxes whose volume exceeds this are rejected with -1: they are almost
 * certainly a caller bug, and the sealed set would be huge. Raw
 * islet_iter_N() carries no such cap — it only allocates for entries found.
 */
#define ISLET_FILL_MAX_VOL 1048576u

/**
 * @brief Initialize the islet subsystem.
 *
 * Registers custom types with libqmap (uint64_t for Morton codes, uint32_t
 * for values) and sets up the Morton code comparator for sorted iteration.
 * Also initializes the internal IDM (ID Manager) for iterator handles.
 *
 * @warning Must be called before any other islet functions. Calling other
 *          functions without initialization results in undefined behavior.
 *
 * @note This function can be called multiple times safely (idempotent if
 *       qmap types are already registered).
 *
 * Example:
 * @code
 * int main() {
 *     islet_init();  // Always call first
 *     uint32_t db = islet_open(NULL, NULL, 0xFF);
 *     // ... use database
 *     return 0;
 * }
 * @endcode
 *
 * @see islet_open
 */
void islet_init(void);

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
 * @return Database handle for use with other islet functions.
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
 * islet_init();
 * uint32_t db = islet_open(NULL, NULL, 0xFF);  // 256-entry initial table
 * @endcode
 *
 * Example (persistent):
 * @code
 * islet_init();
 * uint32_t db = islet_open("world.db", "main", 0xFFF);  // 4096-entry initial table
 * @endcode
 *
 * @see islet_init
 * @see qmap_open
 * @see qmap_save
 * @see qmap_close
 */
uint32_t islet_open(char *filename, char *database, uint32_t mask);

/**
 * @brief Create an iterator for all points within a rectangular bounding box.
 *
 * Queries the spatial database for all entries within the axis-aligned
 * bounding box defined by start point 's' and lengths 'l'. The iterator
 * pre-allocates and collects all matching results using Morton code range
 * queries internally.
 *
 * @param[in] pdb_hd Database handle from islet_open().
 * @param[in] s      Start point (minimum corner of bounding box).
 *                   Array of int16_t with at least the dimension count
 *                   of the function used.
 * @param[in] l      Lengths of bounding box per dimension (unsigned).
 *                   Array of uint16_t with at least the dimension count
 *                   of the function used.
 *                   The box covers s[i]..s[i]+l[i] inclusive per dimension.
 *
 * @note Per-dimension variants islet_iter_1..4; the function name carries
 *       the dimension count (1..4).
 *
 * @return Iterator handle for use with islet_next(). The handle is
 *         automatically freed when islet_next() returns 0.
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
 * uint32_t iter = islet_iter_2(db, start, lengths);
 * @endcode
 *
 * Example (3D region):
 * @code
 * int16_t start[3] = {-10, -10, -10};
 * uint16_t lengths[3] = {20, 20, 20};  // 20x20x20 cube (8000 points)
 * uint32_t iter = islet_iter_3(db, start, lengths);
 * @endcode
 *
 * @see islet_next
 * @see morton_set_1 morton_set_2 morton_set_3 morton_set_4
 */
uint32_t islet_iter_1(uint32_t pdb_hd, int16_t *s, uint16_t *l);
uint32_t islet_iter_2(uint32_t pdb_hd, int16_t *s, uint16_t *l);
uint32_t islet_iter_3(uint32_t pdb_hd, int16_t *s, uint16_t *l);
uint32_t islet_iter_4(uint32_t pdb_hd, int16_t *s, uint16_t *l);

/**
 * @brief 2D x 32-bit box iterator. Same contract as islet_iter_N(),
 *        on int32_t lanes (start AND lengths). Use with islet_next32().
 */
uint32_t islet_iter_2_32(uint32_t pdb_hd, int32_t *s, int32_t *l);

/**
 * @brief Advance iterator and retrieve the next point/value pair.
 *
 * Retrieves the next entry from the iterator created by islet_iter_N().
 * Skips empty cells in the bounding box. When all entries have been
 * returned (or the box was empty), returns 0 and automatically frees
 * the iterator's internal memory.
 *
 * @param[out] p   Output point array. Must have space for the dimension
 *                 count of the matching islet_iter_N() call. Filled with
 *                 the coordinates of the next point.
 * @param[out] ref Output pointer for the stored value (uint32_t).
 *                 Filled with the value stored at this point.
 * @param[in]  cur Iterator handle from islet_iter_N().
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
 *          Do not call islet_next() again with the same handle.
 *
 * @warning The iterator's internal memory is freed when this returns 0.
 *          There is no separate "close" or "free" function.
 *
 * Typical usage pattern:
 * @code
 * int16_t start[3] = {0, 0, 0};
 * uint16_t lengths[3] = {10, 10, 10};
 * uint32_t iter = islet_iter(db, start, lengths, 3);
 * 
 * int16_t point[3];
 * uint32_t value;
 * while (islet_next(point, &value, iter)) {
 *     printf("Point (%d,%d,%d) = %u\n",
 *            point[0], point[1], point[2], value);
 * }
 * // Iterator is automatically freed after loop
 * @endcode
 *
 * @see islet_iter
 */
int islet_next(int16_t *p, uint32_t *ref, uint32_t cur);

/**
 * @brief Advance a 2D x 32-bit iterator (islet_iter_2_32 handle).
 *
 * Same contract as islet_next(), on int32_t lanes. Iterators from
 * islet_iter_1..4 use islet_next(); 2D x 32-bit iterators use this.
 *
 * @param[out] p   Output point array (int32_t[2]).
 * @param[out] ref Output pointer for the stored value (uint32_t).
 * @param[in]  cur Iterator handle from islet_iter_2_32().
 *
 * @return 1 on entry, 0 when complete (handle freed).
 */
int islet_next32(int32_t *p, uint32_t *ref, uint32_t cur);

/**
 * @brief Delete the first value stored at a spatial coordinate.
 *
 * Per-dimension variants (islet_del_1..4); the function name carries the
 * dimension count. Removes the first value stored at the given point
 * from the database. Internally converts the coordinate to a Morton
 * code and calls qmap_del(). If no entry exists at the coordinate, this
 * is a no-op (safe to call). When several values share the cell, only
 * the earliest-inserted one is removed; use islet_del_all_3() to clear
 * the cell.
 *
 * @param[in] pdb_hd Database handle from islet_open().
 * @param[in] p      Point coordinate. Array of int16_t with at least the
 *                   dimension count of the function used.
 *
 * @note This operation invalidates any pointers obtained from islet_get_3()
 *       or islet_next() that refer to this coordinate.
 *
 * @note Safe to call on non-existent coordinates (no error, no effect).
 *
 * Example:
 * @code
 * int16_t pos[3] = {10, 20, 30};
 * islet_del_3(db, pos);  // Remove entry at (10,20,30)
 * @endcode
 *
 * @see islet_get_3
 * @see islet_put_3
 * @see qmap_del
 */
static inline void
islet_del_1(uint32_t pdb_hd, int16_t *p)
{
	uint64_t at = morton_set_1(p);
	qmap_del(pdb_hd, &at);
}

/**
 * @brief 2D delete. See islet_del_1() for the family docs.
 */
static inline void
islet_del_2(uint32_t pdb_hd, int16_t *p)
{
	uint64_t at = morton_set_2(p);
	qmap_del(pdb_hd, &at);
}

/**
 * @brief 3D delete. See islet_del_1() for the family docs.
 */
static inline void
islet_del_3(uint32_t pdb_hd, int16_t *p)
{
	uint64_t at = morton_set_3(p);
	qmap_del(pdb_hd, &at);
}

/**
 * @brief 4D delete. See islet_del_1() for the family docs.
 */
static inline void
islet_del_4(uint32_t pdb_hd, int16_t *p)
{
	uint64_t at = morton_set_4(p);
	qmap_del(pdb_hd, &at);
}

/**
 * @brief 2D x 32-bit delete. See islet_del_1() for the family docs.
 */
static inline void
islet_del_2_32(uint32_t pdb_hd, int32_t *p)
{
	uint64_t at = morton_set_2_32(p);
	qmap_del(pdb_hd, &at);
}

/**
 * @brief Delete every value stored at a spatial coordinate.
 *
 * Per-dimension variants (islet_del_all_1..4). Removes all values stored
 * at the given point from the database. Internally converts the
 * coordinate to a Morton code and calls qmap_del_all(). If no entry
 * exists at the coordinate, this is a no-op (safe to call, returns 0).
 *
 * @param[in] pdb_hd Database handle from islet_open().
 * @param[in] p      Point coordinate. Array of int16_t with at least the
 *                   dimension count of the function used.
 *
 * @return Number of values removed (0 if the cell was empty).
 *
 * Example:
 * @code
 * int16_t pos[3] = {10, 20, 30};
 * uint32_t n = islet_del_all_3(db, pos);  // Clear the cell
 * @endcode
 *
 * @see islet_del_3
 * @see islet_set_3
 * @see qmap_del_all
 */
static inline uint32_t
islet_del_all_1(uint32_t pdb_hd, int16_t *p)
{
	uint64_t at = morton_set_1(p);
	uint32_t n = qmap_count(pdb_hd, &at);

	if (n)
		qmap_del_all(pdb_hd, &at);

	return n;
}

/**
 * @brief 2D delete-all. See islet_del_all_1() for the family docs.
 */
static inline uint32_t
islet_del_all_2(uint32_t pdb_hd, int16_t *p)
{
	uint64_t at = morton_set_2(p);
	uint32_t n = qmap_count(pdb_hd, &at);

	if (n)
		qmap_del_all(pdb_hd, &at);

	return n;
}

/**
 * @brief 3D delete-all. See islet_del_all_1() for the family docs.
 */
static inline uint32_t
islet_del_all_3(uint32_t pdb_hd, int16_t *p)
{
	uint64_t at = morton_set_3(p);
	uint32_t n = qmap_count(pdb_hd, &at);

	if (n)
		qmap_del_all(pdb_hd, &at);

	return n;
}

/**
 * @brief 4D delete-all. See islet_del_all_1() for the family docs.
 */
static inline uint32_t
islet_del_all_4(uint32_t pdb_hd, int16_t *p)
{
	uint64_t at = morton_set_4(p);
	uint32_t n = qmap_count(pdb_hd, &at);

	if (n)
		qmap_del_all(pdb_hd, &at);

	return n;
}

/**
 * @brief 2D x 32-bit delete-all. See islet_del_all_1() for the
 *        family docs.
 */
static inline uint32_t
islet_del_all_2_32(uint32_t pdb_hd, int32_t *p)
{
	uint64_t at = morton_set_2_32(p);
	uint32_t n = qmap_count(pdb_hd, &at);

	if (n)
		qmap_del_all(pdb_hd, &at);

	return n;
}

/**
 * @brief Retrieve the value stored at a spatial coordinate.
 *
 * Per-dimension variants (islet_get_1..4). Looks up the value at the given
 * point in the database. Internally converts the coordinate to a Morton
 * code and calls qmap_get().
 *
 * @param[in] pdb_hd Database handle from islet_open().
 * @param[in] p      Point coordinate. Array of int16_t with at least the
 *                   dimension count of the function used. Coordinates are
 *                   signed 16-bit integers ranging from -32768 to 32767.
 *
 * @return The first stored uint32_t value at this coordinate, or ISLET_MISS
 *         (UINT32_MAX) if no entry exists at this point. When several
 *         values share the cell, the earliest-inserted one is returned;
 *         use islet_get_multi_3() to retrieve them all.
 *
 * @note ISLET_MISS equals UINT32_MAX (0xFFFFFFFF), the same as QM_MISS.
 *       This is the standard sentinel value for missing entries.
 *
 * @note The returned value is a copy, not a pointer. Unlike qmap_get()
 *       which returns a pointer, islet_get_3() returns the actual uint32_t
 *       value.
 *
 * Example:
 * @code
 * int16_t pos[3] = {10, 20, 30};
 * uint32_t value = islet_get_3(db, pos);
 * if (value == ISLET_MISS) {
 *     printf("No entry at (10,20,30)\n");
 * } else {
 *     printf("Value at (10,20,30) = %u\n", value);
 * }
 * @endcode
 *
 * @see islet_put_3
 * @see islet_del_3
 * @see ISLET_MISS
 * @see qmap_get
 */
static inline uint32_t
islet_get_1(uint32_t pdb_hd, int16_t *p)
{
	uint64_t at = morton_set_1(p);
	const void *ref = qmap_get(pdb_hd, &at);

	if (ref)
		return * (uint32_t *) ref;

	return ISLET_MISS;
}

/**
 * @brief 2D retrieve. See islet_get_1() for the family docs.
 */
static inline uint32_t
islet_get_2(uint32_t pdb_hd, int16_t *p)
{
	uint64_t at = morton_set_2(p);
	const void *ref = qmap_get(pdb_hd, &at);

	if (ref)
		return * (uint32_t *) ref;

	return ISLET_MISS;
}

/**
 * @brief 3D retrieve. See islet_get_1() for the family docs.
 */
static inline uint32_t
islet_get_3(uint32_t pdb_hd, int16_t *p)
{
	uint64_t at = morton_set_3(p);
	const void *ref = qmap_get(pdb_hd, &at);

	if (ref)
		return * (uint32_t *) ref;

	return ISLET_MISS;
}

/**
 * @brief 4D retrieve. See islet_get_1() for the family docs.
 */
static inline uint32_t
islet_get_4(uint32_t pdb_hd, int16_t *p)
{
	uint64_t at = morton_set_4(p);
	const void *ref = qmap_get(pdb_hd, &at);

	if (ref)
		return * (uint32_t *) ref;

	return ISLET_MISS;
}

/**
 * @brief 2D x 32-bit retrieve. See islet_get_1() for the family docs.
 */
static inline uint32_t
islet_get_2_32(uint32_t pdb_hd, int32_t *p)
{
	uint64_t at = morton_set_2_32(p);
	const void *ref = qmap_get(pdb_hd, &at);

	if (ref)
		return * (uint32_t *) ref;

	return ISLET_MISS;
}

/**
 * @brief Count the values stored at a spatial coordinate.
 *
 * @param[in] pdb_hd Database handle from islet_open().
 * @param[in] p      Point coordinate. Array of int16_t with at least the
 *                   dimension count of the function used.
 *
 * @return Number of values stored at this coordinate (0 if empty).
 *
 * @see islet_get_1 islet_get_2 islet_get_3 islet_get_4
 * @see islet_get_multi_1 islet_get_multi_2 islet_get_multi_3 islet_get_multi_4
 * @see qmap_count
 */
static inline uint32_t
islet_cell_count_1(uint32_t pdb_hd, int16_t *p)
{
	uint64_t at = morton_set_1(p);

	return qmap_count(pdb_hd, &at);
}

/**
 * @brief 2D cell count. See islet_cell_count_1() for the family docs.
 */
static inline uint32_t
islet_cell_count_2(uint32_t pdb_hd, int16_t *p)
{
	uint64_t at = morton_set_2(p);

	return qmap_count(pdb_hd, &at);
}

/**
 * @brief 3D cell count. See islet_cell_count_1() for the family docs.
 */
static inline uint32_t
islet_cell_count_3(uint32_t pdb_hd, int16_t *p)
{
	uint64_t at = morton_set_3(p);

	return qmap_count(pdb_hd, &at);
}

/**
 * @brief 4D cell count. See islet_cell_count_1() for the family docs.
 */
static inline uint32_t
islet_cell_count_4(uint32_t pdb_hd, int16_t *p)
{
	uint64_t at = morton_set_4(p);

	return qmap_count(pdb_hd, &at);
}

/**
 * @brief 2D x 32-bit cell count. See islet_cell_count_1() for the
 *        family docs.
 */
static inline uint32_t
islet_cell_count_2_32(uint32_t pdb_hd, int32_t *p)
{
	uint64_t at = morton_set_2_32(p);

	return qmap_count(pdb_hd, &at);
}

/**
 * @brief Diagnostic: index entries examined by the last box walk.
 *
 * Returns the number of index entries examined by the most recent
 * islet_iter_N()/rec_axis_fill_bbox_N() box walk in this process. Tests and
 * benchmarks use it to prove the Z-interval skip engages (entries
 * examined well below the morton-interval width on sparse boxes).
 *
 * @return Entry-examination count of the last box walk (0 if none ran yet).
 */
uint32_t islet_last_scan_count(void);

/**
 * @brief Create an iterator for all values stored at one coordinate.
 *
 * Opens a chain-aware cursor over every value stored at the given point,
 * in insertion order. Use islet_cell_next() to retrieve the values.
 *
 * @param[in] pdb_hd Database handle from islet_open().
 * @param[in] p      Point coordinate. Array of int16_t with at least the
 *                   dimension count of the function used.
 *
 * @note Per-dimension variants islet_get_multi_1..4.
 *
 * @return Iterator handle for use with islet_cell_next(), or QM_MISS when
 *         the cell holds no values.
 *
 * Example:
 * @code
 * int16_t pos[3] = {10, 20, 30};
 * uint32_t cur = islet_get_multi_3(db, pos);
 * if (cur != QM_MISS) {
 *     uint32_t value;
 *     while (islet_cell_next(&value, cur))
 *         printf("value %u\n", value);
 * }
 * @endcode
 *
 * @see islet_cell_next
 * @see islet_cell_count_1 islet_cell_count_2 islet_cell_count_3 islet_cell_count_4
 */
uint32_t islet_get_multi_1(uint32_t pdb_hd, int16_t *p);
uint32_t islet_get_multi_2(uint32_t pdb_hd, int16_t *p);
uint32_t islet_get_multi_3(uint32_t pdb_hd, int16_t *p);
uint32_t islet_get_multi_4(uint32_t pdb_hd, int16_t *p);

/**
 * @brief 2D x 32-bit cell iterator. Same contract as islet_get_multi_N()
 *        (use the shared islet_cell_next() to read values).
 */
uint32_t islet_get_multi_2_32(uint32_t pdb_hd, int32_t *p);

/**
 * @brief Advance a cell iterator and retrieve the next value.
 *
 * Retrieves the next value from the iterator created by
 * islet_get_multi_N().
 * When all values have been returned, returns 0 and automatically frees
 * the iterator.
 *
 * @param[out] ref Output pointer for the stored value (uint32_t).
 * @param[in]  cur Iterator handle from islet_get_multi_N().
 *
 * @return 1 if a value was retrieved. 0 when iteration is complete; the
 *         iterator is automatically freed on return 0.
 *
 * @warning After this function returns 0, the iterator handle becomes invalid.
 *
 * @see islet_get_multi_1 islet_get_multi_2 islet_get_multi_3 islet_get_multi_4
 */
int islet_cell_next(uint32_t *ref, uint32_t cur);

/**
 * @brief Store a value at a spatial coordinate (append).
 *
 * Per-dimension variants (islet_put_1..4). Appends the value at the given
 * point in the database. Internally converts the coordinate to a Morton
 * code and calls qmap_put(). If entries already exist at this coordinate,
 * the new value is ADDED alongside them (multi-value cell) - nothing is
 * replaced. Use islet_set_N() for replace semantics, islet_get_multi_N()
 * to read all values back.
 *
 * @param[in] pdb_hd Database handle from islet_open().
 * @param[in] p      Point coordinate. Array of int16_t with at least the
 *                   dimension count of the function used. Coordinates are
 *                   signed 16-bit integers ranging from -32768 to 32767.
 * @param[in] thing  Value to store (uint32_t). Can be any 32-bit value
 *                   including 0. Avoid using QM_MISS (0xFFFFFFFF) as it
 *                   may cause confusion, though it's technically valid.
 *
 * @note If the database is file-backed (filename provided to islet_open()),
 *       changes are automatically saved at process exit. Call qmap_save()
 *       explicitly for mid-execution persistence.
 *
 * Example (store single value):
 * @code
 * int16_t pos[3] = {10, 20, 30};
 * islet_put_3(db, pos, 42);  // Store value 42 at (10,20,30)
 * @endcode
 *
 * Example (update existing value):
 * @code
 * int16_t pos[3] = {10, 20, 30};
 * uint32_t old = islet_get_3(db, pos);
 * if (old != ISLET_MISS) {
 *     islet_set_3(db, pos, old + 1);  // Increment (replace)
 * }
 * @endcode
 *
 * Example (populate a 2D grid):
 * @code
 * for (int16_t x = 0; x < 10; x++) {
 *     for (int16_t y = 0; y < 10; y++) {
 *         int16_t pos[2] = {x, y};
 *         islet_put_2(db, pos, x * 10 + y);
 *     }
 * }
 * @endcode
 *
 * @see islet_get_3
 * @see islet_del_3
 * @see qmap_put
 * @see qmap_save
 */
static inline void
islet_put_1(uint32_t pdb_hd, int16_t *p, uint32_t thing)
{
	uint64_t code = morton_set_1(p);
	qmap_put(pdb_hd, &code, &thing);
}

/**
 * @brief 2D append. See islet_put_1() for the family docs.
 */
static inline void
islet_put_2(uint32_t pdb_hd, int16_t *p, uint32_t thing)
{
	uint64_t code = morton_set_2(p);
	qmap_put(pdb_hd, &code, &thing);
}

/**
 * @brief 3D append. See islet_put_1() for the family docs.
 */
static inline void
islet_put_3(uint32_t pdb_hd, int16_t *p, uint32_t thing)
{
	uint64_t code = morton_set_3(p);
	qmap_put(pdb_hd, &code, &thing);
}

/**
 * @brief 4D append. See islet_put_1() for the family docs.
 */
static inline void
islet_put_4(uint32_t pdb_hd, int16_t *p, uint32_t thing)
{
	uint64_t code = morton_set_4(p);
	qmap_put(pdb_hd, &code, &thing);
}

/**
 * @brief 2D x 32-bit append. See islet_put_1() for the family docs.
 */
static inline void
islet_put_2_32(uint32_t pdb_hd, int32_t *p, uint32_t thing)
{
	uint64_t code = morton_set_2_32(p);
	qmap_put(pdb_hd, &code, &thing);
}

/**
 * @brief Store a value at a spatial coordinate (replace).
 *
 * Per-dimension variants (islet_set_1..4). Replaces every value at the
 * given point with a single new value: islet_del_all_N() followed by
 * islet_put_N(). This preserves the historical single-value overwrite
 * convenience on top of multi-value cells.
 *
 * @param[in] pdb_hd Database handle from islet_open().
 * @param[in] p      Point coordinate. Array of int16_t with at least the
 *                   dimension count of the function used.
 * @param[in] thing  Value to store (uint32_t).
 *
 * Example:
 * @code
 * int16_t pos[3] = {10, 20, 30};
 * islet_put_3(db, pos, 42);
 * islet_set_3(db, pos, 99);  // Cell now holds exactly {99}
 * @endcode
 *
 * @see islet_put_3
 * @see islet_del_all_3
 * @see islet_get_3
 */
static inline void
islet_set_1(uint32_t pdb_hd, int16_t *p, uint32_t thing)
{
	uint64_t code = morton_set_1(p);
	qmap_del_all(pdb_hd, &code);
	qmap_put(pdb_hd, &code, &thing);
}

/**
 * @brief 2D replace. See islet_set_1() for the family docs.
 */
static inline void
islet_set_2(uint32_t pdb_hd, int16_t *p, uint32_t thing)
{
	uint64_t code = morton_set_2(p);
	qmap_del_all(pdb_hd, &code);
	qmap_put(pdb_hd, &code, &thing);
}

/**
 * @brief 3D replace. See islet_set_1() for the family docs.
 */
static inline void
islet_set_3(uint32_t pdb_hd, int16_t *p, uint32_t thing)
{
	uint64_t code = morton_set_3(p);
	qmap_del_all(pdb_hd, &code);
	qmap_put(pdb_hd, &code, &thing);
}

/**
 * @brief 4D replace. See islet_set_1() for the family docs.
 */
static inline void
islet_set_4(uint32_t pdb_hd, int16_t *p, uint32_t thing)
{
	uint64_t code = morton_set_4(p);
	qmap_del_all(pdb_hd, &code);
	qmap_put(pdb_hd, &code, &thing);
}

/**
 * @brief 2D x 32-bit replace. See islet_set_1() for the family docs.
 */
static inline void
islet_set_2_32(uint32_t pdb_hd, int32_t *p, uint32_t thing)
{
	uint64_t code = morton_set_2_32(p);
	qmap_del_all(pdb_hd, &code);
	qmap_put(pdb_hd, &code, &thing);
}

/**
 * @brief Fill a recall candidate set with the values in a bounding box.
 *
 * Per-dimension variants (rec_axis_fill_bbox_1..4). Kernel-form
 * space-axis adapter: streams every value stored in the axis-aligned
 * box [s, s+l] into a recall candidate set, then seals it (sort +
 * dedup). Multiple values sharing one cell each enter the set; sealing
 * collapses exact duplicates. The walk never materializes the box
 * volume — sparse queries over large boxes stay cheap.
 *
 * @param[in] pdb_hd Database handle from islet_open().
 * @param[in] s      Start point (minimum corner). Array of int16_t.
 * @param[in] l      Lengths per dimension. Array of uint16_t.
 * @param[out] out   Caller-owned candidate set; appended to, then sealed.
 *                   May already hold refs (result is the union, sealed).
 *
 * @return 0 on success (out sealed). -1 when out is NULL, or the box
 *         volume exceeds ISLET_FILL_MAX_VOL (note: 4D volumes are 4-way
 *         products, so keep each side small).
 *
 * Example:
 * @code
 * rec_set_t *cands = rec_set_new();
 * int16_t s[3] = {0, 0, 0};
 * uint16_t l[3] = {16, 16, 16};
 * if (rec_axis_fill_bbox_3(db, s, l, cands) == 0) {
 *     // rec_set_count(cands) distinct refs, sorted
 * }
 * rec_set_free(cands);
 * @endcode
 *
 * @see ISLET_FILL_MAX_VOL
 * @see islet_iter_1 islet_iter_2 islet_iter_3 islet_iter_4
 */
int rec_axis_fill_bbox_1(uint32_t pdb_hd, int16_t *s,
		uint16_t *l, rec_set_t *out);
int rec_axis_fill_bbox_2(uint32_t pdb_hd, int16_t *s,
		uint16_t *l, rec_set_t *out);
int rec_axis_fill_bbox_3(uint32_t pdb_hd, int16_t *s,
		uint16_t *l, rec_set_t *out);
int rec_axis_fill_bbox_4(uint32_t pdb_hd, int16_t *s,
		uint16_t *l, rec_set_t *out);

/**
 * @brief 2D x 32-bit box fill. Same contract as
 *        rec_axis_fill_bbox_N(), on int32_t lanes (start AND
 *        lengths). ISLET_FILL_MAX_VOL caps every config alike.
 */
int rec_axis_fill_bbox_2_32(uint32_t pdb_hd, int32_t *s,
		int32_t *l, rec_set_t *out);

/**
 * @brief Per-dimension operation table.
 *
 * A single global, indexed by dimension count (1..4), giving the
 * per-dimension implementation of each core op. Each function pointer
 * takes NO dimension argument — the index into islet_ops[] is the dim.
 * This is the runtime-dim mechanism: when `dim` is a variable, call
 * `islet_ops[dim].morton_set(p)` instead of switching on the value
 * yourself. The type-safe alternative is to call the per-dimension
 * function directly (morton_set_3(), etc.), which the compiler fully
 * unrolls.
 *
 * islet_ops[0] is all-NULL (dimension 0 is invalid); islet_ops[1..4] are
 * initialized by the library.
 *
 * Example (runtime dimension):
 * @code
 * islet_ops[dim].morton_set(p);   // dim in 1..4
 * islet_ops[dim].iter(db, s, l);  // box-iterate the runtime dim
 * @endcode
 *
 * @see morton_set_1 morton_set_2 morton_set_3 morton_set_4
 * @see islet_iter_1 islet_iter_2 islet_iter_3 islet_iter_4
 */
typedef struct {
	uint64_t (*morton_set)(int16_t *p);
	void     (*morton_get)(int16_t *p, uint64_t code);
	void     (*point_add)(int16_t *tar, int16_t *a, int16_t *b);
	void     (*point_copy)(int16_t *tar, int16_t *src);
	void     (*put)(uint32_t pdb_hd, int16_t *p, uint32_t thing);
	uint32_t (*get)(uint32_t pdb_hd, int16_t *p);
	void     (*set)(uint32_t pdb_hd, int16_t *p, uint32_t thing);
	void     (*del)(uint32_t pdb_hd, int16_t *p);
	uint32_t (*del_all)(uint32_t pdb_hd, int16_t *p);
	uint32_t (*cell_count)(uint32_t pdb_hd, int16_t *p);
	uint32_t (*iter)(uint32_t pdb_hd, int16_t *s, uint16_t *l);
	uint32_t (*get_multi)(uint32_t pdb_hd, int16_t *p);
	int      (*fill)(uint32_t pdb_hd, int16_t *s, uint16_t *l,
			 rec_set_t *out);
} islet_ops_t;

/**
 * @brief The per-dimension operation table (islet_ops[0] = all NULL).
 *
 * NOTE: literal 5 (MAX_DIM + 1) here — MAX_DIM lives in libislet.c and
 * this header stays standalone. Keep the two in sync.
 */
extern const islet_ops_t islet_ops[5];

#if ISLET_SIMD_MORTON
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
 * Output codes match morton_set_4() exactly.
 *
 * @param[out] out     Output Morton codes. Array of uint64_t with 'n' elements.
 * @param[in]  points  Input points. Array of int16_t[4] with 'n' entries.
 * @param[in]  n       Number of points to encode (0..UINT32_MAX).
 *
 * @return Number of points encoded (always == n).
 */
uint32_t morton_set_bulk4(uint64_t *out, int16_t points[][4], uint32_t n);

/**
 * @brief Batch-decode multiple 3D Morton codes to points using SIMD.
 *
 * Inverse of morton_set_bulk(): decodes 4 codes at a time under AVX2
 * (dense compact, scalar tail + fallback otherwise). Output points
 * match morton_get_3() exactly.
 *
 * @param[out] points  Output points. Array of int16_t[3] with 'n' entries.
 * @param[in]  codes   Input Morton codes. Array of uint64_t with 'n' elements.
 * @param[in]  n       Number of codes to decode (0..UINT32_MAX).
 *
 * @return Number of points decoded (always == n).
 */
uint32_t morton_get_bulk(int16_t points[][3], const uint64_t *codes, uint32_t n);

/**
 * @brief Batch-decode multiple 4D Morton codes to points using SIMD.
 *
 * 4D analogue of morton_get_bulk(): decodes 4 codes at a time under
 * AVX2 (dense stride-4 compact), scalar tail + fallback otherwise.
 * Output points match morton_get_4() exactly.
 *
 * @param[out] points  Output points. Array of int16_t[4] with 'n' entries.
 * @param[in]  codes   Input Morton codes. Array of uint64_t with 'n' elements.
 * @param[in]  n       Number of codes to decode (0..UINT32_MAX).
 *
 * @return Number of points decoded (always == n).
 */
uint32_t morton_get_bulk4(int16_t points[][4], const uint64_t *codes, uint32_t n);
#endif

/** @} */

#endif
