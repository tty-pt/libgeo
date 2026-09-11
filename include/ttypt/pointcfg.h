#ifndef POINTCFG_H
#define POINTCFG_H

/**
 * @file pointcfg.h
 * @brief Configuration objects: one "static-method class" per point config.
 *
 * Each spatial config — width (bytes per lane) x dimension count — is one
 * exported const object (a struct of function pointers, "a class with
 * all-static methods"). Pick the object matching your coordinates and call
 * its members:
 *
 * @code
 * #include <ttypt/pointcfg.h>
 *
 * islet_init();
 * uint32_t db = islet_open(NULL, NULL, 0xFF);
 *
 * int16_t p[3] = { 10, 20, 30 };
 * Point3_2.put(db, p, 42);            // store 42 at the cell
 * uint32_t v = Point3_2.get(db, p);   // read it back
 *
 * int16_t s[3] = { 0, 0, 0 };
 * uint16_t l[3] = { 4, 4, 4 };
 * uint32_t cur = Point3_2.iter(db, s, l);
 * while (Point3_2.next(p, &v, cur))   // advance the box iterator
 *     ;   // ... use (p, v)
 * @endcode
 *
 * A language server offers autocomplete for every member of a config off a
 * single symbol, so the whole operations set for that config is discoverable
 * without leaving the type.
 *
 * Objects:
 *   Point1_2 Point2_2 Point3_2 Point4_2   int16 lanes (2 bytes), 1..4 dims
 *   Point2_4                              int32 lanes (4 bytes), 2 dims
 *
 * @note The flat per-dimension functions (morton_set_N, islet_put_N,
 *       point_add_N, the islet_*_2_32 family, islet_ops[]) remain for ABI
 *       compatibility, tight codec/vector loops, and genuinely runtime
 *       dims; the config objects are the recommended ergonomic surface and
 *       point at the same implementations.
 *
 * @note One config per database: all configs share the uint64 key space
 *       with different layouts, and the file carries no config tag —
 *       never reopen a database with another config's functions.
 *
 * @see islet_core
 */

#include <stdint.h>

#include <ttypt/islet.h>
#include <ttypt/rec.h>

/** @defgroup islet_pointcfg Config objects
 *  @brief Static-method interfaces, one object per point config.
 *  @{
 */

/**
 * @brief 2-byte-lane (int16) point config object type.
 *
 * One shared type for the 1..4 dimensional int16 configs; the dimension
 * count lives in the object name (Point1_2..Point4_2) and is documented
 * per member below. Box lengths are uint16_t; iterators advance via
 * islet_next().
 *
 * Member-to-flat-function mapping (dimension N fixed per object):
 *   morton_set  = morton_set_N       morton_get = morton_get_N
 *   add/sub/min/max/copy/set/idx/debug/vol = point_*_N
 *   put/get/del/del_all/cell_count = islet_put/get/del/del_all/cell_count_N
 *   replace     = islet_set_N (renamed: .set is the point broadcast)
 *   get_multi   = islet_get_multi_N    iter = islet_iter_N
 *   fill_bbox   = rec_axis_fill_bbox_N
 *   next        = islet_next
 */
typedef struct islet_point2b {
	/** Encode an N-lane point to its Morton code. m = morton_set_N(p). */
	uint64_t (*morton_set)(int16_t *p);
	/** Decode a Morton code back into an N-lane point. m = morton_get_N. */
	void (*morton_get)(int16_t *p, uint64_t code);
	/** Component-wise addition. m = point_add_N(tar, a, b). */
	void (*add)(int16_t *tar, int16_t *a, int16_t *b);
	/** Component-wise subtraction. m = point_sub_N(tar, a, b). */
	void (*sub)(int16_t *tar, int16_t *a, int16_t *b);
	/** Component-wise minimum. m = point_min_N(tar, a, b). */
	void (*min)(int16_t *tar, int16_t *a, int16_t *b);
	/** Component-wise maximum. m = point_max_N(tar, a, b). */
	void (*max)(int16_t *tar, int16_t *a, int16_t *b);
	/** Copy one point into another. m = point_copy_N(tar, src). */
	void (*copy)(int16_t *tar, int16_t *src);
	/** Volume (signed product of all lanes). m = point_vol_N(p). */
	int32_t (*vol)(int16_t *p);
	/** Broadcast one value to all N lanes. m = point_set_N(tar, v). */
	void (*set)(int16_t *tar, int16_t v);
	/** Print "label(p0, ..., pN-1)" to stderr. m = point_debug_N. */
	void (*debug)(char *label, int16_t *p);
	/** Row-major linear index of p within box [s, e). m = point_idx_N. */
	uint64_t (*idx)(int16_t *p, int16_t *s, int16_t *e);

	/** Append ref at a coordinate (multi-value cell). m = islet_put_N. */
	void (*put)(uint32_t pdb_hd, int16_t *p, uint32_t ref);
	/** First value at a coordinate, or ISLET_MISS. m = islet_get_N. */
	uint32_t (*get)(uint32_t pdb_hd, int16_t *p);
	/** Replace every value at a coordinate with ref. m = islet_set_N. */
	void (*replace)(uint32_t pdb_hd, int16_t *p, uint32_t ref);
	/** Delete the first value at a coordinate. m = islet_del_N. */
	void (*del)(uint32_t pdb_hd, int16_t *p);
	/** Delete every value, return count removed. m = islet_del_all_N. */
	uint32_t (*del_all)(uint32_t pdb_hd, int16_t *p);
	/** Count values at a coordinate. m = islet_cell_count_N. */
	uint32_t (*cell_count)(uint32_t pdb_hd, int16_t *p);
	/** Cursor over all values at a cell; islet_cell_next() drains
	 *  it. m = islet_get_multi_N. */
	uint32_t (*get_multi)(uint32_t pdb_hd, int16_t *p);
	/** Region iterator over box [s, s+l). m = islet_iter_N. */
	uint32_t (*iter)(uint32_t pdb_hd, int16_t *s, uint16_t *l);
	/** Stream the box into a sealed recall set. m = rec_axis_fill_bbox_N. */
	int (*fill_bbox)(uint32_t pdb_hd, int16_t *s, uint16_t *l,
			rec_set_t *out);
	/** Advance a box iterator. m = islet_next. */
	int (*next)(int16_t *p, uint32_t *ref, uint32_t cur);
} islet_point2b_t;

/**
 * @brief 4-byte-lane (int32) point config object type.
 *
 * The 2D x 32-bit dense config: 2 lanes of int32_t. Box lengths are also
 * int32_t; iterators advance via islet_next32(). Everything else mirrors
 * islet_point2b_t member-for-member on int32 lanes.
 *
 * Member-to-flat-function mapping:
 *   morton_set/morton_get = morton_set_2_32/morton_get_2_32
 *   add/sub/min/max/copy/set/idx/debug/vol = point_*_2_32
 *   put/get/del/del_all/cell_count = islet_put/get/del/del_all/cell_count_2_32
 *   replace = islet_set_2_32 (renamed: .set is the point broadcast)
 *   get_multi = islet_get_multi_2_32, iter = islet_iter_2_32
 *   fill_bbox = rec_axis_fill_bbox_2_32, next = islet_next32
 */
typedef struct islet_point4b {
	/** Encode the 2-lane point to its dense Morton code. */
	uint64_t (*morton_set)(int32_t *p);
	/** Decode a Morton code back into the 2-lane point. */
	void (*morton_get)(int32_t *p, uint64_t code);
	/** Component-wise addition (mod 2^32). */
	void (*add)(int32_t *tar, int32_t *a, int32_t *b);
	/** Component-wise subtraction (mod 2^32). */
	void (*sub)(int32_t *tar, int32_t *a, int32_t *b);
	/** Component-wise minimum. */
	void (*min)(int32_t *tar, int32_t *a, int32_t *b);
	/** Component-wise maximum. */
	void (*max)(int32_t *tar, int32_t *a, int32_t *b);
	/** Copy one point into another (exact 8-byte move). */
	void (*copy)(int32_t *tar, int32_t *src);
	/** Area (product of both lanes, unsigned, fits uint64_t). */
	uint64_t (*vol)(int32_t *p);
	/** Broadcast one value to both lanes. */
	void (*set)(int32_t *tar, int32_t v);
	/** Print "label(p0, p1)" to stderr. */
	void (*debug)(char *label, int32_t *p);
	/** Row-major linear index of p within box [s, e). */
	uint64_t (*idx)(int32_t *p, int32_t *s, int32_t *e);

	/** Append ref at a coordinate (multi-value cell). */
	void (*put)(uint32_t pdb_hd, int32_t *p, uint32_t ref);
	/** First value at a coordinate, or ISLET_MISS. */
	uint32_t (*get)(uint32_t pdb_hd, int32_t *p);
	/** Replace every value at a coordinate with ref. */
	void (*replace)(uint32_t pdb_hd, int32_t *p, uint32_t ref);
	/** Delete the first value at a coordinate. */
	void (*del)(uint32_t pdb_hd, int32_t *p);
	/** Delete every value, return count removed. */
	uint32_t (*del_all)(uint32_t pdb_hd, int32_t *p);
	/** Count values at a coordinate. */
	uint32_t (*cell_count)(uint32_t pdb_hd, int32_t *p);
	/** Cursor over all values at a cell; islet_cell_next() drains it. */
	uint32_t (*get_multi)(uint32_t pdb_hd, int32_t *p);
	/** Region iterator over box [s, s+l). */
	uint32_t (*iter)(uint32_t pdb_hd, int32_t *s, int32_t *l);
	/** Stream the box into a sealed recall set. */
	int (*fill_bbox)(uint32_t pdb_hd, int32_t *s, int32_t *l,
			rec_set_t *out);
	/** Advance a box iterator. */
	int (*next)(int32_t *p, uint32_t *ref, uint32_t cur);
} islet_point4b_t;

/** @brief 1D int16 config object (2-byte lanes). */
extern const islet_point2b_t Point1_2;
/** @brief 2D int16 config object (2-byte lanes). */
extern const islet_point2b_t Point2_2;
/** @brief 3D int16 config object (2-byte lanes). */
extern const islet_point2b_t Point3_2;
/** @brief 4D int16 config object (2-byte lanes). */
extern const islet_point2b_t Point4_2;
/** @brief 2D int32 config object (4-byte lanes, dense full-key codec). */
extern const islet_point4b_t Point2_4;

/** @} */

#endif