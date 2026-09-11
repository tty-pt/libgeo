#ifndef POINT_H
#define POINT_H

/**
 * @file point.h
 * @brief Point arithmetic and utility functions for spatial operations.
 *
 * Provides helper functions for manipulating multi-dimensional points
 * represented as int16_t arrays. These are building blocks for spatial
 * operations in libislet.
 *
 * All functions are per-dimension specializations (point_copy_1..4, etc.).
 * There is no runtime dimension argument: pick the function matching the
 * dimension count, so the compiler sees fully unrolled bodies.
 *
 * The 2D x 32-bit config (int32_t lanes, morton_set_2_32 / islet_*_2_32)
 * has its own same-shaped family (point_add_2_32, etc.).
 *
 * Prefer the config objects (pointcfg.h: Point1_2..Point4_2, Point2_4)
 * in application code; these inlines remain the tight-loop fast path.
 *
 * @note No bounds checking is performed for efficiency. Caller must ensure
 *       arrays have sufficient space for the dimension count used.
 *
 * @see islet_core
 * @see islet_morton
 */

#include <stdint.h>
#include <stdio.h>

/** @defgroup islet_point Point utilities
 *  @brief Arithmetic and helper functions for int16_t coordinate points.
 *
 *  These functions operate on points represented as arrays of int16_t
 *  coordinates. They provide common vector operations (add, subtract,
 *  min, max), utility functions (copy, set, volume), and spatial
 *  indexing helpers.
 *
 *  Memory Safety:
 *  - All arrays must have at least as many elements as the dimension
 *    count implied by the function name (1..4)
 *  - No bounds checking is performed for performance
 *  - Caller is responsible for ensuring valid memory access
 *  - Output arrays may alias input arrays unless noted otherwise
 *
 *  Overflow Behavior:
 *  - All arithmetic follows standard int16_t overflow semantics (wrapping)
 *  - Large multiplications (e.g., point_vol) may overflow
 *  - No overflow detection or saturation is provided
 *
 *  @see islet_core
 *  @{
 */

/**
 * @brief Add two points component-wise.
 *
 * Performs vector addition: tar[i] = orig[i] + tr[i] for each dimension.
 * Overflow follows standard int16_t arithmetic (wrapping at -32768/32767).
 *
 * @param[out] tar  Output point. May alias orig or tr for in-place operations.
 * @param[in]  orig First operand (addend). Array of int16_t.
 * @param[in]  tr   Second operand (addend). Array of int16_t.
 *
 * @note No overflow detection. 32767 + 1 = -32768 (standard int16_t wrapping).
 *
 * Example:
 * @code
 * int16_t a[3] = {10, 20, 30};
 * int16_t b[3] = {1, 2, 3};
 * int16_t result[3];
 * point_add_3(result, a, b);  // result = {11, 22, 33}
 * @endcode
 *
 * @see point_sub_3
 */
static inline void
point_add_1(int16_t *tar, int16_t *orig, int16_t *tr)
{
	tar[0] = (int16_t)(orig[0] + tr[0]);
}

/**
 * @brief 2D vector addition. See point_add_1() for the family docs.
 */
static inline void
point_add_2(int16_t *tar, int16_t *orig, int16_t *tr)
{
	tar[0] = (int16_t)(orig[0] + tr[0]);
	tar[1] = (int16_t)(orig[1] + tr[1]);
}

/**
 * @brief 3D vector addition. See point_add_1() for the family docs.
 */
static inline void
point_add_3(int16_t *tar, int16_t *orig, int16_t *tr)
{
	tar[0] = (int16_t)(orig[0] + tr[0]);
	tar[1] = (int16_t)(orig[1] + tr[1]);
	tar[2] = (int16_t)(orig[2] + tr[2]);
}

/**
 * @brief 4D vector addition. See point_add_1() for the family docs.
 */
static inline void
point_add_4(int16_t *tar, int16_t *orig, int16_t *tr)
{
	tar[0] = (int16_t)(orig[0] + tr[0]);
	tar[1] = (int16_t)(orig[1] + tr[1]);
	tar[2] = (int16_t)(orig[2] + tr[2]);
	tar[3] = (int16_t)(orig[3] + tr[3]);
}

/**
 * @brief Subtract two points component-wise.
 *
 * Performs vector subtraction: tar[i] = orig[i] - tr[i] for each dimension.
 * Overflow follows standard int16_t arithmetic (wrapping at -32768/32767).
 *
 * @param[out] tar  Output point. May alias orig or tr for in-place operations.
 * @param[in]  orig Minuend (value to subtract from). Array of int16_t.
 * @param[in]  tr   Subtrahend (value to subtract). Array of int16_t.
 *
 * @note No overflow detection. -32768 - 1 = 32767 (standard int16_t wrapping).
 *
 * Example (compute delta):
 * @code
 * int16_t end[3] = {100, 200, 300};
 * int16_t start[3] = {50, 60, 70};
 * int16_t delta[3];
 * point_sub_3(delta, end, start);  // delta = {50, 140, 230}
 * @endcode
 *
 * @see point_add_3
 */
static inline void
point_sub_1(int16_t *tar, int16_t *orig, int16_t *tr)
{
	tar[0] = (int16_t)(orig[0] - tr[0]);
}

/**
 * @brief 2D vector subtraction. See point_sub_1() for the family docs.
 */
static inline void
point_sub_2(int16_t *tar, int16_t *orig, int16_t *tr)
{
	tar[0] = (int16_t)(orig[0] - tr[0]);
	tar[1] = (int16_t)(orig[1] - tr[1]);
}

/**
 * @brief 3D vector subtraction. See point_sub_1() for the family docs.
 */
static inline void
point_sub_3(int16_t *tar, int16_t *orig, int16_t *tr)
{
	tar[0] = (int16_t)(orig[0] - tr[0]);
	tar[1] = (int16_t)(orig[1] - tr[1]);
	tar[2] = (int16_t)(orig[2] - tr[2]);
}

/**
 * @brief 4D vector subtraction. See point_sub_1() for the family docs.
 */
static inline void
point_sub_4(int16_t *tar, int16_t *orig, int16_t *tr)
{
	tar[0] = (int16_t)(orig[0] - tr[0]);
	tar[1] = (int16_t)(orig[1] - tr[1]);
	tar[2] = (int16_t)(orig[2] - tr[2]);
	tar[3] = (int16_t)(orig[3] - tr[3]);
}

/**
 * @brief Compute component-wise minimum of two points.
 *
 * For each dimension i, sets tar[i] to the smaller of a[i] and b[i].
 * Useful for computing bounding box minimum corners or clamping operations.
 *
 * @param[out] tar Output point (minimum per component). May alias a or b.
 * @param[in]  a   First operand. Array of int16_t.
 * @param[in]  b   Second operand. Array of int16_t.
 *
 * @note Comparison uses signed int16_t ordering (-32768 is minimum value).
 *
 * @see point_max_3
 */
static inline void
point_min_1(int16_t *tar, int16_t *a, int16_t *b)
{
	tar[0] = a[0] < b[0] ? a[0] : b[0];
}

/**
 * @brief 2D component-wise minimum. See point_min_1() for the family docs.
 */
static inline void
point_min_2(int16_t *tar, int16_t *a, int16_t *b)
{
	tar[0] = a[0] < b[0] ? a[0] : b[0];
	tar[1] = a[1] < b[1] ? a[1] : b[1];
}

/**
 * @brief 3D component-wise minimum. See point_min_1() for the family docs.
 */
static inline void
point_min_3(int16_t *tar, int16_t *a, int16_t *b)
{
	tar[0] = a[0] < b[0] ? a[0] : b[0];
	tar[1] = a[1] < b[1] ? a[1] : b[1];
	tar[2] = a[2] < b[2] ? a[2] : b[2];
}

/**
 * @brief 4D component-wise minimum. See point_min_1() for the family docs.
 */
static inline void
point_min_4(int16_t *tar, int16_t *a, int16_t *b)
{
	tar[0] = a[0] < b[0] ? a[0] : b[0];
	tar[1] = a[1] < b[1] ? a[1] : b[1];
	tar[2] = a[2] < b[2] ? a[2] : b[2];
	tar[3] = a[3] < b[3] ? a[3] : b[3];
}

/**
 * @brief Compute component-wise maximum of two points.
 *
 * For each dimension i, sets tar[i] to the larger of a[i] and b[i].
 * Useful for computing bounding box maximum corners or clamping operations.
 *
 * @param[out] tar Output point (maximum per component). May alias a or b.
 * @param[in]  a   First operand. Array of int16_t.
 * @param[in]  b   Second operand. Array of int16_t.
 *
 * @note Comparison uses signed int16_t ordering (32767 is maximum value).
 *
 * @see point_min_3
 */
static inline void
point_max_1(int16_t *tar, int16_t *a, int16_t *b)
{
	tar[0] = a[0] > b[0] ? a[0] : b[0];
}

/**
 * @brief 2D component-wise maximum. See point_max_1() for the family docs.
 */
static inline void
point_max_2(int16_t *tar, int16_t *a, int16_t *b)
{
	tar[0] = a[0] > b[0] ? a[0] : b[0];
	tar[1] = a[1] > b[1] ? a[1] : b[1];
}

/**
 * @brief 3D component-wise maximum. See point_max_1() for the family docs.
 */
static inline void
point_max_3(int16_t *tar, int16_t *a, int16_t *b)
{
	tar[0] = a[0] > b[0] ? a[0] : b[0];
	tar[1] = a[1] > b[1] ? a[1] : b[1];
	tar[2] = a[2] > b[2] ? a[2] : b[2];
}

/**
 * @brief 4D component-wise maximum. See point_max_1() for the family docs.
 */
static inline void
point_max_4(int16_t *tar, int16_t *a, int16_t *b)
{
	tar[0] = a[0] > b[0] ? a[0] : b[0];
	tar[1] = a[1] > b[1] ? a[1] : b[1];
	tar[2] = a[2] > b[2] ? a[2] : b[2];
	tar[3] = a[3] > b[3] ? a[3] : b[3];
}

/**
 * @brief Copy a point (component-wise shallow copy).
 *
 * Copies the coordinates of one point into another. The destination must
 * have space for the dimension count of the function used. 1D/2D/4D use
 * exact scalar or multi-byte moves; nothing is ever read or written past
 * the point's own elements.
 *
 * @param[out] tar  Destination point. Must have space for the dim count.
 * @param[in]  orig Source point. Array of int16_t.
 *
 * @warning Do not use a wider variant than the arrays actually hold: e.g.
 *          copying a bare int16_t[3] with point_copy_4() would read/write
 *          past the end.
 *
 * @see point_set_3
 */
static inline void
point_copy_1(int16_t *tar, int16_t *orig)
{
	*tar = *orig;
}

/**
 * @brief 2D point copy (exact 4-byte move). See point_copy_1() for docs.
 */
static inline void
point_copy_2(int16_t *tar, int16_t *orig)
{
	*(int32_t *)tar = *(int32_t *)orig;
}

/**
 * @brief 3D point copy. See point_copy_1() for the family docs.
 */
static inline void
point_copy_3(int16_t *tar, int16_t *orig)
{
	tar[0] = orig[0];
	tar[1] = orig[1];
	tar[2] = orig[2];
}

/**
 * @brief 4D point copy (exact 8-byte move). See point_copy_1() for docs.
 */
static inline void
point_copy_4(int16_t *tar, int16_t *orig)
{
	/* Exact 8 bytes (4 x int16): no overrun possible. */
	*(int64_t *)tar = *(int64_t *)orig;
}

/**
 * @brief Compute the volume (product of all components) of a point.
 *
 * Multiplies all components together to get the volume of a box with
 * dimensions specified by the point. Useful for computing bounding box
 * volumes or array sizes for spatial grids.
 *
 * @param[in] p Point with dimensions. Array of int16_t.
 *
 * @return Product of all components: p[0] * p[1] * ... . For 3D a point
 *         {200, 200, 200} yields 8,000,000 (fits int32_t), but a 4D
 *         product of large sides may overflow int32_t.
 *
 * @note Negative components will produce negative or unexpected results.
 *       This function is designed for positive dimensions (lengths).
 *
 * Example (3D volume):
 * @code
 * int16_t size[3] = {10, 20, 30};
 * int32_t volume = point_vol_3(size);  // volume = 6000
 * @endcode
 */
static inline int32_t
point_vol_1(int16_t *p)
{
	return (int32_t)p[0];
}

/**
 * @brief 2D area. See point_vol_1() for the family docs.
 */
static inline int32_t
point_vol_2(int16_t *p)
{
	return (int32_t)p[0] * (int32_t)p[1];
}

/**
 * @brief 3D volume. See point_vol_1() for the family docs.
 */
static inline int32_t
point_vol_3(int16_t *p)
{
	return (int32_t)p[0] * (int32_t)p[1] * (int32_t)p[2];
}

/**
 * @brief 4D volume. See point_vol_1() for the family docs.
 */
static inline int32_t
point_vol_4(int16_t *p)
{
	return (int32_t)p[0] * (int32_t)p[1]
		* (int32_t)p[2] * (int32_t)p[3];
}

/**
 * @brief Set all components of a point to the same value.
 *
 * Broadcasts a single value to all dimensions: tar[i] = value for all i.
 * Useful for initializing points to zero, setting uniform bounds, or
 * creating uniform scaling factors.
 *
 * @param[out] tar   Output point. Must have space for the dim count.
 * @param[in]  value Value to assign to all components.
 *
 * Example (zero initialization):
 * @code
 * int16_t pos[3];
 * point_set_3(pos, 0);  // pos = {0, 0, 0}
 * @endcode
 *
 * Example (uniform bounds):
 * @code
 * int16_t min_bounds[3];
 * point_set_3(min_bounds, -100);  // {-100, -100, -100}
 * @endcode
 *
 * @see point_copy_3
 */
static inline void
point_set_1(int16_t *tar, int16_t value)
{
	tar[0] = value;
}

/**
 * @brief 2D broadcast set. See point_set_1() for the family docs.
 */
static inline void
point_set_2(int16_t *tar, int16_t value)
{
	tar[0] = value;
	tar[1] = value;
}

/**
 * @brief 3D broadcast set. See point_set_1() for the family docs.
 */
static inline void
point_set_3(int16_t *tar, int16_t value)
{
	tar[0] = value;
	tar[1] = value;
	tar[2] = value;
}

/**
 * @brief 4D broadcast set. See point_set_1() for the family docs.
 */
static inline void
point_set_4(int16_t *tar, int16_t value)
{
	tar[0] = value;
	tar[1] = value;
	tar[2] = value;
	tar[3] = value;
}

/**
 * @brief Print a point to stderr for debugging.
 *
 * Outputs the point in format: "label(p[0], p[1], ..., p[dim-1])\n"
 * to stderr. Useful for debugging spatial algorithms and visualizing
 * coordinate values during development.
 *
 * @param[in] label String label to print before the point (e.g., "pos").
 * @param[in] p     Point to print. Array of int16_t.
 *
 * @note Output goes to stderr, not stdout.
 *
 * Example output:
 * @code
 * int16_t pos[3] = {10, -20, 30};
 * point_debug_3("position", pos);
 * // Output to stderr: "position(10, -20, 30)\n"
 * @endcode
 */
static inline void
point_debug_1(char *label, int16_t *p)
{
	fprintf(stderr, "%s(%d)\n", label, p[0]);
}

/**
 * @brief Print a 2D point. See point_debug_1() for the family docs.
 */
static inline void
point_debug_2(char *label, int16_t *p)
{
	fprintf(stderr, "%s(%d, %d)\n", label, p[0], p[1]);
}

/**
 * @brief Print a 3D point. See point_debug_1() for the family docs.
 */
static inline void
point_debug_3(char *label, int16_t *p)
{
	fprintf(stderr, "%s(%d, %d, %d)\n", label, p[0], p[1], p[2]);
}

/**
 * @brief Print a 4D point. See point_debug_1() for the family docs.
 */
static inline void
point_debug_4(char *label, int16_t *p)
{
	fprintf(stderr, "%s(%d, %d, %d, %d)\n",
			label, p[0], p[1], p[2], p[3]);
}

/**
 * @brief Compute a linear index for a point within a bounding box.
 *
 * Converts a multi-dimensional coordinate to a linear array index using
 * row-major ordering. This is useful for mapping spatial coordinates to
 * flat array indices when storing spatial data in contiguous memory.
 *
 * Row-major means the first dimension (p[0]) varies fastest, last
 * dimension (p[dim-1]) varies slowest. This matches C array layout.
 *
 * @param[in] p   Point to index. Must be within the box [s, e).
 *                Array of int16_t.
 * @param[in] s   Box start (minimum corner, inclusive). Array of int16_t.
 * @param[in] e   Box end (maximum corner, exclusive). Array of int16_t.
 *
 * @return Linear index in range [0, volume-1] where volume is the product
 *         of (e[i] - s[i]) for all dimensions. Returns uint64_t to handle
 *         large volumes, but overflow is still possible.
 *
 * @warning Undefined behavior if p is outside the box [s, e). The function
 *          does not check bounds. Results may be incorrect or overflow.
 *
 * @note The point p[i] must satisfy: s[i] <= p[i] < e[i] for all dimensions.
 *
 * Example (2D grid):
 * @code
 * int16_t start[2] = {0, 0};
 * int16_t end[2] = {10, 10};  // 10x10 grid
 * int16_t point[2] = {3, 5};
 * uint64_t idx = point_idx_2(point, start, end);
 * // idx = 5 * 10 + 3 = 53
 * @endcode
 *
 * @see point_vol_3
 * @see islet_iter
 */
static inline uint64_t
point_idx_1(int16_t *p, int16_t *s, int16_t *e)
{
	(void)e;
	return (uint64_t)(p[0] - s[0]);
}

/**
 * @brief 2D row-major index. See point_idx_1() for the family docs.
 */
static inline uint64_t
point_idx_2(int16_t *p, int16_t *s, int16_t *e)
{
	return (uint64_t)(p[0] - s[0])
		+ (uint64_t)(p[1] - s[1]) * (uint64_t)(e[0] - s[0]);
}

/**
 * @brief 3D row-major index. See point_idx_1() for the family docs.
 */
static inline uint64_t
point_idx_3(int16_t *p, int16_t *s, int16_t *e)
{
	uint64_t w0 = (uint64_t)(e[0] - s[0]);
	uint64_t w1 = (uint64_t)(e[1] - s[1]);

	return (uint64_t)(p[0] - s[0])
		+ (uint64_t)(p[1] - s[1]) * w0
		+ (uint64_t)(p[2] - s[2]) * (w0 * w1);
}

/**
 * @brief 4D row-major index. See point_idx_1() for the family docs.
 */
static inline uint64_t
point_idx_4(int16_t *p, int16_t *s, int16_t *e)
{
	uint64_t w0 = (uint64_t)(e[0] - s[0]);
	uint64_t w1 = (uint64_t)(e[1] - s[1]);
	uint64_t w2 = (uint64_t)(e[2] - s[2]);

	return (uint64_t)(p[0] - s[0])
		+ (uint64_t)(p[1] - s[1]) * w0
		+ (uint64_t)(p[2] - s[2]) * (w0 * w1)
		+ (uint64_t)(p[3] - s[3]) * (w0 * w1 * w2);
}

/* 2D x 32-bit-lane family (morton_set_2_32 / islet_*_2_32 config).
 * Same operations as the int16 families above, on int32_t lanes
 * (-2147483648..2147483647). Arithmetic wraps on overflow; nothing
 * is ever read or written past the point's own 2 elements. */

/**
 * @brief 2D x 32-bit vector addition. See point_add_1() for the
 *        family docs.
 */
static inline void
point_add_2_32(int32_t *tar, int32_t *orig, int32_t *tr)
{
	tar[0] = orig[0] + tr[0];
	tar[1] = orig[1] + tr[1];
}

/**
 * @brief 2D x 32-bit vector subtraction. See point_sub_1() for the
 *        family docs.
 */
static inline void
point_sub_2_32(int32_t *tar, int32_t *orig, int32_t *tr)
{
	tar[0] = orig[0] - tr[0];
	tar[1] = orig[1] - tr[1];
}

/**
 * @brief 2D x 32-bit component-wise minimum. See point_min_1() for
 *        the family docs.
 */
static inline void
point_min_2_32(int32_t *tar, int32_t *a, int32_t *b)
{
	tar[0] = a[0] < b[0] ? a[0] : b[0];
	tar[1] = a[1] < b[1] ? a[1] : b[1];
}

/**
 * @brief 2D x 32-bit component-wise maximum. See point_max_1() for
 *        the family docs.
 */
static inline void
point_max_2_32(int32_t *tar, int32_t *a, int32_t *b)
{
	tar[0] = a[0] > b[0] ? a[0] : b[0];
	tar[1] = a[1] > b[1] ? a[1] : b[1];
}

/**
 * @brief 2D x 32-bit point copy (exact 8-byte move). See
 *        point_copy_1() for the family docs.
 */
static inline void
point_copy_2_32(int32_t *tar, int32_t *orig)
{
	/* Exact 8 bytes (2 x int32): no overrun possible. */
	*(int64_t *)tar = *(int64_t *)orig;
}

/**
 * @brief 2D x 32-bit area. The product of two lanes (each at most
 *        2^32-1 wide) always fits in uint64_t.
 */
static inline uint64_t
point_vol_2_32(int32_t *p)
{
	return (uint64_t)(uint32_t)p[0] * (uint64_t)(uint32_t)p[1];
}

/**
 * @brief 2D x 32-bit broadcast set. See point_set_1() for the
 *        family docs.
 */
static inline void
point_set_2_32(int32_t *tar, int32_t value)
{
	tar[0] = value;
	tar[1] = value;
}

/**
 * @brief Print a 2D x 32-bit point. See point_debug_1() for the
 *        family docs.
 */
static inline void
point_debug_2_32(char *label, int32_t *p)
{
	fprintf(stderr, "%s(%d, %d)\n", label, p[0], p[1]);
}

/**
 * @brief 2D x 32-bit row-major index. See point_idx_1() for the
 *        family docs.
 */
static inline uint64_t
point_idx_2_32(int32_t *p, int32_t *s, int32_t *e)
{
	return (uint64_t)(p[0] - s[0])
		+ (uint64_t)(p[1] - s[1]) * (uint64_t)(e[0] - s[0]);
}

/** @} */

#endif