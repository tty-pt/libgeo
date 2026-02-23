#ifndef POINT_H
#define POINT_H

/**
 * @file point.h
 * @brief Point arithmetic and utility functions for spatial operations.
 *
 * Provides helper functions for manipulating multi-dimensional points
 * represented as int16_t arrays. These are building blocks for spatial
 * operations in libgeo and support arbitrary dimensions.
 *
 * All functions work on coordinate arrays with a runtime dimension parameter,
 * making them flexible for 2D, 3D, or higher-dimensional spaces.
 *
 * @note No bounds checking is performed for efficiency. Caller must ensure
 *       arrays have sufficient space for the specified dimensions.
 *
 * @see geo_core
 * @see geo_morton
 */

#include <stdint.h>
#include <stdio.h>

/** @defgroup geo_point Point utilities
 *  @brief Arithmetic and helper functions for int16_t coordinate points.
 *
 *  These functions operate on points represented as arrays of int16_t
 *  coordinates. They provide common vector operations (add, subtract, min, max),
 *  utility functions (copy, set, volume), and spatial indexing helpers.
 *
 *  Memory Safety:
 *  - All arrays must have at least 'dim' elements allocated
 *  - No bounds checking is performed for performance
 *  - Caller is responsible for ensuring valid memory access
 *  - Output arrays may alias input arrays unless noted otherwise
 *
 *  Overflow Behavior:
 *  - All arithmetic follows standard int16_t overflow semantics (wrapping)
 *  - Large multiplications (e.g., point_vol) may overflow
 *  - No overflow detection or saturation is provided
 *
 *  @see geo_core
 *  @{
 */

/**
 * @brief Add two points component-wise.
 *
 * Performs vector addition: tar[i] = orig[i] + tr[i] for each dimension.
 * Overflow follows standard int16_t arithmetic (wrapping at -32768/32767).
 *
 * @param[out] tar  Output point. May alias orig or tr for in-place operations.
 *                  Must have space for at least 'dim' elements.
 * @param[in]  orig First operand (addend). Array of int16_t with 'dim' elements.
 * @param[in]  tr   Second operand (addend). Array of int16_t with 'dim' elements.
 * @param[in]  dim  Number of dimensions to process.
 *
 * @note All arrays must have at least 'dim' elements allocated.
 * @note No overflow detection. 32767 + 1 = -32768 (standard int16_t wrapping).
 *
 * @warning Caller must ensure all pointers are valid and arrays are sized correctly.
 *
 * Example:
 * @code
 * int16_t a[3] = {10, 20, 30};
 * int16_t b[3] = {1, 2, 3};
 * int16_t result[3];
 * point_add(result, a, b, 3);  // result = {11, 22, 33}
 * @endcode
 *
 * @see point_sub
 */
static inline void
point_add(int16_t *tar, int16_t *orig,
		int16_t *tr, uint8_t dim)
{
	for (uint8_t i = 0; i < dim; i++)
		tar[i] = orig[i] + tr[i];
}

/**
 * @brief Subtract two points component-wise.
 *
 * Performs vector subtraction: tar[i] = orig[i] - tr[i] for each dimension.
 * Overflow follows standard int16_t arithmetic (wrapping at -32768/32767).
 *
 * @param[out] tar  Output point. May alias orig or tr for in-place operations.
 *                  Must have space for at least 'dim' elements.
 * @param[in]  orig Minuend (value to subtract from). Array of int16_t.
 * @param[in]  tr   Subtrahend (value to subtract). Array of int16_t.
 * @param[in]  dim  Number of dimensions to process.
 *
 * @note All arrays must have at least 'dim' elements allocated.
 * @note No overflow detection. -32768 - 1 = 32767 (standard int16_t wrapping).
 *
 * @warning Caller must ensure all pointers are valid and arrays are sized correctly.
 *
 * Example (compute delta):
 * @code
 * int16_t end[3] = {100, 200, 300};
 * int16_t start[3] = {50, 60, 70};
 * int16_t delta[3];
 * point_sub(delta, end, start, 3);  // delta = {50, 140, 230}
 * @endcode
 *
 * @see point_add
 */
static inline void
point_sub(int16_t *tar, int16_t *orig,
		int16_t *tr, uint8_t dim)
{
	for (uint8_t i = 0; i < dim; i++)
		tar[i] = orig[i] - tr[i];
}

/**
 * @brief Compute component-wise minimum of two points.
 *
 * For each dimension i, sets tar[i] to the smaller of a[i] and b[i].
 * Useful for computing bounding box minimum corners or clamping operations.
 *
 * @param[out] tar Output point (minimum per component). May alias a or b.
 *                 Must have space for at least 'dim' elements.
 * @param[in]  a   First operand. Array of int16_t with 'dim' elements.
 * @param[in]  b   Second operand. Array of int16_t with 'dim' elements.
 * @param[in]  dim Number of dimensions to process.
 *
 * @note Comparison uses signed int16_t ordering (-32768 is minimum value).
 *
 * Example (bounding box minimum):
 * @code
 * int16_t p1[3] = {10, 50, 30};
 * int16_t p2[3] = {20, 40, 35};
 * int16_t min[3];
 * point_min(min, p1, p2, 3);  // min = {10, 40, 30}
 * @endcode
 *
 * @see point_max
 */
static inline void
point_min(int16_t *tar, int16_t *a,
		int16_t *b, uint8_t dim)
{
	for (uint8_t i = 0; i < dim; i++)
		tar[i] = a[i] < b[i] ? a[i] : b[i];
}

/**
 * @brief Compute component-wise maximum of two points.
 *
 * For each dimension i, sets tar[i] to the larger of a[i] and b[i].
 * Useful for computing bounding box maximum corners or clamping operations.
 *
 * @param[out] tar Output point (maximum per component). May alias a or b.
 *                 Must have space for at least 'dim' elements.
 * @param[in]  a   First operand. Array of int16_t with 'dim' elements.
 * @param[in]  b   Second operand. Array of int16_t with 'dim' elements.
 * @param[in]  dim Number of dimensions to process.
 *
 * @note Comparison uses signed int16_t ordering (32767 is maximum value).
 *
 * Example (bounding box maximum):
 * @code
 * int16_t p1[3] = {10, 50, 30};
 * int16_t p2[3] = {20, 40, 35};
 * int16_t max[3];
 * point_max(max, p1, p2, 3);  // max = {20, 50, 35}
 * @endcode
 *
 * @see point_min
 */
static inline void
point_max(int16_t *tar, int16_t *a,
		int16_t *b, uint8_t dim)
{
	for (uint8_t i = 0; i < dim; i++)
		tar[i] = a[i] > b[i] ? a[i] : b[i];
}

/**
 * @brief Copy a point from one array to another.
 *
 * Copies 'dim' elements from orig to tar. Functionally equivalent to
 * memcpy(tar, orig, dim * sizeof(int16_t)) but implemented as a loop.
 *
 * @param[out] tar  Destination array. Must have space for 'dim' elements.
 * @param[in]  orig Source array. Must have at least 'dim' elements.
 * @param[in]  dim  Number of dimensions to copy.
 *
 * @note Arrays must not overlap. For overlapping regions, use memmove().
 *
 * @note This is NOT safe for aliasing (tar == orig is wasteful but harmless).
 *
 * Example:
 * @code
 * int16_t original[3] = {10, 20, 30};
 * int16_t copy[3];
 * point_copy(copy, original, 3);  // copy = {10, 20, 30}
 * @endcode
 *
 * @see point_set
 */
static inline void
point_copy(int16_t *tar, int16_t *orig, uint8_t dim)
{
	for (uint8_t i = 0; i < dim; i++)
		tar[i] = orig[i];
}

/**
 * @brief Compute the volume (product of all components).
 *
 * Multiplies all components together to get the volume of a box with
 * dimensions specified by the point. Useful for computing bounding box
 * volumes or array sizes for spatial grids.
 *
 * @param[in] p   Point with dimensions. Array of int16_t with 'dim' elements.
 * @param[in] dim Number of dimensions.
 *
 * @return Product of all components: p[0] * p[1] * ... * p[dim-1].
 *         Result is int32_t to accommodate larger products, but overflow
 *         is still possible for large values.
 *
 * @warning Overflow can occur even with int32_t result. For example,
 *          a 3D box of {1000, 1000, 1000} = 1,000,000,000 (fits in int32_t),
 *          but {2000, 2000, 2000} = 8,000,000,000 (overflows int32_t).
 *
 * @note Negative components will produce negative or unexpected results.
 *       This function is designed for positive dimensions (lengths).
 *
 * Example (2D area):
 * @code
 * int16_t size[2] = {100, 50};
 * int32_t area = point_vol(size, 2);  // area = 5000
 * @endcode
 *
 * Example (3D volume):
 * @code
 * int16_t size[3] = {10, 20, 30};
 * int32_t volume = point_vol(size, 3);  // volume = 6000
 * @endcode
 */
static inline int32_t
point_vol(int16_t *p, uint8_t dim)
{
	int32_t res = 1;

	for (uint8_t i = 0; i < dim; i++)
		res *= p[i];

	return res;
}

/**
 * @brief Set all components of a point to the same value.
 *
 * Broadcasts a single value to all dimensions: tar[i] = value for all i.
 * Useful for initializing points to zero, setting uniform bounds, or
 * creating uniform scaling factors.
 *
 * @param[out] tar   Output point. Must have space for 'dim' elements.
 * @param[in]  value Value to assign to all components.
 * @param[in]  dim   Number of dimensions.
 *
 * Example (zero initialization):
 * @code
 * int16_t pos[3];
 * point_set(pos, 0, 3);  // pos = {0, 0, 0}
 * @endcode
 *
 * Example (uniform bounds):
 * @code
 * int16_t min_bounds[3];
 * int16_t max_bounds[3];
 * point_set(min_bounds, -100, 3);  // {-100, -100, -100}
 * point_set(max_bounds, 100, 3);   // {100, 100, 100}
 * @endcode
 *
 * @see point_copy
 */
static inline void
point_set(int16_t *tar, int16_t value, uint8_t dim)
{
	for (uint8_t i = 0; i < dim; i++)
		tar[i] = value;
}

/**
 * @brief Print a point to stderr for debugging.
 *
 * Outputs the point in format: "label(p[0], p[1], ..., p[dim-1])\n"
 * to stderr. Useful for debugging spatial algorithms and visualizing
 * coordinate values during development.
 *
 * @param[in] label String label to print before the point (e.g., "pos", "min").
 * @param[in] p     Point to print. Array of int16_t with 'dim' elements.
 * @param[in] dim   Number of dimensions to print.
 *
 * @note Output goes to stderr, not stdout. This avoids interfering with
 *       normal program output.
 *
 * Example output:
 * @code
 * int16_t pos[3] = {10, -20, 30};
 * point_debug("position", pos, 3);
 * // Output to stderr: "position(10, -20, 30)\n"
 * @endcode
 */
static inline void
point_debug(char *label, int16_t *p, uint8_t dim)
{
	fprintf(stderr, "%s(", label);
	for (uint8_t i = 0; i < dim; i++)
		fprintf(stderr, "%s%d", i ? ", " : "", p[i]);
	fprintf(stderr, ")\n");
}

/**
 * @brief Compute a linear index for a point within a bounding box.
 *
 * Converts a multi-dimensional coordinate to a linear array index using
 * row-major ordering. This is useful for mapping spatial coordinates to
 * flat array indices when storing spatial data in contiguous memory.
 *
 * The algorithm computes:
 * - For 2D: idx = (p[1] - s[1]) * (e[0] - s[0]) + (p[0] - s[0])
 * - For 3D: idx = ((p[2] - s[2]) * (e[1] - s[1]) + (p[1] - s[1]))
 *                  * (e[0] - s[0]) + (p[0] - s[0])
 * - For N-D: generalizes the pattern (row-major order)
 *
 * @param[in] p   Point to index. Must be within the box [s, e).
 *                Array of int16_t with 'dim' elements.
 * @param[in] s   Box start (minimum corner, inclusive).
 *                Array of int16_t with 'dim' elements.
 * @param[in] e   Box end (maximum corner, exclusive).
 *                Array of int16_t with 'dim' elements.
 * @param[in] dim Number of dimensions.
 *
 * @return Linear index in range [0, volume-1] where volume is the product
 *         of (e[i] - s[i]) for all dimensions. Returns uint64_t to handle
 *         large volumes, but overflow is still possible.
 *
 * @warning Undefined behavior if p is outside the box [s, e). The function
 *          does not check bounds. Results may be incorrect or overflow.
 *
 * @warning For very large boxes, the index may overflow uint64_t.
 *          Keep bounding boxes reasonably sized.
 *
 * @note The point p[i] must satisfy: s[i] <= p[i] < e[i] for all dimensions.
 *
 * @note Row-major ordering means the first dimension (p[0]) varies fastest,
 *       last dimension (p[dim-1]) varies slowest. This matches C array layout.
 *
 * Example (2D grid):
 * @code
 * int16_t start[2] = {0, 0};
 * int16_t end[2] = {10, 10};  // 10x10 grid
 * int16_t point[2] = {3, 5};
 * uint64_t idx = point_idx(point, start, end, 2);
 * // idx = 5 * 10 + 3 = 53
 * @endcode
 *
 * Example (3D array access):
 * @code
 * int16_t start[3] = {0, 0, 0};
 * int16_t end[3] = {100, 100, 100};
 * uint32_t *voxels = malloc(100 * 100 * 100 * sizeof(uint32_t));
 * 
 * int16_t pos[3] = {10, 20, 30};
 * uint64_t idx = point_idx(pos, start, end, 3);
 * voxels[idx] = 42;  // Set voxel at (10,20,30)
 * @endcode
 *
 * @see point_vol
 * @see geo_iter
 */
static inline uint64_t
point_idx(int16_t *p, int16_t *s, int16_t *e, uint8_t dim)
{
	uint64_t res = 0;
	uint64_t mult = 1;

	for (uint8_t i = 0; i < dim; i++) {
		res += (p[i] - s[i]) * mult;
		mult *= (e[i] - s[i]);
	}

	return res;
}

/** @} */

#endif
