#ifndef MORTON_H
#define MORTON_H

/**
 * @file morton.h
 * @brief Morton code (Z-order) encoding and decoding utilities.
 *
 * Morton codes, also known as Z-order codes, are a way to map multi-dimensional
 * coordinates to a single dimension while preserving spatial locality. They
 * interleave the binary representations of each coordinate dimension, creating
 * a space-filling curve that visits all points in a recursive Z pattern.
 *
 * Benefits:
 * - Preserves spatial locality: nearby points in N-D space have nearby Morton codes
 * - Enables efficient range queries on 1D data structures (hash tables, B-trees)
 * - Simple bitwise operations for encoding/decoding
 *
 * This implementation is optimized for 3D coordinates (int16_t per dimension)
 * and produces 64-bit Morton codes (48 bits used for 3D, 16 bits reserved).
 * 4D coordinates use a dense stride-4 packing that fills all 64 bits.
 * The 2D x 32-bit config (int32_t per dimension) is the dense
 * stride-2 partner: 2 x 32 = 64 bits, no reserved bits.
 *
 * Prefer the config objects (pointcfg.h: Point1_2..Point4_2, Point2_4)
 * in application code; these inlines remain the tight-loop fast path.
 *
 * References:
 * - Morton, G.M. (1966). "A computer Oriented Geodetic Data Base"
 * - http://www.vision-tools.com/h-tropf/multidimensionalrangequery.pdf
 * - https://en.wikipedia.org/wiki/Z-order_curve
 *
 * @see geo_core
 */

#include <stdint.h>
#include <limits.h>

/** @defgroup geo_morton Morton code helpers
 *  @brief Encode and decode Morton (Z-order) codes for spatial indexing.
 *
 *  Morton codes interleave the bits of multi-dimensional coordinates to create
 *  a single linear ordering that preserves spatial locality. For example, in 2D:
 *
 *  Coordinate (x=5, y=3):
 *  - x = 5 = 0b0101 (binary)
 *  - y = 3 = 0b0011 (binary)
 *  - Interleaved: 0b00110011 = 51 (y bits in odd positions, x in even)
 *
 *  In 3D, the pattern extends with z bits in every third position.
 *
 *  Coordinate System:
 *  - Input: int16_t signed coordinates (-32768 to 32767)
 *  - Internal: uint16_t unsigned coordinates (0 to 65535) via offset mapping
 *  - Output: uint64_t Morton code (48 bits for 3D with the upper 16 bits
 *    reserved; all 64 bits for 4D)
 *
 *  The conversion from signed to unsigned preserves ordering:
 *  -32768 → 0, -32767 → 1, ..., 0 → 32768, ..., 32767 → 65535
 *
 *  @note Currently optimized for 3D/4D with fast bit-manipulation
 *        algorithms (stride-3 packing for dims 1-3, dense stride-4
 *        packing for dim 4, dense stride-2 packing for the 2D x 32-bit
 *        config). 1D/2D/3D codes are bit-identical to
 *        v0.5.0; 4D codes are new.
 *
 *  @note All dimensions share the uint64 key space: use one dimension
 *        per database.
 *
 *  @see geo_core
 *  @see geo_point
 *  @{
 */

/* When libgeo.c builds the external ABI wrappers it defines
 * GEO_MORTON_RENAME_FOR_WRAPPERS before including this header, so the
 * static inline versions take the _il names and don't clash with the
 * extern definitions in that TU. Consumers never define it. */
#ifdef GEO_MORTON_RENAME_FOR_WRAPPERS
#define morton_set_1  morton_set_1_il
#define morton_set_2  morton_set_2_il
#define morton_set_3  morton_set_3_il
#define morton_set_4  morton_set_4_il
#define morton_set_2_32  morton_set_2_32_il
#define morton_get_1  morton_get_1_il
#define morton_get_2  morton_get_2_il
#define morton_get_3  morton_get_3_il
#define morton_get_4  morton_get_4_il
#define morton_get_2_32  morton_get_2_32_il
#endif

static inline uint16_t
geo_unsign(int16_t n)
{
	return (uint16_t)(n + SHRT_MAX + 1);
}

static inline int16_t
geo_sign(uint16_t n)
{
	return (int16_t)(n - SHRT_MAX - 1);
}

/* 32-bit-lane bias helpers for the 2D x 32-bit dense config. The
 * unsigned lane is the full uint32 range: -2^31 -> 0, ..., 0 ->
 * 0x80000000, ..., 2^31-1 -> 0xFFFFFFFF. All arithmetic is modulo
 * 2^32, so the add/sub wrap exactly like the 16-bit pair above. */
static inline uint32_t
geo_unsign32(int32_t n)
{
	return (uint32_t)n + 0x80000000u;
}

static inline int32_t
geo_sign32(uint32_t n)
{
	return (int32_t)(n - 0x80000000u);
}

static inline uint64_t
geo_spread3(uint32_t x)
{
	uint64_t v = x & 0xFFFFu;

	v = (v | (v << 32)) & 0x1F00000000FFFFULL;
	v = (v | (v << 16)) & 0x1F0000FF0000FFULL;
	v = (v | (v << 8)) & 0x100F00F00F00F00FULL;
	v = (v | (v << 4)) & 0x10C30C30C30C30C3ULL;
	v = (v | (v << 2)) & 0x1249249249249249ULL;

	return v;
}

static inline uint64_t
geo_compact_axis(uint64_t code, uint32_t shift)
{
	code >>= shift;
	code &= 0x1249249249249249ULL;
	code = (code ^ (code >> 2))  & 0x10C30C30C30C30C3ULL;
	code = (code ^ (code >> 4))  & 0x100F00F00F00F00FULL;
	code = (code ^ (code >> 8))  & 0x1F0000FF0000FFULL;
	code = (code ^ (code >> 16)) & 0x1F00000000FFFFULL;
	code = (code ^ (code >> 32)) & 0x00000000001FFFFFULL;
	return (uint32_t) code;
}

/* spread2(x):
 *   Take x in [0..0xFFFFFFFF] and produce a 64-bit word where its
 *   bit-i goes to bit-(2*i) in the result (even positions). 2D x 32-bit
 *   needs all 64 bits (2 x 32), so there are no reserved top bits.
 *
 * Part of a Morton-2D32 encode:  code = spread2(x) | spread2(y)<<1
 */
static inline uint64_t
geo_spread2(uint32_t x)
{
	uint64_t v = x;

	v = (v | (v << 16)) & 0x0000FFFF0000FFFFULL;
	v = (v | (v << 8))  & 0x00FF00FF00FF00FFULL;
	v = (v | (v << 4))  & 0x0F0F0F0F0F0F0F0FULL;
	v = (v | (v << 2))  & 0x3333333333333333ULL;
	v = (v | (v << 1))  & 0x5555555555555555ULL;

	return v;
}

/* compact_axis2(): collect one out of every 2 bits from 'code',
 * starting at 'shift' (0 = x, 1 = y).
 * Returns the low-order 32 bits containing that coordinate.
 */
static inline uint64_t
geo_compact_axis2(uint64_t code, uint32_t shift)
{
	code >>= shift;
	code &= 0x5555555555555555ULL;
	code = (code ^ (code >> 1))  & 0x3333333333333333ULL;
	code = (code ^ (code >> 2))  & 0x0F0F0F0F0F0F0F0FULL;
	code = (code ^ (code >> 4))  & 0x00FF00FF00FF00FFULL;
	code = (code ^ (code >> 8))  & 0x0000FFFF0000FFFFULL;
	code = (code ^ (code >> 16)) & 0x00000000FFFFFFFFULL;
	return (uint32_t) code;
}

/* spread4(x):
 *   Take x ∈ [0..0xFFFF] and produce a 64-bit word where
 *   its bit-i goes to bit-(4*i) in the result. 4D needs all 64 bits
 *   (4 × 16), so unlike the 3D path there are no reserved top bits.
 *
 * Part of a Morton-4D encode:  code = spread4(x)
 *                                   | spread4(y)<<1
 *                                   | spread4(z)<<2
 *                                   | spread4(w)<<3
 */
static inline uint64_t
geo_spread4(uint32_t x)
{
	uint64_t v = x & 0xFFFFu;

	/* separate high/low bytes with a 24-bit gap */
	v = (v | (v << 24)) & 0xFF000000FFULL;
	/* down to nibbles */
	v = (v | (v << 12)) & 0x000F000F000F000FULL;
	/* down to bit pairs */
	v = (v | (v << 6)) & 0x0303030303030303ULL;
	/* final 4-bit interleave */
	v = (v | (v << 3)) & 0x1111111111111111ULL;

	return v;
}

/* compact_axis4(): collect one out of every 4 bits from 'code',
 * starting at 'shift' (0 = x, 1 = y, 2 = z, 3 = w).
 * Returns the low-order 16 bits containing that coordinate.
 */
static inline uint64_t
geo_compact_axis4(uint64_t code, uint32_t shift)
{
	code >>= shift;
	code &= 0x1111111111111111ULL;
	code = (code ^ (code >> 3))  & 0x0303030303030303ULL;
	code = (code ^ (code >> 6))  & 0x000F000F000F000FULL;
	code = (code ^ (code >> 12)) & 0x000000FF000000FFULL;
	code = (code ^ (code >> 24)) & 0x000000000000FFFFULL;
	return (uint16_t) code;
}

static inline void
geo_decode3(uint64_t code, uint32_t *x, uint32_t *y, uint32_t *z)
{
	static const uint64_t mask_off = 0x0000FFFFFFFFFFFFULL;
	code &= mask_off;
	*x = geo_compact_axis(code, 0);
	*y = geo_compact_axis(code, 1);
	*z = geo_compact_axis(code, 2);
}

static inline void
geo_decode4(uint64_t code,
	    uint32_t *x, uint32_t *y, uint32_t *z, uint32_t *w)
{
	*x = (uint32_t)geo_compact_axis4(code, 0);
	*y = (uint32_t)geo_compact_axis4(code, 1);
	*z = (uint32_t)geo_compact_axis4(code, 2);
	*w = (uint32_t)geo_compact_axis4(code, 3);
}

/**
 * @brief Encode multi-dimensional coordinates into Morton codes.
 *
 * These per-dimension specializations are the encode API. Each encodes
 * exactly one dimension count; there is no runtime `dim` argument to
 * branch on, so the compiler sees a fully unrolled expression.
 *
 * Algorithm (dims 1-3):
 * 1. Convert signed int16_t to unsigned uint16_t (add 32768)
 * 2. Spread each 16-bit coordinate across 48 bits (every 3rd bit)
 * 3. Combine: code = spread(x) | spread(y)<<1 | spread(z)<<2
 *
 * Dim 4 uses dense stride-4 packing that fills all 64 bits. Dims 1-3
 * code values are bit-identical to v0.5.0; 4D codes are new.
 *
 * @param[in] p Input point. Array of int16_t with at least the
 *              dimension count of the function used. Coordinates range
 *              from -32768 to 32767.
 *
 * @return 64-bit Morton code. For 1-3D, uses 48 bits (16 bits per
 *         dimension) with the upper 16 bits reserved; 4D fills all 64
 *         bits.
 *
 * Example (3D):
 * @code
 * int16_t point[3] = {10, -5, 100};
 * uint64_t code = morton_set_3(point);
 * // code now contains the interleaved bit representation
 * @endcode
 *
 * Example (round-trip verification):
 * @code
 * int16_t original[3] = {123, -456, 789};
 * uint64_t code = morton_set_3(original);
 * int16_t decoded[3];
 * morton_get_3(decoded, code);
 * // decoded[0]==123, decoded[1]==-456, decoded[2]==789
 * @endcode
 *
 * @see morton_set_1 morton_set_2 morton_set_4
 * @see morton_get_1 morton_get_2 morton_get_3 morton_get_4
 * @see geo_put
 * @see geo_get
 */
static inline uint64_t
morton_set_1(int16_t *p)
{
	uint16_t up0 = geo_unsign(p[0]);

	return geo_spread3(up0);
}

/**
 * @brief Encode a 2D point. See morton_set_1() for the family docs.
 */
static inline uint64_t
morton_set_2(int16_t *p)
{
	uint16_t up0 = geo_unsign(p[0]);
	uint16_t up1 = geo_unsign(p[1]);

	return geo_spread3(up0) | (geo_spread3(up1) << 1);
}

/**
 * @brief Encode a 3D point. See morton_set_1() for the family docs.
 */
static inline uint64_t
morton_set_3(int16_t *p)
{
	uint16_t up0 = geo_unsign(p[0]);
	uint16_t up1 = geo_unsign(p[1]);
	uint16_t up2 = geo_unsign(p[2]);

	return geo_spread3(up0)
		| (geo_spread3(up1) << 1)
		| (geo_spread3(up2) << 2);
}

/**
 * @brief Encode a 4D point. See morton_set_1() for the family docs.
 */
static inline uint64_t
morton_set_4(int16_t *p)
{
	uint16_t up0 = geo_unsign(p[0]);
	uint16_t up1 = geo_unsign(p[1]);
	uint16_t up2 = geo_unsign(p[2]);
	uint16_t up3 = geo_unsign(p[3]);

	return geo_spread4(up0)
		| (geo_spread4(up1) << 1)
		| (geo_spread4(up2) << 2)
		| (geo_spread4(up3) << 3);
}

/**
 * @brief Encode a 2D x 32-bit point. Dense full-64-bit layout
 *        (2 x 32 bits); part of the dense config family, not the
 *        legacy 16-bit-lane sparse layout that morton_set_2() uses.
 *
 * Algorithm:
 * 1. Convert signed int32_t to unsigned uint32_t (add 2^31)
 * 2. Spread each 32-bit coordinate across 64 bits (every 2nd bit)
 * 3. Combine: code = spread2(x) | spread2(y)<<1
 *
 * @param[in] p Input point. Array of int32_t with at least 2 elements.
 *              Coordinates range from -2147483648 to 2147483647.
 *
 * @return 64-bit Morton code filling all 64 bits (no reserved bits).
 */
static inline uint64_t
morton_set_2_32(int32_t *p)
{
	uint32_t up0 = geo_unsign32(p[0]);
	uint32_t up1 = geo_unsign32(p[1]);

	return geo_spread2(up0) | (geo_spread2(up1) << 1);
}

/**
 * @brief Decode a Morton code into a multi-dimensional coordinate.
 *
 * These per-dimension specializations are the decode API. Each decodes
 * exactly one dimension count; there is no runtime `dim` argument.
 * Inverse of the matching morton_set_N().
 *
 * @param[out] pos  Output point. Array of int16_t with space for the
 *                  dimension count of the function used.
 * @param[in]  code Morton code to decode (64-bit).
 *
 * @note The output coordinates will be in range -32768 to 32767 (int16_t).
 *
 * @warning Decoding a Morton code that wasn't created by morton_set_N()
 *          with valid coordinates may produce unexpected results (garbage
 *          coordinates).
 *
 * Example:
 * @code
 * uint64_t code = 0x123456789ABCULL;  // Some Morton code
 * int16_t point[3];
 * morton_get_3(point, code);
 * // point now contains the decoded coordinates
 * @endcode
 *
 * @see morton_get_1 morton_get_2 morton_get_4
 * @see morton_set_1 morton_set_2 morton_set_3 morton_set_4
 * @see geo_iter
 */
static inline void
morton_get_1(int16_t *pos, uint64_t code)
{
	pos[0] = geo_sign((uint16_t)geo_compact_axis(code, 0));
}

/**
 * @brief Decode a 2D Morton code. See morton_get_1() for the family docs.
 */
static inline void
morton_get_2(int16_t *pos, uint64_t code)
{
	pos[0] = geo_sign((uint16_t)geo_compact_axis(code, 0));
	pos[1] = geo_sign((uint16_t)geo_compact_axis(code, 1));
}

/**
 * @brief Decode a 3D Morton code. See morton_get_1() for the family docs.
 */
static inline void
morton_get_3(int16_t *pos, uint64_t code)
{
	uint32_t uup[] = { 0, 0, 0 };

	geo_decode3(code, &uup[0], &uup[1], &uup[2]);
	pos[0] = geo_sign((uint16_t)uup[0]);
	pos[1] = geo_sign((uint16_t)uup[1]);
	pos[2] = geo_sign((uint16_t)uup[2]);
}

/**
 * @brief Decode a 4D Morton code. See morton_get_1() for the family docs.
 */
static inline void
morton_get_4(int16_t *pos, uint64_t code)
{
	uint32_t uup[] = { 0, 0, 0, 0 };

	geo_decode4(code, &uup[0], &uup[1], &uup[2], &uup[3]);
	pos[0] = geo_sign((uint16_t)uup[0]);
	pos[1] = geo_sign((uint16_t)uup[1]);
	pos[2] = geo_sign((uint16_t)uup[2]);
	pos[3] = geo_sign((uint16_t)uup[3]);
}

/**
 * @brief Decode a 2D x 32-bit Morton code. Inverse of morton_set_2_32().
 *
 * @param[out] pos  Output point. Array of int32_t with space for
 *                  2 elements.
 * @param[in]  code Morton code to decode (64-bit, dense layout).
 */
static inline void
morton_get_2_32(int32_t *pos, uint64_t code)
{
	pos[0] = geo_sign32((uint32_t)geo_compact_axis2(code, 0));
	pos[1] = geo_sign32((uint32_t)geo_compact_axis2(code, 1));
}

/** @} */

#endif
