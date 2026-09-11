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

#ifndef FAST_MORTON
#define FAST_MORTON 1
#endif

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
 *        packing for dim 4). 1D/2D/3D codes are bit-identical to
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
#define morton_set  morton_set_il
#define morton_get  morton_get_il
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
 * @brief Encode a multi-dimensional coordinate into a Morton code.
 *
 * Converts a signed coordinate point to a single 64-bit Morton code by
 * interleaving the bits of each dimension. The encoding preserves spatial
 * locality: points that are close in N-dimensional space will have similar
 * Morton codes.
 *
 * Algorithm (3D):
 * 1. Convert signed int16_t to unsigned uint16_t (add 32768)
 * 2. Spread each 16-bit coordinate across 48 bits (every 3rd bit)
 * 3. Combine: code = spread(x) | spread(y)<<1 | spread(z)<<2
 *
 * Dim 4 uses dense stride-4 packing over all 64 bits. Dims 1-3 use the
 * stride-3 packing above and are bit-identical to v0.5.0.
 *
 * @param[in] p   Input point. Array of int16_t with at least 'dim' elements.
 *                Coordinates range from -32768 to 32767.
 * @param[in] dim Number of dimensions (1..4).
 *
 * @return 64-bit Morton code. For 3D, uses 48 bits (16 bits per dimension)
 *         with the upper 16 bits reserved; 4D fills all 64 bits.
 *
 * Example (3D):
 * @code
 * int16_t point[3] = {10, -5, 100};
 * uint64_t code = morton_set(point, 3);
 * // code now contains the interleaved bit representation
 * @endcode
 *
 * Example (round-trip verification):
 * @code
 * int16_t original[3] = {123, -456, 789};
 * uint64_t code = morton_set(original, 3);
 * int16_t decoded[3];
 * morton_get(decoded, code, 3);
 * // decoded[0]==123, decoded[1]==-456, decoded[2]==789
 * @endcode
 *
 * @see morton_get
 * @see geo_put
 * @see geo_get
 */
static inline uint64_t
morton_set(int16_t *p, uint8_t dim)
{
	/* NOTE: literal 4 here — MAX_DIM lives in libgeo.c, and this
	 * header must stay standalone. Keep the two in sync. */
	uint16_t up[4] = {0, 0, 0, 0};

	for (uint8_t i = 0; i < dim && i < 4; i++)
		up[i] = geo_unsign(p[i]);

#if FAST_MORTON
	switch (dim) {
	case 1:
		return geo_spread3(up[0]);
	case 2:
		return geo_spread3(up[0]) | (geo_spread3(up[1]) << 1);
	case 3:
		return geo_spread3(up[0])
			| (geo_spread3(up[1]) << 1)
			| (geo_spread3(up[2]) << 2);
	default:
		/* dim == 4 (dim == 0 encodes all-zero coords to 0, same as
		 * the old 3D-default path did). dim > 4 is invalid input. */
		return geo_spread4(up[0])
			| (geo_spread4(up[1]) << 1)
			| (geo_spread4(up[2]) << 2)
			| (geo_spread4(up[3]) << 3);
	}
#else
	{
		uint64_t mask = 0x1;
		uint64_t result = 0;
		/* Stride must match the FAST paths above: stride-3 packing
		 * for dims 1-3, stride-4 for dim 4. */
		uint8_t dd = (dim == 4) ? 4 : 3;

		for (uint8_t b = 0; b < 16; b++, mask <<= 1)
			for (uint8_t i = 0; i < dd; i++)
				result |= (up[i] & mask) >> b
					<< ((b * dd) + i);
		return result;
	}
#endif
}

/**
 * @brief Decode a Morton code into a multi-dimensional coordinate.
 *
 * Converts a 64-bit Morton code back to the original N-dimensional coordinate
 * by de-interleaving the bits. This is the inverse operation of morton_set().
 *
 * Algorithm (3D):
 * 1. Extract every 3rd bit for each dimension (compact operation)
 * 2. Convert unsigned uint16_t to signed int16_t (subtract 32768)
 * 3. Store in output array
 *
 * @param[out] pos    Output point. Array of int16_t with space for at least
 *                    'dim' elements. Filled with decoded coordinates.
 * @param[in]  code   Morton code to decode (64-bit).
 * @param[in]  dim    Number of dimensions to decode (1..4).
 *
 * @note The output coordinates will be in range -32768 to 32767 (int16_t).
 *
 * @warning Decoding a Morton code that wasn't created by morton_set() with
 *          valid coordinates may produce unexpected results (garbage coordinates).
 *
 * Example:
 * @code
 * uint64_t code = 0x123456789ABCULL;  // Some Morton code
 * int16_t point[3];
 * morton_get(point, code, 3);
 * // point now contains the decoded coordinates
 * @endcode
 *
 * @see morton_set
 * @see geo_iter
 */
static inline void
morton_get(int16_t *pos, uint64_t code, uint8_t dim)
{
	uint32_t uup[] = { 0, 0, 0, 0 };

#if FAST_MORTON
	if (dim == 4)
		geo_decode4(code, &uup[0], &uup[1], &uup[2], &uup[3]);
	else
		geo_decode3(code, &uup[0], &uup[1], &uup[2]);
#else
	{
		uint8_t dd = (dim == 4) ? 4 : 3;

		for (uint8_t b = 0; b < 16; b++)
			for (uint8_t i = 0; i < dd; i++)
				uup[i] |= ((code >> (b * dd + i)) & 0x1) << b;
	}
#endif

	/* Clamped to 4: dims above that are invalid input; never read
	 * past uup. */
	for (uint8_t i = 0; i < dim && i < 4; i++)
		pos[i] = geo_sign(uup[i]);
}


/** @} */

#endif
