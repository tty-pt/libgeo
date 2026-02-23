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
 *
 * References:
 * - Morton, G.M. (1966). "A computer Oriented Geodetic Data Base"
 * - http://www.vision-tools.com/h-tropf/multidimensionalrangequery.pdf
 * - https://en.wikipedia.org/wiki/Z-order_curve
 *
 * @see geo_core
 */

#include <stdint.h>

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
 *  - Output: uint64_t Morton code (48 bits for 3D, upper 16 bits reserved)
 *
 *  The conversion from signed to unsigned preserves ordering:
 *  -32768 → 0, -32767 → 1, ..., 0 → 32768, ..., 32767 → 65535
 *
 *  @note Currently optimized for 3D with fast bit-manipulation algorithms.
 *        Designed for future expansion to arbitrary dimensions.
 *
 *  @see geo_core
 *  @see geo_point
 *  @{
 */

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
 * @param[in] pos Input point. Array of int16_t with at least 'dim' elements.
 *                Coordinates range from -32768 to 32767.
 * @param[in] dim Number of dimensions. Currently optimized for 3D, designed
 *                for future multi-dimensional support.
 *
 * @return 64-bit Morton code. For 3D, uses 48 bits (16 bits per dimension).
 *         Upper 16 bits are reserved for future use (e.g., world layer).
 *
 * @note The current implementation uses optimized bit-manipulation for 3D.
 *       The dim parameter is accepted for future expansion but currently
 *       assumes dim=3 in the optimized path.
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
uint64_t morton_set(int16_t *pos, uint8_t dim);

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
 * @param[in]  morton Morton code to decode (64-bit).
 * @param[in]  dim    Number of dimensions to decode. Currently optimized for 3D.
 *
 * @note The output coordinates will be in range -32768 to 32767 (int16_t).
 *
 * @note The current implementation uses optimized bit-manipulation for 3D.
 *       The dim parameter is accepted for future expansion but currently
 *       assumes dim=3 in the optimized path.
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
void morton_get(int16_t *pos, uint64_t morton, uint8_t dim);

/** @} */

#endif
