#ifndef MORTON_H
#define MORTON_H

/**
 * @file morton.h
 * @brief Morton code utilities.
 */

#include <stdint.h>

/** @defgroup geo_morton Morton code helpers
 *  @brief Encode and decode Morton (Z-order) codes.
 *  @{
 */

/**
 * Encode a coordinate into a Morton code.
 *
 * @param pos Input point.
 * @param dim Number of dimensions.
 * @return Morton code.
 */
uint64_t morton_set(int16_t *pos, uint8_t dim);

/**
 * Decode a Morton code into a coordinate.
 *
 * @param pos Output point.
 * @param morton Morton code.
 * @param dim Number of dimensions.
 */
void morton_get(int16_t *pos, uint64_t morton, uint8_t dim);

/** @} */

#endif
