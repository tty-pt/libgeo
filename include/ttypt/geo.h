#ifndef GEO_H
#define GEO_H

/**
 * @file geo.h
 * @brief Public API for libgeo.
 */

#include <stdint.h>
#include <ttypt/qmap.h>

#include "morton.h"

/** @defgroup geo_core Geo core API
 *  @brief Spatial map operations.
 *  @{
 */

/**
 * Sentinel value returned by geo_get when no entry exists.
 */
#define GEO_MISS (130056652770671ULL)

/**
 * Initialize the geo subsystem.
 */
void geo_init(void);

/**
 * Open or create a database.
 *
 * filename and database can be NULL for in-memory only.
 *
 * @param filename File path for persistence, or NULL.
 * @param database Database name within the file, or NULL.
 * @param mask Hash mask (table size is mask + 1).
 * @return Database handle.
 */
uint32_t geo_open(char *filename, char *database, uint32_t mask);

/**
 * Create an iterator for a rectangular region.
 *
 * @param pdb_hd Database handle.
 * @param s Start point.
 * @param l Lengths per dimension.
 * @param dim Number of dimensions.
 * @return Iterator handle.
 */
uint32_t geo_iter(uint32_t pdb_hd, int16_t *s,
		uint16_t *l, uint8_t dim);

/**
 * Advance an iterator created by geo_iter.
 * Returns 1 on success, 0 when finished.
 *
 * @param p Output point.
 * @param ref Output reference value.
 * @param cur Iterator handle.
 * @return 1 if advanced, 0 when finished.
 */
int geo_next(int16_t *p, uint32_t *ref, uint32_t cur);

/**
 * Delete the entry at a point.
 *
 * @param pdb_hd Database handle.
 * @param p Point.
 * @param dim Number of dimensions.
 */
static inline void
geo_del(uint32_t pdb_hd, int16_t *p, uint8_t dim)
{
	uint64_t at = morton_set(p, dim);
	qmap_del(pdb_hd, &at);
}

/**
 * Fetch the value at a point.
 *
 * @param pdb_hd Database handle.
 * @param p Point.
 * @param dim Number of dimensions.
 * @return Stored value or GEO_MISS.
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
 * Store a value at a point.
 *
 * @param pdb_hd Database handle.
 * @param p Point.
 * @param thing Value to store.
 * @param dim Number of dimensions.
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
