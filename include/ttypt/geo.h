#ifndef GEO_H
#define GEO_H

#include <stdint.h>
#include <ttypt/qmap.h>

#include "morton.h"

#define GEO_MISS (130056652770671ULL)

void geo_init(void);
uint32_t geo_open(char *filename, char *database, uint32_t mask);

uint32_t geo_iter(uint32_t pdb_hd, int16_t *s,
		uint16_t *l, uint8_t dim);

int geo_next(int16_t *p, uint32_t *ref, uint32_t cur);

static inline void
geo_del(uint32_t pdb_hd, int16_t *p, uint8_t dim)
{
	uint64_t at = morton_set(p, dim);
	qmap_del(pdb_hd, &at);
}

static inline uint32_t
geo_get(uint32_t pdb_hd, int16_t *p, uint8_t dim)
{
	uint64_t at = morton_set(p, dim);
	const void *ref = qmap_get(pdb_hd, &at);

	if (ref)
		return * (uint32_t *) ref;

	return QM_MISS;
}

static inline void
geo_put(uint32_t pdb_hd, int16_t *p,
		uint32_t thing, uint8_t dim)
{
	uint64_t code = morton_set(p, dim);
	qmap_put(pdb_hd, &code, &thing);
}

#endif
