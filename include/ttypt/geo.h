#ifndef GEO_H
#define GEO_H

#include <stdint.h>
#include <ttypt/qmap.h>

#include "morton.h"

#define GEO_MISS (130056652770671ULL)

void geo_init(void);
unsigned geo_open(char *database, unsigned mask);

unsigned geo_iter(unsigned pdb_hd, int16_t *s,
		uint16_t *l, uint8_t dim);

int geo_next(int16_t *p, unsigned *ref, unsigned cur);

static inline void
geo_del(unsigned pdb_hd, int16_t *p, uint8_t dim)
{
	uint64_t at = morton_set(p, dim);
	qmap_del(pdb_hd, &at);
}

static inline unsigned
geo_get(unsigned pdb_hd, int16_t *p, uint8_t dim)
{
	uint64_t at = morton_set(p, dim);
	const void *ref = qmap_get(pdb_hd, &at);

	if (ref)
		return * (unsigned *) ref;

	return QM_MISS;
}

static inline void
geo_put(unsigned pdb_hd, int16_t *p,
		unsigned thing, uint8_t dim)
{
	uint64_t code = morton_set(p, dim);
	qmap_put(pdb_hd, &code, &thing);
}

#endif
