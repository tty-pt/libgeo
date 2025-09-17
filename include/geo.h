#ifndef GEO_H
#define GEO_H

#include <stdint.h>
#include <qmap.h>

#include "morton.h"

#define GEO_MISS (130056652770671ULL)

void geo_init(void);
unsigned geo_open(char *database, unsigned mask);
void geo_search(unsigned pdb_hd, unsigned *mat,
		int16_t *quad_s, int16_t *quad_e, uint8_t dim);

static inline uint64_t
geo_mwhere(unsigned pdb_hd, unsigned where)
{
	const void *data;
	data = qmap_get(pdb_hd, &where);

	if (!data)
		return GEO_MISS;

	return * (uint64_t *) data;
}

static inline int
geo_has(unsigned pdb_hd, unsigned room) {
	return geo_mwhere(pdb_hd, room) != GEO_MISS;
}

static inline void
geo_where(unsigned pdb_hd, int16_t *p,
		unsigned room, uint8_t dim)
{
	uint64_t x = geo_mwhere(pdb_hd, room);
	morton_get(p, x, dim);
}

static inline void
geo_del(unsigned pdb_hd, unsigned what)
{
	uint64_t code = geo_mwhere(pdb_hd, what);
	qmap_del(pdb_hd + 1, &code);
}

static inline unsigned
geo_get(unsigned pdb_hd, int16_t *p, uint8_t dim)
{
	uint64_t at = morton_set(p, dim);
	const void *ref = qmap_get(pdb_hd + 1, &at);

	if (ref)
		return * (unsigned *) ref;

	return QM_MISS;
}

static inline void
geo_put(unsigned pdb_hd, int16_t *p,
		unsigned thing, uint8_t dim)
{
	uint64_t code = morton_set(p, dim);
	qmap_put(pdb_hd, &thing, &code);
}

#endif
