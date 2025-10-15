/* see http://www.vision-tools.com/h-tropf/multidimensionalrangequery.pdf
 */

#include "../include/ttypt/geo.h"
#include "../include/ttypt/point.h"
#include "../include/ttypt/morton.h"

#include <limits.h>

#include <ttypt/qsys.h>
#include <ttypt/idm.h>

typedef struct {
	int16_t p[4];
	unsigned ref;
} geo_curi_t;

typedef struct {
	geo_curi_t *items;
	unsigned m, pos;
	uint8_t dim;
} geo_cur_t;

#define MAX_DIM 3
#define FAST_MORTON 1
#define COMPUTE_BMLM 0

static const uint64_t m1 = 0x924924924924ULL;
static const uint64_t m0 = 0x400000000000ULL;

static unsigned qm_u, qm_u64;

static idm_t geo_idm;

geo_cur_t geo_cursors[1024];

static inline uint16_t
unsign(int16_t n)
{
	return (uint16_t)(n + SHRT_MAX + 1);
}

static inline int16_t
sign(uint16_t n)
{
	return (int16_t)(n - SHRT_MAX - 1);
}

/* spread3(x):
 *   Take x ∈ [0..0xFFFF] and produce a 64-bit word where
 *   its bit-i goes to bit-(3*i) in the result.
 *
 * Part of a Morton-3D encode:  code = spread3(x)
 *                                  | spread3(y)<<1
 *                                  | spread3(z)<<2
 */
static inline uint64_t spread3(uint32_t x)
{
	/* keep only low 16 bits */
	uint64_t v = x & 0xFFFFu;  

	/* make room for high triples */
	v = (v | (v << 32)) & 0x1F00000000FFFFULL;
	/* down to 8-bit chunks */
	v = (v | (v << 16)) & 0x1F0000FF0000FFULL;
	/* down to 4-bit groups */
	v = (v | (v << 8)) & 0x100F00F00F00F00FULL;
	/* down to 2-bit groups */
	v = (v | (v << 4)) & 0x10C30C30C30C30C3ULL;
	/* final 3 bit interleave */
	v = (v | (v << 2)) & 0x1249249249249249ULL;

	return v;
}

static inline uint64_t morton3_pack_u16(
		uint16_t x,
		uint16_t y,
		uint16_t z,
		uint16_t world)
{
	return spread3(x)
		| (spread3(y) << 1)
		| (spread3(z) << 2)
		| ((uint64_t)world << 48);
}

uint64_t
morton_set(int16_t *p, uint8_t dim)
{
	uint16_t up[dim];

	for (uint8_t i = 0; i < dim; i++)
		up[i] = unsign(p[i]);

#if FAST_MORTON
	return morton3_pack_u16(up[0], up[1], up[2], 0);
		/* | ((uint64_t) p[3] << 48); */
#else
	uint64_t mask = 0x1;
	uint64_t result = 0;

	for (uint8_t b = 0; b < 16; b++, mask <<= 1)
		for (uint8_t i = 0; i < MAX_DIM; i++)
			result |= (up[i] & mask) >> b
				<< ((b * MAX_DIM) + i);
	return result;
#endif

}

/* compact_axis(): collect one out of every 3 bits from 'code',
 * starting at 'shift' (0 = x, 1 = y, 2 = z).
 * Returns low-order 21 bits containing that coordinate.
 */
static inline uint32_t
compact_axis(uint64_t code, unsigned shift)
{
    code >>= shift;
    /* align the desired series to LSB */
    /* first keep only 1---1---1
     * pattern → mask 0x1249249249249… */
    code &= 0x1249249249249249ULL;

    /* Now collapse gaps:  3→2 → 2→1 → 1→0 */
    code = (code ^ (code >> 2))  & 0x10C30C30C30C30C3ULL;
    code = (code ^ (code >> 4))  & 0x100F00F00F00F00FULL;
    code = (code ^ (code >> 8))  & 0x1F0000FF0000FFULL;
    code = (code ^ (code >> 16)) & 0x1F00000000FFFFULL;
    code = (code ^ (code >> 32)) & 0x00000000001FFFFFULL;

    /* low 21 bits hold the axis value */
    return (uint32_t) code;
}

static inline void decode3(uint64_t code,
                           uint32_t *x,
			   uint32_t *y,
			   uint32_t *z)
{
    *x = compact_axis(code, 0);   /* bits 0,3,6,…   */
    *y = compact_axis(code, 1);   /* bits 1,4,7,…   */
    *z = compact_axis(code, 2);   /* bits 2,5,8,…   */
}

void
morton_get(int16_t *pos, uint64_t code, uint8_t dim)
{
	static const uint64_t mask_off = 0x0000FFFFFFFFFFFFULL;
	uint32_t uup[] = { 0, 0, 0, 0 };

#if FAST_MORTON
	decode3(code & mask_off, &uup[0], &uup[1], &uup[2]);
#else
	for (uint8_t b = 0; b < 16; b++)
		for (uint8_t i = 0; i < MAX_DIM; i++)
			uup[i] |= ((code >> (b * MAX_DIM + i)) & 0x1) << b;
#endif

	for (uint8_t i = 0; i < dim; i++)
		pos[i] = sign(uup[i]);
}


static inline uint64_t
qload_0(uint64_t c, uint64_t lm0, uint64_t lm1) // LOAD(0111...
{
	return (c & (~ lm0)) | lm1;
}

static inline uint64_t
qload_1(uint64_t c, uint64_t lm0, uint64_t lm1) // LOAD(1000...
{
	return (c | lm0) & (~ lm1);
}

static inline int
inrange_p(int16_t *drp, int16_t *min, int16_t *max, uint8_t dim)
{
	for (uint8_t i = 0; i < dim; i++)
		if (drp[i] < min[i] || drp[i] > max[i])
			return 0;

	return 1;
}

static inline void
compute_bmlm(uint64_t *bm, uint64_t *lm,
	     uint64_t dr, uint64_t min, uint64_t max)
{
	register uint64_t lm0 = m0, lm1 = m1;

	for (; lm0; lm0 >>= 1, lm1 >>= 1)
	{
		register uint64_t a = dr & lm0,
			 b = min & lm0,
			 c = max & lm0;

		if (b) { // ? 1 ?
			if (!a && c) { // 0 1 1
				*bm = min;
				break;
			}

		} else if (a) // 1 0 ?
			if (c) { // 1 0 1
				*lm = qload_0(max, lm0, lm1);
				min = qload_1(min, lm0, lm1);
			} else { // 1 0 0
				*lm = max;
				break;
			}

		else if (c) { // 0 0 1
			// BIGMIN
			*bm = qload_1(min, lm0, lm1);
			max = qload_0(max, lm0, lm1);
		}
	}
}

static inline unsigned
geo_search(geo_curi_t *curi, unsigned pdb_hd,
		int16_t *s, uint16_t *l, uint8_t dim)
{
	uint64_t rmin = morton_set(s, dim),
		 rmax, code, idx;
	int16_t e[dim], p[dim];
	const void *key, *value;
	unsigned cur, n = 0;
	geo_curi_t *ci;

	point_add(e, s, (int16_t *) l, dim);
	rmax = morton_set(e, dim);
	cur = qmap_iter(pdb_hd, &rmin, QM_RANGE);

next:	if (!qmap_next(&key, &value, cur))
		return n;

	code = * (uint64_t *) key;

	if (code < rmin || code > rmax) {
#if COMPUTE_BMLM
		compute_bmlm(&rmin, &rmax, code, rmin, rmax);
#endif
		goto next;
	}

	morton_get(p, code, dim);

	if (!inrange_p(p, s, e, dim))
		goto next;

	idx = point_idx(p, s, e, dim);
	ci = &curi[idx];	

	point_copy(ci->p, p, dim);
	ci->ref = * (unsigned *) value;
	n++;

	goto next;
}

unsigned
geo_iter(unsigned pdb_hd, int16_t *s, uint16_t *l, uint8_t dim)
{
	unsigned cur = idm_new(&geo_idm);
	uint32_t m = point_vol((int16_t *) l, dim);
	geo_cur_t *c = &geo_cursors[cur];

	c->items = malloc(sizeof(geo_curi_t) * m);

	for (unsigned i = 0; i < m; i++)
		c->items[i].ref = QM_MISS;

	geo_search(c->items, pdb_hd, s, l, dim);

	c->m = m;
	c->dim = dim;
	c->pos = 0;
	return cur;
}

int
geo_next(int16_t *p, unsigned *ref, unsigned cur)
{
	geo_cur_t *c = &geo_cursors[cur];
	geo_curi_t *ci;

next:	ci = &c->items[c->pos];	

	if (ci->ref != QM_MISS) {
		point_copy(p, ci->p, c->dim);
		*ref = ci->ref;
		c->pos++;
		return 1;
	}

	c->pos++;

	if (c->pos >= c->m)
		goto out;

	goto next;

out:	free(c->items);
	idm_del(&geo_idm, cur);
	return 0;
}

static int
morton_cmp(const void * const va,
		const void * const vb,
		size_t len UNUSED)
{
	uint64_t a = * (uint64_t *) va,
		 b = * (uint64_t *) vb;

	return b > a ? -1 : (a > b ? 1 : 0);
}

void
geo_init(void) {
	qm_u = qmap_reg(sizeof(unsigned));
	qm_u64 = qmap_reg(sizeof(uint64_t));
	qmap_cmp_set(qm_u64, morton_cmp);
	geo_idm = idm_init();
}


unsigned
geo_open(char *database, unsigned mask) {
	return qmap_open(qm_u64, qm_u, mask, 0);
}
