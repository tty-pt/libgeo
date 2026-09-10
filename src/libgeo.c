/* see http://www.vision-tools.com/h-tropf/multidimensionalrangequery.pdf
 */

#include "../include/ttypt/geo.h"
#include "../include/ttypt/point.h"
#include "../include/ttypt/morton.h"

#include <limits.h>
#include <stdlib.h>

#include <ttypt/qsys.h>
#include <ttypt/idm.h>

typedef struct {
	int16_t p[4];
	uint32_t ref;
} geo_curi_t;

typedef struct {
	geo_curi_t *items;
	uint32_t n, pos;
	uint8_t dim;
} geo_cur_t;

#define MAX_DIM 3
#define FAST_MORTON 1

static uint32_t qm_u, qm_u64;

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

static inline uint64_t morton2_pack_u16(uint16_t x, uint16_t y)
{
	return spread3(x) | (spread3(y) << 1);
}

static inline uint64_t morton1_pack_u16(uint16_t x)
{
	return spread3(x);
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
	uint16_t up[3] = {0, 0, 0};

	for (uint8_t i = 0; i < dim && i < 3; i++)
		up[i] = unsign(p[i]);

#if FAST_MORTON
	switch (dim) {
	case 1:
		return morton1_pack_u16(up[0]);
	case 2:
		return morton2_pack_u16(up[0], up[1]);
	default:
		return morton3_pack_u16(up[0], up[1], up[2], 0);
	}
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
compact_axis(uint64_t code, uint32_t shift)
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


static inline int
inrange_p(int16_t *drp, int16_t *min, int16_t *max, uint8_t dim)
{
	for (uint8_t i = 0; i < dim; i++)
		if (drp[i] < min[i] || drp[i] > max[i])
			return 0;

	return 1;
}

/* Largest forward jump past provably out-of-box Z-space.
 *
 * Soundness proof: an aligned 2^k cube in unsigned-coordinate space
 * occupies one contiguous morton interval [base, base + 8^k). When that
 * cube is disjoint from the query box (in any queried dimension), no
 * address in its interval can decode to an in-box point — a decoded
 * point inside the interval shares the cube's coordinate high bits in
 * every queried lane, hence lies inside the cube, hence outside the
 * box. Every stored key in [code, nlb) is therefore a false positive
 * the walker would discard anyway. k = 0 (the point itself, already
 * known out-of-box) always applies, so the walk strictly progresses.
 */
static uint64_t
geo_jump_over_gap(uint64_t code, int16_t *p,
		int16_t *s, uint16_t *l, uint8_t dim)
{
	uint32_t maxd = 0;
	int kmax;

	/* Cell distance from p to the box (p is out-of-box, so maxd >= 1).
	 * Only cubes with side on the order of maxd can be disjoint; the
	 * slack covers alignment luck. Capping merely shrinks jumps. */
	for (uint8_t d = 0; d < dim; d++) {
		int32_t up = (int32_t)p[d] + 32768;
		int32_t bs = (int32_t)s[d] + 32768;
		int32_t be = bs + l[d];
		uint32_t dist = up < bs ? (uint32_t)(bs - up)
			: up > be ? (uint32_t)(up - be) : 0;

		if (dist > maxd)
			maxd = dist;
	}

	/* Adjacent false positives cannot hide a useful cube: a linear
	 * step is cheaper than the descent. Only far misses pay for it. */
	if (maxd < 4)
		return code + 1;

	kmax = 0;
	while (kmax < 15 && (1u << kmax) <= (maxd << 2))
		kmax++;

	for (int k = kmax; k >= 0; k--) {
		uint32_t side = 1u << k;
		uint32_t mask = side - 1;
		int disjoint = 0;

		for (uint8_t d = 0; d < dim; d++) {
			uint32_t up = (uint32_t)((int32_t)p[d] + 32768);
			uint32_t lo = up & ~mask;
			uint32_t hi = lo + side - 1;
			uint32_t bs = (uint32_t)((int32_t)s[d] + 32768);
			uint32_t be = bs + l[d];

			if (lo > be || hi < bs) {
				disjoint = 1;
				break;
			}
		}

		if (disjoint) {
			uint64_t span = (1ULL << (3 * k)) - 1;
			uint64_t end = code | span;

			/* Wrapped past the last address: no smaller cube
			 * can extend further, step down instead. k = 0
			 * (span 0) always terminates the descent. */
			if (end == UINT64_MAX)
				continue;

			return end + 1;
		}
	}

	/* Unreachable: k = 0 always finds the point itself disjoint. */
	return code + 1;
}

/* Shared box walker: visits every stored (point, value) pair whose point
 * lies in [s, s+l), in morton-discovery order. Multiple values sharing one
 * cell are visited as distinct entries. Returns the visit count.
 *
 * On a QM_MULTIVALUE map a plain QM_RANGE walk would only iterate the
 * duplicates of the starting key, so the walk REQUIRES QM_RANGE_GE.
 * The map is QM_SORTED by morton code, so the walk stops at the first
 * key past rmax.
 *
 * Z-interval skip: a stored key inside [rmin, rmax] but outside the box
 * proves its aligned neighborhood may be empty of matches, so the walk
 * ratchets a skip floor past the largest box-disjoint aligned cube
 * containing the key (geo_jump_over_gap) and skips later keys below the
 * floor without decoding them. Duplicates below the floor are safe to
 * skip: chains are contiguous in sorted order, so a whole chain shares
 * one code and one verdict.
 */
typedef int (*geo_visit_fn)(int16_t *p, uint32_t ref, void *ud);

/* Diagnostic: index entries fully examined (decoded) by the most recent
 * box walk. Tests prove the Z-interval skip engages (decoded well below
 * the morton-interval width on dense boxes). */
static uint32_t geo_scan_count;

uint32_t
geo_last_scan_count(void)
{
	return geo_scan_count;
}

static uint32_t
geo_box_visit(uint32_t pdb_hd, int16_t *s, uint16_t *l, uint8_t dim,
		geo_visit_fn visit, void *ud)
{
	uint64_t rmin, rmax, floor, code;
	int16_t e[4], p[4];
	const void *key, *value;
	uint32_t cur, n = 0;

	if (dim == 0 || dim > MAX_DIM)
		return 0;

	rmin = morton_set(s, dim);
	point_add(e, s, (int16_t *) l, dim);
	rmax = morton_set(e, dim);

	/* Single ordered pass. floor ratchets past proven-empty address
	 * spans; keys below it are false positives by construction and are
	 * skipped without decoding. No cursor is ever reopened. */
	geo_scan_count = 0;
	floor = rmin;
	cur = qmap_iter(pdb_hd, &rmin, QM_RANGE | QM_RANGE_GE);

	while (qmap_next(&key, &value, cur)) {
		code = * (uint64_t *) key;

		if (code > rmax)
			break;

		if (code < floor)
			continue;

		geo_scan_count++;
		morton_get(p, code, dim);

		if (!inrange_p(p, s, e, dim)) {
			/* Past the last possible key: nothing left to jump to. */
			if (code == UINT64_MAX)
				break;

			floor = geo_jump_over_gap(code, p, s, l, dim);
			continue;
		}

		n++;

		if (visit(p, * (uint32_t *) value, ud))
			break;
	}

	qmap_fin(cur);
	return n;
}

typedef struct {
	geo_curi_t *items;
	uint32_t n, cap;
	uint8_t dim;
} geo_collect_t;

static int
geo_collect_visit(int16_t *p, uint32_t ref, void *ud)
{
	geo_collect_t *c = ud;

	if (c->n == c->cap) {
		uint32_t ncap = c->cap ? c->cap * 2 : 64;
		geo_curi_t *ni = realloc(c->items, ncap * sizeof *ni);

		if (!ni)
			return 1;

		c->items = ni;
		c->cap = ncap;
	}

	point_copy(c->items[c->n].p, p, c->dim);
	c->items[c->n].ref = ref;
	c->n++;
	return 0;
}

uint32_t
geo_iter(uint32_t pdb_hd, int16_t *s, uint16_t *l, uint8_t dim)
{
	uint32_t cur = idm_new(&geo_idm);
	geo_cur_t *c = &geo_cursors[cur];
	geo_collect_t col = { NULL, 0, 0, dim };

	geo_box_visit(pdb_hd, s, l, dim, geo_collect_visit, &col);

	c->items = col.items;
	c->n = col.n;
	c->dim = dim;
	c->pos = 0;
	return cur;
}

int
geo_next(int16_t *p, uint32_t *ref, uint32_t cur)
{
	geo_cur_t *c = &geo_cursors[cur];

	if (c->pos >= c->n) {
		free(c->items);
		c->items = NULL;
		idm_del(&geo_idm, cur);
		return 0;
	}

	point_copy(p, c->items[c->pos].p, c->dim);
	*ref = c->items[c->pos].ref;
	c->pos++;
	return 1;
}

/* Per-cell chain cursors for geo_get_multi (indexed by idm handle). */
static uint32_t geo_mcursors[1024];

uint32_t
geo_get_multi(uint32_t pdb_hd, int16_t *p, uint8_t dim)
{
	uint64_t code = morton_set(p, dim);
	uint32_t qcur = qmap_get_multi(pdb_hd, &code);
	uint32_t cur;

	if (qcur == QM_MISS)
		return QM_MISS;

	cur = idm_new(&geo_idm);
	geo_mcursors[cur] = qcur;
	return cur;
}

int
geo_cell_next(uint32_t *ref, uint32_t cur)
{
	const void *key, *value;

	if (!qmap_next(&key, &value, geo_mcursors[cur])) {
		qmap_fin(geo_mcursors[cur]);
		idm_del(&geo_idm, cur);
		return 0;
	}

	*ref = * (uint32_t *) value;
	return 1;
}

static int
geo_fill_visit(int16_t *p, uint32_t ref, void *ud)
{
	rec_set_t *out = ud;

	(void) p;
	rec_set_push(out, (rec_ref_t) ref);
	return 0;
}

int
rec_axis_fill_bbox(uint32_t pdb_hd, int16_t *s,
		uint16_t *l, uint8_t dim, rec_set_t *out)
{
	uint64_t v = 1;

	if (!out || dim == 0 || dim > MAX_DIM)
		return -1;

	for (uint8_t i = 0; i < dim; i++) {
		v *= l[i];

		if (v > GEO_FILL_MAX_VOL)
			return -1;
	}

	geo_box_visit(pdb_hd, s, l, dim, geo_fill_visit, out);
	rec_set_seal(out);
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
	qm_u = qmap_reg(sizeof(uint32_t));
	qm_u64 = qmap_reg(sizeof(uint64_t));
	qmap_cmp_set(qm_u64, morton_cmp);
	geo_idm = idm_init();
}


uint32_t
geo_open(char *filename, char *database, uint32_t mask) {
	return qmap_open(filename, database, qm_u64, qm_u, mask,
			QM_SORTED | QM_MULTIVALUE);
}
