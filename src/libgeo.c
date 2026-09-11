/* see http://www.vision-tools.com/h-tropf/multidimensionalrangequery.pdf
 */

/* Ask morton.h to name its static inline versions *_il so this TU can
 * also emit the external ABI symbols without conflicting. Must be
 * defined before any header include. */
#define GEO_MORTON_RENAME_FOR_WRAPPERS

#ifndef FAST_MORTON
#define FAST_MORTON 1
#endif

#include "../include/ttypt/geo.h"
#include "../include/ttypt/point.h"
#include "../include/ttypt/morton.h"

#include <limits.h>
#include <stdlib.h>

#include <ttypt/qsys.h>
#include <ttypt/idm.h>

#define MAX_DIM 4

typedef struct {
	int16_t p[MAX_DIM];
	uint32_t ref;
} geo_curi_t;

typedef struct {
	geo_curi_t *items;
	uint32_t n, pos;
	uint8_t dim;
} geo_cur_t;

static uint32_t qm_u, qm_u64;

static idm_t geo_idm;

geo_cur_t geo_cursors[1024];

/* Extern ABI wrappers: consumers that link -lgeo call these.
 * The header's static inline versions (renamed *_il above) are used
 * for all internal calls. */
#undef morton_set
#undef morton_get

uint64_t
morton_set(int16_t *p, uint8_t dim)
{
	return morton_set_il(p, dim);
}

void
morton_get(int16_t *pos, uint64_t code, uint8_t dim)
{
	morton_get_il(pos, code, dim);
}


static inline int
inrange_p(int16_t *drp, int16_t *min, int16_t *max, uint8_t dim)
{
	if (dim == 1)
		return drp[0] >= min[0] && drp[0] <= max[0];
	if (dim == 2)
		return drp[0] >= min[0] && drp[0] <= max[0]
			&& drp[1] >= min[1] && drp[1] <= max[1];
	if (dim == 3)
		return drp[0] >= min[0] && drp[0] <= max[0]
			&& drp[1] >= min[1] && drp[1] <= max[1]
			&& drp[2] >= min[2] && drp[2] <= max[2];
	if (dim == 4)
		return drp[0] >= min[0] && drp[0] <= max[0]
			&& drp[1] >= min[1] && drp[1] <= max[1]
			&& drp[2] >= min[2] && drp[2] <= max[2]
			&& drp[3] >= min[3] && drp[3] <= max[3];

	/* Unreachable: the sole caller (geo_box_visit) admits only
	 * dims 1..MAX_DIM. An invalid dim matches nothing. */
	return 0;
}

/* Largest forward jump past provably out-of-box Z-space.
 *
 * Soundness proof: an aligned 2^k cube in unsigned-coordinate space
 * occupies one contiguous morton interval [base, base + 2^(D*k)), where
 * D is the dimension count (8^k for 3D, 16^k for 4D). When that cube is
 * disjoint from the query box (in any queried dimension), no address in
 * its interval can decode to an in-box point — a decoded point inside
 * the interval shares the cube's coordinate high bits in every queried
 * lane, hence lies inside the cube, hence outside the box. Every stored
 * key in [code, nlb) is therefore a false positive the walker would
 * discard anyway. k = 0 (the point itself, already known out-of-box)
 * always applies, so the walk strictly progresses.
 */

typedef struct {
	uint32_t lo[MAX_DIM];
	uint32_t hi[MAX_DIM];
} geo_box_t;

static uint64_t
geo_jump_over_gap(uint64_t code, int16_t *p,
		const geo_box_t *ub, uint8_t dim)
{
	uint32_t maxd = 0;
	int kmax;

	for (uint8_t d = 0; d < dim; d++) {
		uint32_t up = (uint32_t)((int32_t)p[d] + 32768);
		uint32_t dist = up < ub->lo[d] ? ub->lo[d] - up
			: up > ub->hi[d] ? up - ub->hi[d] : 0;

		if (dist > maxd)
			maxd = dist;
	}

	if (maxd < 4)
		return code + 1;

	kmax = 31 - __builtin_clz(maxd << 2);
	if (kmax > 15) kmax = 15;

	for (int k = kmax; k >= 0; k--) {
		uint32_t side = 1u << k;
		uint32_t mask = side - 1;
		int disjoint = 0;

		if (dim == 3) {
			uint32_t up0 = (uint32_t)((int32_t)p[0] + 32768);
			uint32_t up1 = (uint32_t)((int32_t)p[1] + 32768);
			uint32_t up2 = (uint32_t)((int32_t)p[2] + 32768);
			disjoint = ((up0 & ~mask) > ub->hi[0]
				|| (up0 | mask) < ub->lo[0]
				|| (up1 & ~mask) > ub->hi[1]
				|| (up1 | mask) < ub->lo[1]
				|| (up2 & ~mask) > ub->hi[2]
				|| (up2 | mask) < ub->lo[2]);
		} else if (dim == 4) {
			uint32_t up0 = (uint32_t)((int32_t)p[0] + 32768);
			uint32_t up1 = (uint32_t)((int32_t)p[1] + 32768);
			uint32_t up2 = (uint32_t)((int32_t)p[2] + 32768);
			uint32_t up3 = (uint32_t)((int32_t)p[3] + 32768);
			disjoint = ((up0 & ~mask) > ub->hi[0]
				|| (up0 | mask) < ub->lo[0]
				|| (up1 & ~mask) > ub->hi[1]
				|| (up1 | mask) < ub->lo[1]
				|| (up2 & ~mask) > ub->hi[2]
				|| (up2 | mask) < ub->lo[2]
				|| (up3 & ~mask) > ub->hi[3]
				|| (up3 | mask) < ub->lo[3]);
		} else
		{
			for (uint8_t d = 0; d < dim; d++) {
				uint32_t up = (uint32_t)((int32_t)p[d] + 32768);
				uint32_t lo = up & ~mask;
				uint32_t hi = lo + side - 1;

				if (lo > ub->hi[d] || hi < ub->lo[d]) {
					disjoint = 1;
					break;
				}
			}
		}

		if (disjoint) {
			/* D*k address bits for an aligned side-2^k cube in
			 * D dimensions (dim <= MAX_DIM <= 4, k <= 15, so the
			 * shift stays within 64 bits). */
			uint64_t span = (1ULL << (dim * k)) - 1;
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
	geo_box_t ub;

	if (dim == 0 || dim > MAX_DIM)
		return 0;

	rmin = morton_set(s, dim);
	point_add(e, s, (int16_t *) l, dim);
	rmax = morton_set(e, dim);

	for (uint8_t d = 0; d < dim; d++) {
		ub.lo[d] = (uint32_t)((int32_t)s[d] + 32768);
		ub.hi[d] = ub.lo[d] + l[d];
	}

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

			floor = geo_jump_over_gap(code, p, &ub, dim);
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

#if GEO_SIMD_MORTON

#include <string.h>

#ifdef __AVX2__
#include <immintrin.h>

static inline unsigned long
geo_axis_u16(int16_t v)
{
	return (unsigned long)((uint16_t)(v + SHRT_MAX + 1));
}

static inline void
morton_spread3_4x(__m256i ux, __m256i uy, __m256i uz,
		__m256i *out_x, __m256i *out_y, __m256i *out_z)
{
	__m256i v;

	v = _mm256_and_si256(ux,
		_mm256_set1_epi64x(0x000000000000FFFFULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 32));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x001F00000000FFFFULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 16));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x001F0000FF0000FFULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 8));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x100F00F00F00F00FULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 4));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x10C30C30C30C30C3ULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 2));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x1249249249249249ULL));
	*out_x = v;

	v = _mm256_and_si256(uy,
		_mm256_set1_epi64x(0x000000000000FFFFULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 32));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x001F00000000FFFFULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 16));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x001F0000FF0000FFULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 8));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x100F00F00F00F00FULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 4));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x10C30C30C30C30C3ULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 2));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x1249249249249249ULL));
	*out_y = _mm256_slli_epi64(v, 1);

	v = _mm256_and_si256(uz,
		_mm256_set1_epi64x(0x000000000000FFFFULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 32));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x001F00000000FFFFULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 16));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x001F0000FF0000FFULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 8));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x100F00F00F00F00FULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 4));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x10C30C30C30C30C3ULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 2));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x1249249249249249ULL));
	*out_z = _mm256_slli_epi64(v, 2);
}

uint32_t
morton_set_bulk(uint64_t *out, int16_t points[][3], uint32_t n)
{
	uint32_t i = 0;

	/* process 4 points at a time with AVX2 */
	for (; i + 4 <= n; i += 4) {
		/* Row-major (x,y,z)-triples: gather each axis into a
		 * 4-lane int64 vector with the unsigned coordinate offset
		 * applied. Lane j holds point i+j's value for that axis. */
		__m256i ux = _mm256_set_epi64x(geo_axis_u16(points[i+3][0]),
				geo_axis_u16(points[i+2][0]),
				geo_axis_u16(points[i+1][0]),
				geo_axis_u16(points[i+0][0]));
		__m256i uy = _mm256_set_epi64x(geo_axis_u16(points[i+3][1]),
				geo_axis_u16(points[i+2][1]),
				geo_axis_u16(points[i+1][1]),
				geo_axis_u16(points[i+0][1]));
		__m256i uz = _mm256_set_epi64x(geo_axis_u16(points[i+3][2]),
				geo_axis_u16(points[i+2][2]),
				geo_axis_u16(points[i+1][2]),
				geo_axis_u16(points[i+0][2]));

		__m256i sx, sy, sz;
		morton_spread3_4x(ux, uy, uz, &sx, &sy, &sz);

		__m256i codes = _mm256_or_si256(
				_mm256_or_si256(sx, sy), sz);
		_mm256_storeu_si256((__m256i *)(out + i), codes);
	}

	/* scalar tail */
	for (; i < n; i++)
		out[i] = morton_set(points[i], 3);

	return i;
}

static inline void
morton_spread4_4x(__m256i ux, __m256i uy, __m256i uz, __m256i uw,
		__m256i *out_x, __m256i *out_y,
		__m256i *out_z, __m256i *out_w)
{
	__m256i v;

	v = _mm256_and_si256(ux,
		_mm256_set1_epi64x(0x000000000000FFFFULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 24));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x000000FF000000FFULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 12));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x000F000F000F000FULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 6));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x0303030303030303ULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 3));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x1111111111111111ULL));
	*out_x = v;

	v = _mm256_and_si256(uy,
		_mm256_set1_epi64x(0x000000000000FFFFULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 24));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x000000FF000000FFULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 12));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x000F000F000F000FULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 6));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x0303030303030303ULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 3));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x1111111111111111ULL));
	*out_y = _mm256_slli_epi64(v, 1);

	v = _mm256_and_si256(uz,
		_mm256_set1_epi64x(0x000000000000FFFFULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 24));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x000000FF000000FFULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 12));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x000F000F000F000FULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 6));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x0303030303030303ULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 3));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x1111111111111111ULL));
	*out_z = _mm256_slli_epi64(v, 2);

	v = _mm256_and_si256(uw,
		_mm256_set1_epi64x(0x000000000000FFFFULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 24));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x000000FF000000FFULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 12));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x000F000F000F000FULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 6));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x0303030303030303ULL));
	v = _mm256_or_si256(v, _mm256_slli_epi64(v, 3));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x1111111111111111ULL));
	*out_w = _mm256_slli_epi64(v, 3);
}

uint32_t
morton_set_bulk4(uint64_t *out, int16_t points[][4], uint32_t n)
{
	uint32_t i = 0;

	/* process 4 points at a time with AVX2 */
	for (; i + 4 <= n; i += 4) {
		/* Row-major (x,y,z,w)-quads: gather each axis into a
		 * 4-lane int64 vector with the unsigned coordinate offset
		 * applied. Lane j holds point i+j's value for that axis. */
		__m256i ux = _mm256_set_epi64x(geo_axis_u16(points[i+3][0]),
				geo_axis_u16(points[i+2][0]),
				geo_axis_u16(points[i+1][0]),
				geo_axis_u16(points[i+0][0]));
		__m256i uy = _mm256_set_epi64x(geo_axis_u16(points[i+3][1]),
				geo_axis_u16(points[i+2][1]),
				geo_axis_u16(points[i+1][1]),
				geo_axis_u16(points[i+0][1]));
		__m256i uz = _mm256_set_epi64x(geo_axis_u16(points[i+3][2]),
				geo_axis_u16(points[i+2][2]),
				geo_axis_u16(points[i+1][2]),
				geo_axis_u16(points[i+0][2]));
		__m256i uw = _mm256_set_epi64x(geo_axis_u16(points[i+3][3]),
				geo_axis_u16(points[i+2][3]),
				geo_axis_u16(points[i+1][3]),
				geo_axis_u16(points[i+0][3]));

		__m256i sx, sy, sz, sw;
		morton_spread4_4x(ux, uy, uz, uw, &sx, &sy, &sz, &sw);

		__m256i codes = _mm256_or_si256(
				_mm256_or_si256(sx, sy),
				_mm256_or_si256(sz, sw));
		_mm256_storeu_si256((__m256i *)(out + i), codes);
	}

	/* scalar tail */
	for (; i < n; i++)
		out[i] = morton_set(points[i], 4);

	return i;
}

#elif defined(__ARM_NEON)
#include <arm_neon.h>

uint32_t
morton_set_bulk(uint64_t *out, int16_t points[][3], uint32_t n)
{
	uint32_t i = 0;

	/* scalar on ARM NEON — proper NEON spread3 TBD */
	for (; i < n; i++)
		out[i] = morton_set(points[i], 3);

	return i;
}

uint32_t
morton_set_bulk4(uint64_t *out, int16_t points[][4], uint32_t n)
{
	uint32_t i = 0;

	/* scalar on ARM NEON — proper NEON spread4 TBD */
	for (; i < n; i++)
		out[i] = morton_set(points[i], 4);

	return i;
}

#else /* no SIMD */

uint32_t
morton_set_bulk(uint64_t *out, int16_t points[][3], uint32_t n)
{
	for (uint32_t i = 0; i < n; i++)
		out[i] = morton_set(points[i], 3);
	return n;
}

uint32_t
morton_set_bulk4(uint64_t *out, int16_t points[][4], uint32_t n)
{
	for (uint32_t i = 0; i < n; i++)
		out[i] = morton_set(points[i], 4);
	return n;
}

#endif /* __AVX2__ / __ARM_NEON / scalar */

#endif /* GEO_SIMD_MORTON */
