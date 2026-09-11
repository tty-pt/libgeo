/* see http://www.vision-tools.com/h-tropf/multidimensionalrangequery.pdf
 */

/* Ask morton.h to name its static inline versions *_il so this TU can
 * also emit the external ABI symbols without conflicting. Must be
 * defined before any header include. */
#define ISLET_MORTON_RENAME_FOR_WRAPPERS

#include "../include/ttypt/islet.h"
#include "../include/ttypt/point.h"
#include "../include/ttypt/morton.h"
#include "../include/ttypt/pointcfg.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include <ttypt/qsys.h>
#include <ttypt/idm.h>

#define MAX_DIM 4

/* Cursor item slots are 8 bytes: they hold either an int16_t[4]
 * point (configs 1..4) or an int32_t[2] point (config 2_32). The
 * per-cursor copy function moves exactly the width of the config
 * that created the cursor, so one pool and one idm serve both. */
typedef struct {
	int32_t p[2];
	uint32_t ref;
} islet_curi_t;

typedef struct {
	islet_curi_t *items;
	uint32_t n, pos;
	void (*copy)(void *, void *);
} islet_cur_t;

static uint32_t qm_u, qm_u64;

static idm_t islet_idm;

islet_cur_t islet_cursors[1024];

/* Extern ABI wrappers: consumers that link -lislet call these.
 * The header's static inline versions (renamed *_il above) are used
 * for all internal calls. */
#undef morton_set_1
#undef morton_set_2
#undef morton_set_3
#undef morton_set_4
#undef morton_set_2_32
#undef morton_get_1
#undef morton_get_2
#undef morton_get_3
#undef morton_get_4
#undef morton_get_2_32

uint64_t
morton_set_1(int16_t *p)
{
	return morton_set_1_il(p);
}

uint64_t
morton_set_2(int16_t *p)
{
	return morton_set_2_il(p);
}

uint64_t
morton_set_3(int16_t *p)
{
	return morton_set_3_il(p);
}

uint64_t
morton_set_4(int16_t *p)
{
	return morton_set_4_il(p);
}

void
morton_get_1(int16_t *pos, uint64_t code)
{
	morton_get_1_il(pos, code);
}

void
morton_get_2(int16_t *pos, uint64_t code)
{
	morton_get_2_il(pos, code);
}

void
morton_get_3(int16_t *pos, uint64_t code)
{
	morton_get_3_il(pos, code);
}

void
morton_get_4(int16_t *pos, uint64_t code)
{
	morton_get_4_il(pos, code);
}

uint64_t
morton_set_2_32(int32_t *p)
{
	return morton_set_2_32_il(p);
}

void
morton_get_2_32(int32_t *pos, uint64_t code)
{
	morton_get_2_32_il(pos, code);
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

	/* Unreachable: the sole callers (the islet_box_walk_N walk loops)
	 * pass only literal dims 1..4. An invalid dim matches nothing. */
	return 0;
}

/* 2D x 32-bit inrange. Same contract as inrange_p on int32_t lanes;
 * the sole caller (islet_box_walk_2_32) always passes 2 lanes. */
static inline int
inrange_p32(int32_t *drp, int32_t *min, int32_t *max)
{
	return drp[0] >= min[0] && drp[0] <= max[0]
		&& drp[1] >= min[1] && drp[1] <= max[1];
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
} islet_box_t;

/* Forced inline into the islet_box_walk_N monomorphizations below so
 * the literal dimension count reaches this body: the maxd loop then
 * unrolls and the dim==3/4/else chain inside the k-loop folds to a
 * single path. */
static inline __attribute__((always_inline)) uint64_t
islet_jump_over_gap(uint64_t code, int16_t *p,
		const islet_box_t *ub, uint8_t dim)
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

/* 2D x 32-bit gap jump. Same soundness argument as islet_jump_over_gap
 * (an aligned 2^k cube in unsigned-lane space occupies the contiguous
 * morton interval [base, base + 2^(D*k)), here D = 2): disjoint cubes
 * are provably empty of matches, and k = 0 always applies, so the
 * walk strictly progresses.
 *
 * Differences from the 16-bit path are mechanical: the 2^31 lane
 * bias, and a starting guess that never shifts past the lane width.
 * `33 - clz(maxd)` agrees with `31 - clz(maxd << 2)` exactly for
 * maxd < 2^30; larger distances clamp to kmax = 31, which keeps
 * `1u << k` defined and `dim*k` = 62 inside 64 bits. */
static inline __attribute__((always_inline)) uint64_t
islet_jump_over_gap32(uint64_t code, int32_t *p, const islet_box_t *ub)
{
	uint32_t maxd = 0;
	int kmax;

	for (uint8_t d = 0; d < 2; d++) {
		uint32_t up = (uint32_t)p[d] + 0x80000000u;
		uint32_t dist = up < ub->lo[d] ? ub->lo[d] - up
			: up > ub->hi[d] ? up - ub->hi[d] : 0;

		if (dist > maxd)
			maxd = dist;
	}

	if (maxd < 4)
		return code + 1;

	kmax = maxd >= ((uint32_t)1 << 30) ? 31 : 33 - __builtin_clz(maxd);
	if (kmax > 31)
		kmax = 31;

	for (int k = kmax; k >= 0; k--) {
		uint32_t side = 1u << k;
		uint32_t mask = side - 1;
		int disjoint = 0;

		for (uint8_t d = 0; d < 2; d++) {
			uint32_t up = (uint32_t)p[d] + 0x80000000u;
			uint32_t lo = up & ~mask;
			uint32_t hi = lo + side - 1;

			if (lo > ub->hi[d] || hi < ub->lo[d]) {
				disjoint = 1;
				break;
			}
		}

		if (disjoint) {
			/* 2*k address bits for an aligned side-2^k cube in
			 * 2 dimensions (k <= 31, so the shift stays within
			 * 64 bits). */
			uint64_t span = (1ULL << (2 * k)) - 1;
			uint64_t end = code | span;

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
 * containing the key (islet_jump_over_gap) and skips later keys below the
 * floor without decoding them. Duplicates below the floor are safe to
 * skip: chains are contiguous in sorted order, so a whole chain shares
 * one code and one verdict.
 */
/* The visit callback takes the decoded point as void *: walkers are
 * stamped per lane type (int16_t/int32_t) and each visitor casts it
 * back to its own lane type. */
typedef int (*islet_visit_fn)(void *p, uint32_t ref, void *ud);

/* Diagnostic: index entries fully examined (decoded) by the most recent
 * box walk. Tests prove the Z-interval skip engages (decoded well below
 * the morton-interval width on dense boxes). */
static uint32_t islet_scan_count;

uint32_t
islet_last_scan_count(void)
{
	return islet_scan_count;
}

/* Per-config box walkers. ISLET_BOX_WALK_CFG stamps out islet_box_walk_NAME
 * with the lane type, length type, unsigned-lane bias, codec and point
 * ops as compile-time parameters, so every op in the hot loop folds to
 * its exact path with no runtime dim dispatch. NAME selects the
 * matching inrange/jump pair: configs 1..4 use the dim-taking int16
 * functions (folded by the literal D), config 2_32 uses its own
 * lane-typed, dim-free pair. The int16 instantiations below emit the
 * same expressions as the former ISLET_BOX_WALK(N) macro; see the 2_32
 * instantiation for the 32-bit-lane config. */
#define ISLET_INR_1(p, s, e, D)    inrange_p(p, s, e, D)
#define ISLET_INR_2(p, s, e, D)    inrange_p(p, s, e, D)
#define ISLET_INR_3(p, s, e, D)    inrange_p(p, s, e, D)
#define ISLET_INR_4(p, s, e, D)    inrange_p(p, s, e, D)
#define ISLET_INR_2_32(p, s, e, D) inrange_p32(p, s, e)
#define ISLET_JUMP_1(c, p, ub, D)    islet_jump_over_gap(c, p, ub, D)
#define ISLET_JUMP_2(c, p, ub, D)    islet_jump_over_gap(c, p, ub, D)
#define ISLET_JUMP_3(c, p, ub, D)    islet_jump_over_gap(c, p, ub, D)
#define ISLET_JUMP_4(c, p, ub, D)    islet_jump_over_gap(c, p, ub, D)
#define ISLET_JUMP_2_32(c, p, ub, D) islet_jump_over_gap32(c, p, ub)
#define ISLET_BOX_WALK_CFG(NAME, D, PT, LT, BIAS, MSET, PADD, MGET) \
static uint32_t \
islet_box_walk_##NAME(uint32_t pdb_hd, PT *s, LT *l, \
		islet_visit_fn visit, void *ud) \
{ \
	uint64_t rmin, rmax, floor, code; \
	PT e[MAX_DIM], p[MAX_DIM]; \
	const void *key, *value; \
	uint32_t cur, n = 0; \
	islet_box_t ub; \
 \
	rmin = MSET(s); \
	PADD(e, s, (PT *) l); \
	rmax = MSET(e); \
 \
	for (uint8_t d = 0; d < D; d++) { \
		ub.lo[d] = (uint32_t)((int64_t)s[d] + (int64_t)BIAS); \
		ub.hi[d] = ub.lo[d] + (uint32_t)l[d]; \
	} \
 \
	/* Single ordered pass. floor ratchets past proven-empty address \
	 * spans; keys below it are false positives by construction and are \
	 * skipped without decoding. No cursor is ever reopened. */ \
	islet_scan_count = 0; \
	floor = rmin; \
	cur = qmap_iter(pdb_hd, &rmin, QM_RANGE | QM_RANGE_GE); \
 \
	while (qmap_next(&key, &value, cur)) { \
		code = * (uint64_t *) key; \
 \
		if (code > rmax) \
			break; \
 \
		if (code < floor) \
			continue; \
 \
		islet_scan_count++; \
		MGET(p, code); \
 \
		if (!ISLET_INR_##NAME(p, s, e, D)) { \
			/* Past the last possible key: nothing left to jump to. */ \
			if (code == UINT64_MAX) \
				break; \
 \
			floor = ISLET_JUMP_##NAME(code, p, &ub, D); \
			continue; \
		} \
 \
		n++; \
 \
		if (visit(p, * (uint32_t *) value, ud)) \
			break; \
	} \
 \
	qmap_fin(cur); \
	return n; \
}

ISLET_BOX_WALK_CFG(1, 1, int16_t, uint16_t, 32768,
	morton_set_1_il, point_add_1, morton_get_1_il)
ISLET_BOX_WALK_CFG(2, 2, int16_t, uint16_t, 32768,
	morton_set_2_il, point_add_2, morton_get_2_il)
ISLET_BOX_WALK_CFG(3, 3, int16_t, uint16_t, 32768,
	morton_set_3_il, point_add_3, morton_get_3_il)
ISLET_BOX_WALK_CFG(4, 4, int16_t, uint16_t, 32768,
	morton_set_4_il, point_add_4, morton_get_4_il)
ISLET_BOX_WALK_CFG(2_32, 2, int32_t, int32_t, 0x80000000u,
	morton_set_2_32_il, point_add_2_32, morton_get_2_32_il)
#undef ISLET_BOX_WALK_CFG
#undef ISLET_INR_1
#undef ISLET_INR_2
#undef ISLET_INR_3
#undef ISLET_INR_4
#undef ISLET_INR_2_32
#undef ISLET_JUMP_1
#undef ISLET_JUMP_2
#undef ISLET_JUMP_3
#undef ISLET_JUMP_4
#undef ISLET_JUMP_2_32

typedef struct {
	islet_curi_t *items;
	uint32_t n, cap;
	void (*copy)(void *, void *);
} islet_collect_t;

static int
islet_collect_visit(void *vp, uint32_t ref, void *ud)
{
	int16_t *p = vp;
	islet_collect_t *c = ud;

	if (c->n == c->cap) {
		uint32_t ncap = c->cap ? c->cap * 2 : 64;
		islet_curi_t *ni = realloc(c->items, ncap * sizeof *ni);

		if (!ni)
			return 1;

		c->items = ni;
		c->cap = ncap;
	}

	c->copy(c->items[c->n].p, p);
	c->items[c->n].ref = ref;
	c->n++;
	return 0;
}

/* 2D x 32-bit collector. Same shape as islet_collect_visit on
 * int32_t lanes; the shared islet_cursors[] pool and idm serve both. */
static int
islet_collect_visit32(void *vp, uint32_t ref, void *ud)
{
	int32_t *p = vp;
	islet_collect_t *c = ud;

	if (c->n == c->cap) {
		uint32_t ncap = c->cap ? c->cap * 2 : 64;
		islet_curi_t *ni = realloc(c->items, ncap * sizeof *ni);

		if (!ni)
			return 1;

		c->items = ni;
		c->cap = ncap;
	}

	c->copy(c->items[c->n].p, p);
	c->items[c->n].ref = ref;
	c->n++;
	return 0;
}

/* Per-config iterators. ISLET_ITER_CFG stamps islet_iter_NAME with the
 * point/length types, box walker, copy op and collector of that
 * config; one cursor pool and idm handle both lane widths. */
#define ISLET_ITER_CFG(NAME, PT, LT, WALK, COPY, CVIS) \
uint32_t \
islet_iter_##NAME(uint32_t pdb_hd, PT *s, LT *l) \
{ \
	uint32_t cur = idm_new(&islet_idm); \
	islet_cur_t *c = &islet_cursors[cur]; \
	islet_collect_t col = { NULL, 0, 0, (void (*)(void *, void *))COPY }; \
 \
	WALK(pdb_hd, s, l, CVIS, &col); \
 \
	c->items = col.items; \
	c->n = col.n; \
	c->copy = (void (*)(void *, void *))COPY; \
	c->pos = 0; \
	return cur; \
}

ISLET_ITER_CFG(1, int16_t, uint16_t,
	islet_box_walk_1, point_copy_1, islet_collect_visit)
ISLET_ITER_CFG(2, int16_t, uint16_t,
	islet_box_walk_2, point_copy_2, islet_collect_visit)
ISLET_ITER_CFG(3, int16_t, uint16_t,
	islet_box_walk_3, point_copy_3, islet_collect_visit)
ISLET_ITER_CFG(4, int16_t, uint16_t,
	islet_box_walk_4, point_copy_4, islet_collect_visit)
ISLET_ITER_CFG(2_32, int32_t, int32_t,
	islet_box_walk_2_32, point_copy_2_32, islet_collect_visit32)
#undef ISLET_ITER_CFG

static int
islet_next_impl(void *p, uint32_t *ref, uint32_t cur)
{
	islet_cur_t *c = &islet_cursors[cur];

	if (c->pos >= c->n) {
		free(c->items);
		c->items = NULL;
		idm_del(&islet_idm, cur);
		return 0;
	}

	c->copy(p, c->items[c->pos].p);
	*ref = c->items[c->pos].ref;
	c->pos++;
	return 1;
}

int
islet_next(int16_t *p, uint32_t *ref, uint32_t cur)
{
	return islet_next_impl(p, ref, cur);
}

int
islet_next32(int32_t *p, uint32_t *ref, uint32_t cur)
{
	return islet_next_impl(p, ref, cur);
}

/* Per-cell chain cursors for islet_get_multi (indexed by idm handle). */
static uint32_t islet_mcursors[1024];

/* Per-config cell-chain cursors. The qmap chain handles are
 * lane-agnostic, so one pool serves every config; islet_cell_next
 * stays the shared value-only advance for all of them. */
#define ISLET_GET_MULTI_CFG(NAME, PT, MSET) \
uint32_t \
islet_get_multi_##NAME(uint32_t pdb_hd, PT *p) \
{ \
	uint64_t code = MSET(p); \
	uint32_t qcur = qmap_get_multi(pdb_hd, &code); \
	uint32_t cur; \
 \
	if (qcur == QM_MISS) \
		return QM_MISS; \
 \
	cur = idm_new(&islet_idm); \
	islet_mcursors[cur] = qcur; \
	return cur; \
}

ISLET_GET_MULTI_CFG(1, int16_t, morton_set_1_il)
ISLET_GET_MULTI_CFG(2, int16_t, morton_set_2_il)
ISLET_GET_MULTI_CFG(3, int16_t, morton_set_3_il)
ISLET_GET_MULTI_CFG(4, int16_t, morton_set_4_il)
ISLET_GET_MULTI_CFG(2_32, int32_t, morton_set_2_32_il)
#undef ISLET_GET_MULTI_CFG

int
islet_cell_next(uint32_t *ref, uint32_t cur)
{
	const void *key, *value;

	if (!qmap_next(&key, &value, islet_mcursors[cur])) {
		qmap_fin(islet_mcursors[cur]);
		idm_del(&islet_idm, cur);
		return 0;
	}

	*ref = * (uint32_t *) value;
	return 1;
}

static int
islet_fill_visit(void *vp, uint32_t ref, void *ud)
{
	rec_set_t *out = ud;

	(void) vp;
	rec_set_push(out, (rec_ref_t) ref);
	return 0;
}

/* Per-config box fills. ISLET_FILL_CFG stamps rec_axis_fill_bbox_NAME
 * with the point/length types and box walker; the volume product is
 * uint64 so wide lanes cannot overflow it, and ISLET_FILL_MAX_VOL caps
 * every config identically. */
#define ISLET_FILL_CFG(NAME, D, PT, LT, WALK) \
int \
rec_axis_fill_bbox_##NAME(uint32_t pdb_hd, PT *s, \
		LT *l, rec_set_t *out) \
{ \
	uint64_t v = 1; \
 \
	if (!out) \
		return -1; \
 \
	for (uint8_t i = 0; i < D; i++) { \
		v *= (uint64_t)l[i]; \
 \
		if (v > ISLET_FILL_MAX_VOL) \
			return -1; \
	} \
 \
	WALK(pdb_hd, s, l, islet_fill_visit, out); \
	rec_set_seal(out); \
	return 0; \
}

ISLET_FILL_CFG(1, 1, int16_t, uint16_t, islet_box_walk_1)
ISLET_FILL_CFG(2, 2, int16_t, uint16_t, islet_box_walk_2)
ISLET_FILL_CFG(3, 3, int16_t, uint16_t, islet_box_walk_3)
ISLET_FILL_CFG(4, 4, int16_t, uint16_t, islet_box_walk_4)
ISLET_FILL_CFG(2_32, 2, int32_t, int32_t, islet_box_walk_2_32)
#undef ISLET_FILL_CFG

/* The per-dimension operation table: islet_ops[N] points at the
 * N-dimensional implementations (no dim argument on the calls; the
 * index is the dim). islet_ops[0] is all NULL. The walker never goes
 * through this table — it calls the islet_box_walk_N loops directly. */
const islet_ops_t islet_ops[5] = {
	[0] = { NULL },
	[1] = { morton_set_1, morton_get_1, point_add_1, point_copy_1,
		islet_put_1, islet_get_1, islet_set_1, islet_del_1,
		islet_del_all_1, islet_cell_count_1,
		islet_iter_1, islet_get_multi_1, rec_axis_fill_bbox_1 },
	[2] = { morton_set_2, morton_get_2, point_add_2, point_copy_2,
		islet_put_2, islet_get_2, islet_set_2, islet_del_2,
		islet_del_all_2, islet_cell_count_2,
		islet_iter_2, islet_get_multi_2, rec_axis_fill_bbox_2 },
	[3] = { morton_set_3, morton_get_3, point_add_3, point_copy_3,
		islet_put_3, islet_get_3, islet_set_3, islet_del_3,
		islet_del_all_3, islet_cell_count_3,
		islet_iter_3, islet_get_multi_3, rec_axis_fill_bbox_3 },
	[4] = { morton_set_4, morton_get_4, point_add_4, point_copy_4,
		islet_put_4, islet_get_4, islet_set_4, islet_del_4,
		islet_del_all_4, islet_cell_count_4,
		islet_iter_4, islet_get_multi_4, rec_axis_fill_bbox_4 },
};

/* =====================================================================
 * Config objects (pointcfg.h): the public "static-method" interface.
 * Each config instance is backed by tiny functions mirroring the flat
 * inline families (point_*_N / islet_*_N / morton_set_N / ...) so the
 * members are addressable. The header inlines stay the inlinable fast
 * path for tight loops; these exist so Point<D>_<B>.member is one
 * direct call through an exported object.
 * ===================================================================== */

/* Vector / point-utility backers for the int16 configs (dim D literal).
 * Bodies mirror point_add_N/sub_N/min_N/max_N/copy_N/set_N/debug_N/
 * vol_N/idx_N exactly. */
#define ISLET_P2B_VEC(D) \
static void isletc_p2_add_##D(int16_t *tar, int16_t *a, int16_t *b) \
{ \
	for (int i = 0; i < D; i++) \
		tar[i] = (int16_t)(a[i] + b[i]); \
} \
static void isletc_p2_sub_##D(int16_t *tar, int16_t *a, int16_t *b) \
{ \
	for (int i = 0; i < D; i++) \
		tar[i] = (int16_t)(a[i] - b[i]); \
} \
static void isletc_p2_min_##D(int16_t *tar, int16_t *a, int16_t *b) \
{ \
	for (int i = 0; i < D; i++) \
		tar[i] = a[i] < b[i] ? a[i] : b[i]; \
} \
static void isletc_p2_max_##D(int16_t *tar, int16_t *a, int16_t *b) \
{ \
	for (int i = 0; i < D; i++) \
		tar[i] = a[i] > b[i] ? a[i] : b[i]; \
} \
static void isletc_p2_copy_##D(int16_t *tar, int16_t *src) \
{ \
	for (int i = 0; i < D; i++) \
		tar[i] = src[i]; \
} \
static void isletc_p2_set_##D(int16_t *tar, int16_t v) \
{ \
	for (int i = 0; i < D; i++) \
		tar[i] = v; \
} \
static void isletc_p2_debug_##D(char *label, int16_t *p) \
{ \
	fprintf(stderr, "%s(", label); \
	for (int i = 0; i < D; i++) \
		fprintf(stderr, "%s%d", i ? ", " : "", p[i]); \
	fprintf(stderr, ")\n"); \
} \
static int32_t isletc_p2_vol_##D(int16_t *p) \
{ \
	int32_t acc = 1; \
	for (int i = 0; i < D; i++) \
		acc *= (int32_t)p[i]; \
	return acc; \
} \
static uint64_t isletc_p2_idx_##D(int16_t *p, int16_t *s, int16_t *e) \
{ \
	uint64_t acc = 0, stride = 1; \
	for (int i = 0; i < D; i++) { \
		acc += (uint64_t)(p[i] - s[i]) * stride; \
		stride *= (uint64_t)(e[i] - s[i]); \
	} \
	return acc; \
}

/* Database-op backers for the int16 configs. Bodies mirror the
 * islet_put_N/get_N/set_N/del_N/del_all_N/cell_count_N inlines. */
#define ISLET_P2B_DB(D) \
static void isletc_p2_put_##D(uint32_t db, int16_t *p, uint32_t ref) \
{ \
	uint64_t c = morton_set_##D(p); \
	qmap_put(db, &c, &ref); \
} \
static uint32_t isletc_p2_get_##D(uint32_t db, int16_t *p) \
{ \
	uint64_t c = morton_set_##D(p); \
	const void *v = qmap_get(db, &c); \
	return v ? *(uint32_t *)v : ISLET_MISS; \
} \
static void isletc_p2_replace_##D(uint32_t db, int16_t *p, uint32_t ref) \
{ \
	uint64_t c = morton_set_##D(p); \
	qmap_del_all(db, &c); \
	qmap_put(db, &c, &ref); \
} \
static void isletc_p2_del_##D(uint32_t db, int16_t *p) \
{ \
	uint64_t c = morton_set_##D(p); \
	qmap_del(db, &c); \
} \
static uint32_t isletc_p2_del_all_##D(uint32_t db, int16_t *p) \
{ \
	uint64_t c = morton_set_##D(p); \
	uint32_t n = qmap_count(db, &c); \
	if (n) \
		qmap_del_all(db, &c); \
	return n; \
} \
static uint32_t isletc_p2_cell_count_##D(uint32_t db, int16_t *p) \
{ \
	uint64_t c = morton_set_##D(p); \
	return qmap_count(db, &c); \
}

ISLET_P2B_VEC(1)
ISLET_P2B_VEC(2)
ISLET_P2B_VEC(3)
ISLET_P2B_VEC(4)
ISLET_P2B_DB(1)
ISLET_P2B_DB(2)
ISLET_P2B_DB(3)
ISLET_P2B_DB(4)

#undef ISLET_P2B_VEC
#undef ISLET_P2B_DB

/* 2D x 32-bit backers: int32 lanes, both vector and DB ops. */
#define ISLET_P4B_BACKERS \
static void isletc_p4_add(int32_t *tar, int32_t *a, int32_t *b) \
{ \
	tar[0] = a[0] + b[0]; \
	tar[1] = a[1] + b[1]; \
} \
static void isletc_p4_sub(int32_t *tar, int32_t *a, int32_t *b) \
{ \
	tar[0] = a[0] - b[0]; \
	tar[1] = a[1] - b[1]; \
} \
static void isletc_p4_min(int32_t *tar, int32_t *a, int32_t *b) \
{ \
	tar[0] = a[0] < b[0] ? a[0] : b[0]; \
	tar[1] = a[1] < b[1] ? a[1] : b[1]; \
} \
static void isletc_p4_max(int32_t *tar, int32_t *a, int32_t *b) \
{ \
	tar[0] = a[0] > b[0] ? a[0] : b[0]; \
	tar[1] = a[1] > b[1] ? a[1] : b[1]; \
} \
static void isletc_p4_copy(int32_t *tar, int32_t *src) \
{ \
	*(int64_t *)tar = *(int64_t *)src; \
} \
static void isletc_p4_set(int32_t *tar, int32_t v) \
{ \
	tar[0] = v; \
	tar[1] = v; \
} \
static void isletc_p4_debug(char *label, int32_t *p) \
{ \
	fprintf(stderr, "%s(%d, %d)\n", label, p[0], p[1]); \
} \
static uint64_t isletc_p4_vol(int32_t *p) \
{ \
	return (uint64_t)(uint32_t)p[0] * (uint64_t)(uint32_t)p[1]; \
} \
static uint64_t isletc_p4_idx(int32_t *p, int32_t *s, int32_t *e) \
{ \
	return (uint64_t)(p[0] - s[0]) \
		+ (uint64_t)(p[1] - s[1]) * (uint64_t)(e[0] - s[0]); \
} \
static void isletc_p4_put(uint32_t db, int32_t *p, uint32_t ref) \
{ \
	uint64_t c = morton_set_2_32(p); \
	qmap_put(db, &c, &ref); \
} \
static uint32_t isletc_p4_get(uint32_t db, int32_t *p) \
{ \
	uint64_t c = morton_set_2_32(p); \
	const void *v = qmap_get(db, &c); \
	return v ? *(uint32_t *)v : ISLET_MISS; \
} \
static void isletc_p4_replace(uint32_t db, int32_t *p, uint32_t ref) \
{ \
	uint64_t c = morton_set_2_32(p); \
	qmap_del_all(db, &c); \
	qmap_put(db, &c, &ref); \
} \
static void isletc_p4_del(uint32_t db, int32_t *p) \
{ \
	uint64_t c = morton_set_2_32(p); \
	qmap_del(db, &c); \
} \
static uint32_t isletc_p4_del_all(uint32_t db, int32_t *p) \
{ \
	uint64_t c = morton_set_2_32(p); \
	uint32_t n = qmap_count(db, &c); \
	if (n) \
		qmap_del_all(db, &c); \
	return n; \
} \
static uint32_t isletc_p4_cell_count(uint32_t db, int32_t *p) \
{ \
	uint64_t c = morton_set_2_32(p); \
	return qmap_count(db, &c); \
}

ISLET_P4B_BACKERS

#undef ISLET_P4B_BACKERS

/* The public config objects. Point1_2..Point4_2 share the int16 type;
 * Point2_4 is the int32 2D config. iter/get_multi/fill_bbox/next reuse
 * the exported islet_iter_*, islet_get_multi_*, rec_axis_fill_bbox_* and
 * islet_next / islet_next32 / islet_cell_next symbols directly. */
const islet_point2b_t Point1_2 = {
	.morton_set = morton_set_1,
	.morton_get = morton_get_1,
	.add = isletc_p2_add_1,
	.sub = isletc_p2_sub_1,
	.min = isletc_p2_min_1,
	.max = isletc_p2_max_1,
	.copy = isletc_p2_copy_1,
	.vol = isletc_p2_vol_1,
	.set = isletc_p2_set_1,
	.debug = isletc_p2_debug_1,
	.idx = isletc_p2_idx_1,
	.put = isletc_p2_put_1,
	.get = isletc_p2_get_1,
	.replace = isletc_p2_replace_1,
	.del = isletc_p2_del_1,
	.del_all = isletc_p2_del_all_1,
	.cell_count = isletc_p2_cell_count_1,
	.get_multi = islet_get_multi_1,
	.iter = islet_iter_1,
	.fill_bbox = rec_axis_fill_bbox_1,
	.next = islet_next,
};

const islet_point2b_t Point2_2 = {
	.morton_set = morton_set_2,
	.morton_get = morton_get_2,
	.add = isletc_p2_add_2,
	.sub = isletc_p2_sub_2,
	.min = isletc_p2_min_2,
	.max = isletc_p2_max_2,
	.copy = isletc_p2_copy_2,
	.vol = isletc_p2_vol_2,
	.set = isletc_p2_set_2,
	.debug = isletc_p2_debug_2,
	.idx = isletc_p2_idx_2,
	.put = isletc_p2_put_2,
	.get = isletc_p2_get_2,
	.replace = isletc_p2_replace_2,
	.del = isletc_p2_del_2,
	.del_all = isletc_p2_del_all_2,
	.cell_count = isletc_p2_cell_count_2,
	.get_multi = islet_get_multi_2,
	.iter = islet_iter_2,
	.fill_bbox = rec_axis_fill_bbox_2,
	.next = islet_next,
};

const islet_point2b_t Point3_2 = {
	.morton_set = morton_set_3,
	.morton_get = morton_get_3,
	.add = isletc_p2_add_3,
	.sub = isletc_p2_sub_3,
	.min = isletc_p2_min_3,
	.max = isletc_p2_max_3,
	.copy = isletc_p2_copy_3,
	.vol = isletc_p2_vol_3,
	.set = isletc_p2_set_3,
	.debug = isletc_p2_debug_3,
	.idx = isletc_p2_idx_3,
	.put = isletc_p2_put_3,
	.get = isletc_p2_get_3,
	.replace = isletc_p2_replace_3,
	.del = isletc_p2_del_3,
	.del_all = isletc_p2_del_all_3,
	.cell_count = isletc_p2_cell_count_3,
	.get_multi = islet_get_multi_3,
	.iter = islet_iter_3,
	.fill_bbox = rec_axis_fill_bbox_3,
	.next = islet_next,
};

const islet_point2b_t Point4_2 = {
	.morton_set = morton_set_4,
	.morton_get = morton_get_4,
	.add = isletc_p2_add_4,
	.sub = isletc_p2_sub_4,
	.min = isletc_p2_min_4,
	.max = isletc_p2_max_4,
	.copy = isletc_p2_copy_4,
	.vol = isletc_p2_vol_4,
	.set = isletc_p2_set_4,
	.debug = isletc_p2_debug_4,
	.idx = isletc_p2_idx_4,
	.put = isletc_p2_put_4,
	.get = isletc_p2_get_4,
	.replace = isletc_p2_replace_4,
	.del = isletc_p2_del_4,
	.del_all = isletc_p2_del_all_4,
	.cell_count = isletc_p2_cell_count_4,
	.get_multi = islet_get_multi_4,
	.iter = islet_iter_4,
	.fill_bbox = rec_axis_fill_bbox_4,
	.next = islet_next,
};

const islet_point4b_t Point2_4 = {
	.morton_set = morton_set_2_32,
	.morton_get = morton_get_2_32,
	.add = isletc_p4_add,
	.sub = isletc_p4_sub,
	.min = isletc_p4_min,
	.max = isletc_p4_max,
	.copy = isletc_p4_copy,
	.vol = isletc_p4_vol,
	.set = isletc_p4_set,
	.debug = isletc_p4_debug,
	.idx = isletc_p4_idx,
	.put = isletc_p4_put,
	.get = isletc_p4_get,
	.replace = isletc_p4_replace,
	.del = isletc_p4_del,
	.del_all = isletc_p4_del_all,
	.cell_count = isletc_p4_cell_count,
	.get_multi = islet_get_multi_2_32,
	.iter = islet_iter_2_32,
	.fill_bbox = rec_axis_fill_bbox_2_32,
	.next = islet_next32,
};

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
islet_init(void) {
	static int inited = 0;

	if (inited)
		return;
	inited = 1;
	qm_u = qmap_reg(sizeof(uint32_t));
	qm_u64 = qmap_reg(sizeof(uint64_t));
	qmap_cmp_set(qm_u64, morton_cmp);
	islet_idm = idm_init();
}


uint32_t
islet_open(char *filename, char *database, uint32_t mask) {
	return qmap_open(filename, database, qm_u64, qm_u, mask,
			QM_SORTED | QM_MULTIVALUE);
}

#if ISLET_SIMD_MORTON

#include <string.h>

#ifdef __AVX2__
#include <immintrin.h>

static inline unsigned long
islet_axis_u16(int16_t v)
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
		__m256i ux = _mm256_set_epi64x(islet_axis_u16(points[i+3][0]),
				islet_axis_u16(points[i+2][0]),
				islet_axis_u16(points[i+1][0]),
				islet_axis_u16(points[i+0][0]));
		__m256i uy = _mm256_set_epi64x(islet_axis_u16(points[i+3][1]),
				islet_axis_u16(points[i+2][1]),
				islet_axis_u16(points[i+1][1]),
				islet_axis_u16(points[i+0][1]));
		__m256i uz = _mm256_set_epi64x(islet_axis_u16(points[i+3][2]),
				islet_axis_u16(points[i+2][2]),
				islet_axis_u16(points[i+1][2]),
				islet_axis_u16(points[i+0][2]));

		__m256i sx, sy, sz;
		morton_spread3_4x(ux, uy, uz, &sx, &sy, &sz);

		__m256i codes = _mm256_or_si256(
				_mm256_or_si256(sx, sy), sz);
		_mm256_storeu_si256((__m256i *)(out + i), codes);
	}

	/* scalar tail */
	for (; i < n; i++)
		out[i] = morton_set_3_il(points[i]);

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
		__m256i ux = _mm256_set_epi64x(islet_axis_u16(points[i+3][0]),
				islet_axis_u16(points[i+2][0]),
				islet_axis_u16(points[i+1][0]),
				islet_axis_u16(points[i+0][0]));
		__m256i uy = _mm256_set_epi64x(islet_axis_u16(points[i+3][1]),
				islet_axis_u16(points[i+2][1]),
				islet_axis_u16(points[i+1][1]),
				islet_axis_u16(points[i+0][1]));
		__m256i uz = _mm256_set_epi64x(islet_axis_u16(points[i+3][2]),
				islet_axis_u16(points[i+2][2]),
				islet_axis_u16(points[i+1][2]),
				islet_axis_u16(points[i+0][2]));
		__m256i uw = _mm256_set_epi64x(islet_axis_u16(points[i+3][3]),
				islet_axis_u16(points[i+2][3]),
				islet_axis_u16(points[i+1][3]),
				islet_axis_u16(points[i+0][3]));

		__m256i sx, sy, sz, sw;
		morton_spread4_4x(ux, uy, uz, uw, &sx, &sy, &sz, &sw);

		__m256i codes = _mm256_or_si256(
				_mm256_or_si256(sx, sy),
				_mm256_or_si256(sz, sw));
		_mm256_storeu_si256((__m256i *)(out + i), codes);
	}

	/* scalar tail */
	for (; i < n; i++)
		out[i] = morton_set_4_il(points[i]);

	return i;
}

/* Bulk decode: compact axis from 4 codes simultaneously. Each lane of
 * 'codes' holds a 64-bit Morton code.  Returns __m256i with the
 * compacted unsigned coordinate in the low 16 bits of each lane. */
static inline __m256i
morton_compact_axis_4x(__m256i codes, uint32_t shift)
{
	__m256i v;

	if (shift == 1)
		v = _mm256_srli_epi64(codes, 1);
	else if (shift == 2)
		v = _mm256_srli_epi64(codes, 2);
	else
		v = codes;

	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x1249249249249249ULL));
	v = _mm256_xor_si256(v, _mm256_srli_epi64(v, 2));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x10C30C30C30C30C3ULL));
	v = _mm256_xor_si256(v, _mm256_srli_epi64(v, 4));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x100F00F00F00F00FULL));
	v = _mm256_xor_si256(v, _mm256_srli_epi64(v, 8));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x1F0000FF0000FFULL));
	v = _mm256_xor_si256(v, _mm256_srli_epi64(v, 16));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x00000000001FFFFFULL));

	return v;
}

uint32_t
morton_get_bulk(int16_t points[][3], const uint64_t *codes, uint32_t n)
{
	uint32_t i = 0;
	__m256i bias = _mm256_set1_epi16(SHRT_MAX + 1);

	/* process 4 points at a time with AVX2 */
	for (; i + 4 <= n; i += 4) {
		__m256i c = _mm256_loadu_si256((const __m256i *)(codes + i));
		__m256i sx = morton_compact_axis_4x(c, 0);
		__m256i sy = morton_compact_axis_4x(c, 1);
		__m256i sz = morton_compact_axis_4x(c, 2);

		sx = _mm256_sub_epi16(sx, bias);
		sy = _mm256_sub_epi16(sy, bias);
		sz = _mm256_sub_epi16(sz, bias);

		/* Each axis vector holds one value per 64-bit lane (low 16
		 * bits only), so as a uint16_t[16] the 4 point values sit
		 * at stride 4 (indices 0,4,8,12). Scatter to the
		 * interleaved output via temp aligned storage — the
		 * compact is the expensive part this parallelizes. */
		__attribute__((aligned(32))) uint16_t tx[16], ty[16], tz[16];
		_mm256_store_si256((__m256i *)tx, sx);
		_mm256_store_si256((__m256i *)ty, sy);
		_mm256_store_si256((__m256i *)tz, sz);
		for (int k = 0; k < 4; k++) {
			points[i + k][0] = (int16_t)tx[k * 4];
			points[i + k][1] = (int16_t)ty[k * 4];
			points[i + k][2] = (int16_t)tz[k * 4];
		}
	}

	/* scalar tail */
	for (; i < n; i++)
		morton_get_3_il(points[i], codes[i]);

	return i;
}

/* 4D compact_axis4 4-wide: collect every 4th bit starting at 'shift',
 * one code per 64-bit lane, low 16 bits of each lane hold the result. */
static inline __m256i
morton_compact_axis4_4x(__m256i codes, uint32_t shift)
{
	__m256i v;

	if (shift == 1)
		v = _mm256_srli_epi64(codes, 1);
	else if (shift == 2)
		v = _mm256_srli_epi64(codes, 2);
	else if (shift == 3)
		v = _mm256_srli_epi64(codes, 3);
	else
		v = codes;

	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x1111111111111111ULL));
	v = _mm256_xor_si256(v, _mm256_srli_epi64(v, 3));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x0303030303030303ULL));
	v = _mm256_xor_si256(v, _mm256_srli_epi64(v, 6));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x000F000F000F000FULL));
	v = _mm256_xor_si256(v, _mm256_srli_epi64(v, 12));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x000000FF000000FFULL));
	v = _mm256_xor_si256(v, _mm256_srli_epi64(v, 24));
	v = _mm256_and_si256(v,
		_mm256_set1_epi64x(0x000000000000FFFFULL));

	return v;
}

uint32_t
morton_get_bulk4(int16_t points[][4], const uint64_t *codes, uint32_t n)
{
	uint32_t i = 0;
	__m256i bias = _mm256_set1_epi16(SHRT_MAX + 1);

	/* process 4 points at a time with AVX2 */
	for (; i + 4 <= n; i += 4) {
		__m256i c = _mm256_loadu_si256((const __m256i *)(codes + i));
		__m256i sx = morton_compact_axis4_4x(c, 0);
		__m256i sy = morton_compact_axis4_4x(c, 1);
		__m256i sz = morton_compact_axis4_4x(c, 2);
		__m256i sw = morton_compact_axis4_4x(c, 3);

		sx = _mm256_sub_epi16(sx, bias);
		sy = _mm256_sub_epi16(sy, bias);
		sz = _mm256_sub_epi16(sz, bias);
		sw = _mm256_sub_epi16(sw, bias);

		/* Same stride-4 layout as morton_get_bulk(); scatter via
		 * temp aligned storage. */
		__attribute__((aligned(32))) uint16_t tx[16], ty[16], tz[16], tw[16];
		_mm256_store_si256((__m256i *)tx, sx);
		_mm256_store_si256((__m256i *)ty, sy);
		_mm256_store_si256((__m256i *)tz, sz);
		_mm256_store_si256((__m256i *)tw, sw);
		for (int k = 0; k < 4; k++) {
			points[i + k][0] = (int16_t)tx[k * 4];
			points[i + k][1] = (int16_t)ty[k * 4];
			points[i + k][2] = (int16_t)tz[k * 4];
			points[i + k][3] = (int16_t)tw[k * 4];
		}
	}

	/* scalar tail */
	for (; i < n; i++)
		morton_get_4_il(points[i], codes[i]);

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
		out[i] = morton_set_3_il(points[i]);

	return i;
}

uint32_t
morton_set_bulk4(uint64_t *out, int16_t points[][4], uint32_t n)
{
	uint32_t i = 0;

	/* scalar on ARM NEON — proper NEON spread4 TBD */
	for (; i < n; i++)
		out[i] = morton_set_4_il(points[i]);

	return i;
}

uint32_t
morton_get_bulk(int16_t points[][3], const uint64_t *codes, uint32_t n)
{
	for (uint32_t i = 0; i < n; i++)
		morton_get_3_il(points[i], codes[i]);
	return n;
}

uint32_t
morton_get_bulk4(int16_t points[][4], const uint64_t *codes, uint32_t n)
{
	for (uint32_t i = 0; i < n; i++)
		morton_get_4_il(points[i], codes[i]);
	return n;
}

#else /* no SIMD */

uint32_t
morton_set_bulk(uint64_t *out, int16_t points[][3], uint32_t n)
{
	for (uint32_t i = 0; i < n; i++)
		out[i] = morton_set_3_il(points[i]);
	return n;
}

uint32_t
morton_set_bulk4(uint64_t *out, int16_t points[][4], uint32_t n)
{
	for (uint32_t i = 0; i < n; i++)
		out[i] = morton_set_4_il(points[i]);
	return n;
}

uint32_t
morton_get_bulk(int16_t points[][3], const uint64_t *codes, uint32_t n)
{
	for (uint32_t i = 0; i < n; i++)
		morton_get_3_il(points[i], codes[i]);
	return n;
}

uint32_t
morton_get_bulk4(int16_t points[][4], const uint64_t *codes, uint32_t n)
{
	for (uint32_t i = 0; i < n; i++)
		morton_get_4_il(points[i], codes[i]);
	return n;
}

#endif /* __AVX2__ / __ARM_NEON / scalar */

#endif /* ISLET_SIMD_MORTON */

/* ---- rec_query axis registration (islet / space) ---- */

struct rec_islet_params {
	int16_t  s[4];
	uint16_t l[4];
	int      dim;
};

static int islet_fill(void *ctx, void *params, rec_set_t *out)
{
	uint32_t pdb_hd = (uint32_t)(uintptr_t)ctx;
	struct rec_islet_params *p = params;

	if (!p || p->dim < 1 || p->dim > 4)
		return -1;
	return islet_ops[p->dim].fill(pdb_hd, p->s, p->l, out);
}

/*
 * Decode "dim=N s=x,y,... l=dx,dy,..." into a heap-owned rec_islet_params
 * (freed never — one-shot CLI process lifetime, matches the other axis
 * decode fns). `dim` (1..4) is required; `s` and `l` must each list
 * exactly `dim` comma-separated integers. NULL on malformed/missing
 * fields or OOM.
 */
static void *islet_decode(const char *str)
{
	struct rec_islet_params *p;
	char *buf, *cur;
	const char *sval = NULL, *lval = NULL;
	int dim = 0;

	if (!str)
		return NULL;
	buf = malloc(strlen(str) + 1);
	if (!buf)
		return NULL;
	strcpy(buf, str);
	cur = buf;
	while (*cur) {
		char *key, *val;

		while (*cur == ' ')
			cur++;
		if (!*cur)
			break;
		key = cur;
		while (*cur && *cur != '=' && *cur != ' ')
			cur++;
		if (*cur != '=') {
			if (*cur)
				cur++;
			continue;
		}
		*cur++ = '\0';
		val = cur;
		while (*cur && *cur != ' ')
			cur++;
		if (*cur)
			*cur++ = '\0';
		if (!strcmp(key, "dim"))
			dim = atoi(val);
		else if (!strcmp(key, "s"))
			sval = val;
		else if (!strcmp(key, "l"))
			lval = val;
	}
	if (dim < 1 || dim > 4 || !sval || !lval) {
		free(buf);
		return NULL;
	}

	p = calloc(1, sizeof(*p));
	if (!p) {
		free(buf);
		return NULL;
	}
	p->dim = dim;

	{
		const char *c = sval;
		int i;

		for (i = 0; i < dim; i++) {
			if (!*c) {
				free(p);
				free(buf);
				return NULL;
			}
			p->s[i] = (int16_t)strtol(c, (char **)&c, 10);
			if (i < dim - 1) {
				if (*c != ',') {
					free(p);
					free(buf);
					return NULL;
				}
				c++;
			}
		}
	}
	{
		const char *c = lval;
		int i;

		for (i = 0; i < dim; i++) {
			if (!*c) {
				free(p);
				free(buf);
				return NULL;
			}
			p->l[i] = (uint16_t)strtol(c, (char **)&c, 10);
			if (i < dim - 1) {
				if (*c != ',') {
					free(p);
					free(buf);
					return NULL;
				}
				c++;
			}
		}
	}

	free(buf);
	return p;
}

__attribute__((constructor)) static void islet_rec_axis_init(void)
{
	static const rec_axis_t islet_axis = {
		"islet", islet_fill, NULL, NULL, islet_decode
	};

	islet_init();
	rec_axis_register(&islet_axis);
}

/*
 * rec_axis_open convention (PLAN-REC-QUERY.md §4.3, optional CLI-open
 * convention, not part of libqmap's core rec_query registry API): spec
 * is "filename:database:mask" (`:`-separated, any/all fields may be
 * empty for islet_open()'s NULL/0 defaults). Returns the uint32_t db
 * handle widened to a pointer via uintptr_t, same cast the constructor's
 * own rec_axis_set_ctx() callers already use.
 */
void *rec_axis_open(const char *spec)
{
	char *buf, *cur, *fname, *dbname, *maskstr;
	uint32_t mask;
	uint32_t db;

	if (!spec)
		spec = "";
	buf = strdup(spec);
	if (!buf)
		return NULL;

	cur = buf;
	fname = cur;
	cur = strchr(cur, ':');
	if (cur)
		*cur++ = '\0';
	else
		cur = buf + strlen(buf);
	dbname = cur;
	cur = strchr(cur, ':');
	if (cur)
		*cur++ = '\0';
	maskstr = cur;

	mask = (maskstr && *maskstr) ? (uint32_t)strtoul(maskstr, NULL, 10) : 0;
	db = islet_open(*fname ? fname : NULL, *dbname ? dbname : NULL, mask);

	free(buf);
	return (void *)(uintptr_t)db;
}
