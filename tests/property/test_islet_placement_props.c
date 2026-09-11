/**
 * @file test_islet_placement_props.c
 * @brief Property tests: every placed (point, value) comes back exactly at
 *        its own coordinate, and nothing else does.
 *
 * True oracle (independent of the walker and of morton round-trips): the
 * test keeps its OWN model of every islet_put/islet_del_all/islet_set call.
 * For each random box the model is filtered geometrically (pure coordinate
 * arithmetic, no morton) to build the expected pair multiset; the islet_iter
 * walk must match it EXACTLY — same points, same values, same multiplicity.
 * Single-cell walks and islet_get_multi chains pin the in-cell ordering.
 * Deterministic PRNG; edit churn exercises the dirty-rebuild path.
 */

#include "../test_common.h"
#include "../../include/ttypt/islet.h"
#include "../../include/ttypt/point.h"
#include "../../include/ttypt/morton.h"
#include <stdlib.h>
#include <string.h>

#define PLACE_CLOUDS 12
#define PLACE_BOXES 30
#define PLACE_POINTS 150
#define MODEL_CAP (PLACE_POINTS * 2 + 16)

typedef struct {
	int16_t p[3];
	uint32_t ref;
} pair_t;

static uint8_t g_pdim;

static int cmp_pair(const void *va, const void *vb)
{
	const pair_t *a = va, *b = vb;
	uint64_t ca = islet_ops[g_pdim].morton_set((int16_t *)(void *)a->p);
	uint64_t cb = islet_ops[g_pdim].morton_set((int16_t *)(void *)b->p);

	if (ca != cb)
		return ca < cb ? -1 : 1;
	if (a->ref != b->ref)
		return a->ref < b->ref ? -1 : 1;
	return 0;
}

static int cell_eq(const pair_t *a, const pair_t *b, uint8_t dim)
{
	for (uint8_t i = 0; i < dim; i++)
		if (a->p[i] != b->p[i])
			return 0;
	return 1;
}

/* --- model bookkeeping (the ONLY record of placement) --- */
static void model_put(pair_t *m, size_t *nm, int16_t *p, uint32_t ref,
		uint8_t dim)
{
	m[*nm].ref = ref;
	for (uint8_t i = 0; i < dim; i++)
		m[*nm].p[i] = p[i];
	(*nm)++;
}

static void model_del_cell(pair_t *m, size_t *nm, int16_t *p, uint8_t dim)
{
	pair_t probe;
	uint8_t i;

	memset(&probe, 0, sizeof probe);
	for (i = 0; i < dim; i++)
		probe.p[i] = p[i];

	/* Order-preserving removal: the model is the ONLY record of insertion
	 * order, and assert_cell_matches_model() compares islet_get_multi chains
	 * against it. Swap-compaction would jump the last row into the removed
	 * slot and silently reorder every surviving row's array position, so
	 * the model would no longer reflect the true insertion order of the
	 * surviving puts (the DB chain tail-appends and is preserved by
	 * qmap_rebuild_map). */
	for (size_t j = 0; j < *nm;) {
		if (cell_eq(&m[j], &probe, dim)) {
			memmove(&m[j], &m[j + 1], (*nm - j - 1) * sizeof *m);
			(*nm)--;
		} else {
			j++;
		}
	}
}

static void model_set_cell(pair_t *m, size_t *nm, int16_t *p, uint32_t ref,
		uint8_t dim)
{
	model_del_cell(m, nm, p, dim);
	model_put(m, nm, p, ref, dim);
}

/* --- box walk collecting (point, value) pairs --- */
static size_t walk_collect(uint32_t db, int16_t *s, uint16_t *l, uint8_t dim,
		pair_t *out, size_t cap)
{
	int16_t e[4];
	islet_ops[dim].point_add(e, s, (int16_t *)l);
	uint32_t it = islet_ops[dim].iter(db, s, l);
	size_t n = 0;

	while (n < cap && islet_next(out[n].p, &out[n].ref, it)) {
		for (uint8_t i = 0; i < dim; i++)
			ASSERT(out[n].p[i] >= s[i] && out[n].p[i] <= e[i]);
		n++;
	}
	return n;
}

/* --- model-filtered expected pairs for a box (pure coordinate math) --- */
static size_t model_in_box(pair_t *m, size_t nm, int16_t *s, uint16_t *l,
		uint8_t dim, pair_t *out)
{
	int16_t e[4];
	islet_ops[dim].point_add(e, s, (int16_t *)l);
	size_t n = 0;

	for (size_t i = 0; i < nm; i++) {
		int inside = 1;
		for (uint8_t d = 0; d < dim; d++)
			if (m[i].p[d] < s[d] || m[i].p[d] > e[d]) {
				inside = 0;
				break;
			}
		if (inside)
			out[n++] = m[i];
	}
	return n;
}

static void assert_box_matches_model(uint32_t db, int16_t *s, uint16_t *l,
		uint8_t dim, pair_t *m, size_t nm)
{
	pair_t exp[PLACE_POINTS * 2 + 16];
	pair_t walk[PLACE_POINTS * 2 + 16];

	size_t ne = model_in_box(m, nm, s, l, dim, exp);
	size_t nw = walk_collect(db, s, l, dim, walk, PLACE_POINTS * 2 + 16);

	ASSERT_EQ(nw, ne);
	if (ne == 0)
		return;

	g_pdim = dim;
	pair_t *ea = malloc(ne * sizeof *ea);
	pair_t *wa = malloc(nw * sizeof *wa);
	memcpy(ea, exp, ne * sizeof *ea);
	memcpy(wa, walk, nw * sizeof *wa);
	qsort(ea, ne, sizeof *ea, cmp_pair);
	qsort(wa, nw, sizeof *wa, cmp_pair);
	for (size_t i = 0; i < ne; i++) {
		ASSERT_EQ(ea[i].ref, wa[i].ref);
		ASSERT_POINT_EQ(ea[i].p, wa[i].p, dim);
	}
	free(ea);
	free(wa);
}

/* Single cell: get_multi must yield the model entries in insertion order. */
static void assert_cell_matches_model(uint32_t db, int16_t *p, uint8_t dim,
		pair_t *m, size_t nm)
{
	/* Per-cell bound (1D packs ~200 rows into 49 cells; a cell can hold a
	 * dozen). Overflow asserts loudly instead of silently truncating nwant,
	 * which would compare the DB's true count against a clipped expected. */
#define CELL_CAP 64
	pair_t want[CELL_CAP];
	size_t nwant = 0;

	for (size_t i = 0; i < nm; i++) {
		int same = 1;
		for (uint8_t d = 0; d < dim; d++)
			if (m[i].p[d] != p[d]) {
				same = 0;
				break;
			}
		if (same) {
			ASSERT(nwant < CELL_CAP);
			want[nwant++] = m[i];
		}
	}

	ASSERT_EQ(islet_ops[dim].cell_count(db, p), nwant);

	uint32_t cur = islet_ops[dim].get_multi(db, p);
	uint32_t ref;

	if (nwant == 0) {
		ASSERT_EQ(cur, QM_MISS);
		return;
	}
	for (size_t i = 0; i < nwant; i++) {
		ASSERT_EQ(islet_cell_next(&ref, cur), 1);
		ASSERT_EQ(ref, want[i].ref);
	}
	ASSERT_EQ(islet_cell_next(&ref, cur), 0);

	/* single-cell box walk: exactly the model pairs, one per value */
	pair_t walk[CELL_CAP];
	uint16_t z[3] = {0, 0, 0};
	size_t nw = walk_collect(db, p, z, dim, walk, CELL_CAP);
	ASSERT(nw < CELL_CAP);
	ASSERT_EQ(nw, nwant);

	g_pdim = dim;
	pair_t ea[CELL_CAP];
	for (size_t i = 0; i < nwant; i++)
		ea[i] = want[i];
	qsort(ea, nwant, sizeof *ea, cmp_pair);
	qsort(walk, nw, sizeof *walk, cmp_pair);
	for (size_t i = 0; i < nwant; i++) {
		ASSERT_EQ(ea[i].ref, walk[i].ref);
		ASSERT_POINT_EQ(ea[i].p, walk[i].p, dim);
	}
#undef CELL_CAP
}

static void run_clouds(uint8_t dim, int span, int npoints, uint64_t seed,
		const char *tag)
{
	static char names[PLACE_CLOUDS][64];

	test_seed_rng(seed);
	for (int c = 0; c < PLACE_CLOUDS; c++) {
		pair_t model[MODEL_CAP];
		size_t nm = 0;
		int16_t p[4];
		uint32_t ref;

		snprintf(names[c], sizeof names[c], "plc_prop_%s_%d", tag, c);
		uint32_t db = islet_open(NULL, names[c], 4095);

		for (int i = 0; i < npoints; i++) {
			for (uint8_t d = 0; d < dim; d++)
				p[d] = test_rand_coord_range(0, span);
			ref = (uint32_t)(test_rand64() % 5000);
			islet_ops[dim].put(db, p, ref);
			model_put(model, &nm, p, ref, dim);
			/* every 4th point: a second value at the SAME cell */
			if (i % 4 == 0) {
				pair_t ccell = model[nm - 1];
				uint32_t ref2 = (uint32_t)(test_rand64() % 5000);
				islet_ops[dim].put(db, ccell.p, ref2);
				model_put(model, &nm, ccell.p, ref2, dim);
			}
		}

		/* edit churn on even clouds: del_all + set (dirty rebuild) */
		if (c % 2 == 0) {
			for (int i = 0; i < npoints / 6; i++) {
				for (uint8_t d = 0; d < dim; d++)
					p[d] = test_rand_coord_range(0, span);
				islet_ops[dim].del_all(db, p);
				model_del_cell(model, &nm, p, dim);
			}
			for (int i = 0; i < npoints / 8; i++) {
				for (uint8_t d = 0; d < dim; d++)
					p[d] = test_rand_coord_range(0, span);
				ref = (uint32_t)(test_rand64() % 5000);
				islet_ops[dim].set(db, p, ref);
				model_set_cell(model, &nm, p, ref, dim);
			}
		}

		for (int b = 0; b < PLACE_BOXES; b++) {
			int16_t s[4] = {0, 0, 0, 0};
			uint16_t l[4] = {1, 1, 1, 1};
			if (nm > 0 && b % 2 == 0) {
				/* box anchored at a real placed point */
				int16_t *sp = model[test_rand64() % nm].p;
				for (uint8_t d = 0; d < dim; d++) {
					s[d] = sp[d];
					l[d] = (uint16_t)(test_rand64() % 8);
				}
			} else {
				for (uint8_t d = 0; d < dim; d++) {
					s[d] = test_rand_coord_range(-4, span);
					l[d] = (uint16_t)(test_rand64() % (span / 3 + 1) + 1);
				}
			}
			assert_box_matches_model(db, s, l, dim, model, nm);
		}

		/* single-cell exactness at a handful of model cells */
		for (int i = 0; i < 4; i++) {
			if (nm == 0)
				break;
			pair_t pick = model[test_rand64() % nm];
			assert_cell_matches_model(db, pick.p, dim, model, nm);
		}
	}
}

TEST(placement_oracle_3d) {
	islet_init();
	run_clouds(3, 48, PLACE_POINTS, 4041, "3d");
}

TEST(placement_oracle_2d) {
	islet_init();
	run_clouds(2, 48, PLACE_POINTS, 4042, "2d");
}

TEST(placement_oracle_1d) {
	islet_init();
	run_clouds(1, 48, PLACE_POINTS, 4043, "1d");
}

int main(void)
{
	test_suite_begin("Islet Placement Oracle Property Tests");
	RUN_TEST(placement_oracle_3d);
	RUN_TEST(placement_oracle_2d);
	RUN_TEST(placement_oracle_1d);
	return test_suite_end();
}