/*
 * Integration tests: persisted geometry is exact — every value that went to
 * disk at coordinate P reads back from P, and box walks after reopen return
 * exactly the placed (point, value) pairs (nothing relocated, nothing lost,
 * nothing invented).
 */

#include "../test_common.h"
#include "../../include/ttypt/geo.h"
#include "../../include/ttypt/point.h"
#include "../../include/ttypt/morton.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define PLC_F "/tmp/test_geo_placement_persist.db"

typedef struct {
	int16_t p[3];
	uint32_t ref;
} pair_t;

static int cmp_pair(const void *va, const void *vb)
{
	const pair_t *a = va, *b = vb;
	uint64_t ca = morton_set((int16_t *)(void *)a->p, 3);
	uint64_t cb = morton_set((int16_t *)(void *)b->p, 3);

	if (ca != cb)
		return ca < cb ? -1 : 1;
	return a->ref < b->ref ? -1 : (a->ref > b->ref ? 1 : 0);
}

static void setup_once(void)
{
	static int initialized = 0;

	if (!initialized) {
		geo_init();
		initialized = 1;
	}
}

/* Cloud + MV cells survive save/close/reopen with exact coordinates */
TEST(placement_persist_exact_cloud) {
	setup_once();
	unlink(PLC_F);

	uint32_t db = geo_open(PLC_F, "plc", 4095);

	pair_t placed[24];
	/* unique cells: latticed, negative + positive + extremes */
	int idx = 0;
	for (int i = -3; i < 3 && idx < 21; i++)
		for (int j = -2; j < 2 && idx < 21; j++) {
			placed[idx].p[0] = (int16_t)(i * 7 - 100);
			placed[idx].p[1] = (int16_t)(j * 11 + 50);
			placed[idx].p[2] = (int16_t)((i + j) * 5 + 32000);
			placed[idx].ref = 1024u + (uint32_t)idx;
			geo_put(db, placed[idx].p, placed[idx].ref, 3);
			idx++;
		}
	/* MV cell: three values at ONE coordinate (off the lattice) */
	int16_t mv[3] = {-100, 50, 32001};
	size_t nmv = 3;
	placed[idx++] = (pair_t){{-100, 50, 32001}, 7001};
	geo_put(db, mv, 7001, 3);
	placed[idx++] = (pair_t){{-100, 50, 32001}, 7002};
	geo_put(db, mv, 7002, 3);
	placed[idx++] = (pair_t){{-100, 50, 32001}, 7003};
	geo_put(db, mv, 7003, 3);

	ASSERT_EQ(idx, 24);
	qmap_save();
	qmap_close(db);

	db = geo_open(PLC_F, "plc", 4095);

	/* exact-point reads on every placed coordinate */
	for (int i = 0; i < 24; i++)
		if (placed[i].p[0] == -100 && placed[i].p[1] == 50 &&
				placed[i].p[2] == 32001)
			continue; /* MV: covered below */
		else
			ASSERT_EQ(geo_get(db, placed[i].p, 3), placed[i].ref);

	ASSERT_EQ(geo_cell_count(db, mv, 3), nmv);
	uint32_t cur = geo_get_multi(db, mv, 3);
	ASSERT(cur != QM_MISS);
	uint32_t ref;
	ASSERT_EQ(geo_cell_next(&ref, cur), 1);
	ASSERT_EQ(ref, 7001);
	ASSERT_EQ(geo_cell_next(&ref, cur), 1);
	ASSERT_EQ(ref, 7002);
	ASSERT_EQ(geo_cell_next(&ref, cur), 1);
	ASSERT_EQ(ref, 7003);
	ASSERT_EQ(geo_cell_next(&ref, cur), 0);

	/* enclosing box: the exact placed pair multiset (pure box filter) */
	int16_t s[3] = {-130, -70, 31900};
	int16_t e[3] = {-130 + 160, -70 + 120, 31900 + 140};
	uint16_t l[3] = {160, 120, 140};
	pair_t walk[32];
	uint32_t it = geo_iter(db, s, l, 3);
	size_t nw = 0;
	while (nw < 32 && geo_next(walk[nw].p, &walk[nw].ref, it))
		nw++;

	size_t ne = 0;
	for (int i = 0; i < 24; i++) {
		int inside = 1;
		for (int d = 0; d < 3; d++)
			if (placed[i].p[d] < s[d] || placed[i].p[d] > e[d]) {
				inside = 0;
				break;
			}
		if (inside)
			ne++;
	}
	ASSERT_EQ(ne, 19); /* 5 lattice points at y=61 lie outside y<=50 */
	ASSERT_EQ(nw, ne);

	pair_t *exp = malloc(ne * sizeof *exp);
	size_t k = 0;
	for (int i = 0; i < 24; i++) {
		int inside = 1;
		for (int d = 0; d < 3; d++)
			if (placed[i].p[d] < s[d] || placed[i].p[d] > e[d]) {
				inside = 0;
				break;
			}
		if (inside)
			exp[k++] = placed[i];
	}
	qsort(exp, ne, sizeof *exp, cmp_pair);
	qsort(walk, nw, sizeof *walk, cmp_pair);
	for (size_t i = 0; i < ne; i++) {
		ASSERT_EQ(exp[i].ref, walk[i].ref);
		ASSERT_POINT_EQ(exp[i].p, walk[i].p, 3);
	}
	free(exp);

	qmap_close(db);
	unlink(PLC_F);
}

int main(void)
{
	test_suite_begin("Geo Placement Persistence Integration Tests");
	RUN_TEST(placement_persist_exact_cloud);
	return test_suite_end();
}