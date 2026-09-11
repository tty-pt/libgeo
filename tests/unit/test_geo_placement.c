/**
 * @file test_geo_placement.c
 * @brief Unit tests: every value is found exactly where it was placed.
 *
 * The load-bearing property here is geographic IDENTITY, not just counts:
 * * geo_get_3(P) returns the value put at P, and every immediate neighbor of P
 *   misses;
 * * a box walk returns exactly the (point, value) pairs placed inside the
 *   box — nothing from outside, nothing relocated to another cell, no
 *   invented pairs;
 * * bounding-box endpoints and the int16 coordinate extremes round-trip;
 * * the storage dimension is part of the geography (dim 2 vs dim 3 at the
 *   same physical prefix are different cells).
 *
 * The pair comparison canonicalizes on (morton code, value), so a value
 * stored at one coordinate can never satisfy a query for another.
 */

#include "../test_common.h"
#include "../../include/ttypt/geo.h"
#include "../../include/ttypt/point.h"
#include "../../include/ttypt/morton.h"
#include <stdlib.h>
#include <string.h>

typedef struct {
	int16_t p[3];
	uint32_t ref;
} pair_t;

static uint8_t g_pdim;

static int cmp_pair(const void *va, const void *vb)
{
	const pair_t *a = va, *b = vb;
	uint64_t ca = geo_ops[g_pdim].morton_set((int16_t *)(void *)a->p);
	uint64_t cb = geo_ops[g_pdim].morton_set((int16_t *)(void *)b->p);

	if (ca != cb)
		return ca < cb ? -1 : 1;
	if (a->ref != b->ref)
		return a->ref < b->ref ? -1 : 1;
	return 0;
}

static void assert_pair_multiset(const pair_t *exp, size_t ne,
		const pair_t *walk, size_t nw, uint8_t dim)
{
	ASSERT_EQ(nw, ne);
	if (ne == 0)
		return;

	g_pdim = dim;
	pair_t *a = malloc(ne * sizeof *a);
	pair_t *b = malloc(nw * sizeof *b);

	memcpy(a, exp, ne * sizeof *a);
	memcpy(b, walk, nw * sizeof *b);
	qsort(a, ne, sizeof *a, cmp_pair);
	qsort(b, nw, sizeof *b, cmp_pair);
	for (size_t i = 0; i < ne; i++) {
		ASSERT_EQ(a[i].ref, b[i].ref);
		ASSERT_POINT_EQ(a[i].p, b[i].p, dim);
	}
	free(a);
	free(b);
}

/* Walk a box, asserting every returned point is geometrically inside it and
 * that codes come out non-decreasing (morton-discovery order). */
static size_t walk_collect(uint32_t db, int16_t *s, uint16_t *l,
		uint8_t dim, pair_t *out, size_t cap)
{
	int16_t e[4];
	geo_ops[dim].point_add(e, s, (int16_t *)l);
	uint32_t it = geo_ops[dim].iter(db, s, l);
	size_t n = 0;
	uint64_t prev = 0;
	int first = 1;

	while (n < cap && geo_next(out[n].p, &out[n].ref, it)) {
		for (uint8_t i = 0; i < dim; i++)
			ASSERT(out[n].p[i] >= s[i] && out[n].p[i] <= e[i]);
		uint64_t code = geo_ops[dim].morton_set(out[n].p);
		if (!first)
			ASSERT_GE(code, prev);
		prev = code;
		first = 0;
		n++;
	}
	return n;
}

static void setup_once(void)
{
	static int initialized = 0;

	if (!initialized) {
		geo_init();
		initialized = 1;
	}
}

/* --- exact point retrieval; immediate neighbors must MISS --- */
TEST(placement_exact_get_neighbors) {
	setup_once();
	uint32_t db = geo_open(NULL, "plc_exact", 1023);

	int16_t p[3] = {12, -7, 300};
	geo_put_3(db, p, 0xDEADBEEFu);
	ASSERT_EQ(geo_get_3(db, p), 0xDEADBEEFu);

	const int delta[6][3] = {
		{-1, 0, 0}, {1, 0, 0}, {0, -1, 0},
		{0, 1, 0}, {0, 0, -1}, {0, 0, 1},
	};
	for (int i = 0; i < 6; i++) {
		int16_t n[3] = {p[0] + delta[i][0], p[1] + delta[i][1],
			p[2] + delta[i][2]};
		ASSERT_EQ(geo_get_3(db, n), GEO_MISS);
	}

	/* single-cell box at P yields exactly the placed pair */
	int16_t s[3] = {p[0], p[1], p[2]};
	uint16_t l[3] = {0, 0, 0};
	pair_t walk[8];
	size_t nw = walk_collect(db, s, l, 3, walk, 8);
	ASSERT_EQ(nw, 1);
	ASSERT_EQ(walk[0].ref, 0xDEADBEEFu);
	ASSERT_POINT_EQ(walk[0].p, p, 3);
}

/* --- bounding-box endpoints are inclusive: s and s+l are both found --- */
TEST(placement_inclusive_corners) {
	setup_once();
	uint32_t db = geo_open(NULL, "plc_corners", 1023);

	int16_t s[3] = {-2, -3, -4};
	uint16_t l[3] = {5, 7, 9};
	int16_t e[3] = {s[0] + l[0], s[1] + l[1], s[2] + l[2]};
	int16_t mid[3] = {0, 1, 0};

	pair_t exp[3] = {
		{{s[0], s[1], s[2]}, 1001}, {{e[0], e[1], e[2]}, 1002},
		{{mid[0], mid[1], mid[2]}, 1003},
	};
	for (int i = 0; i < 3; i++)
		geo_put_3(db, exp[i].p, exp[i].ref);

	/* one cell outside each face must stay excluded */
	geo_put_3(db, (int16_t[3]){e[0] + 1, e[1], e[2]}, 2001);
	geo_put_3(db, (int16_t[3]){s[0] - 1, s[1], s[2]}, 2002);
	geo_put_3(db, (int16_t[3]){e[0], e[1] + 1, e[2]}, 2003);
	geo_put_3(db, (int16_t[3]){e[0], e[1], e[2] + 1}, 2004);

	pair_t walk[16];
	size_t nw = walk_collect(db, s, l, 3, walk, 16);
	assert_pair_multiset(exp, 3, walk, nw, 3);
}

/* --- a cloud of distinct (point -> value): nothing lost, nothing extra,
 *      nothing relocated --- */
TEST(placement_value_identity_cloud) {
	setup_once();
	uint32_t db = geo_open(NULL, "plc_cloud", 4095);

	pair_t placed[40];
	for (int i = 0; i < 40; i++) {
		placed[i].p[0] = (int16_t)((i * 7) % 51 - 13);
		placed[i].p[1] = (int16_t)((i * 13) % 49 - 17);
		placed[i].p[2] = (int16_t)((i * 17) % 47 - 19);
		placed[i].ref = 50000u + (uint32_t)i * 3u;
		geo_put_3(db, placed[i].p, placed[i].ref);
	}

	/* enclosing box bound from the model (pad by 1) */
	int16_t s[3] = {placed[0].p[0], placed[0].p[1], placed[0].p[2]};
	int16_t e[3] = {s[0], s[1], s[2]};
	for (int i = 1; i < 40; i++) {
		for (int d = 0; d < 3; d++) {
			if (placed[i].p[d] < s[d])
				s[d] = placed[i].p[d];
			if (placed[i].p[d] > e[d])
				e[d] = placed[i].p[d];
		}
	}
	uint16_t l[3] = { (uint16_t)(e[0] - s[0] + 1),
		(uint16_t)(e[1] - s[1] + 1), (uint16_t)(e[2] - s[2] + 1) };

	pair_t walk[64];
	size_t nw = walk_collect(db, s, l, 3, walk, 64);
	assert_pair_multiset(placed, 40, walk, nw, 3);

	/* every point is retrievable at exactly its own coordinate */
	for (int i = 0; i < 40; i++)
		ASSERT_EQ(geo_get_3(db, placed[i].p), placed[i].ref);
}

/* --- int16 extremes round-trip at exactly the placed coordinate --- */
TEST(placement_boundary_int16) {
	setup_once();
	uint32_t db = geo_open(NULL, "plc_bound", 1023);

	pair_t placed[5] = {
		{{-32768, -32768, -32768}, 1},
		{{-32767, 0, 5}, 2},
		{{0, 0, 0}, 3},
		{{32766, 32767, 32766}, 4},
		{{32767, 32767, 32767}, 5},
	};
	for (int i = 0; i < 5; i++)
		geo_put_3(db, placed[i].p, placed[i].ref);

	for (int i = 0; i < 5; i++)
		ASSERT_EQ(geo_get_3(db, placed[i].p), placed[i].ref);

	/* exact single-cell boxes at both extremes */
	uint16_t z[3] = {0, 0, 0};
	pair_t walk[4];
	ASSERT_EQ(walk_collect(db, placed[0].p, z, 3, walk, 4), 1);
	ASSERT_EQ(walk[0].ref, placed[0].ref);
	ASSERT_EQ(walk_collect(db, placed[4].p, z, 3, walk, 4), 1);
	ASSERT_EQ(walk[0].ref, placed[4].ref);

	/* small boxes rooted in deep-negative / deep-positive territory */
	int16_t s1[3] = {-32768, -32768, -32768};
	uint16_t l1[3] = {0, 1, 5}; /* only placed[0] (y == -32768) */
	ASSERT_EQ(walk_collect(db, s1, l1, 3, walk, 4), 1);
	ASSERT_EQ(walk[0].ref, placed[0].ref);

	int16_t s2[3] = {32766, 32767, 32766};
	uint16_t l2[3] = {1, 0, 1}; /* placed[3] and placed[4] */
	pair_t exp[2] = {placed[3], placed[4]};
	ASSERT_EQ(walk_collect(db, s2, l2, 3, walk, 4), 2);
	assert_pair_multiset(exp, 2, walk, 2, 3);
}

/* --- storage dimension is part of the geography --- */
TEST(placement_dim_isolation) {
	setup_once();
	uint32_t d2 = geo_open(NULL, "plc_dim2", 1023);
	uint32_t d3 = geo_open(NULL, "plc_dim3", 1023);

	int16_t p2[2] = {5, 5};
	int16_t p3[3] = {5, 5, 0};

	geo_put_2(d2, p2, 1);
	geo_put_3(d3, p3, 2);
	geo_put_3(d3, p3, 3);

	ASSERT_EQ(geo_get_2(d2, p2), 1);
	ASSERT_EQ(geo_get_3(d3, p3), 2);
	ASSERT_EQ(geo_cell_count_2(d2, p2), 1);
	ASSERT_EQ(geo_cell_count_3(d3, p3), 2);

	/* the other-dim code is a different cell: MISS */
	ASSERT_EQ(geo_get_3(d2, p3), GEO_MISS);
	ASSERT_EQ(geo_get_2(d3, p2), GEO_MISS);

	/* 2D box finds only the 2D value; 3D box only the 3D siblings */
	int16_t s2[2] = {4, 4}, s3[3] = {3, 3, -1};
	uint16_t l2[2] = {3, 3}, l3[3] = {5, 5, 3};
	pair_t w2[8], w3[8], e3[2] = {{{5, 5, 0}, 2}, {{5, 5, 0}, 3}};

	ASSERT_EQ(walk_collect(d2, s2, l2, 2, w2, 8), 1);
	ASSERT_EQ(w2[0].ref, 1);
	ASSERT_EQ(w2[0].p[0], 5);
	ASSERT_EQ(w2[0].p[1], 5);

	ASSERT_EQ(walk_collect(d3, s3, l3, 3, w3, 8), 2);
	assert_pair_multiset(e3, 2, w3, 2, 3);
}

/* --- sparse diagonal across the full int16 range stays exact --- */
TEST(placement_diagonal_full_range) {
	setup_once();
	uint32_t db = geo_open(NULL, "plc_diag", 1023);

	pair_t placed[9];
	for (int i = 0; i < 9; i++) {
		int16_t c = (int16_t)(-32768 + i * (65535 / 8));
		placed[i].p[0] = c;
		placed[i].p[1] = (int16_t)(-c);
		placed[i].p[2] = (int16_t)(c / 2 - 100);
		placed[i].ref = 9000u + (uint32_t)i;
		geo_put_3(db, placed[i].p, placed[i].ref);
		ASSERT_EQ(geo_get_3(db, placed[i].p), placed[i].ref);
	}

	uint16_t z[3] = {0, 0, 0};
	pair_t walk[4];
	for (int i = 0; i < 9; i++) {
		ASSERT_EQ(walk_collect(db, placed[i].p, z, 3, walk, 4), 1);
		ASSERT_EQ(walk[0].ref, placed[i].ref);
		ASSERT_POINT_EQ(walk[0].p, placed[i].p, 3);
	}
}

int main(void)
{
	test_suite_begin("Geo Placement Unit Tests");
	RUN_TEST(placement_exact_get_neighbors);
	RUN_TEST(placement_inclusive_corners);
	RUN_TEST(placement_value_identity_cloud);
	RUN_TEST(placement_boundary_int16);
	RUN_TEST(placement_dim_isolation);
	RUN_TEST(placement_diagonal_full_range);
	return test_suite_end();
}