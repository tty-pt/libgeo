/*
 * Integration tests: rec_axis_fill_bbox agrees with the raw geo_iter path.
 * Seeded clouds (distinct + MV + duplicate values); fill (sealed) must equal
 * the deduplicated, sorted raw multiset for every box.
 */

#include "../test_common.h"
#include "../../include/ttypt/geo.h"
#include "../../include/ttypt/point.h"
#include "../../include/ttypt/morton.h"
#include <stdlib.h>
#include <string.h>

static void setup_once(void) {
    static int initialized = 0;
    if (!initialized) {
        geo_init();
        initialized = 1;
    }
}

static int cmp_u32(const void *a, const void *b) {
    uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;
    return x > y ? 1 : (x < y ? -1 : 0);
}

/* Collect the raw path into a malloc'd array; caller frees. */
static uint32_t *raw_collect(uint32_t db, int16_t *s, uint16_t *l, uint8_t dim,
                             size_t *n_out) {
    size_t cap = 64, n = 0;
    uint32_t *vals = malloc(cap * sizeof *vals);
    uint32_t iter = geo_iter(db, s, l, dim);
    int16_t p[4];
    uint32_t ref;

    while (geo_next(p, &ref, iter)) {
        if (n == cap) {
            cap *= 2;
            vals = realloc(vals, cap * sizeof *vals);
        }
        vals[n++] = ref;
    }
    *n_out = n;
    return vals;
}

/* Sort + dedup in place; returns the deduped length. */
static size_t sort_dedup(uint32_t *vals, size_t n) {
    size_t w = 0;

    if (n == 0)
        return 0;
    qsort(vals, n, sizeof *vals, cmp_u32);
    for (size_t r = 1; r < n; r++)
        if (vals[r] != vals[w])
            vals[++w] = vals[r];
    return w + 1;
}

/* Assert fill(db, box) == dedup(sort(raw(db, box))). */
static void assert_parity(uint32_t db, int16_t *s, uint16_t *l, uint8_t dim) {
    size_t nraw = 0;
    uint32_t *raw = raw_collect(db, s, l, dim, &nraw);
    size_t nunion = sort_dedup(raw, nraw);

    rec_set_t *out = rec_set_new();
    ASSERT_EQ(rec_axis_fill_bbox(db, s, l, dim, out), 0);
    ASSERT_EQ(rec_set_count(out), nunion);
    for (size_t i = 0; i < nunion; i++)
        ASSERT_EQ(rec_set_at(out)[i], (rec_ref_t)raw[i]);

    rec_set_free(out);
    free(raw);
}

/* Distinct-values cloud: exact count equality too */
TEST(fill_parity_distinct) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_par_distinct", 4095);

    test_seed_rng(1001);
    for (int i = 0; i < 300; i++) {
        int16_t p[3] = {
            test_rand_coord_range(0, 64),
            test_rand_coord_range(0, 64),
            test_rand_coord_range(0, 64),
        };
        geo_set(db, p, 1000 + i, 3); /* distinct values, one per cell */
    }

    int16_t boxes[][3] = {{0, 0, 0}, {10, 10, 10}, {0, 32, 0}, {50, 50, 50}};
    uint16_t lens[][3] = {{64, 64, 64}, {20, 20, 20}, {64, 32, 64}, {4, 4, 4}};
    for (int b = 0; b < 4; b++)
        assert_parity(db, boxes[b], lens[b], 3);
}

/* MV cloud with duplicate values across cells and siblings */
TEST(fill_parity_mv_dupvals) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_par_mv", 4095);

    test_seed_rng(2002);
    for (int i = 0; i < 200; i++) {
        int16_t p[3] = {
            test_rand_coord_range(0, 32),
            test_rand_coord_range(0, 32),
            test_rand_coord_range(0, 32),
        };
        geo_put(db, p, (uint32_t)(i % 37), 3); /* heavy value overlap + cell collisions */
    }

    int16_t boxes[][3] = {{0, 0, 0}, {5, 5, 5}, {16, 0, 16}, {0, 0, 0}};
    uint16_t lens[][3] = {{32, 32, 32}, {10, 10, 10}, {16, 32, 16}, {1, 1, 1}};
    for (int b = 0; b < 4; b++)
        assert_parity(db, boxes[b], lens[b], 3);
}

/* Empty boxes agree (zero on both paths) */
TEST(fill_parity_empty) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_par_empty", 1023);

    int16_t p[3] = {90, 90, 90};
    geo_put(db, p, 1, 3);

    int16_t s[3] = {0, 0, 0};
    uint16_t l[3] = {10, 10, 10};
    assert_parity(db, s, l, 3);
}

int main(void) {
    test_suite_begin("Geo Fill Parity Integration Tests");
    RUN_TEST(fill_parity_distinct);
    RUN_TEST(fill_parity_mv_dupvals);
    RUN_TEST(fill_parity_empty);
    return test_suite_end();
}
