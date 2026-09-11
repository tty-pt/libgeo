/**
 * @file test_fill_props.c
 * @brief Property tests: fill == raw on randomized clouds and boxes.
 *
 * True cross-path oracle: the sealed fill must equal the deduplicated raw
 * multiset for every random (cloud, box) pair, in 2D and 3D, with MV
 * collisions and duplicate values. Deterministic PRNG.
 */

#include "../test_common.h"
#include "../../include/ttypt/islet.h"
#include "../../include/ttypt/point.h"
#include "../../include/ttypt/morton.h"
#include <stdlib.h>
#include <string.h>

#define FILL_PROP_CLOUDS 12
#define FILL_PROP_BOXES 25
#define FILL_PROP_POINTS 200

static int cmp_u32(const void *a, const void *b) {
    uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;
    return x > y ? 1 : (x < y ? -1 : 0);
}

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

static void check_box(uint32_t db, int16_t *s, uint16_t *l, uint8_t dim) {
    size_t cap = 64, n = 0;
    uint32_t *raw = malloc(cap * sizeof *raw);
    uint32_t iter = islet_ops[dim].iter(db, s, l);
    int16_t p[4];
    uint32_t ref;

    while (islet_next(p, &ref, iter)) {
        if (n == cap) {
            cap *= 2;
            raw = realloc(raw, cap * sizeof *raw);
        }
        raw[n++] = ref;
    }
    size_t nunion = sort_dedup(raw, n);

    rec_set_t *out = rec_set_new();
    ASSERT_EQ(islet_ops[dim].fill(db, s, l, out), 0);
    ASSERT_EQ(rec_set_count(out), nunion);
    for (size_t i = 0; i < nunion; i++)
        ASSERT_EQ(rec_set_at(out)[i], (rec_ref_t)raw[i]);

    rec_set_free(out);
    free(raw);
}

TEST(property_fill_matches_raw_3d) {
    static char names[FILL_PROP_CLOUDS][64];
    islet_init();
    test_seed_rng(777);

    for (int c = 0; c < FILL_PROP_CLOUDS; c++) {
        snprintf(names[c], sizeof names[c], "prop_fill3d_%d", c);
        uint32_t db = islet_open(NULL, names[c], 4095);

        for (int i = 0; i < FILL_PROP_POINTS; i++) {
            int16_t p[3] = {
                test_rand_coord_range(0, 48),
                test_rand_coord_range(0, 48),
                test_rand_coord_range(0, 48),
            };
            islet_put_3(db, p, (uint32_t)(test_rand64() % 61));
        }

        for (int b = 0; b < FILL_PROP_BOXES; b++) {
            int16_t s[3] = {
                test_rand_coord_range(0, 40),
                test_rand_coord_range(0, 40),
                test_rand_coord_range(0, 40),
            };
            uint16_t l[3] = {
                (uint16_t)(test_rand64() % 24 + 1),
                (uint16_t)(test_rand64() % 24 + 1),
                (uint16_t)(test_rand64() % 24 + 1),
            };
            check_box(db, s, l, 3);
        }
    }
}

TEST(property_fill_matches_raw_2d) {
    static char names[FILL_PROP_CLOUDS][64];
    islet_init();
    test_seed_rng(4242);

    for (int c = 0; c < FILL_PROP_CLOUDS; c++) {
        snprintf(names[c], sizeof names[c], "prop_fill2d_%d", c);
        uint32_t db = islet_open(NULL, names[c], 4095);

        for (int i = 0; i < FILL_PROP_POINTS; i++) {
            int16_t p[2] = {
                test_rand_coord_range(0, 48),
                test_rand_coord_range(0, 48),
            };
            islet_put_2(db, p, (uint32_t)(test_rand64() % 61));
        }

        for (int b = 0; b < FILL_PROP_BOXES; b++) {
            int16_t s[2] = {
                test_rand_coord_range(0, 40),
                test_rand_coord_range(0, 40),
            };
            uint16_t l[2] = {
                (uint16_t)(test_rand64() % 24 + 1),
                (uint16_t)(test_rand64() % 24 + 1),
            };
            check_box(db, s, l, 2);
        }
    }
}

int main(void) {
    test_suite_begin("Islet Fill Property Tests");
    RUN_TEST(property_fill_matches_raw_3d);
    RUN_TEST(property_fill_matches_raw_2d);
    return test_suite_end();
}
