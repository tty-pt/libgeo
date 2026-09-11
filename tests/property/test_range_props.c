/**
 * @file test_range_props.c
 * @brief Property tests: the skipping walker always agrees with a
 *        brute-force full-map-scan oracle.
 *
 * True oracle (not self-referential): the reference path scans the whole
 * map with qmap_iter NULL and filters by box membership, independently of
 * islet_box_visit's GE seeks, restarts, and jumps. Randomized clouds with MV
 * collisions and duplicate values, random boxes, 1D/2D/3D, plus edit churn
 * between queries (dirty sorted-index rebuild path). Deterministic PRNG.
 */

#include "../test_common.h"
#include "../../include/ttypt/islet.h"
#include "../../include/ttypt/point.h"
#include "../../include/ttypt/morton.h"
#include "../../include/ttypt/qmap.h"
#include <stdlib.h>
#include <string.h>

#define RANGE_PROP_CLOUDS 10
#define RANGE_PROP_BOXES 30
#define RANGE_PROP_POINTS 180

static int cmp_u32(const void *a, const void *b) {
    uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;
    return x > y ? 1 : (x < y ? -1 : 0);
}

static void check_box(uint32_t db, int16_t *s, uint16_t *l, uint8_t dim) {
    int16_t e[4], p[4];
    const void *key, *value;
    size_t bcap = 64, bn = 0, interval = 0;
    uint32_t *brute = malloc(bcap * sizeof *brute);

    islet_ops[dim].point_add(e, s, (int16_t *)l);
    uint64_t rmin = islet_ops[dim].morton_set(s);
    uint64_t rmax = islet_ops[dim].morton_set(e);

    uint32_t bcur = qmap_iter(db, NULL, 0);
    while (qmap_next(&key, &value, bcur)) {
        uint64_t code = *(uint64_t *)key;

        if (code < rmin || code > rmax)
            continue;
        interval++;

        islet_ops[dim].morton_get(p, code);
        int inside = 1;
        for (uint8_t i = 0; i < dim; i++)
            if (p[i] < s[i] || p[i] > e[i]) {
                inside = 0;
                break;
            }
        if (!inside)
            continue;

        if (bn == bcap) {
            bcap *= 2;
            brute = realloc(brute, bcap * sizeof *brute);
        }
        brute[bn++] = *(uint32_t *)value;
    }
    qmap_fin(bcur);

    size_t wcap = 64, wn = 0;
    uint32_t *walk = malloc(wcap * sizeof *walk);

    /* Collect via the walker under test. */
    {
        uint32_t witer = islet_ops[dim].iter(db, s, l);
        int16_t wp[4];
        uint32_t wref;
        while (islet_next(wp, &wref, witer)) {
            if (wn == wcap) {
                wcap *= 2;
                walk = realloc(walk, wcap * sizeof *walk);
            }
            walk[wn++] = wref;
        }
    }
    uint32_t scanned = islet_last_scan_count();

    ASSERT_EQ(wn, bn);
    qsort(brute, bn, sizeof *brute, cmp_u32);
    qsort(walk, wn, sizeof *walk, cmp_u32);
    for (size_t i = 0; i < bn; i++)
        ASSERT_EQ(walk[i], brute[i]);
    ASSERT(scanned <= interval);

    free(brute);
    free(walk);
}

static void run_clouds(uint8_t dim, int span, int npoints, uint64_t seed,
                       const char *tag) {
    static char names[RANGE_PROP_CLOUDS][64];

    test_seed_rng(seed);
    for (int c = 0; c < RANGE_PROP_CLOUDS; c++) {
        snprintf(names[c], sizeof names[c], "prop_range_%s_%d", tag, c);
        uint32_t db = islet_open(NULL, names[c], 4095);

        for (int i = 0; i < npoints; i++) {
            int16_t p[4] = {0, 0, 0, 0};
            for (uint8_t d = 0; d < dim; d++)
                p[d] = test_rand_coord_range(0, span);
            /* Small value domain => duplicate values + cell collisions. */
            islet_ops[dim].put(db, p, (uint32_t)(test_rand64() % 41));
        }

        /* Edit churn on even clouds: delete + re-add (dirty rebuild). */
        if (c % 2 == 0) {
            for (int i = 0; i < npoints / 6; i++) {
                int16_t p[4] = {0, 0, 0, 0};
                for (uint8_t d = 0; d < dim; d++)
                    p[d] = test_rand_coord_range(0, span);
                islet_ops[dim].del_all(db, p);
            }
            for (int i = 0; i < npoints / 6; i++) {
                int16_t p[4] = {0, 0, 0, 0};
                for (uint8_t d = 0; d < dim; d++)
                    p[d] = test_rand_coord_range(0, span);
                islet_ops[dim].put(db, p, (uint32_t)(test_rand64() % 41));
            }
        }

        for (int b = 0; b < RANGE_PROP_BOXES; b++) {
            int16_t s[4] = {0, 0, 0, 0};
            uint16_t l[4] = {1, 1, 1, 1};
            for (uint8_t d = 0; d < dim; d++) {
                s[d] = test_rand_coord_range(-4, span);
                l[d] = (uint16_t)(test_rand64() % (span / 3) + 1);
            }
            check_box(db, s, l, dim);
        }
    }
}

TEST(property_range_oracle_3d) {
    islet_init();
    run_clouds(3, 48, RANGE_PROP_POINTS, 9001, "3d");
}

TEST(property_range_oracle_2d) {
    islet_init();
    run_clouds(2, 48, RANGE_PROP_POINTS, 9002, "2d");
}

TEST(property_range_oracle_1d) {
    islet_init();
    run_clouds(1, 48, RANGE_PROP_POINTS, 9003, "1d");
}

int main(void) {
    test_suite_begin("Islet Range Oracle Property Tests");
    RUN_TEST(property_range_oracle_3d);
    RUN_TEST(property_range_oracle_2d);
    RUN_TEST(property_range_oracle_1d);
    return test_suite_end();
}
