/**
 * @file test_geo_fill.c
 * @brief Unit tests for rec_axis_fill_bbox() (kernel-form space adapter).
 */

#include "../test_common.h"
#include "../../include/ttypt/geo.h"
#include "../../include/ttypt/point.h"
#include "../../include/ttypt/morton.h"
#include <string.h>

static void setup_once(void) {
    static int initialized = 0;
    if (!initialized) {
        geo_init();
        initialized = 1;
    }
}

/* Empty database, small box: success, zero refs */
TEST(fill_empty_box) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_fill_empty", 1023);

    int16_t s[3] = {0, 0, 0};
    uint16_t l[3] = {10, 10, 10};
    rec_set_t *out = rec_set_new();

    ASSERT_EQ(rec_axis_fill_bbox(db, s, l, 3, out), 0);
    ASSERT_EQ(rec_set_count(out), (size_t)0);

    rec_set_free(out);
}

/* Single cell, single value */
TEST(fill_single_cell) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_fill_single", 1023);

    int16_t at[3] = {5, 5, 5};
    geo_put(db, at, 77, 3);

    int16_t s[3] = {0, 0, 0};
    uint16_t l[3] = {10, 10, 10};
    rec_set_t *out = rec_set_new();

    ASSERT_EQ(rec_axis_fill_bbox(db, s, l, 3, out), 0);
    ASSERT_EQ(rec_set_count(out), (size_t)1);
    ASSERT_EQ(rec_set_at(out)[0], (rec_ref_t)77);

    rec_set_free(out);
}

/* Multi-value cell: every sibling enters the set */
TEST(fill_mv_cell) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_fill_mv", 1023);

    int16_t at[3] = {5, 5, 5};
    geo_put(db, at, 11, 3);
    geo_put(db, at, 22, 3);

    int16_t s[3] = {0, 0, 0};
    uint16_t l[3] = {10, 10, 10};
    rec_set_t *out = rec_set_new();

    ASSERT_EQ(rec_axis_fill_bbox(db, s, l, 3, out), 0);
    ASSERT_EQ(rec_set_count(out), (size_t)2);
    ASSERT_EQ(rec_set_at(out)[0], (rec_ref_t)11);
    ASSERT_EQ(rec_set_at(out)[1], (rec_ref_t)22);

    rec_set_free(out);
}

/* Several cells: all values collected */
TEST(fill_multi_cell) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_fill_multi", 1023);

    for (int i = 0; i < 5; i++) {
        int16_t p[3] = {i, i, i};
        geo_put(db, p, 100 + i, 3);
    }

    int16_t s[3] = {0, 0, 0};
    uint16_t l[3] = {10, 10, 10};
    rec_set_t *out = rec_set_new();

    ASSERT_EQ(rec_axis_fill_bbox(db, s, l, 3, out), 0);
    ASSERT_EQ(rec_set_count(out), (size_t)5);
    for (int i = 0; i < 5; i++)
        ASSERT_EQ(rec_set_at(out)[i], (rec_ref_t)(100 + i));

    rec_set_free(out);
}

/* 2D boxes work */
TEST(fill_dim2) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_fill_dim2", 1023);

    int16_t a[2] = {3, 4};
    int16_t b[2] = {7, 1};
    geo_put(db, a, 5, 2);
    geo_put(db, b, 6, 2);

    int16_t s[2] = {0, 0};
    uint16_t l[2] = {10, 10};
    rec_set_t *out = rec_set_new();

    ASSERT_EQ(rec_axis_fill_bbox(db, s, l, 2, out), 0);
    ASSERT_EQ(rec_set_count(out), (size_t)2);

    rec_set_free(out);
}

/* Points inside the morton interval but outside the box are excluded */
TEST(fill_excludes_morton_false_positives) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_fill_fp", 1023);

    int16_t s[3] = {0, 0, 0};
    uint16_t l[3] = {4, 4, 4};
    int16_t e[3];
    point_add(e, s, (int16_t *)l, 3);
    uint64_t rmin = morton_set(s, 3);
    uint64_t rmax = morton_set(e, 3);

    /* Find witnesses: outside the box, inside the morton interval */
    int16_t wit[8][3];
    int nwit = 0;
    for (int16_t x = -8; x <= 8 && nwit < 8; x++)
        for (int16_t y = -8; y <= 8 && nwit < 8; y++)
            for (int16_t z = -8; z <= 8 && nwit < 8; z++) {
                int16_t p[3] = {x, y, z};
                int inside = (x >= s[0] && x <= e[0] &&
                              y >= s[1] && y <= e[1] &&
                              z >= s[2] && z <= e[2]);
                if (inside)
                    continue;
                uint64_t code = morton_set(p, 3);
                if (code >= rmin && code <= rmax) {
                    wit[nwit][0] = x;
                    wit[nwit][1] = y;
                    wit[nwit][2] = z;
                    nwit++;
                }
            }
    ASSERT_GT(nwit, 0); /* the trap must exist for the test to mean anything */

    /* One true in-box point plus the false-positive witnesses */
    int16_t good[3] = {1, 1, 1};
    geo_put(db, good, 1, 3);
    for (int i = 0; i < nwit; i++)
        geo_put(db, wit[i], 500 + i, 3);

    rec_set_t *out = rec_set_new();
    ASSERT_EQ(rec_axis_fill_bbox(db, s, l, 3, out), 0);
    ASSERT_EQ(rec_set_count(out), (size_t)1);
    ASSERT_EQ(rec_set_at(out)[0], (rec_ref_t)1);

    rec_set_free(out);
}

/* Sealed output: sorted, deduplicated across cells and siblings */
TEST(fill_sealed_sorted_unique) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_fill_seal", 1023);

    int16_t a[3] = {1, 1, 1};
    int16_t b[3] = {2, 2, 2};
    int16_t c[3] = {3, 3, 3};
    geo_put(db, a, 30, 3);
    geo_put(db, a, 10, 3);
    geo_put(db, b, 20, 3);
    geo_put(db, b, 10, 3); /* duplicate value across cells */
    geo_put(db, c, 20, 3); /* duplicate value across cells */

    int16_t s[3] = {0, 0, 0};
    uint16_t l[3] = {10, 10, 10};
    rec_set_t *out = rec_set_new();

    ASSERT_EQ(rec_axis_fill_bbox(db, s, l, 3, out), 0);
    ASSERT_EQ(rec_set_count(out), (size_t)3);
    ASSERT_EQ(rec_set_at(out)[0], (rec_ref_t)10);
    ASSERT_EQ(rec_set_at(out)[1], (rec_ref_t)20);
    ASSERT_EQ(rec_set_at(out)[2], (rec_ref_t)30);

    rec_set_free(out);
}

/* Oversize boxes are rejected, set untouched */
TEST(fill_rejects_oversize) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_fill_big", 1023);

    int16_t s[3] = {0, 0, 0};
    uint16_t l[3] = {1024, 1024, 2}; /* 2M cells > cap */
    rec_set_t *out = rec_set_new();

    ASSERT_EQ(rec_axis_fill_bbox(db, s, l, 3, out), -1);
    ASSERT_EQ(rec_set_count(out), (size_t)0);

    rec_set_free(out);
}

/* Exactly GEO_FILL_MAX_VOL cells is still accepted */
TEST(fill_allows_exact_cap) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_fill_cap", 1023);

    int16_t s[3] = {0, 0, 0};
    uint16_t l[3] = {100, 100, 100}; /* exactly 1M cells */
    rec_set_t *out = rec_set_new();

    ASSERT_EQ(rec_axis_fill_bbox(db, s, l, 3, out), 0);
    ASSERT_EQ(rec_set_count(out), (size_t)0);

    rec_set_free(out);
}

/* Bad dimensions and NULL set are rejected */
TEST(fill_rejects_bad_args) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_fill_args", 1023);

    int16_t s[3] = {0, 0, 0};
    uint16_t l[3] = {4, 4, 4};
    rec_set_t *out = rec_set_new();

    ASSERT_EQ(rec_axis_fill_bbox(db, s, l, 0, out), -1);
    ASSERT_EQ(rec_axis_fill_bbox(db, s, l, 4, out), -1);
    ASSERT_EQ(rec_axis_fill_bbox(db, s, l, 3, NULL), -1);
    ASSERT_EQ(rec_set_count(out), (size_t)0);

    rec_set_free(out);
}

/* Reusing a set unions (then re-seals) */
TEST(fill_reuses_set) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_fill_reuse", 1023);

    int16_t a[3] = {1, 1, 1};
    int16_t b[3] = {8, 8, 8};
    geo_put(db, a, 5, 3);
    geo_put(db, b, 9, 3);

    rec_set_t *out = rec_set_new();

    int16_t s1[3] = {0, 0, 0};
    uint16_t l1[3] = {2, 2, 2};
    ASSERT_EQ(rec_axis_fill_bbox(db, s1, l1, 3, out), 0);
    ASSERT_EQ(rec_set_count(out), (size_t)1);

    int16_t s2[3] = {7, 7, 7};
    uint16_t l2[3] = {2, 2, 2};
    ASSERT_EQ(rec_axis_fill_bbox(db, s2, l2, 3, out), 0);
    ASSERT_EQ(rec_set_count(out), (size_t)2);
    ASSERT_EQ(rec_set_at(out)[0], (rec_ref_t)5);
    ASSERT_EQ(rec_set_at(out)[1], (rec_ref_t)9);

    rec_set_free(out);
}

int main(void) {
    test_suite_begin("Geo Fill (rec_axis_fill_bbox) Unit Tests");

    RUN_TEST(fill_empty_box);
    RUN_TEST(fill_single_cell);
    RUN_TEST(fill_mv_cell);
    RUN_TEST(fill_multi_cell);
    RUN_TEST(fill_dim2);
    RUN_TEST(fill_excludes_morton_false_positives);
    RUN_TEST(fill_sealed_sorted_unique);
    RUN_TEST(fill_rejects_oversize);
    RUN_TEST(fill_allows_exact_cap);
    RUN_TEST(fill_rejects_bad_args);
    RUN_TEST(fill_reuses_set);

    return test_suite_end();
}
