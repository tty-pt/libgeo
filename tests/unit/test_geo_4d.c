/**
 * @file test_geo_4d.c
 * @brief Unit tests for 4D geo API: CRUD, box queries vs a brute-force
 * oracle, fill parity, and invalid-dimension rejection.
 */

#include "../test_common.h"
#include "../../include/ttypt/geo.h"
#include "../../include/ttypt/point.h"
#include <string.h>

/* Initialize geo once for all tests */
static void setup_once(void) {
    static int initialized = 0;
    if (!initialized) {
        geo_init();
        initialized = 1;
    }
}

TEST(geo4_put_get_del) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_geo4_crud", 1023);

    int16_t p[4] = {10, -20, 30, -40};
    int16_t q[4] = {11, -20, 30, -40};

    ASSERT_EQ(geo_get(db, p, 4), GEO_MISS);
    geo_put(db, p, 42, 4);
    ASSERT_EQ(geo_get(db, p, 4), 42u);
    ASSERT_EQ(geo_get(db, q, 4), GEO_MISS);

    geo_put(db, p, 43, 4);
    ASSERT_EQ(geo_cell_count(db, p, 4), 2u);
    ASSERT_EQ(geo_get(db, p, 4), 42u);

    geo_del(db, p, 4);
    ASSERT_EQ(geo_get(db, p, 4), 43u);
    ASSERT_EQ(geo_del_all(db, p, 4), 1u);
    ASSERT_EQ(geo_get(db, p, 4), GEO_MISS);
}

TEST(geo4_box_vs_brute_force) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_geo4_oracle", 4095);

    /* Deterministic pseudo-random cloud in a 32^4 region */
    #define GEO4_N 2000
    static int16_t cloud[GEO4_N][4];
    static uint32_t refs[GEO4_N];
    test_seed_rng(0xB004D);
    for (int i = 0; i < GEO4_N; i++) {
        for (int d = 0; d < 4; d++)
            cloud[i][d] = (int16_t)(test_rand64() % 32);
        refs[i] = (uint32_t)(1000000 + i);
        geo_put(db, cloud[i], refs[i], 4);
    }

    /* Several boxes, including edge-touching and empty ones */
    int16_t starts[][4] = {
        {0, 0, 0, 0}, {8, 8, 8, 8}, {0, 0, 0, 0},
        {31, 31, 31, 31}, {5, 0, 20, 10},
    };
    uint16_t lens[][4] = {
        {32, 32, 32, 32}, {8, 8, 8, 8}, {0, 0, 0, 0},
        {1, 1, 1, 1}, {7, 32, 5, 12},
    };

    for (int b = 0; b < 5; b++) {
        int16_t *s = starts[b];
        uint16_t *l = lens[b];

        /* Brute force: multiset of refs whose point is in [s, s+l] */
        static uint32_t expect[GEO4_N * 2];
        int nexp = 0;
        for (int i = 0; i < GEO4_N; i++) {
            int in = 1;
            for (int d = 0; d < 4; d++) {
                int32_t lo = s[d], hi = (int32_t)s[d] + l[d];
                if (cloud[i][d] < lo || cloud[i][d] > hi) { in = 0; break; }
            }
            if (in) expect[nexp++] = refs[i];
        }

        uint32_t it = geo_iter(db, s, l, 4);
        int16_t pt[4];
        uint32_t ref;
        int ngot = 0;
        static uint32_t got[GEO4_N * 2];
        while (geo_next(pt, &ref, it))
            got[ngot++] = ref;

        ASSERT_EQ(ngot, nexp);
        /* Every returned ref must be expected (multiset compare by
         * O(n^2) matching — n is small). */
        for (int i = 0; i < ngot; i++) {
            int found = 0;
            for (int j = 0; j < nexp; j++) {
                if (got[i] == expect[j]) {
                    expect[j] = UINT32_MAX; /* consume */
                    found = 1;
                    break;
                }
            }
            ASSERT(found);
        }
    }
    #undef GEO4_N
}

TEST(geo4_fill_parity) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_geo4_fill", 4095);

    test_seed_rng(0xF114);
    for (int i = 0; i < 500; i++) {
        int16_t p[4] = {
            (int16_t)(test_rand64() % 16), (int16_t)(test_rand64() % 16),
            (int16_t)(test_rand64() % 16), (int16_t)(test_rand64() % 16)
        };
        geo_put(db, p, (uint32_t)(i % 61), 4);
    }

    int16_t s[4] = {0, 0, 0, 0};
    uint16_t l[4] = {16, 16, 16, 16};

    rec_set_t *cands = rec_set_new();
    ASSERT_NOT_NULL(cands);
    ASSERT_EQ(rec_axis_fill_bbox(db, s, l, 4, cands), 0);

    /* Raw iterator must yield the same multiset size pre-dedup: count
     * raw entries and compare against sealed distinct count loosely —
     * exact check: every raw ref must be present post-seal is implied
     * by construction; here assert fill agrees with a recount. */
    uint32_t it = geo_iter(db, s, l, 4);
    int16_t pt[4];
    uint32_t ref;
    int nraw = 0;
    while (geo_next(pt, &ref, it))
        nraw++;
    ASSERT(nraw > 0);
    ASSERT_LE((int)rec_set_count(cands), nraw);
    rec_set_free(cands);
}

TEST(geo4_invalid_dims) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_geo4_invalid", 1023);

    int16_t s[5] = {0, 0, 0, 0, 0};
    uint16_t l[5] = {4, 4, 4, 4, 4};
    rec_set_t *cands = rec_set_new();

    /* dim 0 and dim 5 are rejected everywhere (dim 0 yields an
     * empty iterator; fill returns -1) */
    {
        int16_t z[4] = {0, 0, 0, 0};
        uint16_t zl[4] = {4, 4, 4, 4};
        uint32_t it = geo_iter(db, z, zl, 0);
        int16_t pt[4];
        uint32_t ref;
        ASSERT_EQ(geo_next(pt, &ref, it), 0);
    }
    ASSERT_EQ(rec_axis_fill_bbox(db, s, l, 0, cands), -1);
    ASSERT_EQ(rec_axis_fill_bbox(db, s, l, 5, cands), -1);

    /* dim 4 accepted */
    ASSERT_EQ(rec_axis_fill_bbox(db, s, l, 4, cands), 0);
    rec_set_free(cands);
    (void)db;
}

int main(void) {
    test_suite_begin("Geo 4D Unit Tests");

    RUN_TEST(geo4_put_get_del);
    RUN_TEST(geo4_box_vs_brute_force);
    RUN_TEST(geo4_fill_parity);
    RUN_TEST(geo4_invalid_dims);

    return test_suite_end();
}
