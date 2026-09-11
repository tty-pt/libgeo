/**
 * @file test_geo_2d32.c
 * @brief Unit tests for the 2D x 32-bit dense config (geo_*_2_32).
 */

#include "../test_common.h"
#include "../../include/ttypt/geo.h"
#include "../../include/ttypt/point.h"
#include "../../include/ttypt/morton.h"
#include "../../include/ttypt/qmap.h"
#include <string.h>
#include <stdlib.h>
#include <limits.h>

static void setup_once(void) {
    static int initialized = 0;
    if (!initialized) {
        geo_init();
        initialized = 1;
    }
}

static int32_t rand32(void) {
    return (int32_t)test_rand64();
}

/* Edge lanes must all round-trip exactly. */
TEST(codec_roundtrip_edges) {
    static const int32_t lanes[] = {
        INT32_MIN, INT32_MIN + 1, -1, 0, 1, INT32_MAX - 1, INT32_MAX
    };
    size_t n = sizeof lanes / sizeof lanes[0];

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            int32_t p[2] = { lanes[i], lanes[j] };
            uint64_t code = morton_set_2_32(p);
            int32_t back[2] = { 0, 0 };

            morton_get_2_32(back, code);
            ASSERT_EQ(back[0], p[0]);
            ASSERT_EQ(back[1], p[1]);
        }
    }
}

/* Random lanes round-trip; spread preserves order per lane. */
TEST(codec_roundtrip_random) {
    test_seed_rng(0x2d32);

    for (int i = 0; i < 2000; i++) {
        int32_t p[2] = { rand32(), rand32() };
        uint64_t code = morton_set_2_32(p);
        int32_t back[2];

        morton_get_2_32(back, code);
        ASSERT_EQ(back[0], p[0]);
        ASSERT_EQ(back[1], p[1]);
    }
}

/* Known answers + lane monotonicity (spread2 is order-preserving). */
TEST(codec_known_answers) {
    int32_t lo[2] = { INT32_MIN, INT32_MIN };
    int32_t hi[2] = { INT32_MAX, INT32_MAX };

    /* Bias maps -2^31 -> lane 0, 2^31-1 -> lane 0xFFFFFFFF. */
    ASSERT_EQ(morton_set_2_32(lo), (uint64_t)0);
    ASSERT_EQ(morton_set_2_32(hi), UINT64_MAX);

    int32_t a[2] = { 100, -5 };
    int32_t b[2] = { 101, -5 };
    int32_t c[2] = { 100, -4 };
    ASSERT(morton_set_2_32(a) < morton_set_2_32(b));
    ASSERT(morton_set_2_32(a) < morton_set_2_32(c));
}

/* 32-bit-lane point ops smoke test. */
TEST(point_ops_32) {
    int32_t a[2] = { 10, -20 };
    int32_t b[2] = { 3, 7 };
    int32_t t[2];

    point_add_2_32(t, a, b);
    ASSERT_EQ(t[0], 13);
    ASSERT_EQ(t[1], -13);

    point_sub_2_32(t, a, b);
    ASSERT_EQ(t[0], 7);
    ASSERT_EQ(t[1], -27);

    point_min_2_32(t, a, b);
    ASSERT_EQ(t[0], 3);
    ASSERT_EQ(t[1], -20);

    point_max_2_32(t, a, b);
    ASSERT_EQ(t[0], 10);
    ASSERT_EQ(t[1], 7);

    point_copy_2_32(t, a);
    ASSERT_EQ(t[0], 10);
    ASSERT_EQ(t[1], -20);

    point_set_2_32(t, -99);
    ASSERT_EQ(t[0], -99);
    ASSERT_EQ(t[1], -99);

    int32_t l[2] = { 100000, 300 };
    ASSERT_EQ(point_vol_2_32(l), (uint64_t)30000000);

    int32_t s[2] = { 0, 0 };
    int32_t e[2] = { 10, 10 };
    int32_t p[2] = { 3, 5 };
    ASSERT_EQ(point_idx_2_32(p, s, e), (uint64_t)53);
}

/* Basic put/get, including coordinates outside the int16 world. */
TEST(put_get_basic) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_2d32_basic", 1023);

    int32_t p1[2] = { 10, 20 };
    int32_t p2[2] = { 100000, -200000 };
    int32_t p3[2] = { INT32_MIN + 5, INT32_MAX - 5 };

    geo_put_2_32(db, p1, 1);
    geo_put_2_32(db, p2, 2);
    geo_put_2_32(db, p3, 3);

    ASSERT_EQ(geo_get_2_32(db, p1), 1);
    ASSERT_EQ(geo_get_2_32(db, p2), 2);
    ASSERT_EQ(geo_get_2_32(db, p3), 3);

    int32_t miss[2] = { 11, 20 };
    ASSERT_EQ(geo_get_2_32(db, miss), GEO_MISS);
}

/* Replace semantics. */
TEST(set_replace) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_2d32_set", 1023);

    int32_t p[2] = { -500000, 700000 };
    geo_put_2_32(db, p, 1);
    geo_put_2_32(db, p, 2);
    ASSERT_EQ(geo_cell_count_2_32(db, p), 2);

    geo_set_2_32(db, p, 99);
    ASSERT_EQ(geo_get_2_32(db, p), 99);
    ASSERT_EQ(geo_cell_count_2_32(db, p), 1);
}

/* Delete paths. */
TEST(del_paths) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_2d32_del", 1023);

    int32_t p[2] = { 123456789, -123456789 };
    geo_put_2_32(db, p, 10);
    geo_put_2_32(db, p, 20);

    geo_del_2_32(db, p); /* removes earliest only */
    ASSERT_EQ(geo_get_2_32(db, p), 20);
    ASSERT_EQ(geo_cell_count_2_32(db, p), 1);

    ASSERT_EQ(geo_del_all_2_32(db, p), 1);
    ASSERT_EQ(geo_get_2_32(db, p), GEO_MISS);
    ASSERT_EQ(geo_cell_count_2_32(db, p), 0);

    geo_del_2_32(db, p); /* no-op on empty cell */
    ASSERT_EQ(geo_del_all_2_32(db, p), 0);
}

/* Multi-value chain reads back in insertion order. */
TEST(get_multi_chain) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_2d32_mv", 1023);

    int32_t p[2] = { INT32_MIN, INT32_MAX };
    geo_put_2_32(db, p, 111);
    geo_put_2_32(db, p, 222);
    geo_put_2_32(db, p, 333);

    uint32_t cur = geo_get_multi_2_32(db, p);
    ASSERT(cur != QM_MISS);
    uint32_t v;
    ASSERT(geo_cell_next(&v, cur));
    ASSERT_EQ(v, 111);
    ASSERT(geo_cell_next(&v, cur));
    ASSERT_EQ(v, 222);
    ASSERT(geo_cell_next(&v, cur));
    ASSERT_EQ(v, 333);
    ASSERT(!geo_cell_next(&v, cur));

    int32_t empty[2] = { 0, 0 };
    ASSERT_EQ(geo_get_multi_2_32(db, empty), QM_MISS);
}

/* Iterator returns every stored pair; decoded points re-encode. */
TEST(iter_collect_roundtrip) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_2d32_iter", 1023);

    test_seed_rng(1234);
    int32_t pts[64][2];
    for (int i = 0; i < 64; i++) {
        pts[i][0] = rand32();
        pts[i][1] = rand32();
        geo_put_2_32(db, pts[i], 1000 + i);
    }

    int32_t s[2] = { INT32_MIN, INT32_MIN };
    int32_t l[2] = { INT32_MAX, INT32_MAX }; /* e = {-1,-1}: full low half */
    uint32_t iter = geo_iter_2_32(db, s, l);
    int32_t p[2];
    uint32_t ref;
    int n = 0;
    while (geo_next32(p, &ref, iter)) {
        /* Decoded point re-encodes and the value maps back. */
        ASSERT_LT(ref - 1000, 64);
        ASSERT_EQ(pts[ref - 1000][0], p[0]);
        ASSERT_EQ(pts[ref - 1000][1], p[1]);
        ASSERT_EQ(morton_set_2_32(pts[ref - 1000]), morton_set_2_32(p));
        n++;
    }

    /* Only the points with both lanes < 0 fall in the box. */
    int want = 0;
    for (int i = 0; i < 64; i++)
        if (pts[i][0] < 0 && pts[i][1] < 0)
            want++;
    ASSERT_EQ(n, want);
}

static int cmp_u32(const void *a, const void *b) {
    uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;
    return x > y ? 1 : (x < y ? -1 : 0);
}

/* Brute-force oracle: full map scan over [s, s+l], plus the number of
 * stored keys inside the morton interval. Lane-agnostic except for
 * the 32-bit decode. */
static uint32_t *brute_collect32(uint32_t db, int32_t *s, int32_t *l,
                                 size_t *n_out, size_t *interval_out) {
    size_t cap = 64, n = 0, interval = 0;
    uint32_t *vals = malloc(cap * sizeof *vals);
    int32_t e[2], p[2];
    const void *key, *value;

    point_add_2_32(e, s, l);
    uint64_t rmin = morton_set_2_32(s);
    uint64_t rmax = morton_set_2_32(e);

    uint32_t cur = qmap_iter(db, NULL, 0);
    while (qmap_next(&key, &value, cur)) {
        uint64_t code = *(uint64_t *)key;

        if (code < rmin || code > rmax)
            continue;
        interval++;

        morton_get_2_32(p, code);
        if (p[0] < s[0] || p[0] > e[0] || p[1] < s[1] || p[1] > e[1])
            continue;

        if (n == cap) {
            cap *= 2;
            vals = realloc(vals, cap * sizeof *vals);
        }
        vals[n++] = *(uint32_t *)value;
    }
    qmap_fin(cur);

    *n_out = n;
    *interval_out = interval;
    return vals;
}

static uint32_t *walk_collect32(uint32_t db, int32_t *s, int32_t *l,
                                size_t *n_out) {
    size_t cap = 64, n = 0;
    uint32_t *vals = malloc(cap * sizeof *vals);
    uint32_t iter = geo_iter_2_32(db, s, l);
    int32_t p[2];
    uint32_t ref;

    while (geo_next32(p, &ref, iter)) {
        if (n == cap) {
            cap *= 2;
            vals = realloc(vals, cap * sizeof *vals);
        }
        vals[n++] = ref;
    }
    *n_out = n;
    return vals;
}

static void assert_oracle32(uint32_t db, int32_t *s, int32_t *l,
                            int require_skip) {
    size_t nshall = 0, interval = 0, nwalk = 0;
    uint32_t *expected = brute_collect32(db, s, l, &nshall, &interval);
    uint32_t *got = walk_collect32(db, s, l, &nwalk);
    uint32_t scanned = geo_last_scan_count();

    ASSERT_EQ(nwalk, nshall);
    qsort(expected, nshall, sizeof *expected, cmp_u32);
    qsort(got, nwalk, sizeof *got, cmp_u32);
    for (size_t i = 0; i < nshall; i++)
        ASSERT_EQ(got[i], expected[i]);

    ASSERT(scanned <= interval);
    if (require_skip)
        ASSERT(scanned < interval);

    free(expected);
    free(got);
}

/* Mixed-sign database: decoded lanes far below the box drive maxd to
 * ~2^31, exercising the 32-bit kmax clamp, 1u<<31 masks and 2*31-bit
 * spans. Multiset equality proves the clamped jumps stay sound. */
TEST(walk_oracle_mixed) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_2d32_oracle", 4095);

    test_seed_rng(777);
    for (int i = 0; i < 300; i++) {
        int32_t p[2] = { rand32() % 2000 - 1000, rand32() % 2000 - 1000 };
        geo_put_2_32(db, p, 5000 + i);
    }
    /* Far negative points: in-code-interval, far out-of-box. */
    int32_t far[2] = { INT32_MIN + 7, INT32_MIN + 9 };
    geo_put_2_32(db, far, 9999);

    int32_t s[2] = { -100, -100 };
    int32_t l[2] = { 200, 200 };
    assert_oracle32(db, s, l, 0);

    int32_t s2[2] = { 0, 0 };
    int32_t l2[2] = { 50, 50 };
    assert_oracle32(db, s2, l2, 0);
}

/* Dense grid with a queried sub-box: the Z-spill must be jumped over. */
TEST(walk_skip_engages) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_2d32_skip", 8191);

    for (int32_t x = 0; x < 64; x++)
        for (int32_t y = 0; y < 64; y++) {
            int32_t p[2] = { x, y };
            geo_put_2_32(db, p, (uint32_t)(x * 64 + y));
        }

    int32_t s[2] = { 8, 8 };
    int32_t l[2] = { 16, 16 };
    assert_oracle32(db, s, l, 1);
}

/* Fill candidates equal the (deduped) iterator multiset on one box. */
TEST(fill_parity) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_2d32_fillpar", 4095);

    test_seed_rng(4242);
    for (int i = 0; i < 200; i++) {
        int32_t p[2] = { rand32() % 500, rand32() % 500 };
        geo_put_2_32(db, p, 20000 + i); /* unique values: seal is exact */
    }

    int32_t s[2] = { -50, -50 };
    int32_t l[2] = { 200, 200 };

    size_t nwalk = 0;
    uint32_t *got = walk_collect32(db, s, l, &nwalk);

    rec_set_t *out = rec_set_new();
    ASSERT_EQ(rec_axis_fill_bbox_2_32(db, s, l, out), 0);
    size_t nfill = rec_set_count(out);
    ASSERT_EQ(nfill, nwalk);

    qsort(got, nwalk, sizeof *got, cmp_u32);
    const rec_ref_t *refs = rec_set_at(out);
    for (size_t i = 0; i < nfill; i++)
        ASSERT_EQ((uint32_t)refs[i], got[i]);

    free(got);
    rec_set_free(out);
}

/* Fill argument edges. */
TEST(fill_edges) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_2d32_filledge", 1023);

    int32_t s[2] = { 0, 0 };
    int32_t l[2] = { 10, 10 };
    rec_set_t *out = rec_set_new();

    ASSERT_EQ(rec_axis_fill_bbox_2_32(db, s, l, out), 0);
    ASSERT_EQ(rec_set_count(out), (size_t)0);

    ASSERT_EQ(rec_axis_fill_bbox_2_32(db, s, l, NULL), -1);

    int32_t big[2] = { 2000, 2000 }; /* 4M cells > cap */
    ASSERT_EQ(rec_axis_fill_bbox_2_32(db, s, big, out), -1);

    int32_t cap[2] = { 1024, 1024 }; /* exactly GEO_FILL_MAX_VOL */
    ASSERT_EQ(rec_axis_fill_bbox_2_32(db, s, cap, out), 0);

    rec_set_free(out);
}

/* Box ending exactly at INT32_MAX (no lane wrap): codes adjacent to
 * UINT64_MAX still walk. */
TEST(top_edge_box) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_2d32_topedge", 1023);

    int32_t p[2] = { INT32_MAX - 1, INT32_MAX - 1 };
    geo_put_2_32(db, p, 31337);

    int32_t s[2] = { INT32_MAX - 50, INT32_MAX - 50 };
    int32_t l[2] = { 50, 50 };
    size_t nwalk = 0;
    uint32_t *got = walk_collect32(db, s, l, &nwalk);
    ASSERT_EQ(nwalk, 1);
    ASSERT_EQ(got[0], 31337);
    free(got);

    int32_t slow[2] = { INT32_MIN, INT32_MIN };
    int32_t lfull[2] = { INT32_MAX, INT32_MAX }; /* e = {-1,-1} */
    size_t n2 = 0;
    uint32_t *got2 = walk_collect32(db, slow, lfull, &n2);
    ASSERT_EQ(n2, 0); /* the stored point is above the box end */
    free(got2);
}

/* The dense 32-bit codes live in a different region than the legacy
 * sparse 16-bit codes for the same numerals: configs never alias. */
TEST(config_isolation) {
    int16_t p16[2] = { 5, 7 };
    int32_t p32[2] = { 5, 7 };

    ASSERT(morton_set_2_32(p32) != morton_set_2(p16));
}

int main(void) {
    test_suite_begin("Geo 2D x 32-bit Config Unit Tests");

    RUN_TEST(codec_roundtrip_edges);
    RUN_TEST(codec_roundtrip_random);
    RUN_TEST(codec_known_answers);
    RUN_TEST(point_ops_32);
    RUN_TEST(put_get_basic);
    RUN_TEST(set_replace);
    RUN_TEST(del_paths);
    RUN_TEST(get_multi_chain);
    RUN_TEST(iter_collect_roundtrip);
    RUN_TEST(walk_oracle_mixed);
    RUN_TEST(walk_skip_engages);
    RUN_TEST(fill_parity);
    RUN_TEST(fill_edges);
    RUN_TEST(top_edge_box);
    RUN_TEST(config_isolation);

    return test_suite_end();
}
