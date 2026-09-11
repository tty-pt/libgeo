/**
 * @file test_islet_range.c
 * @brief Unit tests for the Z-interval skip in the box walker.
 *
 * Strategy: a brute-force full-map-scan oracle (independent code path:
 * qmap_iter NULL scan + own box check) is the truth. The walker must agree
 * with it exactly on every configuration, including MV chains and dirty
 * (post-edit) maps. A dense adversarial configuration (every morton-interval
 * address occupied) locks the soundness proof: jumps may skip false
 * positives but never an in-box key. islet_last_scan_count() proves the skip
 * engages (examined < interval width).
 */

#include "../test_common.h"
#include "../../include/ttypt/islet.h"
#include "../../include/ttypt/point.h"
#include "../../include/ttypt/morton.h"
#include "../../include/ttypt/qmap.h"
#include <stdlib.h>
#include <string.h>

static void setup_once(void) {
    static int initialized = 0;
    if (!initialized) {
        islet_init();
        initialized = 1;
    }
}

static int cmp_u32(const void *a, const void *b) {
    uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;
    return x > y ? 1 : (x < y ? -1 : 0);
}

/* Brute-force oracle: full map scan. Returns malloc'd multiset of values
 * whose point lies in [s, s+l] (inclusive end, matching inrange_p), plus
 * the number of stored keys inside the morton interval. Caller frees. */
static uint32_t *brute_collect(uint32_t db, int16_t *s, uint16_t *l, uint8_t dim,
                               size_t *n_out, size_t *interval_out) {
    size_t cap = 64, n = 0, interval = 0;
    uint32_t *vals = malloc(cap * sizeof *vals);
    int16_t e[4], p[4];
    const void *key, *value;

    islet_ops[dim].point_add(e, s, (int16_t *)l);
    uint64_t rmin = islet_ops[dim].morton_set(s);
    uint64_t rmax = islet_ops[dim].morton_set(e);

    uint32_t cur = qmap_iter(db, NULL, 0);
    while (qmap_next(&key, &value, cur)) {
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

/* Raw walker collect into a malloc'd multiset. Caller frees. */
static uint32_t *walk_collect(uint32_t db, int16_t *s, uint16_t *l, uint8_t dim,
                              size_t *n_out) {
    size_t cap = 64, n = 0;
    uint32_t *vals = malloc(cap * sizeof *vals);
    uint32_t iter = islet_ops[dim].iter(db, s, l);
    int16_t p[4];
    uint32_t ref;

    while (islet_next(p, &ref, iter)) {
        if (n == cap) {
            cap *= 2;
            vals = realloc(vals, cap * sizeof *vals);
        }
        vals[n++] = ref;
    }
    *n_out = n;
    return vals;
}

/* Assert walker multiset == oracle multiset, and the skip engaged
 * (examined entries strictly below the interval width) when required. */
static void assert_oracle(uint32_t db, int16_t *s, uint16_t *l, uint8_t dim,
                          int require_skip) {
    size_t nshall = 0, interval = 0, nwalk = 0;
    uint32_t *expected = brute_collect(db, s, l, dim, &nshall, &interval);
    uint32_t *got = walk_collect(db, s, l, dim, &nwalk);
    uint32_t scanned = islet_last_scan_count();

    ASSERT_EQ(nwalk, nshall);
    qsort(expected, nshall, sizeof *expected, cmp_u32);
    qsort(got, nwalk, sizeof *got, cmp_u32);
    for (size_t i = 0; i < nshall; i++)
        ASSERT_EQ(got[i], expected[i]);

    /* Examined (decoded) entries are a subset of the interval keys ... */
    ASSERT(scanned <= interval);
    /* ... and on dense intervals the skip must fire. */
    if (require_skip)
        ASSERT(scanned < interval);

    free(expected);
    free(got);
}

/* Dense adversarial configuration: every address of the morton interval
 * holds a stored key (125 in-box + 324 false positives). A wrong jump
 * skips in-box keys here with near certainty; the exact count locks it. */
TEST(range_dense_interval_exact) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_range_dense", 4095);

    int16_t s[3] = {0, 0, 0};
    uint16_t l[3] = {4, 4, 4};
    int16_t e[3];
    point_add_3(e, s, (int16_t *)l);
    uint64_t rmin = morton_set_3(s);
    uint64_t rmax = morton_set_3(e);

    int in = 0, fp = 0;
    for (int16_t x = -16; x <= 16; x++)
        for (int16_t y = -16; y <= 16; y++)
            for (int16_t z = -16; z <= 16; z++) {
                int16_t p[3] = {x, y, z};
                int inside = (x >= 0 && x <= 4 &&
                              y >= 0 && y <= 4 &&
                              z >= 0 && z <= 4);
                uint64_t code = morton_set_3(p);
                if (inside) {
                    islet_put_3(db, p, 1000000u + in);
                    in++;
                } else if (code >= rmin && code <= rmax) {
                    islet_put_3(db, p, 2000000u + fp);
                    fp++;
                }
            }
    ASSERT_EQ(in, 125);
    ASSERT(fp > 0);

    /* Exactness is the lock here (this layout caught an unsound BIGMIN);
     * engagement is proven by range_far_slab_engages below. */
    assert_oracle(db, s, l, 3, 0);
}

/* Far-miss slab: false positives well clear of the box sit in large
 * box-disjoint cubes, so jumps fire in bulk and the decoded count drops
 * strictly below the interval width. Layout verified by simulation
 * (729 in-box + in-interval slab, ~300 examinations saved). */
TEST(range_far_slab_engages) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_range_slab", 8191);

    int16_t s[3] = {0, 0, 0};
    uint16_t l[3] = {8, 8, 8};

    int in = 0;
    for (int16_t x = 0; x <= 8; x++)
        for (int16_t y = 0; y <= 8; y++)
            for (int16_t z = 0; z <= 8; z++) {
                int16_t p[3] = {x, y, z};
                islet_put_3(db, p, 1000000u + in);
                in++;
            }
    ASSERT_EQ(in, 729);

    for (int16_t x = 0; x <= 8; x++)
        for (int16_t y = 0; y <= 8; y++)
            for (int16_t z = 12; z <= 20; z++) {
                int16_t p[3] = {x, y, z};
                islet_put_3(db, p, 2000000u + x * 100 + y * 10 + z);
            }

    assert_oracle(db, s, l, 3, 1);
}

/* Empty box amid a populated interval: zero results, still consistent. */
TEST(range_empty_box_consistent) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_range_empty", 1023);

    islet_put_3(db, (int16_t[]){0, 0, 0}, 1);
    islet_put_3(db, (int16_t[]){50, 50, 50}, 2);
    islet_put_3(db, (int16_t[]){-50, -50, -50}, 3);

    int16_t s[3] = {10, 10, 10};
    uint16_t l[3] = {5, 5, 5};
    assert_oracle(db, s, l, 3, 0);
}

/* MV chains on both sides of a jump: outside chains contribute nothing,
 * inside chains yield every sibling exactly once. */
TEST(range_mv_chains_across_jumps) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_range_mvjump", 4095);

    int16_t s[3] = {0, 0, 0};
    uint16_t l[3] = {4, 4, 4};
    int16_t e[3];
    point_add_3(e, s, (int16_t *)l);
    uint64_t rmin = morton_set_3(s);
    uint64_t rmax = morton_set_3(e);

    /* Inside cell with 3 siblings + outside false-positive cell with 3. */
    int16_t inside[3] = {1, 1, 1};
    islet_put_3(db, inside, 11);
    islet_put_3(db, inside, 22);
    islet_put_3(db, inside, 33);

    int found_fp = 0;
    for (int16_t x = -16; x <= 16 && !found_fp; x++)
        for (int16_t y = -16; y <= 16 && !found_fp; y++)
            for (int16_t z = -16; z <= 16 && !found_fp; z++) {
                int16_t p[3] = {x, y, z};
                int is_inside = (x >= 0 && x <= 4 &&
                                 y >= 0 && y <= 4 &&
                                 z >= 0 && z <= 4);
                if (is_inside)
                    continue;
                uint64_t code = morton_set_3(p);
                if (code >= rmin && code <= rmax) {
                    islet_put_3(db, p, 101);
                    islet_put_3(db, p, 102);
                    islet_put_3(db, p, 103);
                    found_fp = 1;
                }
            }
    ASSERT(found_fp);

    assert_oracle(db, s, l, 3, 0);

    /* And the inside chain reads back whole via both paths. */
    ASSERT_EQ(islet_cell_count_3(db, inside), 3);
    uint32_t cur = islet_get_multi_3(db, inside);
    uint32_t ref;
    ASSERT_EQ(islet_cell_next(&ref, cur), 1);
    ASSERT_EQ(ref, 11);
    ASSERT_EQ(islet_cell_next(&ref, cur), 1);
    ASSERT_EQ(ref, 22);
    ASSERT_EQ(islet_cell_next(&ref, cur), 1);
    ASSERT_EQ(ref, 33);
    ASSERT_EQ(islet_cell_next(&ref, cur), 0);
}

/* GE landing exactly on a duplicate chain head yields every sibling. */
TEST(range_landing_on_dup_chain) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_range_landing", 1023);

    int16_t cell[3] = {7, 7, 7};
    islet_put_3(db, cell, 1);
    islet_put_3(db, cell, 2);
    islet_put_3(db, cell, 3);
    islet_put_3(db, cell, 4);

    /* Box corner == the cell: the walk GE-seeks straight onto the chain. */
    int16_t s[3] = {7, 7, 7};
    uint16_t l[3] = {2, 2, 2};
    assert_oracle(db, s, l, 3, 0);
}

/* Edits between queries (dirty sorted-index rebuild) stay consistent. */
TEST(range_after_edits) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_range_edits", 4095);

    int16_t s[3] = {0, 0, 0};
    uint16_t l[3] = {8, 8, 8};

    test_seed_rng(5150);
    for (int i = 0; i < 150; i++) {
        int16_t p[3] = {
            test_rand_coord_range(0, 24),
            test_rand_coord_range(0, 24),
            test_rand_coord_range(0, 24),
        };
        islet_put_3(db, p, (uint32_t)i);
    }
    assert_oracle(db, s, l, 3, 0);

    /* Churn: deletes + overwrites + new points, then re-verify. */
    for (int i = 0; i < 50; i++) {
        int16_t p[3] = {
            test_rand_coord_range(0, 24),
            test_rand_coord_range(0, 24),
            test_rand_coord_range(0, 24),
        };
        islet_del_all_3(db, p);
    }
    for (int i = 0; i < 50; i++) {
        int16_t p[3] = {
            test_rand_coord_range(0, 24),
            test_rand_coord_range(0, 24),
            test_rand_coord_range(0, 24),
        };
        islet_put_3(db, p, 1000 + (uint32_t)i);
    }
    assert_oracle(db, s, l, 3, 0);
}

/* 1D and 2D boxes obey the same oracle. */
TEST(range_dims_1d_2d) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_range_dims", 4095);

    test_seed_rng(808);
    for (int i = 0; i < 100; i++) {
        int16_t p[2] = {
            test_rand_coord_range(-20, 20),
            test_rand_coord_range(-20, 20),
        };
        islet_put_2(db, p, (uint32_t)i);
    }
    int16_t s2[2] = {-10, -10};
    uint16_t l2[2] = {15, 15};
    assert_oracle(db, s2, l2, 2, 0);

    uint32_t db1 = islet_open(NULL, "test_range_dims1", 4095);
    for (int i = 0; i < 100; i++) {
        int16_t p[1] = {test_rand_coord_range(-20, 20)};
        islet_put_1(db1, p, (uint32_t)i);
    }
    int16_t s1[1] = {-10};
    uint16_t l1[1] = {15};
    assert_oracle(db1, s1, l1, 1, 0);
}

int main(void) {
    test_suite_begin("Islet Range (Z-interval skip) Unit Tests");

    RUN_TEST(range_dense_interval_exact);
    RUN_TEST(range_far_slab_engages);
    RUN_TEST(range_empty_box_consistent);
    RUN_TEST(range_mv_chains_across_jumps);
    RUN_TEST(range_landing_on_dup_chain);
    RUN_TEST(range_after_edits);
    RUN_TEST(range_dims_1d_2d);

    return test_suite_end();
}
