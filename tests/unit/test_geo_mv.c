/**
 * @file test_geo_mv.c
 * @brief Unit tests for multi-value cells (kernel-form scope B).
 *
 * Contract under QM_SORTED|QM_MULTIVALUE:
 * - geo_put APPENDS a duplicate at the cell (insertion order preserved).
 * - geo_get returns the FIRST value; GEO_MISS when the cell is empty.
 * - geo_get_multi / geo_cell_next iterate ALL values at the cell.
 * - geo_del removes the FIRST value; geo_del_all removes all.
 * - geo_set replaces (del_all + put).
 * - geo_iter / geo_next yield every stored (point, value) pair, including
 *   MV siblings, in morton-discovery order.
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

/* geo_put appends; geo_get returns the first value */
TEST(mv_put_appends_get_returns_first) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_mv_append", 1023);

    int16_t pos[3] = {10, 20, 30};
    geo_put_3(db, pos, 100);
    ASSERT_EQ(geo_get_3(db, pos), 100);

    geo_put_3(db, pos, 200);
    ASSERT_EQ(geo_get_3(db, pos), 100); /* first, not last */
}

/* geo_get_multi iterates all values in insertion order */
TEST(mv_get_multi_yields_all) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_mv_multi", 1023);

    int16_t pos[3] = {1, 2, 3};
    geo_put_3(db, pos, 11);
    geo_put_3(db, pos, 22);
    geo_put_3(db, pos, 33);

    ASSERT_EQ(geo_cell_count_3(db, pos), 3);

    uint32_t cur = geo_get_multi_3(db, pos);
    ASSERT(cur != QM_MISS);

    uint32_t ref;
    ASSERT_EQ(geo_cell_next(&ref, cur), 1);
    ASSERT_EQ(ref, 11);
    ASSERT_EQ(geo_cell_next(&ref, cur), 1);
    ASSERT_EQ(ref, 22);
    ASSERT_EQ(geo_cell_next(&ref, cur), 1);
    ASSERT_EQ(ref, 33);
    ASSERT_EQ(geo_cell_next(&ref, cur), 0); /* exhausted, auto-freed */
}

/* get_multi on an empty cell misses */
TEST(mv_get_multi_empty_misses) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_mv_multi_empty", 1023);

    int16_t pos[3] = {7, 7, 7};
    ASSERT_EQ(geo_cell_count_3(db, pos), 0);
    ASSERT_EQ(geo_get_multi_3(db, pos), QM_MISS);
    ASSERT_EQ(geo_get_3(db, pos), GEO_MISS);
}

/* geo_del removes only the first value */
TEST(mv_del_removes_first) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_mv_del", 1023);

    int16_t pos[3] = {4, 5, 6};
    geo_put_3(db, pos, 11);
    geo_put_3(db, pos, 22);

    geo_del_3(db, pos);
    ASSERT_EQ(geo_get_3(db, pos), 22);
    ASSERT_EQ(geo_cell_count_3(db, pos), 1);

    geo_del_3(db, pos);
    ASSERT_EQ(geo_get_3(db, pos), GEO_MISS);
    ASSERT_EQ(geo_cell_count_3(db, pos), 0);
}

/* geo_del on an empty cell is a safe no-op */
TEST(mv_del_empty_noop) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_mv_del_empty", 1023);

    int16_t pos[3] = {9, 9, 9};
    geo_del_3(db, pos); /* must not crash */
    ASSERT_EQ(geo_get_3(db, pos), GEO_MISS);
}

/* geo_del_all removes every value at the cell */
TEST(mv_del_all_empties_cell) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_mv_del_all", 1023);

    int16_t pos[3] = {3, 3, 3};
    geo_put_3(db, pos, 11);
    geo_put_3(db, pos, 22);
    geo_put_3(db, pos, 33);

    ASSERT_EQ(geo_del_all_3(db, pos), 3);
    ASSERT_EQ(geo_get_3(db, pos), GEO_MISS);
    ASSERT_EQ(geo_cell_count_3(db, pos), 0);
    ASSERT_EQ(geo_del_all_3(db, pos), 0); /* already empty */
}

/* geo_set replaces all values at the cell with one */
TEST(mv_set_replaces) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_mv_set", 1023);

    int16_t pos[3] = {10, 20, 30};
    geo_put_3(db, pos, 100);
    geo_put_3(db, pos, 200);
    ASSERT_EQ(geo_cell_count_3(db, pos), 2);

    geo_set_3(db, pos, 300);
    ASSERT_EQ(geo_get_3(db, pos), 300);
    ASSERT_EQ(geo_cell_count_3(db, pos), 1);
}

/* raw iterator yields MV siblings, each exactly once */
TEST(mv_iter_yields_siblings) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_mv_iter", 1023);

    int16_t a[3] = {5, 5, 5};
    int16_t b[3] = {6, 6, 6};
    geo_put_3(db, a, 11);
    geo_put_3(db, a, 22);
    geo_put_3(db, b, 33);

    int16_t start[3] = {0, 0, 0};
    uint16_t len[3] = {10, 10, 10};
    uint32_t iter = geo_iter_3(db, start, len);

    int16_t p[3];
    uint32_t ref;
    int n11 = 0, n22 = 0, n33 = 0, total = 0;
    while (geo_next(p, &ref, iter)) {
        total++;
        if (ref == 11) { n11++; ASSERT_POINT_EQ(p, a, 3); }
        else if (ref == 22) { n22++; ASSERT_POINT_EQ(p, a, 3); }
        else if (ref == 33) { n33++; ASSERT_POINT_EQ(p, b, 3); }
        else ASSERT(0); /* unexpected value */
    }
    ASSERT_EQ(total, 3);
    ASSERT_EQ(n11, 1);
    ASSERT_EQ(n22, 1);
    ASSERT_EQ(n33, 1);
}

/* iterator results come out in morton-discovery (non-decreasing code) order */
TEST(mv_iter_morton_order) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_mv_order", 1023);

    int16_t pts[8][3] = {
        {9, 1, 4}, {2, 7, 3}, {5, 5, 5}, {0, 0, 1},
        {8, 8, 0}, {1, 9, 2}, {4, 3, 8}, {7, 0, 6},
    };
    for (int i = 0; i < 8; i++)
        geo_put_3(db, pts[i], 100 + i);

    int16_t start[3] = {0, 0, 0};
    uint16_t len[3] = {10, 10, 10};
    uint32_t iter = geo_iter_3(db, start, len);

    int16_t p[3];
    uint32_t ref;
    uint64_t prev = 0;
    int first = 1, count = 0;
    while (geo_next(p, &ref, iter)) {
        uint64_t code = morton_set_3(p);
        if (!first)
            ASSERT_GE(code, prev);
        prev = code;
        first = 0;
        count++;
    }
    ASSERT_EQ(count, 8);
}

/* duplicate identical values are stored as distinct chain entries */
TEST(mv_duplicate_values_kept) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_mv_dupvals", 1023);

    int16_t pos[3] = {2, 2, 2};
    geo_put_3(db, pos, 42);
    geo_put_3(db, pos, 42);

    ASSERT_EQ(geo_cell_count_3(db, pos), 2);
    ASSERT_EQ(geo_del_all_3(db, pos), 2);
    ASSERT_EQ(geo_get_3(db, pos), GEO_MISS);
}

int main(void) {
    test_suite_begin("Geo Multi-Value Cell Unit Tests");

    RUN_TEST(mv_put_appends_get_returns_first);
    RUN_TEST(mv_get_multi_yields_all);
    RUN_TEST(mv_get_multi_empty_misses);
    RUN_TEST(mv_del_removes_first);
    RUN_TEST(mv_del_empty_noop);
    RUN_TEST(mv_del_all_empties_cell);
    RUN_TEST(mv_set_replaces);
    RUN_TEST(mv_iter_yields_siblings);
    RUN_TEST(mv_iter_morton_order);
    RUN_TEST(mv_duplicate_values_kept);

    return test_suite_end();
}
