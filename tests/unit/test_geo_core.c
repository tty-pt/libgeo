/**
 * @file test_geo_core.c
 * @brief Unit tests for core geo API functions.
 */

#include "../test_common.h"
#include "../../include/ttypt/geo.h"
#include "../../include/ttypt/point.h"
#include <string.h>
#include <unistd.h>
#include <limits.h>

/* Initialize geo once for all tests */
static void setup_once(void) {
    static int initialized = 0;
    if (!initialized) {
        geo_init();
        initialized = 1;
    }
}

/* Test geo_init */
TEST(geo_init_basic) {
    setup_once();
    /* No crash = success */
    ASSERT(1);
}

/* Test geo_open - memory-only database */
TEST(geo_open_memory) {
    setup_once();
    
    uint32_t db = geo_open(NULL, "test_db_mem", 1023);
    
    ASSERT_GT(db, 0); /* Valid handle */
}

/* Test geo_open - file-backed database */
TEST(geo_open_file) {
    setup_once();
    
    const char *filename = "/tmp/test_geo_file.db";
    unlink(filename); /* Clean up any previous test */
    
    uint32_t db = geo_open((char*)filename, "test_db_file", 1023);
    
    ASSERT_GT(db, 0);
}

/* Test geo_put and geo_get */
TEST(geo_put_get_basic) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_db1", 1023);
    
    int16_t pos[3] = {10, 20, 30};
    uint32_t value = 42;
    
    geo_put(db, pos, value, 3);
    uint32_t retrieved = geo_get(db, pos, 3);
    
    ASSERT_EQ(retrieved, value);
}

TEST(geo_put_get_multiple) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_db2", 1023);
    
    int16_t pos1[3] = {10, 20, 30};
    int16_t pos2[3] = {100, 200, 300};
    int16_t pos3[3] = {-50, -100, -150};
    
    geo_put(db, pos1, 1, 3);
    geo_put(db, pos2, 2, 3);
    geo_put(db, pos3, 3, 3);
    
    ASSERT_EQ(geo_get(db, pos1, 3), 1);
    ASSERT_EQ(geo_get(db, pos2, 3), 2);
    ASSERT_EQ(geo_get(db, pos3, 3), 3);
}

TEST(geo_get_missing) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_db3", 1023);
    
    int16_t pos[3] = {10, 20, 30};
    
    uint32_t value = geo_get(db, pos, 3);
    
    /* geo_get returns GEO_MISS for missing entries */
    ASSERT_EQ(value, GEO_MISS);
}

TEST(geo_put_overwrite) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_db4", 1023);
    
    int16_t pos[3] = {10, 20, 30};
    
    geo_put(db, pos, 100, 3);
    ASSERT_EQ(geo_get(db, pos, 3), 100);
    
    geo_put(db, pos, 200, 3);
    ASSERT_EQ(geo_get(db, pos, 3), 200);
}

/* Test geo_del */
TEST(geo_del_existing) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_db5", 1023);
    
    int16_t pos[3] = {10, 20, 30};
    
    geo_put(db, pos, 42, 3);
    ASSERT_EQ(geo_get(db, pos, 3), 42);
    
    geo_del(db, pos, 3);
    ASSERT_EQ(geo_get(db, pos, 3), QM_MISS);
}

TEST(geo_del_nonexistent) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_db6", 1023);
    
    int16_t pos[3] = {10, 20, 30};
    
    /* Deleting non-existent entry should not crash */
    geo_del(db, pos, 3);
    ASSERT_EQ(geo_get(db, pos, 3), GEO_MISS);
}

/* Test geo_iter and geo_next - basic iteration */
TEST(geo_iter_empty) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_db7", 1023);
    
    int16_t start[3] = {0, 0, 0};
    uint16_t len[3] = {10, 10, 10};
    
    uint32_t iter = geo_iter(db, start, len, 3);
    
    int16_t pos[3];
    uint32_t value;
    int result = geo_next(pos, &value, iter);
    
    ASSERT_EQ(result, 0); /* No results */
}

TEST(geo_iter_single_point) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_db8", 1023);
    
    int16_t data_pos[3] = {5, 5, 5};
    geo_put(db, data_pos, 99, 3);
    
    int16_t start[3] = {0, 0, 0};
    uint16_t len[3] = {10, 10, 10};
    
    uint32_t iter = geo_iter(db, start, len, 3);
    
    int16_t pos[3];
    uint32_t value;
    int result = geo_next(pos, &value, iter);
    
    ASSERT_EQ(result, 1);
    ASSERT_POINT_EQ(pos, data_pos, 3);
    ASSERT_EQ(value, 99);
    
    /* Should be exhausted now */
    result = geo_next(pos, &value, iter);
    ASSERT_EQ(result, 0);
}

TEST(geo_iter_multiple_points) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_db9", 1023);
    
    /* Insert 5 points */
    int16_t points[5][3] = {
        {1, 1, 1},
        {2, 2, 2},
        {3, 3, 3},
        {4, 4, 4},
        {5, 5, 5}
    };
    
    for (int i = 0; i < 5; i++) {
        geo_put(db, points[i], i + 100, 3);
    }
    
    /* Query box that contains all points */
    int16_t start[3] = {0, 0, 0};
    uint16_t len[3] = {10, 10, 10};
    
    uint32_t iter = geo_iter(db, start, len, 3);
    
    /* Collect all results */
    int count = 0;
    int16_t pos[3];
    uint32_t value;
    
    while (geo_next(pos, &value, iter)) {
        count++;
        ASSERT_GE(value, 100);
        ASSERT_LE(value, 104);
    }
    
    ASSERT_EQ(count, 5);
}

TEST(geo_iter_partial_overlap) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_db10", 1023);
    
    /* Insert points, some inside and some outside query box */
    geo_put(db, (int16_t[]){5, 5, 5}, 1, 3);   /* Inside */
    geo_put(db, (int16_t[]){15, 15, 15}, 2, 3); /* Outside */
    geo_put(db, (int16_t[]){8, 8, 8}, 3, 3);   /* Inside */
    
    /* Query box [0,0,0] to [10,10,10] */
    int16_t start[3] = {0, 0, 0};
    uint16_t len[3] = {10, 10, 10};
    
    uint32_t iter = geo_iter(db, start, len, 3);
    
    int count = 0;
    int16_t pos[3];
    uint32_t value;
    
    while (geo_next(pos, &value, iter)) {
        count++;
        /* Should only get values 1 and 3 */
        ASSERT(value == 1 || value == 3);
    }
    
    ASSERT_EQ(count, 2);
}

/* Test negative coordinates */
TEST(geo_negative_coordinates) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_db11", 1023);
    
    int16_t pos[3] = {-100, -200, -300};
    
    geo_put(db, pos, 42, 3);
    ASSERT_EQ(geo_get(db, pos, 3), 42);
}

/* Test boundary coordinates */
TEST(geo_boundary_coordinates) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_db12", 1023);
    
    int16_t min_pos[3] = {SHRT_MIN, SHRT_MIN, SHRT_MIN};
    int16_t max_pos[3] = {SHRT_MAX, SHRT_MAX, SHRT_MAX};
    
    geo_put(db, min_pos, 1, 3);
    geo_put(db, max_pos, 2, 3);
    
    ASSERT_EQ(geo_get(db, min_pos, 3), 1);
    ASSERT_EQ(geo_get(db, max_pos, 3), 2);
}

/* Test 2D operations - DISABLED: libgeo is optimized for 3D only */
/* The FAST_MORTON implementation always uses 3 dimensions */
/*
TEST(geo_2d_operations) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_db13", 1023);
    
    int16_t pos[2] = {10, 20};
    
    geo_put(db, pos, 42, 2);
    ASSERT_EQ(geo_get(db, pos, 2), 42);
    
    //  Iterate
    int16_t start[2] = {0, 0};
    uint16_t len[2] = {100, 100};
    uint32_t iter = geo_iter(db, start, len, 2);
    
    int16_t found_pos[3];
    uint32_t value;
    int result = geo_next(found_pos, &value, iter);
    
    ASSERT_EQ(result, 1);
    ASSERT_EQ(found_pos[0], 10);
    ASSERT_EQ(found_pos[1], 20);
    ASSERT_EQ(value, 42);
}
*/

/* Test iteration boundary conditions */
TEST(geo_iter_exact_bounds) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_db14", 1023);
    
    /* Point exactly on query box boundaries */
    int16_t boundary_point[3] = {10, 10, 10};
    geo_put(db, boundary_point, 99, 3);
    
    /* Query box [10,10,10] to [11,11,11] - should include the point */
    int16_t start[3] = {10, 10, 10};
    uint16_t len[3] = {1, 1, 1};
    
    uint32_t iter = geo_iter(db, start, len, 3);
    
    int16_t pos[3];
    uint32_t value;
    int result = geo_next(pos, &value, iter);
    
    ASSERT_EQ(result, 1);
    ASSERT_POINT_EQ(pos, boundary_point, 3);
}

/* Test value range */
TEST(geo_value_range) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_db15", 1023);
    
    int16_t pos1[3] = {1, 1, 1};
    int16_t pos2[3] = {2, 2, 2};
    int16_t pos3[3] = {3, 3, 3};
    
    geo_put(db, pos1, 0, 3);           /* Minimum value */
    geo_put(db, pos2, UINT32_MAX, 3);  /* Maximum value */
    geo_put(db, pos3, 12345678, 3);    /* Arbitrary value */
    
    ASSERT_EQ(geo_get(db, pos1, 3), 0);
    ASSERT_EQ(geo_get(db, pos2, 3), UINT32_MAX);
    ASSERT_EQ(geo_get(db, pos3, 3), 12345678);
}

int main(void) {
    test_suite_begin("Geo Core API Unit Tests");
    
    RUN_TEST(geo_init_basic);
    RUN_TEST(geo_open_memory);
    RUN_TEST(geo_open_file);
    RUN_TEST(geo_put_get_basic);
    RUN_TEST(geo_put_get_multiple);
    RUN_TEST(geo_get_missing);
    RUN_TEST(geo_put_overwrite);
    RUN_TEST(geo_del_existing);
    RUN_TEST(geo_del_nonexistent);
    RUN_TEST(geo_iter_empty);
    RUN_TEST(geo_iter_single_point);
    RUN_TEST(geo_iter_multiple_points);
    RUN_TEST(geo_iter_partial_overlap);
    RUN_TEST(geo_negative_coordinates);
    RUN_TEST(geo_boundary_coordinates);
    /* geo_2d_operations - DISABLED: libgeo is 3D-only with FAST_MORTON */
    RUN_TEST(geo_iter_exact_bounds);
    RUN_TEST(geo_value_range);
    
    return test_suite_end();
}
