/**
 * @file test_islet_core.c
 * @brief Unit tests for core islet API functions.
 */

#include "../test_common.h"
#include "../../include/ttypt/islet.h"
#include "../../include/ttypt/point.h"
#include <string.h>
#include <unistd.h>
#include <limits.h>

/* Initialize islet once for all tests */
static void setup_once(void) {
    static int initialized = 0;
    if (!initialized) {
        islet_init();
        initialized = 1;
    }
}

/* Test islet_init */
TEST(islet_init_basic) {
    setup_once();
    /* No crash = success */
    ASSERT(1);
}

/* Test islet_open - memory-only database */
TEST(islet_open_memory) {
    setup_once();
    
    uint32_t db = islet_open(NULL, "test_db_mem", 1023);
    
    ASSERT_GT(db, 0); /* Valid handle */
}

/* Test islet_open - file-backed database */
TEST(islet_open_file) {
    setup_once();
    
    const char *filename = "/tmp/test_islet_file.db";
    unlink(filename); /* Clean up any previous test */
    
    uint32_t db = islet_open((char*)filename, "test_db_file", 1023);
    
    ASSERT_GT(db, 0);
}

/* Test islet_put and islet_get */
TEST(islet_put_get_basic) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_db1", 1023);
    
    int16_t pos[3] = {10, 20, 30};
    uint32_t value = 42;
    
    islet_put_3(db, pos, value);
    uint32_t retrieved = islet_get_3(db, pos);
    
    ASSERT_EQ(retrieved, value);
}

TEST(islet_put_get_multiple) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_db2", 1023);
    
    int16_t pos1[3] = {10, 20, 30};
    int16_t pos2[3] = {100, 200, 300};
    int16_t pos3[3] = {-50, -100, -150};
    
    islet_put_3(db, pos1, 1);
    islet_put_3(db, pos2, 2);
    islet_put_3(db, pos3, 3);
    
    ASSERT_EQ(islet_get_3(db, pos1), 1);
    ASSERT_EQ(islet_get_3(db, pos2), 2);
    ASSERT_EQ(islet_get_3(db, pos3), 3);
}

TEST(islet_get_missing) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_db3", 1023);
    
    int16_t pos[3] = {10, 20, 30};
    
    uint32_t value = islet_get_3(db, pos);
    
    /* islet_get returns ISLET_MISS for missing entries */
    ASSERT_EQ(value, ISLET_MISS);
}

TEST(islet_put_overwrite) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_db4", 1023);

    int16_t pos[3] = {10, 20, 30};

    islet_put_3(db, pos, 100);
    ASSERT_EQ(islet_get_3(db, pos), 100);

    /* islet_put appends under multi-value cells: first value still wins */
    islet_put_3(db, pos, 200);
    ASSERT_EQ(islet_get_3(db, pos), 100);
    ASSERT_EQ(islet_cell_count_3(db, pos), 2);

    /* islet_set replaces: the cell holds exactly the new value */
    islet_set_3(db, pos, 200);
    ASSERT_EQ(islet_get_3(db, pos), 200);
    ASSERT_EQ(islet_cell_count_3(db, pos), 1);
}

/* Test islet_del */
TEST(islet_del_existing) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_db5", 1023);
    
    int16_t pos[3] = {10, 20, 30};
    
    islet_put_3(db, pos, 42);
    ASSERT_EQ(islet_get_3(db, pos), 42);
    
    islet_del_3(db, pos);
    ASSERT_EQ(islet_get_3(db, pos), QM_MISS);
}

TEST(islet_del_nonexistent) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_db6", 1023);
    
    int16_t pos[3] = {10, 20, 30};
    
    /* Deleting non-existent entry should not crash */
    islet_del_3(db, pos);
    ASSERT_EQ(islet_get_3(db, pos), ISLET_MISS);
}

/* Test islet_iter and islet_next - basic iteration */
TEST(islet_iter_empty) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_db7", 1023);
    
    int16_t start[3] = {0, 0, 0};
    uint16_t len[3] = {10, 10, 10};
    
    uint32_t iter = islet_iter_3(db, start, len);
    
    int16_t pos[3];
    uint32_t value;
    int result = islet_next(pos, &value, iter);
    
    ASSERT_EQ(result, 0); /* No results */
}

TEST(islet_iter_single_point) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_db8", 1023);
    
    int16_t data_pos[3] = {5, 5, 5};
    islet_put_3(db, data_pos, 99);
    
    int16_t start[3] = {0, 0, 0};
    uint16_t len[3] = {10, 10, 10};
    
    uint32_t iter = islet_iter_3(db, start, len);
    
    int16_t pos[3];
    uint32_t value;
    int result = islet_next(pos, &value, iter);
    
    ASSERT_EQ(result, 1);
    ASSERT_POINT_EQ(pos, data_pos, 3);
    ASSERT_EQ(value, 99);
    
    /* Should be exhausted now */
    result = islet_next(pos, &value, iter);
    ASSERT_EQ(result, 0);
}

TEST(islet_iter_multiple_points) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_db9", 1023);
    
    /* Insert 5 points */
    int16_t points[5][3] = {
        {1, 1, 1},
        {2, 2, 2},
        {3, 3, 3},
        {4, 4, 4},
        {5, 5, 5}
    };
    
    for (int i = 0; i < 5; i++) {
        islet_put_3(db, points[i], i + 100);
    }
    
    /* Query box that contains all points */
    int16_t start[3] = {0, 0, 0};
    uint16_t len[3] = {10, 10, 10};
    
    uint32_t iter = islet_iter_3(db, start, len);
    
    /* Collect all results */
    int count = 0;
    int16_t pos[3];
    uint32_t value;
    
    while (islet_next(pos, &value, iter)) {
        count++;
        ASSERT_GE(value, 100);
        ASSERT_LE(value, 104);
    }
    
    ASSERT_EQ(count, 5);
}

TEST(islet_iter_partial_overlap) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_db10", 1023);
    
    /* Insert points, some inside and some outside query box */
    islet_put_3(db, (int16_t[]){5, 5, 5}, 1);   /* Inside */
    islet_put_3(db, (int16_t[]){15, 15, 15}, 2); /* Outside */
    islet_put_3(db, (int16_t[]){8, 8, 8}, 3);   /* Inside */
    
    /* Query box [0,0,0] to [10,10,10] */
    int16_t start[3] = {0, 0, 0};
    uint16_t len[3] = {10, 10, 10};
    
    uint32_t iter = islet_iter_3(db, start, len);
    
    int count = 0;
    int16_t pos[3];
    uint32_t value;
    
    while (islet_next(pos, &value, iter)) {
        count++;
        /* Should only get values 1 and 3 */
        ASSERT(value == 1 || value == 3);
    }
    
    ASSERT_EQ(count, 2);
}

/* Test negative coordinates */
TEST(islet_negative_coordinates) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_db11", 1023);
    
    int16_t pos[3] = {-100, -200, -300};
    
    islet_put_3(db, pos, 42);
    ASSERT_EQ(islet_get_3(db, pos), 42);
}

/* Test boundary coordinates */
TEST(islet_boundary_coordinates) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_db12", 1023);
    
    int16_t min_pos[3] = {SHRT_MIN, SHRT_MIN, SHRT_MIN};
    int16_t max_pos[3] = {SHRT_MAX, SHRT_MAX, SHRT_MAX};
    
    islet_put_3(db, min_pos, 1);
    islet_put_3(db, max_pos, 2);
    
    ASSERT_EQ(islet_get_3(db, min_pos), 1);
    ASSERT_EQ(islet_get_3(db, max_pos), 2);
}

/* Test 2D operations - DISABLED: libislet is optimized for 3D only */
/* The FAST_MORTON implementation always uses 3 dimensions */
/*
TEST(islet_2d_operations) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_db13", 1023);
    
    int16_t pos[2] = {10, 20};
    
    islet_put_2(db, pos, 42);
    ASSERT_EQ(islet_get_2(db, pos), 42);
    
    //  Iterate
    int16_t start[2] = {0, 0};
    uint16_t len[2] = {100, 100};
    uint32_t iter = islet_iter_2(db, start, len);
    
    int16_t found_pos[3];
    uint32_t value;
    int result = islet_next(found_pos, &value, iter);
    
    ASSERT_EQ(result, 1);
    ASSERT_EQ(found_pos[0], 10);
    ASSERT_EQ(found_pos[1], 20);
    ASSERT_EQ(value, 42);
}
*/

/* Test iteration boundary conditions */
TEST(islet_iter_exact_bounds) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_db14", 1023);
    
    /* Point exactly on query box boundaries */
    int16_t boundary_point[3] = {10, 10, 10};
    islet_put_3(db, boundary_point, 99);
    
    /* Query box [10,10,10] to [11,11,11] - should include the point */
    int16_t start[3] = {10, 10, 10};
    uint16_t len[3] = {1, 1, 1};
    
    uint32_t iter = islet_iter_3(db, start, len);
    
    int16_t pos[3];
    uint32_t value;
    int result = islet_next(pos, &value, iter);
    
    ASSERT_EQ(result, 1);
    ASSERT_POINT_EQ(pos, boundary_point, 3);
}

/* Test value range */
TEST(islet_value_range) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_db15", 1023);
    
    int16_t pos1[3] = {1, 1, 1};
    int16_t pos2[3] = {2, 2, 2};
    int16_t pos3[3] = {3, 3, 3};
    
    islet_put_3(db, pos1, 0);           /* Minimum value */
    islet_put_3(db, pos2, UINT32_MAX);  /* Maximum value */
    islet_put_3(db, pos3, 12345678);    /* Arbitrary value */
    
    ASSERT_EQ(islet_get_3(db, pos1), 0);
    ASSERT_EQ(islet_get_3(db, pos2), UINT32_MAX);
    ASSERT_EQ(islet_get_3(db, pos3), 12345678);
}

int main(void) {
    test_suite_begin("Islet Core API Unit Tests");
    
    RUN_TEST(islet_init_basic);
    RUN_TEST(islet_open_memory);
    RUN_TEST(islet_open_file);
    RUN_TEST(islet_put_get_basic);
    RUN_TEST(islet_put_get_multiple);
    RUN_TEST(islet_get_missing);
    RUN_TEST(islet_put_overwrite);
    RUN_TEST(islet_del_existing);
    RUN_TEST(islet_del_nonexistent);
    RUN_TEST(islet_iter_empty);
    RUN_TEST(islet_iter_single_point);
    RUN_TEST(islet_iter_multiple_points);
    RUN_TEST(islet_iter_partial_overlap);
    RUN_TEST(islet_negative_coordinates);
    RUN_TEST(islet_boundary_coordinates);
    /* islet_2d_operations - DISABLED: libislet is 3D-only with FAST_MORTON */
    RUN_TEST(islet_iter_exact_bounds);
    RUN_TEST(islet_value_range);
    
    return test_suite_end();
}
