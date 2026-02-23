/**
 * @file test_point.c
 * @brief Unit tests for point arithmetic and utility functions.
 */

#include "../test_common.h"
#include "../../include/ttypt/point.h"
#include <limits.h>

/* Test point_add */
TEST(point_add_basic) {
    int16_t a[3] = {10, 20, 30};
    int16_t b[3] = {1, 2, 3};
    int16_t result[3];
    
    point_add(result, a, b, 3);
    
    ASSERT_EQ(result[0], 11);
    ASSERT_EQ(result[1], 22);
    ASSERT_EQ(result[2], 33);
}

TEST(point_add_negative) {
    int16_t a[3] = {10, 20, 30};
    int16_t b[3] = {-5, -10, -15};
    int16_t result[3];
    
    point_add(result, a, b, 3);
    
    ASSERT_EQ(result[0], 5);
    ASSERT_EQ(result[1], 10);
    ASSERT_EQ(result[2], 15);
}

TEST(point_add_overflow) {
    int16_t a[3] = {SHRT_MAX, SHRT_MAX, 0};
    int16_t b[3] = {1, 2, 0};
    int16_t result[3];
    
    point_add(result, a, b, 3);
    
    /* Overflow wraps around */
    ASSERT_EQ(result[0], SHRT_MIN);
    ASSERT_EQ(result[1], (int16_t)(SHRT_MIN + 1));
    ASSERT_EQ(result[2], 0);
}

TEST(point_add_in_place) {
    int16_t a[3] = {10, 20, 30};
    int16_t b[3] = {1, 2, 3};
    
    point_add(a, a, b, 3); /* In-place: a = a + b */
    
    ASSERT_EQ(a[0], 11);
    ASSERT_EQ(a[1], 22);
    ASSERT_EQ(a[2], 33);
}

/* Test point_sub */
TEST(point_sub_basic) {
    int16_t a[3] = {100, 200, 300};
    int16_t b[3] = {50, 60, 70};
    int16_t result[3];
    
    point_sub(result, a, b, 3);
    
    ASSERT_EQ(result[0], 50);
    ASSERT_EQ(result[1], 140);
    ASSERT_EQ(result[2], 230);
}

TEST(point_sub_negative_result) {
    int16_t a[3] = {10, 20, 30};
    int16_t b[3] = {20, 30, 40};
    int16_t result[3];
    
    point_sub(result, a, b, 3);
    
    ASSERT_EQ(result[0], -10);
    ASSERT_EQ(result[1], -10);
    ASSERT_EQ(result[2], -10);
}

TEST(point_sub_underflow) {
    int16_t a[3] = {SHRT_MIN, SHRT_MIN, 0};
    int16_t b[3] = {1, 2, 0};
    int16_t result[3];
    
    point_sub(result, a, b, 3);
    
    /* Underflow wraps around */
    ASSERT_EQ(result[0], SHRT_MAX);
    ASSERT_EQ(result[1], (int16_t)(SHRT_MAX - 1));
    ASSERT_EQ(result[2], 0);
}

/* Test point_min */
TEST(point_min_basic) {
    int16_t a[3] = {10, 50, 30};
    int16_t b[3] = {20, 40, 35};
    int16_t result[3];
    
    point_min(result, a, b, 3);
    
    ASSERT_EQ(result[0], 10);
    ASSERT_EQ(result[1], 40);
    ASSERT_EQ(result[2], 30);
}

TEST(point_min_negative) {
    int16_t a[3] = {-10, -50, 0};
    int16_t b[3] = {-20, -40, -5};
    int16_t result[3];
    
    point_min(result, a, b, 3);
    
    ASSERT_EQ(result[0], -20);
    ASSERT_EQ(result[1], -50);
    ASSERT_EQ(result[2], -5);
}

TEST(point_min_equal) {
    int16_t a[3] = {10, 20, 30};
    int16_t b[3] = {10, 20, 30};
    int16_t result[3];
    
    point_min(result, a, b, 3);
    
    ASSERT_POINT_EQ(result, a, 3);
}

/* Test point_max */
TEST(point_max_basic) {
    int16_t a[3] = {10, 50, 30};
    int16_t b[3] = {20, 40, 35};
    int16_t result[3];
    
    point_max(result, a, b, 3);
    
    ASSERT_EQ(result[0], 20);
    ASSERT_EQ(result[1], 50);
    ASSERT_EQ(result[2], 35);
}

TEST(point_max_negative) {
    int16_t a[3] = {-10, -50, 0};
    int16_t b[3] = {-20, -40, -5};
    int16_t result[3];
    
    point_max(result, a, b, 3);
    
    ASSERT_EQ(result[0], -10);
    ASSERT_EQ(result[1], -40);
    ASSERT_EQ(result[2], 0);
}

/* Test point_copy */
TEST(point_copy_basic) {
    int16_t src[3] = {10, 20, 30};
    int16_t dst[3] = {0, 0, 0};
    
    point_copy(dst, src, 3);
    
    ASSERT_POINT_EQ(dst, src, 3);
}

TEST(point_copy_different_dimensions) {
    int16_t src[4] = {10, 20, 30, 40};
    int16_t dst[4] = {0, 0, 0, 0};
    
    point_copy(dst, src, 2);
    
    ASSERT_EQ(dst[0], 10);
    ASSERT_EQ(dst[1], 20);
    ASSERT_EQ(dst[2], 0); /* Unchanged */
    ASSERT_EQ(dst[3], 0); /* Unchanged */
}

/* Test point_vol */
TEST(point_vol_2d) {
    int16_t size[2] = {100, 50};
    
    int32_t vol = point_vol(size, 2);
    
    ASSERT_EQ(vol, 5000);
}

TEST(point_vol_3d) {
    int16_t size[3] = {10, 20, 30};
    
    int32_t vol = point_vol(size, 3);
    
    ASSERT_EQ(vol, 6000);
}

TEST(point_vol_unit_cube) {
    int16_t size[3] = {1, 1, 1};
    
    int32_t vol = point_vol(size, 3);
    
    ASSERT_EQ(vol, 1);
}

TEST(point_vol_large) {
    int16_t size[3] = {100, 100, 100};
    
    int32_t vol = point_vol(size, 3);
    
    ASSERT_EQ(vol, 1000000);
}

TEST(point_vol_zero) {
    int16_t size[3] = {10, 0, 20};
    
    int32_t vol = point_vol(size, 3);
    
    ASSERT_EQ(vol, 0);
}

/* Test point_set */
TEST(point_set_zero) {
    int16_t p[3] = {99, 99, 99};
    
    point_set(p, 0, 3);
    
    ASSERT_EQ(p[0], 0);
    ASSERT_EQ(p[1], 0);
    ASSERT_EQ(p[2], 0);
}

TEST(point_set_uniform) {
    int16_t p[3];
    
    point_set(p, 42, 3);
    
    ASSERT_EQ(p[0], 42);
    ASSERT_EQ(p[1], 42);
    ASSERT_EQ(p[2], 42);
}

TEST(point_set_negative) {
    int16_t p[3];
    
    point_set(p, -100, 3);
    
    ASSERT_EQ(p[0], -100);
    ASSERT_EQ(p[1], -100);
    ASSERT_EQ(p[2], -100);
}

/* Test point_idx */
TEST(point_idx_2d_origin) {
    int16_t start[2] = {0, 0};
    int16_t end[2] = {10, 10};
    int16_t point[2] = {0, 0};
    
    uint64_t idx = point_idx(point, start, end, 2);
    
    ASSERT_EQ(idx, 0);
}

TEST(point_idx_2d_corner) {
    int16_t start[2] = {0, 0};
    int16_t end[2] = {10, 10};
    int16_t point[2] = {9, 9};
    
    uint64_t idx = point_idx(point, start, end, 2);
    
    /* Row-major: idx = y * width + x = 9 * 10 + 9 = 99 */
    ASSERT_EQ(idx, 99);
}

TEST(point_idx_2d_middle) {
    int16_t start[2] = {0, 0};
    int16_t end[2] = {10, 10};
    int16_t point[2] = {3, 5};
    
    uint64_t idx = point_idx(point, start, end, 2);
    
    /* idx = 5 * 10 + 3 = 53 */
    ASSERT_EQ(idx, 53);
}

TEST(point_idx_3d_origin) {
    int16_t start[3] = {0, 0, 0};
    int16_t end[3] = {10, 10, 10};
    int16_t point[3] = {0, 0, 0};
    
    uint64_t idx = point_idx(point, start, end, 3);
    
    ASSERT_EQ(idx, 0);
}

TEST(point_idx_3d_simple) {
    int16_t start[3] = {0, 0, 0};
    int16_t end[3] = {10, 10, 10};
    int16_t point[3] = {1, 0, 0};
    
    uint64_t idx = point_idx(point, start, end, 3);
    
    /* idx = 1 (x varies fastest) */
    ASSERT_EQ(idx, 1);
}

TEST(point_idx_3d_y_varies) {
    int16_t start[3] = {0, 0, 0};
    int16_t end[3] = {10, 10, 10};
    int16_t point[3] = {0, 1, 0};
    
    uint64_t idx = point_idx(point, start, end, 3);
    
    /* idx = 1 * 10 = 10 */
    ASSERT_EQ(idx, 10);
}

TEST(point_idx_3d_z_varies) {
    int16_t start[3] = {0, 0, 0};
    int16_t end[3] = {10, 10, 10};
    int16_t point[3] = {0, 0, 1};
    
    uint64_t idx = point_idx(point, start, end, 3);
    
    /* idx = 1 * 10 * 10 = 100 */
    ASSERT_EQ(idx, 100);
}

TEST(point_idx_offset_box) {
    int16_t start[3] = {10, 20, 30};
    int16_t end[3] = {20, 30, 40};
    int16_t point[3] = {15, 25, 35};
    
    uint64_t idx = point_idx(point, start, end, 3);
    
    /* Normalized: (5, 5, 5) in a 10x10x10 box */
    /* idx = 5 + 5*10 + 5*10*10 = 5 + 50 + 500 = 555 */
    ASSERT_EQ(idx, 555);
}

/* Test dimension variations */
TEST(point_operations_1d) {
    int16_t a[1] = {10};
    int16_t b[1] = {5};
    int16_t result[1];
    
    point_add(result, a, b, 1);
    ASSERT_EQ(result[0], 15);
    
    point_sub(result, a, b, 1);
    ASSERT_EQ(result[0], 5);
    
    point_min(result, a, b, 1);
    ASSERT_EQ(result[0], 5);
    
    point_max(result, a, b, 1);
    ASSERT_EQ(result[0], 10);
}

TEST(point_operations_4d) {
    int16_t a[4] = {10, 20, 30, 40};
    int16_t b[4] = {5, 15, 25, 35};
    int16_t result[4];
    
    point_add(result, a, b, 4);
    int16_t expected_add[4] = {15, 35, 55, 75};
    ASSERT_POINT_EQ(result, expected_add, 4);
    
    point_sub(result, a, b, 4);
    int16_t expected_sub[4] = {5, 5, 5, 5};
    ASSERT_POINT_EQ(result, expected_sub, 4);
}

int main(void) {
    test_suite_begin("Point Utilities Unit Tests");
    
    /* point_add tests */
    RUN_TEST(point_add_basic);
    RUN_TEST(point_add_negative);
    RUN_TEST(point_add_overflow);
    RUN_TEST(point_add_in_place);
    
    /* point_sub tests */
    RUN_TEST(point_sub_basic);
    RUN_TEST(point_sub_negative_result);
    RUN_TEST(point_sub_underflow);
    
    /* point_min tests */
    RUN_TEST(point_min_basic);
    RUN_TEST(point_min_negative);
    RUN_TEST(point_min_equal);
    
    /* point_max tests */
    RUN_TEST(point_max_basic);
    RUN_TEST(point_max_negative);
    
    /* point_copy tests */
    RUN_TEST(point_copy_basic);
    RUN_TEST(point_copy_different_dimensions);
    
    /* point_vol tests */
    RUN_TEST(point_vol_2d);
    RUN_TEST(point_vol_3d);
    RUN_TEST(point_vol_unit_cube);
    RUN_TEST(point_vol_large);
    RUN_TEST(point_vol_zero);
    
    /* point_set tests */
    RUN_TEST(point_set_zero);
    RUN_TEST(point_set_uniform);
    RUN_TEST(point_set_negative);
    
    /* point_idx tests */
    RUN_TEST(point_idx_2d_origin);
    RUN_TEST(point_idx_2d_corner);
    RUN_TEST(point_idx_2d_middle);
    RUN_TEST(point_idx_3d_origin);
    RUN_TEST(point_idx_3d_simple);
    RUN_TEST(point_idx_3d_y_varies);
    RUN_TEST(point_idx_3d_z_varies);
    RUN_TEST(point_idx_offset_box);
    
    /* Dimension variation tests */
    RUN_TEST(point_operations_1d);
    RUN_TEST(point_operations_4d);
    
    return test_suite_end();
}
