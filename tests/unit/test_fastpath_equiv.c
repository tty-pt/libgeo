/**
 * @file test_fastpath_equiv.c
 * @brief Correctness of the per-dimension point_copy_N implementations.
 *
 * Each point_copy_N() is the only implementation for its dimension
 * count (no generic fallback), so these tests validate them directly.
 */

#include "../test_common.h"
#include "../../include/ttypt/point.h"
#include <limits.h>

TEST(fastpath_copy_dim1) {
    int16_t src[1] = {SHRT_MIN};
    int16_t dst[1] = {0};
    point_copy_1(dst, src);
    ASSERT_EQ(dst[0], SHRT_MIN);

    src[0] = SHRT_MAX;
    point_copy_1(dst, src);
    ASSERT_EQ(dst[0], SHRT_MAX);
}

TEST(fastpath_copy_dim2) {
    int16_t src[2] = {SHRT_MIN, SHRT_MAX};
    int16_t dst[2] = {0, 0};
    point_copy_2(dst, src);
    ASSERT_EQ(dst[0], SHRT_MIN);
    ASSERT_EQ(dst[1], SHRT_MAX);

    int16_t z[2] = {-1, 1};
    int16_t zd[2];
    point_copy_2(zd, z);
    ASSERT_EQ(zd[0], -1);
    ASSERT_EQ(zd[1], 1);
}

TEST(fastpath_copy_dim3) {
    int16_t src[3] = {SHRT_MIN, 0, SHRT_MAX};
    int16_t dst[3] = {0, 0, 0};
    point_copy_3(dst, src);
    ASSERT_POINT_EQ(dst, src, 3);
}

TEST(fastpath_copy_dim4) {
    int16_t src[4] = {SHRT_MIN, -1, 0, SHRT_MAX};
    int16_t dst[4] = {0, 0, 0, 0};
    point_copy_4(dst, src);
    ASSERT_POINT_EQ(dst, src, 4);
}

TEST(fastpath_copy_random_all_dims) {
    test_seed_rng(0xC97);
    for (int i = 0; i < 2000; i++) {
        int16_t src[4] = {
            test_rand_coord(), test_rand_coord(),
            test_rand_coord(), test_rand_coord()
        };
        int16_t d1[1], d2[2], d3[3], d4[4];
        point_copy_1(d1, src);
        point_copy_2(d2, src);
        point_copy_3(d3, src);
        point_copy_4(d4, src);
        ASSERT_EQ(d1[0], src[0]);
        ASSERT_POINT_EQ(d2, src, 2);
        ASSERT_POINT_EQ(d3, src, 3);
        ASSERT_POINT_EQ(d4, src, 4);
    }
}

int main(void) {
    test_suite_begin("Fast-path Equivalence Tests");

    RUN_TEST(fastpath_copy_dim1);
    RUN_TEST(fastpath_copy_dim2);
    RUN_TEST(fastpath_copy_dim3);
    RUN_TEST(fastpath_copy_dim4);
    RUN_TEST(fastpath_copy_random_all_dims);

    return test_suite_end();
}
