/**
 * @file test_morton_4d.c
 * @brief Unit tests for 4D Morton code encoding and decoding.
 *
 * 4D uses dense stride-4 packing over the full 64-bit key space
 * (4 x 16 bits = 64). 1D/2D/3D codes are bit-identical to v0.5.0.
 */

#include "../test_common.h"
#include "../../include/ttypt/morton.h"
#include "../../include/ttypt/point.h"
#include <limits.h>

TEST(morton4_round_trip_origin) {
    int16_t p[4] = {0, 0, 0, 0};
    int16_t decoded[4];

    uint64_t code = morton_set(p, 4);
    morton_get(decoded, code, 4);

    ASSERT_POINT_EQ(p, decoded, 4);
}

TEST(morton4_round_trip_positive) {
    int16_t p[4] = {100, 200, 300, 400};
    int16_t decoded[4];

    uint64_t code = morton_set(p, 4);
    morton_get(decoded, code, 4);

    ASSERT_POINT_EQ(p, decoded, 4);
}

TEST(morton4_round_trip_negative) {
    int16_t p[4] = {-100, -200, -300, -400};
    int16_t decoded[4];

    uint64_t code = morton_set(p, 4);
    morton_get(decoded, code, 4);

    ASSERT_POINT_EQ(p, decoded, 4);
}

TEST(morton4_round_trip_mixed) {
    int16_t p[4] = {-100, 0, 200, -300};
    int16_t decoded[4];

    uint64_t code = morton_set(p, 4);
    morton_get(decoded, code, 4);

    ASSERT_POINT_EQ(p, decoded, 4);
}

TEST(morton4_round_trip_min_values) {
    int16_t p[4] = {SHRT_MIN, SHRT_MIN, SHRT_MIN, SHRT_MIN};
    int16_t decoded[4];

    uint64_t code = morton_set(p, 4);
    morton_get(decoded, code, 4);

    ASSERT_POINT_EQ(p, decoded, 4);
}

TEST(morton4_round_trip_max_values) {
    int16_t p[4] = {SHRT_MAX, SHRT_MAX, SHRT_MAX, SHRT_MAX};
    int16_t decoded[4];

    uint64_t code = morton_set(p, 4);
    morton_get(decoded, code, 4);

    ASSERT_POINT_EQ(p, decoded, 4);
}

TEST(morton4_round_trip_boundaries) {
    int16_t test_cases[][4] = {
        {SHRT_MIN, 0, 0, 0},
        {0, SHRT_MIN, 0, 0},
        {0, 0, SHRT_MIN, 0},
        {0, 0, 0, SHRT_MIN},
        {SHRT_MAX, 0, 0, 0},
        {0, SHRT_MAX, 0, 0},
        {0, 0, SHRT_MAX, 0},
        {0, 0, 0, SHRT_MAX},
        {SHRT_MIN, SHRT_MAX, 0, 0},
        {0, 0, SHRT_MIN, SHRT_MAX},
        {SHRT_MIN, SHRT_MIN, SHRT_MAX, SHRT_MAX},
        {-1, -1, -1, -1},
        {1, 2, 3, 4},
    };
    int num_cases = sizeof(test_cases) / sizeof(test_cases[0]);

    for (int i = 0; i < num_cases; i++) {
        int16_t decoded[4];
        uint64_t code = morton_set(test_cases[i], 4);
        morton_get(decoded, code, 4);
        ASSERT_POINT_EQ(test_cases[i], decoded, 4);
    }
}

/* Minimum corner maps to code 0 in every dimension */
TEST(morton4_min_code_is_zero) {
    int16_t p1[1] = {SHRT_MIN};
    int16_t p2[2] = {SHRT_MIN, SHRT_MIN};
    int16_t p3[3] = {SHRT_MIN, SHRT_MIN, SHRT_MIN};
    int16_t p4[4] = {SHRT_MIN, SHRT_MIN, SHRT_MIN, SHRT_MIN};

    ASSERT_EQ(morton_set(p1, 1), 0u);
    ASSERT_EQ(morton_set(p2, 2), 0u);
    ASSERT_EQ(morton_set(p3, 3), 0u);
    ASSERT_EQ(morton_set(p4, 4), 0u);
}

/* Maximum 4D corner fills all 64 key bits (dense packing, no reserve) */
TEST(morton4_max_code_is_u64max) {
    int16_t p[4] = {SHRT_MAX, SHRT_MAX, SHRT_MAX, SHRT_MAX};

    ASSERT_EQ(morton_set(p, 4), 0xFFFFFFFFFFFFFFFFULL);
}

/* Determinism: same input always yields the same code */
TEST(morton4_determinism) {
    int16_t p[4] = {1234, -5678, 9012, -3456};

    uint64_t c1 = morton_set(p, 4);
    uint64_t c2 = morton_set(p, 4);

    ASSERT_EQ(c1, c2);
}

/* Uniqueness over a small 4D grid, with round-trip on each point */
TEST(morton4_grid_uniqueness) {
    uint64_t seen[625];
    int n = 0;

    for (int16_t x = 0; x < 5; x++) {
        for (int16_t y = 0; y < 5; y++) {
            for (int16_t z = 0; z < 5; z++) {
                for (int16_t w = 0; w < 5; w++) {
                    int16_t p[4] = {x, y, z, w};
                    uint64_t code = morton_set(p, 4);

                    for (int i = 0; i < n; i++)
                        ASSERT_NEQ(code, seen[i]);
                    seen[n++] = code;

                    int16_t decoded[4];
                    morton_get(decoded, code, 4);
                    ASSERT_POINT_EQ(p, decoded, 4);
                }
            }
        }
    }
    ASSERT_EQ(n, 625);
}

/* Randomized round-trip across the full coordinate range */
TEST(morton4_random_round_trip) {
    test_seed_rng(0x4D4D4D4D);

    for (int i = 0; i < 5000; i++) {
        int16_t p[4] = {
            test_rand_coord(), test_rand_coord(),
            test_rand_coord(), test_rand_coord()
        };
        int16_t decoded[4];

        uint64_t code = morton_set(p, 4);
        morton_get(decoded, code, 4);
        ASSERT_POINT_EQ(p, decoded, 4);
    }
}

int main(void) {
    test_suite_begin("Morton 4D Unit Tests");

    RUN_TEST(morton4_round_trip_origin);
    RUN_TEST(morton4_round_trip_positive);
    RUN_TEST(morton4_round_trip_negative);
    RUN_TEST(morton4_round_trip_mixed);
    RUN_TEST(morton4_round_trip_min_values);
    RUN_TEST(morton4_round_trip_max_values);
    RUN_TEST(morton4_round_trip_boundaries);
    RUN_TEST(morton4_min_code_is_zero);
    RUN_TEST(morton4_max_code_is_u64max);
    RUN_TEST(morton4_determinism);
    RUN_TEST(morton4_grid_uniqueness);
    RUN_TEST(morton4_random_round_trip);

    return test_suite_end();
}
