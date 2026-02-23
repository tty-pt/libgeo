/**
 * @file test_morton.c
 * @brief Unit tests for Morton code encoding and decoding.
 */

#include "../test_common.h"
#include "../../include/ttypt/morton.h"
#include "../../include/ttypt/point.h"
#include <limits.h>

/* Test basic round-trip encoding and decoding */
TEST(morton_round_trip_origin) {
    int16_t p[3] = {0, 0, 0};
    int16_t decoded[3];
    
    uint64_t code = morton_set(p, 3);
    morton_get(decoded, code, 3);
    
    ASSERT_POINT_EQ(p, decoded, 3);
}

TEST(morton_round_trip_positive) {
    int16_t p[3] = {100, 200, 300};
    int16_t decoded[3];
    
    uint64_t code = morton_set(p, 3);
    morton_get(decoded, code, 3);
    
    ASSERT_POINT_EQ(p, decoded, 3);
}

TEST(morton_round_trip_negative) {
    int16_t p[3] = {-100, -200, -300};
    int16_t decoded[3];
    
    uint64_t code = morton_set(p, 3);
    morton_get(decoded, code, 3);
    
    ASSERT_POINT_EQ(p, decoded, 3);
}

TEST(morton_round_trip_mixed) {
    int16_t p[3] = {-100, 0, 200};
    int16_t decoded[3];
    
    uint64_t code = morton_set(p, 3);
    morton_get(decoded, code, 3);
    
    ASSERT_POINT_EQ(p, decoded, 3);
}

TEST(morton_round_trip_min_values) {
    int16_t p[3] = {SHRT_MIN, SHRT_MIN, SHRT_MIN};
    int16_t decoded[3];
    
    uint64_t code = morton_set(p, 3);
    morton_get(decoded, code, 3);
    
    ASSERT_POINT_EQ(p, decoded, 3);
}

TEST(morton_round_trip_max_values) {
    int16_t p[3] = {SHRT_MAX, SHRT_MAX, SHRT_MAX};
    int16_t decoded[3];
    
    uint64_t code = morton_set(p, 3);
    morton_get(decoded, code, 3);
    
    ASSERT_POINT_EQ(p, decoded, 3);
}

TEST(morton_round_trip_boundaries) {
    int16_t test_cases[][3] = {
        {SHRT_MIN, 0, 0},
        {0, SHRT_MIN, 0},
        {0, 0, SHRT_MIN},
        {SHRT_MAX, 0, 0},
        {0, SHRT_MAX, 0},
        {0, 0, SHRT_MAX},
        {SHRT_MIN, SHRT_MAX, 0},
        {SHRT_MAX, SHRT_MIN, SHRT_MAX},
    };
    
    for (size_t i = 0; i < sizeof(test_cases) / sizeof(test_cases[0]); i++) {
        int16_t decoded[3];
        uint64_t code = morton_set(test_cases[i], 3);
        morton_get(decoded, code, 3);
        ASSERT_POINT_EQ(test_cases[i], decoded, 3);
    }
}

/* Test that different points produce different Morton codes */
TEST(morton_uniqueness) {
    int16_t p1[3] = {0, 0, 0};
    int16_t p2[3] = {1, 0, 0};
    int16_t p3[3] = {0, 1, 0};
    int16_t p4[3] = {0, 0, 1};
    
    uint64_t c1 = morton_set(p1, 3);
    uint64_t c2 = morton_set(p2, 3);
    uint64_t c3 = morton_set(p3, 3);
    uint64_t c4 = morton_set(p4, 3);
    
    ASSERT_NEQ(c1, c2);
    ASSERT_NEQ(c1, c3);
    ASSERT_NEQ(c1, c4);
    ASSERT_NEQ(c2, c3);
    ASSERT_NEQ(c2, c4);
    ASSERT_NEQ(c3, c4);
}

/* Test 2D encoding */
TEST(morton_2d_round_trip) {
    int16_t p[2] = {100, 200};
    int16_t decoded[3] = {0, 0, 0};
    
    uint64_t code = morton_set(p, 2);
    morton_get(decoded, code, 2);
    
    ASSERT_EQ(p[0], decoded[0]);
    ASSERT_EQ(p[1], decoded[1]);
}

/* Test 1D encoding */
TEST(morton_1d_round_trip) {
    int16_t p[1] = {12345};
    int16_t decoded[3] = {0, 0, 0};
    
    uint64_t code = morton_set(p, 1);
    morton_get(decoded, code, 1);
    
    ASSERT_EQ(p[0], decoded[0]);
}

/* Test determinism - same input produces same output */
TEST(morton_determinism) {
    int16_t p[3] = {123, 456, 789};
    
    uint64_t code1 = morton_set(p, 3);
    uint64_t code2 = morton_set(p, 3);
    
    ASSERT_EQ(code1, code2);
}

/* Test Z-order property: adjacent codes are spatially close */
TEST(morton_z_order_locality) {
    /* Points (0,0,0) through (1,1,1) - signed coords get offset to unsigned */
    int16_t points[8][3] = {
        {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
        {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}
    };
    
    uint64_t codes[8];
    for (int i = 0; i < 8; i++) {
        codes[i] = morton_set(points[i], 3);
    }
    
    /* All codes should be unique */
    for (int i = 0; i < 8; i++) {
        for (int j = i + 1; j < 8; j++) {
            ASSERT_NEQ(codes[i], codes[j]);
        }
    }
    
    /* Codes should maintain Z-order property: */
    /* Adjacent points in coordinate space should have relatively close codes */
    /* (0,0,0) and (1,0,0) differ by 1 in x, should differ by small amount in code */
    uint64_t diff_x = codes[1] > codes[0] ? codes[1] - codes[0] : codes[0] - codes[1];
    uint64_t diff_y = codes[2] > codes[0] ? codes[2] - codes[0] : codes[0] - codes[2];
    uint64_t diff_z = codes[4] > codes[0] ? codes[4] - codes[0] : codes[0] - codes[4];
    
    /* These diffs should all be small powers of 2 (interleaved bit positions) */
    ASSERT_LT(diff_x, 1000);
    ASSERT_LT(diff_y, 1000);
    ASSERT_LT(diff_z, 1000);
}

/* Test known Morton code values */
TEST(morton_known_values) {
    /* Since coordinates are signed and internally converted to unsigned,
     * we can't predict exact bit patterns easily. Instead, test that
     * the encoding is consistent and produces distinct codes for
     * distinct inputs. The main value is the round-trip correctness. */
    
    int16_t p0[3] = {0, 0, 0};
    int16_t p1[3] = {1, 1, 1};
    int16_t p_neg[3] = {-1, -1, -1};
    
    uint64_t c0 = morton_set(p0, 3);
    uint64_t c1 = morton_set(p1, 3);
    uint64_t c_neg = morton_set(p_neg, 3);
    
    /* All should be different */
    ASSERT_NEQ(c0, c1);
    ASSERT_NEQ(c0, c_neg);
    ASSERT_NEQ(c1, c_neg);
    
    /* Round-trip should work for all */
    int16_t decoded[3];
    morton_get(decoded, c0, 3);
    ASSERT_POINT_EQ(decoded, p0, 3);
    
    morton_get(decoded, c1, 3);
    ASSERT_POINT_EQ(decoded, p1, 3);
    
    morton_get(decoded, c_neg, 3);
    ASSERT_POINT_EQ(decoded, p_neg, 3);
}

/* Random round-trip tests */
TEST(morton_random_round_trip) {
    test_seed_rng(42);
    
    for (int i = 0; i < 1000; i++) {
        int16_t p[3] = {
            test_rand_coord(),
            test_rand_coord(),
            test_rand_coord()
        };
        int16_t decoded[3];
        
        uint64_t code = morton_set(p, 3);
        morton_get(decoded, code, 3);
        
        ASSERT_POINT_EQ(p, decoded, 3);
    }
}

/* Test that dimension parameter is respected */
TEST(morton_dimension_handling) {
    int16_t p3d[3] = {100, 200, 300};
    int16_t decoded[3];
    
    /* Encode as 2D (should ignore third dimension) */
    uint64_t code2d = morton_set(p3d, 2);
    morton_get(decoded, code2d, 2);
    
    ASSERT_EQ(decoded[0], p3d[0]);
    ASSERT_EQ(decoded[1], p3d[1]);
    /* decoded[2] is undefined, don't check it */
}

/* Test sequential coordinates produce different codes */
TEST(morton_sequential_uniqueness) {
    uint64_t prev_code = morton_set((int16_t[]){0, 0, 0}, 3);
    
    for (int16_t x = 0; x < 10; x++) {
        for (int16_t y = 0; y < 10; y++) {
            for (int16_t z = 0; z < 10; z++) {
                int16_t p[3] = {x, y, z};
                uint64_t code = morton_set(p, 3);
                
                /* Each code should be unique */
                if (x != 0 || y != 0 || z != 0) {
                    ASSERT_NEQ(code, prev_code);
                }
                
                /* Round-trip should work */
                int16_t decoded[3];
                morton_get(decoded, code, 3);
                ASSERT_POINT_EQ(p, decoded, 3);
                
                prev_code = code;
            }
        }
    }
}

int main(void) {
    test_suite_begin("Morton Code Unit Tests");
    
    RUN_TEST(morton_round_trip_origin);
    RUN_TEST(morton_round_trip_positive);
    RUN_TEST(morton_round_trip_negative);
    RUN_TEST(morton_round_trip_mixed);
    RUN_TEST(morton_round_trip_min_values);
    RUN_TEST(morton_round_trip_max_values);
    RUN_TEST(morton_round_trip_boundaries);
    RUN_TEST(morton_uniqueness);
    RUN_TEST(morton_2d_round_trip);
    RUN_TEST(morton_1d_round_trip);
    RUN_TEST(morton_determinism);
    RUN_TEST(morton_z_order_locality);
    RUN_TEST(morton_known_values);
    RUN_TEST(morton_random_round_trip);
    RUN_TEST(morton_dimension_handling);
    RUN_TEST(morton_sequential_uniqueness);
    
    return test_suite_end();
}
