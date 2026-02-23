/**
 * @file test_morton_props.c
 * @brief Property-based tests for Morton code invariants.
 */

#include "../test_common.h"
#include "../../include/ttypt/morton.h"

#define PROPERTY_TEST_ITERATIONS 10000

/* Property: Morton encoding is reversible */
TEST(property_morton_reversible) {
    test_seed_rng(12345);
    
    for (int i = 0; i < PROPERTY_TEST_ITERATIONS; i++) {
        int16_t original[3] = {
            test_rand_coord(),
            test_rand_coord(),
            test_rand_coord()
        };
        int16_t decoded[3];
        
        uint64_t code = morton_set(original, 3);
        morton_get(decoded, code, 3);
        
        ASSERT_POINT_EQ(original, decoded, 3);
    }
}

/* Property: Different points produce different codes */
TEST(property_morton_unique) {
    test_seed_rng(54321);
    
    /* Test a sample of point pairs */
    for (int i = 0; i < 1000; i++) {
        int16_t p1[3] = {
            test_rand_coord(),
            test_rand_coord(),
            test_rand_coord()
        };
        int16_t p2[3] = {
            test_rand_coord(),
            test_rand_coord(),
            test_rand_coord()
        };
        
        uint64_t c1 = morton_set(p1, 3);
        uint64_t c2 = morton_set(p2, 3);
        
        /* If points are different, codes must be different */
        int points_equal = (p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2]);
        int codes_equal = (c1 == c2);
        
        ASSERT_EQ(points_equal, codes_equal);
    }
}

/* Property: Morton encoding is deterministic */
TEST(property_morton_deterministic) {
    test_seed_rng(99999);
    
    for (int i = 0; i < 1000; i++) {
        int16_t p[3] = {
            test_rand_coord(),
            test_rand_coord(),
            test_rand_coord()
        };
        
        uint64_t code1 = morton_set(p, 3);
        uint64_t code2 = morton_set(p, 3);
        
        ASSERT_EQ(code1, code2);
    }
}

int main(void) {
    test_suite_begin("Morton Code Property Tests");
    
    RUN_TEST(property_morton_reversible);
    RUN_TEST(property_morton_unique);
    RUN_TEST(property_morton_deterministic);
    
    return test_suite_end();
}
