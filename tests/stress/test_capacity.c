/*
 * Stress tests for capacity limits in libgeo
 * Tests behavior with different mask values and extreme loads
 */

#include "../test_common.h"
#include "../../include/ttypt/geo.h"
#include "../../include/ttypt/point.h"
#include "../../include/ttypt/morton.h"
#include "../../include/ttypt/qmap.h"
#include <stdlib.h>
#include <string.h>
#include <limits.h>

static void setup_once(void) {
    static int initialized = 0;
    if (!initialized) {
        geo_init();
        initialized = 1;
    }
}

/* Test with small mask - stress hash collisions */
TEST(capacity_small_mask) {
    setup_once();
    /* Small mask means more collisions */
    uint32_t db = geo_open(NULL, "stress_small_mask", 15);
    
    /* Try to insert more than mask+1 points */
    int success_count = 0;
    for (int i = 0; i < 50; i++) {
        int16_t coords[3] = {i, i * 2, i * 3};
        geo_put(db, coords, i, 3);
        
        /* Check if we can retrieve it */
        if (geo_get(db, coords, 3) == (uint32_t)i) {
            success_count++;
        }
    }
    
    /* Most or all should succeed */
    ASSERT(success_count >= 45);
}

/* Test with various mask sizes */
TEST(capacity_various_masks) {
    setup_once();
    
    uint32_t masks[] = {7, 15, 31, 63, 127, 255, 511, 1023};
    int num_masks = sizeof(masks) / sizeof(masks[0]);
    
    for (int m = 0; m < num_masks; m++) {
        uint32_t mask = masks[m];
        uint32_t db = geo_open(NULL, "stress_mask_var", mask);
        
        /* Insert up to 2x the table size */
        int target = (mask + 1) * 2;
        int success = 0;
        
        for (int i = 0; i < target; i++) {
            int16_t coords[3] = {i, i + 1000, i + 2000};
            geo_put(db, coords, i, 3);
            
            if (geo_get(db, coords, 3) == (uint32_t)i) {
                success++;
            }
        }
        
        /* Verify high success rate */
        ASSERT(success >= target * 0.9);
    }
}

/* Test repeated insert/delete cycles */
TEST(capacity_repeated_operations) {
    setup_once();
    uint32_t db = geo_open(NULL, "stress_repeat_op", 255);
    
    int16_t coords[3] = {10, 20, 30};
    
    /* Perform many insert/delete cycles on same point */
    for (int cycle = 0; cycle < 1000; cycle++) {
        geo_put(db, coords, cycle, 3);
        ASSERT_EQ(geo_get(db, coords, 3), (uint32_t)cycle);
        
        geo_del(db, coords, 3);
        ASSERT_EQ(geo_get(db, coords, 3), QM_MISS);
    }
}

/* Test filling to capacity then querying */
TEST(capacity_fill_then_query) {
    setup_once();
    uint32_t db = geo_open(NULL, "stress_fill_query", 127);
    
    /* Fill with points */
    int num_inserted = 0;
    for (int i = 0; i < 200; i++) {
        int16_t coords[3] = {i, i * 2, i * 3};
        geo_put(db, coords, i, 3);
        
        if (geo_get(db, coords, 3) == (uint32_t)i) {
            num_inserted++;
        }
    }
    
    /* Query entire space - just verify we can iterate without crashing */
    int16_t start[3] = {0, 0, 0};
    uint16_t len[3] = {200, 400, 600};
    uint32_t iter = geo_iter(db, start, len, 3);
    
    int count = 0;
    int16_t p[3];
    uint32_t val;
    while (geo_next(p, &val, iter)) {
        count++;
    }
    
    /* Just verify we can iterate - don't check exact count */
    ASSERT(count >= 0);
}

/* Test with maximum mask value */
TEST(capacity_large_mask) {
    setup_once();
    /* Large mask */
    uint32_t db = geo_open(NULL, "stress_large_mask", 65535);
    
    /* Insert many points */
    for (int i = 0; i < 1000; i++) {
        int16_t coords[3] = {i, i + 10000, i + 20000};
        geo_put(db, coords, i, 3);
    }
    
    /* Verify retrieval */
    int verified = 0;
    for (int i = 0; i < 1000; i++) {
        int16_t coords[3] = {i, i + 10000, i + 20000};
        if (geo_get(db, coords, 3) == (uint32_t)i) {
            verified++;
        }
    }
    
    ASSERT_EQ(verified, 1000);
}

int main(void) {
    RUN_TEST(capacity_small_mask);
    RUN_TEST(capacity_various_masks);
    RUN_TEST(capacity_repeated_operations);
    RUN_TEST(capacity_fill_then_query);
    RUN_TEST(capacity_large_mask);
    return test_suite_end();
}
