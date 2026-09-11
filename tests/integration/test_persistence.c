/*
 * Integration tests for persistence in libislet
 * Note: Full file persistence testing is limited due to API constraints
 * These tests verify in-memory database behavior
 */

#include "../test_common.h"
#include "../../include/ttypt/islet.h"
#include "../../include/ttypt/point.h"
#include "../../include/ttypt/morton.h"
#include "../../include/ttypt/qmap.h"
#include <stdlib.h>
#include <string.h>
#include <limits.h>

static void setup_once(void) {
    static int initialized = 0;
    if (!initialized) {
        islet_init();
        initialized = 1;
    }
}

/* Test basic database operations work */
TEST(persistence_basic_ops) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_persist1", 1023);
    
    int16_t coords1[3] = {100, 200, 300};
    int16_t coords2[3] = {-50, -100, -150};
    
    islet_put_3(db, coords1, 42);
    islet_put_3(db, coords2, 99);
    
    /* Verify we can read the data */
    ASSERT_EQ(islet_get_3(db, coords1), 42u);
    ASSERT_EQ(islet_get_3(db, coords2), 99u);
}

/* Test persistence with many points */
TEST(persistence_many_points) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_persist2", 16383);
    
    int num_points = 500;
    
    /* Insert grid of points */
    for (int i = 0; i < num_points; i++) {
        int16_t coords[3] = {i % 50, (i / 50) % 50, i / 2500};
        islet_put_3(db, coords, (uint32_t)(1000 + i));
    }
    
    /* Verify all points */
    int verified = 0;
    for (int i = 0; i < num_points; i++) {
        int16_t coords[3] = {i % 50, (i / 50) % 50, i / 2500};
        uint32_t val = islet_get_3(db, coords);
        if (val == (uint32_t)(1000 + i)) {
            verified++;
        }
    }
    
    ASSERT_EQ(verified, num_points);
}

/* Test persistence after deletions */
TEST(persistence_with_deletions) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_persist3", 1023);
    
    /* Insert 10 points */
    for (int i = 0; i < 10; i++) {
        int16_t coords[3] = {i, i, i};
        islet_put_3(db, coords, (uint32_t)i);
    }
    
    /* Delete every other point */
    for (int i = 0; i < 10; i += 2) {
        int16_t coords[3] = {i, i, i};
        islet_del_3(db, coords);
    }
    
    /* Verify deleted points are gone, others remain */
    for (int i = 0; i < 10; i++) {
        int16_t coords[3] = {i, i, i};
        uint32_t val = islet_get_3(db, coords);
        
        if (i % 2 == 0) {
            ASSERT_EQ(val, QM_MISS);
        } else {
            ASSERT_EQ(val, (uint32_t)i);
        }
    }
}

/* Test persistence after updates */
TEST(persistence_with_updates) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_persist4", 1023);
    
    /* Insert points */
    for (int i = 0; i < 5; i++) {
        int16_t coords[3] = {i, i, i};
        islet_put_3(db, coords, (uint32_t)i);
    }
    
    /* Update values (replace semantics) */
    for (int i = 0; i < 5; i++) {
        int16_t coords[3] = {i, i, i};
        islet_set_3(db, coords, (uint32_t)(1000 + i));
    }
    
    /* Verify updated values */
    for (int i = 0; i < 5; i++) {
        int16_t coords[3] = {i, i, i};
        ASSERT_EQ(islet_get_3(db, coords), (uint32_t)(1000 + i));
    }
}

/* Test incremental updates */
TEST(persistence_incremental_updates) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_persist5", 1023);
    
    /* Insert initial data */
    for (int i = 0; i < 10; i++) {
        int16_t coords[3] = {i, 0, 0};
        islet_put_3(db, coords, (uint32_t)i);
    }
    
    /* Add more data */
    for (int i = 10; i < 20; i++) {
        int16_t coords[3] = {i, 0, 0};
        islet_put_3(db, coords, (uint32_t)i);
    }
    
    /* Verify all data */
    for (int i = 0; i < 20; i++) {
        int16_t coords[3] = {i, 0, 0};
        ASSERT_EQ(islet_get_3(db, coords), (uint32_t)i);
    }
}

int main(void) {
    RUN_TEST(persistence_basic_ops);
    RUN_TEST(persistence_many_points);
    RUN_TEST(persistence_with_deletions);
    RUN_TEST(persistence_with_updates);
    RUN_TEST(persistence_incremental_updates);
    return test_suite_end();
}
