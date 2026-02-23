/*
 * Integration tests for spatial queries in libgeo
 * Tests complex query scenarios with realistic data patterns
 */

#include "../test_common.h"
#include "../../include/ttypt/geo.h"
#include "../../include/ttypt/point.h"
#include "../../include/ttypt/morton.h"
#include "../../include/ttypt/qmap.h"
#include <stdlib.h>
#include <string.h>

static void setup_once(void) {
    static int initialized = 0;
    if (!initialized) {
        geo_init();
        initialized = 1;
    }
}

/* Test query with empty result */
TEST(query_empty_region) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_spatial_3", 1023);
    
    /* Insert some points */
    int16_t coords[3] = {100, 100, 100};
    geo_put(db, coords, 42, 3);
    
    /* Query a different region */
    int16_t start[3] = {0, 0, 0};
    uint16_t len[3] = {10, 10, 10};
    uint32_t iter = geo_iter(db, start, len, 3);
    
    int count = 0;
    int16_t p[3];
    uint32_t val;
    while (geo_next(p, &val, iter)) {
        count++;
    }
    
    ASSERT_EQ(count, 0);
}

/* Test query after deletions */
TEST(query_after_deletions) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_spatial_5", 1023);
    
    /* Insert a few points */
    for (int i = 0; i < 5; i++) {
        int16_t coords[3] = {i, i, i};
        geo_put(db, coords, (uint32_t)i, 3);
    }
    
    /* Delete some */
    geo_del(db, (int16_t[]){0, 0, 0}, 3);
    geo_del(db, (int16_t[]){2, 2, 2}, 3);
    
    /* Verify using geo_get */
    ASSERT_EQ(geo_get(db, (int16_t[]){0, 0, 0}, 3), QM_MISS);
    ASSERT_EQ(geo_get(db, (int16_t[]){2, 2, 2}, 3), QM_MISS);
    ASSERT_EQ(geo_get(db, (int16_t[]){1, 1, 1}, 3), 1u);
}

/* Test query with point updates */
TEST(query_after_updates) {
    setup_once();
    uint32_t db = geo_open(NULL, "test_spatial_6", 1023);
    
    /* Insert initial points */
    for (int i = 0; i < 3; i++) {
        int16_t coords[3] = {i, i, i};
        geo_put(db, coords, (uint32_t)i, 3);
    }
    
    /* Update values */
    for (int i = 0; i < 3; i++) {
        int16_t coords[3] = {i, i, i};
        geo_put(db, coords, (uint32_t)(1000 + i), 3);
    }
    
    /* Verify updated values are retrievable */
    for (int i = 0; i < 3; i++) {
        int16_t coords[3] = {i, i, i};
        ASSERT_EQ(geo_get(db, coords, 3), (uint32_t)(1000 + i));
    }
}

int main(void) {
    RUN_TEST(query_empty_region);
    RUN_TEST(query_after_deletions);
    RUN_TEST(query_after_updates);
    return test_suite_end();
}
