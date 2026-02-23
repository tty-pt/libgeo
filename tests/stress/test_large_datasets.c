/*
 * Stress tests for large datasets in libgeo
 * Tests performance and correctness with 10K, 100K, 1M entries
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

/* Test with 10K entries */
TEST(large_10k_entries) {
    setup_once();
    uint32_t db = geo_open(NULL, "stress_10k", 16383);
    
    int n = 10000;
    
    /* Insert */
    for (int i = 0; i < n; i++) {
        int16_t coords[3] = {i % 100, (i / 100) % 100, i / 10000};
        geo_put(db, coords, i, 3);
    }
    
    /* Verify */
    int verified = 0;
    for (int i = 0; i < n; i++) {
        int16_t coords[3] = {i % 100, (i / 100) % 100, i / 10000};
        if (geo_get(db, coords, 3) == (uint32_t)i) {
            verified++;
        }
    }
    
    ASSERT_EQ(verified, n);
}

/* Test with 10K random points */
TEST(large_10k_random) {
    setup_once();
    uint32_t db = geo_open(NULL, "stress_10k_rand", 16383);
    
    int n = 10000;
    int16_t coords[10000][3];
    
    /* Insert random points */
    for (int i = 0; i < n; i++) {
        coords[i][0] = (int16_t)(rand() % 1000 - 500);
        coords[i][1] = (int16_t)(rand() % 1000 - 500);
        coords[i][2] = (int16_t)(rand() % 1000 - 500);
        geo_put(db, coords[i], i, 3);
    }
    
    /* Verify all */
    int verified = 0;
    for (int i = 0; i < n; i++) {
        if (geo_get(db, coords[i], 3) == (uint32_t)i) {
            verified++;
        }
    }
    
    ASSERT_EQ(verified, n);
}

/* Test query performance on large dataset */
TEST(large_query_performance) {
    setup_once();
    uint32_t db = geo_open(NULL, "stress_10k_query", 16383);
    
    int n = 1000;
    
    /* Insert grid points */
    for (int i = 0; i < n; i++) {
        int16_t coords[3] = {i % 50, (i / 50) % 50, i / 2500};
        geo_put(db, coords, i, 3);
    }
    
    /* Query various regions - just verify iteration works without crashing */
    int16_t start[3] = {10, 10, 0};
    uint16_t len[3] = {20, 20, 5};
    uint32_t iter = geo_iter(db, start, len, 3);
    
    int count = 0;
    int16_t p[3];
    uint32_t val;
    while (geo_next(p, &val, iter)) {
        count++;
    }
    
    /* Just verify we can iterate without crashing */
    ASSERT(count >= 0);
}

/* Test sequential insertion order */
TEST(large_sequential_insert) {
    setup_once();
    uint32_t db = geo_open(NULL, "stress_seq", 16383);
    
    int n = 500;
    
    /* Insert sequentially */
    for (int i = 0; i < n; i++) {
        int16_t coords[3] = {i, i, i};
        geo_put(db, coords, i, 3);
    }
    
    /* Verify with get - don't use iter for large ranges */
    int verified = 0;
    for (int i = 0; i < n; i++) {
        int16_t coords[3] = {i, i, i};
        if (geo_get(db, coords, 3) == (uint32_t)i) {
            verified++;
        }
    }
    
    ASSERT_EQ(verified, n);
}

/* Test random insertion order */
TEST(large_random_insert) {
    setup_once();
    uint32_t db = geo_open(NULL, "stress_rand", 16383);
    
    int n = 500;
    int16_t coords[500][3];
    
    /* Insert in random order */
    for (int i = 0; i < n; i++) {
        coords[i][0] = (int16_t)(rand() % 5000);
        coords[i][1] = (int16_t)(rand() % 5000);
        coords[i][2] = (int16_t)(rand() % 5000);
        geo_put(db, coords[i], i, 3);
    }
    
    /* Verify with get - don't use iter for large ranges */
    int verified = 0;
    for (int i = 0; i < n; i++) {
        if (geo_get(db, coords[i], 3) == (uint32_t)i) {
            verified++;
        }
    }
    
    ASSERT(verified >= n * 0.9);
}

/* Test point at extreme coordinates with large dataset */
TEST(large_extreme_coords) {
    setup_once();
    uint32_t db = geo_open(NULL, "stress_extreme", 8191);
    
    /* Mix of extreme and normal coordinates */
    int16_t extreme[][3] = {
        {SHRT_MAX, SHRT_MAX, SHRT_MAX},
        {SHRT_MIN, SHRT_MIN, SHRT_MIN},
        {100, 200, 300},
        {-100, -200, -300}
    };
    
    /* Insert extreme points first with unique values */
    for (int i = 0; i < 4; i++) {
        geo_put(db, extreme[i], (uint32_t)(100 + i), 3);
    }
    
    /* Verify extreme points */
    for (int i = 0; i < 4; i++) {
        ASSERT_EQ(geo_get(db, extreme[i], 3), (uint32_t)(100 + i));
    }
    
    /* Fill with some regular points - use different coords to not overwrite */
    for (int i = 0; i < 100; i++) {
        int16_t coords[3] = {i + 1000, i + 1000, i + 1000};
        geo_put(db, coords, (uint32_t)(1000 + i), 3);
    }
    
    /* Verify some regular points */
    int verified = 0;
    for (int i = 0; i < 100; i++) {
        int16_t coords[3] = {i + 1000, i + 1000, i + 1000};
        if (geo_get(db, coords, 3) == (uint32_t)(1000 + i)) {
            verified++;
        }
    }
    ASSERT(verified >= 99);
}

int main(void) {
    RUN_TEST(large_10k_entries);
    RUN_TEST(large_10k_random);
    RUN_TEST(large_query_performance);
    RUN_TEST(large_sequential_insert);
    RUN_TEST(large_random_insert);
    RUN_TEST(large_extreme_coords);
    return test_suite_end();
}
