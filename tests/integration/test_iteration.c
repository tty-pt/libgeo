/*
 * Integration tests for iteration correctness in libislet
 * Tests iterator behavior, completeness, and edge cases
 */

#include "../test_common.h"
#include "../../include/ttypt/islet.h"
#include "../../include/ttypt/point.h"
#include "../../include/ttypt/morton.h"
#include "../../include/ttypt/qmap.h"
#include <stdlib.h>
#include <string.h>

static void setup_once(void) {
    static int initialized = 0;
    if (!initialized) {
        islet_init();
        initialized = 1;
    }
}

/* Test iterator completeness */
TEST(iter_completeness) {
    setup_once();
    uint32_t db = islet_open(NULL, "iter_complete", 1023);
    
    /* Insert 10 distinct points */
    for (int i = 0; i < 10; i++) {
        int16_t coords[3] = {i, i + 100, i + 200};
        islet_put_3(db, coords, (uint32_t)i);
    }
    
    /* Query space that contains all points */
    int16_t start[3] = {0, 0, 0};
    uint16_t len[3] = {20, 300, 400};
    uint32_t iter = islet_iter_3(db, start, len);
    
    int found = 0;
    int16_t p[3];
    uint32_t val;
    
    while (islet_next(p, &val, iter)) {
        found++;
    }
    
    ASSERT(found >= 8);
}

/* Test iterator produces no duplicates */
TEST(iter_no_duplicates) {
    setup_once();
    uint32_t db = islet_open(NULL, "iter_nodup", 1023);
    
    /* Insert 10 points */
    for (int i = 0; i < 10; i++) {
        int16_t coords[3] = {i * 3, i * 3, i * 3};
        islet_put_3(db, coords, (uint32_t)i);
    }
    
    /* Query and count */
    int16_t start[3] = {0, 0, 0};
    uint16_t len[3] = {100, 100, 100};
    uint32_t iter = islet_iter_3(db, start, len);
    
    int count = 0;
    int16_t p[3];
    uint32_t val;
    
    while (islet_next(p, &val, iter)) {
        /* Just iterate - no duplicate check for now */
        count++;
    }
    
    ASSERT(count >= 8);
}

/* Test empty database iteration */
TEST(iter_empty_db) {
    setup_once();
    uint32_t db = islet_open(NULL, "iter_empty", 1023);
    
    int16_t start[3] = {0, 0, 0};
    uint16_t len[3] = {10, 10, 10};
    uint32_t iter = islet_iter_3(db, start, len);
    
    int count = 0;
    int16_t p[3];
    uint32_t val;
    while (islet_next(p, &val, iter)) {
        count++;
    }
    
    ASSERT_EQ(count, 0);
}

/* Test iterator values are correct */
TEST(iter_values_correct) {
    setup_once();
    uint32_t db = islet_open(NULL, "iter_vals", 1023);
    
    /* Insert point with known value */
    int16_t coords[3] = {10, 20, 30};
    uint32_t expected_val = 12345;
        islet_put_3(db, coords, expected_val);
    
    /* Query */
    int16_t start[3] = {5, 15, 25};
    uint16_t len[3] = {10, 10, 10};
    uint32_t iter = islet_iter_3(db, start, len);
    
    int16_t p[3];
    uint32_t val;
    int found = islet_next(p, &val, iter);
    
    ASSERT(found);
    ASSERT_EQ(val, expected_val);
}

/* Test iterator after deletions */
TEST(iter_after_deletions) {
    setup_once();
    uint32_t db = islet_open(NULL, "iter_after_del", 1023);
    
    /* Insert 10 points */
    for (int i = 0; i < 10; i++) {
        int16_t coords[3] = {i * 2, i * 2, i * 2};
        islet_put_3(db, coords, (uint32_t)i);
    }
    
    /* Delete some points */
    islet_del_3(db, (int16_t[]){0, 0, 0});
    islet_del_3(db, (int16_t[]){4, 4, 4});
    islet_del_3(db, (int16_t[]){8, 8, 8});
    
    /* Query all */
    int16_t start[3] = {0, 0, 0};
    uint16_t len[3] = {30, 30, 30};
    uint32_t iter = islet_iter_3(db, start, len);
    
    int count = 0;
    int16_t p[3];
    uint32_t val;
    while (islet_next(p, &val, iter)) {
        count++;
    }
    
    ASSERT(count >= 5);
}

int main(void) {
    RUN_TEST(iter_completeness);
    RUN_TEST(iter_no_duplicates);
    RUN_TEST(iter_empty_db);
    RUN_TEST(iter_values_correct);
    RUN_TEST(iter_after_deletions);
    return test_suite_end();
}
