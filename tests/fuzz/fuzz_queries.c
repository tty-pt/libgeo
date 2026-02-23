/*
 * Fuzz test for spatial queries
 * Tests random bounding box queries and verifies results
 */

#include "../../include/ttypt/geo.h"
#include "../../include/ttypt/point.h"
#include "../../include/ttypt/morton.h"
#include "../../include/ttypt/qmap.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef LIBFUZZER
#include <fuzzer/FuzzedDataProvider.h>
#endif

#define MAX_POINTS 100
#define MASK 255

static uint32_t db;
static int16_t stored_coords[MAX_POINTS][3];
static int16_t query_start[3];
static uint16_t query_len[3];
static int stored_count = 0;

static void setup(void) {
    static int initialized = 0;
    if (!initialized) {
        geo_init();
        initialized = 1;
    }
    db = geo_open(NULL, "fuzz_query_db", MASK);
    stored_count = 0;
}

static void insert_points(int16_t start_x, int16_t start_y, int16_t start_z, int count) {
    for (int i = 0; i < count && stored_count < MAX_POINTS; i++) {
        int16_t x = start_x + i;
        int16_t y = start_y + (i / 10);
        int16_t z = start_z + (i / 50);
        
        stored_coords[stored_count][0] = x;
        stored_coords[stored_count][1] = y;
        stored_coords[stored_count][2] = z;
        stored_count++;
        
        int16_t coords[3] = {x, y, z};
        geo_put(db, coords, (uint32_t)i, 3);
    }
}

static int is_in_range(int16_t x, int16_t y, int16_t z) {
    int16_t end_x = query_start[0] + query_len[0];
    int16_t end_y = query_start[1] + query_len[1];
    int16_t end_z = query_start[2] + query_len[2];
    
    return (x >= query_start[0] && x < end_x &&
            y >= query_start[1] && y < end_y &&
            z >= query_start[2] && z < end_z);
}

static void test_query(void) {
    uint32_t iter = geo_iter(db, query_start, query_len, 3);
    
    int16_t p[3];
    uint32_t val;
    int count = 0;
    
    while (geo_next(p, &val, iter)) {
        /* Verify point is within query range */
        if (!is_in_range(p[0], p[1], p[2])) {
            fprintf(stderr, "Query returned point outside range: (%d,%d,%d) not in [%d,%d,%d] + [%d,%d,%d]\n",
                p[0], p[1], p[2],
                query_start[0], query_start[1], query_start[2],
                query_len[0], query_len[1], query_len[2]);
            abort();
        }
        count++;
        
        /* Safety limit */
        if (count > MAX_POINTS * 2) {
            fprintf(stderr, "Query returned too many points\n");
            abort();
        }
    }
}

/* Process fuzz input */
static void process_fuzz(const uint8_t *data, size_t size) {
    setup();
    
    if (size < 12) return;
    
    /* Extract query bounds from fuzz data */
    query_start[0] = (int16_t)(data[0] | (data[1] << 8)) - 64;
    query_start[1] = (int16_t)(data[2] | (data[3] << 8)) - 64;
    query_start[2] = (int16_t)(data[4] | (data[5] << 8)) - 64;
    query_len[0] = (uint16_t)(data[6] | (data[7] << 8)) % 20 + 1;
    query_len[1] = (uint16_t)(data[8] | (data[9] << 8)) % 20 + 1;
    query_len[2] = (uint16_t)(data[10] | (data[11] << 8)) % 20 + 1;
    
    /* Insert some points in and around the query region */
    insert_points(query_start[0] - 10, query_start[1] - 10, query_start[2] - 10, 50);
    
    /* Run query and verify results */
    test_query();
}

/* Standalone test mode */
#ifdef STANDALONE
int main(void) {
    printf("Running spatial query fuzz tests (standalone mode)...\n");
    
    /* Test various query patterns */
    uint8_t test_data[64];
    
    /* Small query regions */
    for (int i = 0; i < 10; i++) {
        test_data[0] = i; test_data[1] = 0;
        test_data[2] = i; test_data[3] = 0;
        test_data[4] = 0; test_data[5] = 0;
        test_data[6] = 10; test_data[7] = 0;
        test_data[8] = 10; test_data[9] = 0;
        test_data[10] = 1; test_data[11] = 0;
        process_fuzz(test_data, sizeof(test_data));
    }
    
    /* Larger regions */
    for (int i = 0; i < 10; i++) {
        memset(test_data, 0, sizeof(test_data));
        test_data[0] = (uint8_t)(-50 & 0xFF);
        test_data[1] = (uint8_t)((-50 >> 8) & 0xFF);
        test_data[6] = 50; test_data[7] = 0;
        test_data[8] = 50; test_data[9] = 0;
        test_data[10] = 10; test_data[11] = 0;
        process_fuzz(test_data, sizeof(test_data));
    }
    
    printf("All standalone tests passed!\n");
    return 0;
}
#endif

/* libFuzzer test mode */
#ifdef LIBFUZZER
extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    if (size < 12) return 0;
    
    FuzzedDataProvider provider(data, size);
    
    setup();
    
    /* Extract query bounds */
    query_start[0] = provider.ConsumeIntegral<int16_t>();
    query_start[1] = provider.ConsumeIntegral<int16_t>();
    query_start[2] = provider.ConsumeIntegral<int16_t>();
    query_len[0] = provider.ConsumeIntegral<uint16_t>() % 50 + 1;
    query_len[1] = provider.ConsumeIntegral<uint16_t>() % 50 + 1;
    query_len[2] = provider.ConsumeIntegral<uint16_t>() % 50 + 1;
    
    /* Insert points */
    insert_points(query_start[0] - 5, query_start[1] - 5, query_start[2] - 5, 30);
    
    /* Run query and verify */
    test_query();
    
    return 0;
}
#endif
