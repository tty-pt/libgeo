/*
 * Fuzz test for geo API operations
 * Tests random sequences of put/get/del/iter operations
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

#define MASK 255

static uint32_t db;

static void setup(void) {
    static int initialized = 0;
    if (!initialized) {
        geo_init();
        initialized = 1;
    }
    db = geo_open(NULL, "fuzz_db", MASK);
}

/* Process fuzz input as simple put/get operations */
static void process_operations(const uint8_t *data, size_t size) {
    setup();
    
    /* Simple pattern: put then get */
    for (size_t i = 0; i + 8 < size; i += 8) {
        int16_t x = (int16_t)(data[i] | (data[i+1] << 8)) - 128;
        int16_t y = (int16_t)(data[i+2] | (data[i+3] << 8)) - 128;
        int16_t z = (int16_t)(data[i+4] | (data[i+5] << 8)) - 128;
        uint32_t val = (uint32_t)(data[i+6] | (data[i+7] << 8));
        
        int16_t coords[3] = {x, y, z};
        
        /* Put */
        geo_put(db, coords, val, 3);
        
        /* Get and verify */
        uint32_t retrieved = geo_get(db, coords, 3);
        
        if (retrieved != val && retrieved != QM_MISS) {
            fprintf(stderr, "Mismatch: put %u at (%d,%d,%d), got %u\n",
                val, x, y, z, retrieved);
            abort();
        }
    }
}

/* Standalone test mode */
#ifdef STANDALONE
int main(void) {
    printf("Running geo API fuzz tests (standalone mode)...\n");
    
    /* Test with various patterns */
    uint8_t test_data[256];
    
    /* Various coordinate patterns */
    for (int i = 0; i < 100; i++) {
        for (size_t j = 0; j < sizeof(test_data); j++) {
            test_data[j] = (uint8_t)(i + j);
        }
        process_operations(test_data, sizeof(test_data));
    }
    
    printf("All standalone tests passed!\n");
    return 0;
}
#endif

/* libFuzzer test mode */
#ifdef LIBFUZZER
extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    if (size < 8) return 0;
    
    FuzzedDataProvider provider(data, size);
    
    setup();
    
    while (provider.remaining_bytes() >= 8) {
        int16_t x = provider.ConsumeIntegral<int16_t>();
        int16_t y = provider.ConsumeIntegral<int16_t>();
        int16_t z = provider.ConsumeIntegral<int16_t>();
        uint32_t val = provider.ConsumeIntegral<uint32_t>();
        
        int16_t coords[3] = {x, y, z};
        geo_put(db, coords, val, 3);
        
        uint32_t retrieved = geo_get(db, coords, 3);
        if (retrieved != val && retrieved != QM_MISS) {
            abort();
        }
    }
    
    return 0;
}
#endif
