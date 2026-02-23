/*
 * Fuzz test for Morton encoding/decoding
 * Uses libFuzzer if available, otherwise can be run standalone
 */

#include "../../include/ttypt/morton.h"
#include "../../include/ttypt/point.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#ifdef LIBFUZZER
#include <fuzzer/FuzzedDataProvider.h>
#endif

/* Morton encode/decode round-trip property */
static int test_roundtrip(int16_t x, int16_t y, int16_t z) {
    int16_t coords[3] = {x, y, z};
    int16_t decoded[3];
    
    uint64_t code = morton_set(coords, 3);
    morton_get(decoded, code, 3);
    
    return (decoded[0] == x && decoded[1] == y && decoded[2] == z);
}

/* Test a single coordinate triple */
static void test_coords(int16_t x, int16_t y, int16_t z) {
    /* Round-trip test */
    if (!test_roundtrip(x, y, z)) {
        fprintf(stderr, "Round-trip failed: (%d, %d, %d)\n", x, y, z);
        abort();
    }
    
    /* Test 2D (uses x, y only) */
    if (!test_roundtrip(x, y, 0)) {
        fprintf(stderr, "2D round-trip failed: (%d, %d)\n", x, y);
        abort();
    }
}

/* Standalone test mode */
#ifdef STANDALONE
int main(void) {
    printf("Running Morton fuzz tests (standalone mode)...\n");
    
    /* Test various coordinate patterns */
    for (int i = -1000; i < 1000; i++) {
        test_coords(i, i, i);
        test_coords(i + 100, i - 50, i + 200);
        test_coords(0, i, -i);
        test_coords(SHRT_MIN, SHRT_MAX, 0);
    }
    
    printf("All standalone tests passed!\n");
    return 0;
}
#endif

/* libFuzzer test mode */
#ifdef LIBFUZZER
extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    FuzzedDataProvider provider(data, size);
    
    /* Extract coordinates from fuzz data */
    int16_t x = provider.ConsumeIntegral<int16_t>();
    int16_t y = provider.ConsumeIntegral<int16_t>();
    int16_t z = provider.ConsumeIntegral<int16_t>();
    
    test_coords(x, y, z);
    
    return 0;
}
#endif
