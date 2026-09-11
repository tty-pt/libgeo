/*
 * Benchmark for query performance in libislet
 */

#include "../test_common.h"
#include "../../include/ttypt/islet.h"
#include "../../include/ttypt/point.h"
#include "../../include/ttypt/morton.h"
#include <stdlib.h>
#include <string.h>

#define BENCH_ITERATIONS 1000

static void setup_once(void) {
    static int initialized = 0;
    if (!initialized) {
        islet_init();
        initialized = 1;
    }
}

int main(void) {
    benchmark_t bench;
    
    setup_once();
    
    printf("\n%s%s=== Point Lookup Benchmarks ===%s\n\n", 
           COLOR_BOLD, COLOR_MAGENTA, COLOR_RESET);
    
    /* Create database with some points */
    uint32_t db = islet_open(NULL, "bench_lookup", 1023);
    for (int i = 0; i < 100; i++) {
        int16_t coords[3] = {i, i, i};
        islet_put_3(db, coords, (uint32_t)i);
    }
    
    /* Benchmark point lookup */
    bench_start(&bench, "Point Lookup (1K lookups)");
    for (int i = 0; i < BENCH_ITERATIONS; i++) {
        int16_t coords[3] = {i % 100, i % 100, i % 100};
        islet_get_3(db, coords);
    }
    bench_end(&bench, BENCH_ITERATIONS);
    
    printf("\n");
    return 0;
}
