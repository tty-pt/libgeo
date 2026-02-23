/*
 * Benchmark for iteration performance in libgeo
 */

#include "../test_common.h"
#include "../../include/ttypt/geo.h"
#include "../../include/ttypt/point.h"
#include "../../include/ttypt/morton.h"
#include <stdlib.h>
#include <string.h>

#define BENCH_ITERATIONS 1000

static void setup_once(void) {
    static int initialized = 0;
    if (!initialized) {
        geo_init();
        initialized = 1;
    }
}

int main(void) {
    benchmark_t bench;
    
    setup_once();
    
    printf("\n%s%s=== Point Insert/Update Benchmarks ===%s\n\n", 
           COLOR_BOLD, COLOR_MAGENTA, COLOR_RESET);
    
    /* Benchmark point insert */
    bench_start(&bench, "Point Insert (1K inserts)");
    for (int i = 0; i < BENCH_ITERATIONS; i++) {
        uint32_t db = geo_open(NULL, "bench_iter", 1023);
        int16_t coords[3] = {i % 100, i % 100, i % 100};
        geo_put(db, coords, (uint32_t)i, 3);
    }
    bench_end(&bench, BENCH_ITERATIONS);
    
    printf("\n");
    return 0;
}
