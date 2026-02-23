/**
 * @file bench_morton.c
 * @brief Benchmark for Morton encoding/decoding performance.
 */

#include "../test_common.h"
#include "../../include/ttypt/morton.h"
#include <stdio.h>

#define BENCH_ITERATIONS 1000000

int main(void) {
    benchmark_t bench;
    int16_t pos[3];
    uint64_t code;
    
    printf("\n%s%s=== Morton Code Benchmarks ===%s\n\n", 
           COLOR_BOLD, COLOR_MAGENTA, COLOR_RESET);
    
    /* Benchmark encoding */
    bench_start(&bench, "Morton Encode (3D)");
    test_seed_rng(42);
    for (int i = 0; i < BENCH_ITERATIONS; i++) {
        pos[0] = test_rand_coord();
        pos[1] = test_rand_coord();
        pos[2] = test_rand_coord();
        code = morton_set(pos, 3);
        (void)code; /* Prevent optimization */
    }
    bench_end(&bench, BENCH_ITERATIONS);
    
    /* Benchmark decoding */
    bench_start(&bench, "Morton Decode (3D)");
    test_seed_rng(42);
    for (int i = 0; i < BENCH_ITERATIONS; i++) {
        code = test_rand64();
        morton_get(pos, code, 3);
    }
    bench_end(&bench, BENCH_ITERATIONS);
    
    /* Benchmark round-trip */
    bench_start(&bench, "Morton Round-Trip (3D)");
    test_seed_rng(42);
    int16_t decoded[3];
    for (int i = 0; i < BENCH_ITERATIONS; i++) {
        pos[0] = test_rand_coord();
        pos[1] = test_rand_coord();
        pos[2] = test_rand_coord();
        code = morton_set(pos, 3);
        morton_get(decoded, code, 3);
    }
    bench_end(&bench, BENCH_ITERATIONS);
    
    printf("\n");
    return 0;
}
