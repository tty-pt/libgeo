/**
 * @file bench_bulk4.c
 * @brief Paired benchmark: scalar morton_set loop vs morton_set_bulk4.
 *
 * Both legs run back-to-back in one process so CPU drift cancels;
 * each leg checksums its output to defeat dead-code elimination.
 * Prints both throughputs plus the paired speedup.
 */

#include "../test_common.h"
#include "../../include/ttypt/geo.h"
#include <stdio.h>

#define BULK4_N 1000000
#define BULK4_ROUNDS 7

static int16_t pts[BULK4_N][4];
static uint64_t out_s[BULK4_N], out_b[BULK4_N];

int main(void) {
    printf("\n%s%s=== Bulk4 Encode Benchmarks ===%s\n\n",
           COLOR_BOLD, COLOR_MAGENTA, COLOR_RESET);

#if GEO_SIMD_MORTON
    test_seed_rng(0xB94);
    for (int i = 0; i < BULK4_N; i++)
        for (int d = 0; d < 4; d++)
            pts[i][d] = test_rand_coord();

    /* Correctness first: bulk must match scalar exactly. */
    morton_set_bulk4(out_b, pts, BULK4_N);
    for (int i = 0; i < BULK4_N; i++)
        out_s[i] = morton_set_4(pts[i]);
    if (memcmp(out_s, out_b, sizeof out_s) != 0) {
        printf("bulk4 MISMATCH vs scalar\n");
        return 1;
    }
    printf("bulk4 correctness: OK\n");
    printf("note: the AVX2 bulk path needs -mavx2 at lib build time; "
           "without it bulk == scalar loop (ratio < 1 vs an -O3 caller "
           "loop is expected). See docs/PERF.md.\n");

    for (int r = 0; r < BULK4_ROUNDS; r++) {
        uint64_t t0 = get_time_usec();
        for (int i = 0; i < BULK4_N; i++)
            out_s[i] = morton_set_4(pts[i]);
        uint64_t t1 = get_time_usec();
        uint32_t n = morton_set_bulk4(out_b, pts, BULK4_N);
        uint64_t t2 = get_time_usec();
        (void)n;

        double s_s = (t1 - t0) / 1000000.0, s_b = (t2 - t1) / 1000000.0;
        uint64_t ck = out_s[BULK4_N - 1] + out_b[BULK4_N - 1];
        printf("[BENCH] bulk4 round %d: scalar %.2f M/s, bulk %.2f M/s, "
               "speedup %.3f (ck=%llu)\n",
               r, BULK4_N / s_s / 1e6, BULK4_N / s_b / 1e6,
               s_s / s_b, (unsigned long long)ck);
    }
#else
    printf("GEO_SIMD_MORTON=0: bulk4 API unavailable, nothing to bench.\n");
#endif

    printf("\n");
    return 0;
}
