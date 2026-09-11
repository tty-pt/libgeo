/**
 * @file bench_bulk_decode.c
 * @brief Paired benchmark: scalar morton_get loop vs morton_get_bulk/bulk4.
 *
 * Both legs run back-to-back in one process so CPU drift cancels;
 * each leg checksums its output to defeat dead-code elimination.
 * Prints both throughputs plus the paired speedup.
 */

#include "../test_common.h"
#include "../../include/ttypt/geo.h"
#include <stdio.h>

#define BULK_N 1000000
#define BULK_ROUNDS 7

static uint64_t codes_3d[BULK_N], codes_4d[BULK_N];
static int16_t pts_s3[BULK_N][3], pts_b3[BULK_N][3];
static int16_t pts_s4[BULK_N][4], pts_b4[BULK_N][4];

int main(void) {
    printf("\n%s%s=== Bulk Decode Benchmarks ===%s\n\n",
           COLOR_BOLD, COLOR_MAGENTA, COLOR_RESET);

#if GEO_SIMD_MORTON
    test_seed_rng(0xA37);
    for (int i = 0; i < BULK_N; i++) {
        for (int d = 0; d < 3; d++)
            pts_s3[i][d] = test_rand_coord();
        codes_3d[i] = morton_set_3(pts_s3[i]);
        for (int d = 0; d < 4; d++)
            pts_s4[i][d] = test_rand_coord();
        codes_4d[i] = morton_set_4(pts_s4[i]);
    }

    /* --- 3D decode --- */

    /* Correctness: bulk must match scalar exactly */
    morton_get_bulk(pts_b3, codes_3d, BULK_N);
    for (int i = 0; i < BULK_N; i++)
        morton_get_3(pts_s3[i], codes_3d[i]);
    if (memcmp(pts_s3, pts_b3, sizeof pts_s3) != 0) {
        printf("bulk3 DECODE MISMATCH\n");
        return 1;
    }
    printf("bulk3 decode correctness: OK\n");

    for (int r = 0; r < BULK_ROUNDS; r++) {
        uint64_t t0 = get_time_usec();
        for (uint32_t i = 0; i < BULK_N; i++)
            morton_get_3(pts_s3[i], codes_3d[i]);
        uint64_t t1 = get_time_usec();
        uint32_t n = morton_get_bulk(pts_b3, codes_3d, BULK_N);
        uint64_t t2 = get_time_usec();
        (void)n;

        double s_s = (t1 - t0) / 1e6, s_b = (t2 - t1) / 1e6;
        uint64_t ck = pts_s3[BULK_N - 1][0] + pts_b3[BULK_N - 1][0];
        printf("[BENCH] bulk3 decode round %d: scalar %.2f M/s, bulk %.2f M/s, "
               "speedup %.3f (ck=%llu)\n",
               r, BULK_N / s_s / 1e6, BULK_N / s_b / 1e6,
               s_s / s_b, (unsigned long long)ck);
    }

    /* --- 4D decode --- */

    morton_get_bulk4(pts_b4, codes_4d, BULK_N);
    for (int i = 0; i < BULK_N; i++)
        morton_get_4(pts_s4[i], codes_4d[i]);
    if (memcmp(pts_s4, pts_b4, sizeof pts_s4) != 0) {
        printf("bulk4 DECODE MISMATCH\n");
        return 1;
    }
    printf("bulk4 decode correctness: OK\n");

    for (int r = 0; r < BULK_ROUNDS; r++) {
        uint64_t t0 = get_time_usec();
        for (uint32_t i = 0; i < BULK_N; i++)
            morton_get_4(pts_s4[i], codes_4d[i]);
        uint64_t t1 = get_time_usec();
        uint32_t n = morton_get_bulk4(pts_b4, codes_4d, BULK_N);
        uint64_t t2 = get_time_usec();
        (void)n;

        double s_s = (t1 - t0) / 1e6, s_b = (t2 - t1) / 1e6;
        uint64_t ck = pts_s4[BULK_N - 1][0] + pts_b4[BULK_N - 1][0];
        printf("[BENCH] bulk4 decode round %d: scalar %.2f M/s, bulk %.2f M/s, "
               "speedup %.3f (ck=%llu)\n",
               r, BULK_N / s_s / 1e6, BULK_N / s_b / 1e6,
               s_s / s_b, (unsigned long long)ck);
    }
#else
    printf("GEO_SIMD_MORTON=0: bulk API unavailable, nothing to bench.\n");
#endif

    printf("\n");
    return 0;
}
