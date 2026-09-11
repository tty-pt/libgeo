/**
 * @file bench_morton.c
 * @brief Benchmark for Morton encoding/decoding performance.
 */

#include "../test_common.h"
#include "../../include/ttypt/morton.h"
#include "../../include/ttypt/pointcfg.h"
#include <stdio.h>

#define BENCH_ITERATIONS 1000000
#define BENCH_PDEP_ROUNDS 7

/* Scalar reference: always spread3/compact_axis path, used to verify
 * bit-identity and measure the PDEP speedup when built with -mbmi2. */
static inline uint64_t
morton_set_3_ref(int16_t *p)
{
	uint16_t up0 = islet_unsign(p[0]);
	uint16_t up1 = islet_unsign(p[1]);
	uint16_t up2 = islet_unsign(p[2]);
	return islet_spread3(up0)
		| (islet_spread3(up1) << 1)
		| (islet_spread3(up2) << 2);
}

static inline void
morton_get_3_ref(int16_t *pos, uint64_t code)
{
	uint32_t uup[] = { 0, 0, 0 };
	islet_decode3(code, &uup[0], &uup[1], &uup[2]);
	pos[0] = islet_sign((uint16_t)uup[0]);
	pos[1] = islet_sign((uint16_t)uup[1]);
	pos[2] = islet_sign((uint16_t)uup[2]);
}

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
        code = morton_set_3(pos);
        (void)code; /* Prevent optimization */
    }
    bench_end(&bench, BENCH_ITERATIONS);
    
    /* Benchmark decoding */
    bench_start(&bench, "Morton Decode (3D)");
    test_seed_rng(42);
    for (int i = 0; i < BENCH_ITERATIONS; i++) {
        code = test_rand64();
        morton_get_3(pos, code);
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
        code = morton_set_3(pos);
        morton_get_3(decoded, code);
    }
    bench_end(&bench, BENCH_ITERATIONS);

    /* Round-trip through the config object: same implementations, but
     * one indirect call per member (the struct-vs-inline tax). */
    bench_start(&bench, "Morton Round-Trip (Point3_2)");
    test_seed_rng(42);
    for (int i = 0; i < BENCH_ITERATIONS; i++) {
        pos[0] = test_rand_coord();
        pos[1] = test_rand_coord();
        pos[2] = test_rand_coord();
        code = Point3_2.morton_set(pos);
        Point3_2.morton_get(decoded, code);
    }
    bench_end(&bench, BENCH_ITERATIONS);

    /* Paired PDEP vs scalar: same TU, same pre-generated inputs (so RNG
     * cost doesn't dilute the codec measurement), bit-identity
     * enforced.  Only active when built with -mbmi2 (the public
     * morton_set_3 uses PDEP; the _ref functions use spread3). */
#if defined(__BMI2__) && ISLET_USE_PDEP
    printf("\n%s%s=== BMI2 PDEP vs Scalar (paired) ===%s\n\n",
           COLOR_BOLD, COLOR_MAGENTA, COLOR_RESET);

    static int16_t pdep_pts[BENCH_ITERATIONS][3];
    static uint64_t pdep_codes[BENCH_ITERATIONS];
    static uint64_t pdep_out_s[BENCH_ITERATIONS], pdep_out_p[BENCH_ITERATIONS];
    static int16_t pdep_dec_s[BENCH_ITERATIONS][3], pdep_dec_p[BENCH_ITERATIONS][3];

    test_seed_rng(0xDAD);
    for (int i = 0; i < BENCH_ITERATIONS; i++) {
        pdep_pts[i][0] = test_rand_coord();
        pdep_pts[i][1] = test_rand_coord();
        pdep_pts[i][2] = test_rand_coord();
        pdep_codes[i] = test_rand64();
    }

    /* Correctness: bit-identity check over the full pre-generated set */
    for (int i = 0; i < BENCH_ITERATIONS; i++) {
        uint64_t c_fast = morton_set_3(pdep_pts[i]);
        uint64_t c_ref  = morton_set_3_ref(pdep_pts[i]);
        if (c_fast != c_ref) {
            printf("MISMATCH at i=%d: pdep=%llu ref=%llu\n",
                   i, (unsigned long long)c_fast,
                   (unsigned long long)c_ref);
            return 1;
        }
        int16_t d_fast[3], d_ref[3];
        morton_get_3(d_fast, c_fast);
        morton_get_3_ref(d_ref, c_ref);
        if (d_fast[0] != d_ref[0] || d_fast[1] != d_ref[1] ||
            d_fast[2] != d_ref[2]) {
            printf("DECODE MISMATCH at i=%d\n", i);
            return 1;
        }
    }
    printf("BMI2 bit-identity: OK\n");

    for (int r = 0; r < BENCH_PDEP_ROUNDS; r++) {
        uint64_t t0, t1, t2;

        /* encode: scalar vs pdep over the same pre-generated points */
        t0 = get_time_usec();
        for (int i = 0; i < BENCH_ITERATIONS; i++)
            pdep_out_s[i] = morton_set_3_ref(pdep_pts[i]);
        t1 = get_time_usec();
        for (int i = 0; i < BENCH_ITERATIONS; i++)
            pdep_out_p[i] = morton_set_3(pdep_pts[i]);
        t2 = get_time_usec();

        double s_s = (t1 - t0) / 1e6, s_p = (t2 - t1) / 1e6;
        uint64_t ck = pdep_out_s[BENCH_ITERATIONS - 1]
                    + pdep_out_p[BENCH_ITERATIONS - 1];
        printf("[BENCH] encode  round %d: scalar %.2f M/s, pdep %.2f M/s, "
               "speedup %.2fx (ck=%llu)\n",
               r, BENCH_ITERATIONS / s_s / 1e6,
               BENCH_ITERATIONS / s_p / 1e6, s_s / s_p,
               (unsigned long long)ck);

        /* decode: scalar vs pdep over the same pre-generated codes */
        t0 = get_time_usec();
        for (int i = 0; i < BENCH_ITERATIONS; i++)
            morton_get_3_ref(pdep_dec_s[i], pdep_codes[i]);
        t1 = get_time_usec();
        for (int i = 0; i < BENCH_ITERATIONS; i++)
            morton_get_3(pdep_dec_p[i], pdep_codes[i]);
        t2 = get_time_usec();

        s_s = (t1 - t0) / 1e6;
        s_p = (t2 - t1) / 1e6;
        ck = (uint64_t)(uint16_t)pdep_dec_s[BENCH_ITERATIONS - 1][0]
           + (uint64_t)(uint16_t)pdep_dec_p[BENCH_ITERATIONS - 1][0];
        printf("[BENCH] decode  round %d: scalar %.2f M/s, pdep %.2f M/s, "
               "speedup %.2fx (ck=%llu)\n",
               r, BENCH_ITERATIONS / s_s / 1e6,
               BENCH_ITERATIONS / s_p / 1e6, s_s / s_p,
               (unsigned long long)ck);

        /* round-trip: scalar vs pdep over the same pre-generated points */
        t0 = get_time_usec();
        for (int i = 0; i < BENCH_ITERATIONS; i++) {
            uint64_t c = morton_set_3_ref(pdep_pts[i]);
            morton_get_3_ref(pdep_dec_s[i], c);
        }
        t1 = get_time_usec();
        for (int i = 0; i < BENCH_ITERATIONS; i++) {
            uint64_t c = morton_set_3(pdep_pts[i]);
            morton_get_3(pdep_dec_p[i], c);
        }
        t2 = get_time_usec();

        s_s = (t1 - t0) / 1e6;
        s_p = (t2 - t1) / 1e6;
        printf("[BENCH] roundtrip round %d: scalar %.2f M/s, pdep %.2f M/s, "
               "speedup %.2fx\n",
               r, BENCH_ITERATIONS / s_s / 1e6,
               BENCH_ITERATIONS / s_p / 1e6, s_s / s_p);
    }
#else
    printf("\n(BMI2 PDEP rows require -mbmi2 at bench compile time)\n");
#endif
    
    printf("\n");
    return 0;
}
