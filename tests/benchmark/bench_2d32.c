/*
 * Benchmark for the 2D x 32-bit dense config (islet_*_2_32).
 */

#include "../test_common.h"
#include "../../include/ttypt/islet.h"
#include "../../include/ttypt/point.h"
#include "../../include/ttypt/morton.h"
#include <stdlib.h>
#include <string.h>

#define BENCH_ITERATIONS 1000000

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

    printf("\n%s%s=== 2D x 32-bit Config Benchmarks ===%s\n\n",
           COLOR_BOLD, COLOR_MAGENTA, COLOR_RESET);

    /* Codec loop over wide lanes */
    int32_t pos[2];
    uint64_t code;
    bench_start(&bench, "Morton Encode+Decode (2D32)");
    test_seed_rng(42);
    for (int i = 0; i < BENCH_ITERATIONS; i++) {
        pos[0] = (int32_t)test_rand64();
        pos[1] = (int32_t)test_rand64();
        code = morton_set_2_32(pos);
        morton_get_2_32(pos, code);
    }
    bench_end(&bench, BENCH_ITERATIONS);

    /* Scatter store over the wide world */
    uint32_t db = islet_open(NULL, "bench_2d32", 8191);
    bench_start(&bench, "Scatter Put (2D32)");
    test_seed_rng(43);
    for (int i = 0; i < 100000; i++) {
        pos[0] = (int32_t)test_rand64();
        pos[1] = (int32_t)test_rand64();
        islet_put_2_32(db, pos, (uint32_t)i);
    }
    bench_end(&bench, 100000);

    /* Scatter lookup over the same keys */
    bench_start(&bench, "Scatter Get (2D32)");
    test_seed_rng(43);
    for (int i = 0; i < 100000; i++) {
        pos[0] = (int32_t)test_rand64();
        pos[1] = (int32_t)test_rand64();
        (void)islet_get_2_32(db, pos);
    }
    bench_end(&bench, 100000);

    /* Box fill over a dense region */
    for (int32_t x = 0; x < 256; x++)
        for (int32_t y = 0; y < 256; y++) {
            int32_t p[2] = { x, y };
            islet_put_2_32(db, p, (uint32_t)(x * 256 + y));
        }
    int32_t s[2] = { 0, 0 };
    int32_t l[2] = { 256, 256 };
    bench_start(&bench, "Box Fill 256x256 (2D32)");
    for (int i = 0; i < 20; i++) {
        rec_set_t *o = rec_set_new();
        rec_axis_fill_bbox_2_32(db, s, l, o);
        rec_set_free(o);
    }
    bench_end(&bench, 20);

    printf("\n");
    return 0;
}
