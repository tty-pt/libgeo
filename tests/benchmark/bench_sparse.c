/*
 * Benchmark: Z-interval skip on a two-region layout.
 *
 * City A fills the query box; region B sits outside the box but inside
 * the box's morton interval (same x/y, z just above the box top), so a
 * linear walk must decode + discard every B key while the skipping walk
 * jumps whole box-disjoint cubes. Reports wall time and examined entries
 * (islet_last_scan_count) for the fill and raw paths.
 */

#include "../test_common.h"
#include "../../include/ttypt/islet.h"
#include "../../include/ttypt/point.h"
#include "../../include/ttypt/morton.h"
#include <stdlib.h>
#include <string.h>

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

    printf("\n%s%s=== Sparse Two-Region Skip Benchmarks ===%s\n\n",
           COLOR_BOLD, COLOR_MAGENTA, COLOR_RESET);

    uint32_t db = islet_open(NULL, "bench_sparse", 65535);

    /* City A: inside the query box [0,32]^3. */
    test_seed_rng(60606);
    for (int i = 0; i < 2000; i++) {
        int16_t p[3] = {
            test_rand_coord_range(0, 32),
            test_rand_coord_range(0, 32),
            test_rand_coord_range(0, 32),
        };
        islet_set_3(db, p, (uint32_t)i);
    }

    /* Region B: outside the box (z = 33..40), inside the morton interval. */
    for (int i = 0; i < 2000; i++) {
        int16_t p[3] = {
            test_rand_coord_range(0, 32),
            test_rand_coord_range(0, 32),
            test_rand_coord_range(33, 41),
        };
        islet_set_3(db, p, 100000u + (uint32_t)i);
    }

    int16_t s[3] = {0, 0, 0};
    uint16_t l[3] = {32, 32, 32};

    bench_start(&bench, "Fill bbox (two-region)");
    size_t nfill = 0;
    for (int r = 0; r < 20; r++) {
        rec_set_t *out = rec_set_new();
        rec_axis_fill_bbox_3(db, s, l, out);
        nfill = rec_set_count(out);
        rec_set_free(out);
    }
    bench_end(&bench, 20);
    printf("    (fill refs: %zu, examined per walk: %u)\n",
           nfill, islet_last_scan_count());

    bench_start(&bench, "Raw collect (two-region)");
    size_t nraw = 0;
    for (int r = 0; r < 20; r++) {
        uint32_t iter = islet_iter_3(db, s, l);
        int16_t p[3];
        uint32_t ref;
        nraw = 0;
        while (islet_next(p, &ref, iter))
            nraw++;
    }
    bench_end(&bench, 20);
    printf("    (raw entries: %zu, examined per walk: %u)\n",
           nraw, islet_last_scan_count());

    printf("\n");
    return 0;
}
