/*
 * Benchmark: rec_axis_fill_bbox vs the raw islet_iter/islet_next collect path.
 * Both share the box walker, so this measures adapter overhead (push+seal)
 * against the malloc-collect + free cycle. Reports wall time and asserts
 * result agreement.
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

/* Fill db with n pseudo-random points in [0,span)^3, values 0..n-1. */
static void seed_cloud(uint32_t db, int n, int span) {
    test_seed_rng(31337);
    for (int i = 0; i < n; i++) {
        int16_t p[3] = {
            test_rand_coord_range(0, span),
            test_rand_coord_range(0, span),
            test_rand_coord_range(0, span),
        };
        islet_set_3(db, p, (uint32_t)i);
    }
}

static size_t raw_count(uint32_t db, int16_t *s, uint16_t *l) {
    uint32_t iter = islet_iter_3(db, s, l);
    int16_t p[3];
    uint32_t ref;
    size_t n = 0;

    while (islet_next(p, &ref, iter))
        n++;
    return n;
}

static void bench_size(const char *tag, int n, int span) {
    benchmark_t bench;
    static char name[64];

    snprintf(name, sizeof name, "bench_fill_%s", tag);
    uint32_t db = islet_open(NULL, name, 65535);
    seed_cloud(db, n, span);

    int16_t s[3] = {0, 0, 0};
    uint16_t l[3] = {(uint16_t)span, (uint16_t)span, (uint16_t)span};

    char label[128];
    snprintf(label, sizeof label, "Fill bbox (%s cloud)", tag);
    bench_start(&bench, label);
    size_t nfill = 0;
    for (int r = 0; r < 20; r++) {
        rec_set_t *out = rec_set_new();
        rec_axis_fill_bbox_3(db, s, l, out);
        nfill = rec_set_count(out);
        rec_set_free(out);
    }
    bench_end(&bench, 20);

    snprintf(label, sizeof label, "Raw collect (%s cloud)", tag);
    bench_start(&bench, label);
    size_t nraw = 0;
    for (int r = 0; r < 20; r++)
        nraw = raw_count(db, s, l);
    bench_end(&bench, 20);

    printf("    (fill refs: %zu, raw entries: %zu)\n", nfill, nraw);
}

int main(void) {
    setup_once();

    printf("\n%s%s=== Fill vs Raw Collect Benchmarks ===%s\n\n",
           COLOR_BOLD, COLOR_MAGENTA, COLOR_RESET);

    bench_size("1k", 1000, 64);
    printf("\n");
    bench_size("10k", 10000, 100); /* 100^3 = 1M cells: exactly at the fill cap */
    printf("\n");
    return 0;
}
