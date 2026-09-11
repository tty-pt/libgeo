/**
 * @file test_dim4_props.c
 * @brief Property-style tests for 4D support: randomized codec
 * round-trips, bulk4-vs-scalar equivalence, and randomized box queries
 * vs a brute-force oracle.
 */

#include "../test_common.h"
#include "../../include/ttypt/islet.h"
#include "../../include/ttypt/morton.h"
#include "../../include/ttypt/point.h"

static void setup_once(void) {
    static int initialized = 0;
    if (!initialized) {
        islet_init();
        initialized = 1;
    }
}

TEST(dim4_props_round_trip) {
    test_seed_rng(0xD14D4);

    for (int i = 0; i < 10000; i++) {
        int16_t p[4] = {
            test_rand_coord(), test_rand_coord(),
            test_rand_coord(), test_rand_coord()
        };
        int16_t d[4];

        morton_get_4(d, morton_set_4(p));
        ASSERT_POINT_EQ(p, d, 4);
    }
}

TEST(dim4_props_codes_cover_full_u64) {
    /* Random 4D codes must routinely use the top 16 bits (unlike 3D,
     * which reserves them): assert at least one sample in 1000 has
     * bit 63 set. Deterministic seed keeps this stable. */
    test_seed_rng(0xF00D);
    int saw_top = 0;

    for (int i = 0; i < 1000; i++) {
        int16_t p[4] = {
            test_rand_coord(), test_rand_coord(),
            test_rand_coord(), test_rand_coord()
        };
        if (morton_set_4(p) >> 63)
            saw_top = 1;
    }
    ASSERT(saw_top);
}

#if ISLET_SIMD_MORTON
TEST(dim4_props_bulk4_matches_scalar) {
    #define N4 1031
    static int16_t pts[N4][4];
    static uint64_t out[N4];

    test_seed_rng(0xB014);
    for (int i = 0; i < N4; i++)
        for (int d = 0; d < 4; d++)
            pts[i][d] = test_rand_coord();

    ASSERT_EQ(morton_set_bulk4(out, pts, N4), (uint32_t)N4);
    for (int i = 0; i < N4; i++)
        ASSERT_EQ(out[i], morton_set_4(pts[i]));
    #undef N4
}
#endif

TEST(dim4_props_box_oracle) {
    setup_once();
    uint32_t db = islet_open(NULL, "test_dim4_props_box", 8191);

    #define C4N 1500
    static int16_t cloud[C4N][4];
    static uint32_t refs[C4N];
    test_seed_rng(0xC10D);
    for (int i = 0; i < C4N; i++) {
        for (int d = 0; d < 4; d++)
            cloud[i][d] = (int16_t)(test_rand64() % 24);
        refs[i] = (uint32_t)i;
        islet_put_4(db, cloud[i], refs[i]);
    }

    test_seed_rng(0xB0B);
    for (int trial = 0; trial < 30; trial++) {
        int16_t s[4];
        uint16_t l[4];
        for (int d = 0; d < 4; d++) {
            s[d] = (int16_t)(test_rand64() % 24);
            l[d] = (uint16_t)(test_rand64() % 10);
        }

        int expect = 0;
        for (int i = 0; i < C4N; i++) {
            int in = 1;
            for (int d = 0; d < 4; d++) {
                int32_t lo = s[d], hi = (int32_t)s[d] + l[d];
                if (cloud[i][d] < lo || cloud[i][d] > hi) { in = 0; break; }
            }
            if (in) expect++;
        }

        uint32_t it = islet_iter_4(db, s, l);
        int16_t pt[4];
        uint32_t ref;
        int got = 0;
        while (islet_next(pt, &ref, it))
            got++;
        ASSERT_EQ(got, expect);
    }
    #undef C4N
}

int main(void) {
    test_suite_begin("Dim4 Property Tests");

    RUN_TEST(dim4_props_round_trip);
    RUN_TEST(dim4_props_codes_cover_full_u64);
#if ISLET_SIMD_MORTON
    RUN_TEST(dim4_props_bulk4_matches_scalar);
#endif
    RUN_TEST(dim4_props_box_oracle);

    return test_suite_end();
}
