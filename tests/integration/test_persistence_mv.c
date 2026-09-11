/*
 * Integration tests for MV persistence in libgeo.
 * - A file written with unique keys reopens and reads fine under the MV build.
 * - An MV-written file (duplicate values at a cell) round-trips: all
 *   siblings survive save/close/reopen and read back via every path.
 */

#include "../test_common.h"
#include "../../include/ttypt/geo.h"
#include "../../include/ttypt/point.h"
#include "../../include/ttypt/morton.h"
#include "../../include/ttypt/qmap.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static void setup_once(void) {
    static int initialized = 0;
    if (!initialized) {
        geo_init();
        initialized = 1;
    }
}

#define MV_F1 "/tmp/test_geo_mv_compat.db"
#define MV_F2 "/tmp/test_geo_mv_roundtrip.db"

/* Unique-key file reopens under the MV build with all cells readable */
TEST(mv_compat_unique_keys_reopen) {
    setup_once();
    unlink(MV_F1);

    uint32_t db = geo_open(MV_F1, "compat", 1023);

    for (int i = 0; i < 50; i++) {
        int16_t p[3] = {i, i * 2, i * 3};
        geo_put_3(db, p, 1000 + i);
    }
    qmap_save();
    qmap_close(db);

    db = geo_open(MV_F1, "compat", 1023);

    int verified = 0;
    for (int i = 0; i < 50; i++) {
        int16_t p[3] = {i, i * 2, i * 3};
        if (geo_get_3(db, p) == (uint32_t)(1000 + i))
            verified++;
        ASSERT_EQ(geo_cell_count_3(db, p), 1);
    }
    ASSERT_EQ(verified, 50);
    qmap_close(db);
    unlink(MV_F1);
}

/* MV duplicates survive save/close/reopen on every read path */
TEST(mv_roundtrip_duplicates) {
    setup_once();
    unlink(MV_F2);

    uint32_t db = geo_open(MV_F2, "mvdata", 1023);

    int16_t a[3] = {5, 5, 5};
    int16_t b[3] = {6, 6, 6};
    geo_put_3(db, a, 11);
    geo_put_3(db, a, 22);
    geo_put_3(db, b, 33);
    qmap_save();
    qmap_close(db);

    db = geo_open(MV_F2, "mvdata", 1023);

    /* get-first + count path */
    ASSERT_EQ(geo_get_3(db, a), 11);
    ASSERT_EQ(geo_cell_count_3(db, a), 2);
    ASSERT_EQ(geo_get_3(db, b), 33);

    /* chain path */
    uint32_t cur = geo_get_multi_3(db, a);
    ASSERT(cur != QM_MISS);
    uint32_t ref;
    ASSERT_EQ(geo_cell_next(&ref, cur), 1);
    ASSERT_EQ(ref, 11);
    ASSERT_EQ(geo_cell_next(&ref, cur), 1);
    ASSERT_EQ(ref, 22);
    ASSERT_EQ(geo_cell_next(&ref, cur), 0);

    /* raw iterator path: both siblings exactly once */
    int16_t start[3] = {0, 0, 0};
    uint16_t len[3] = {10, 10, 10};
    uint32_t iter = geo_iter_3(db, start, len);
    int16_t p[3];
    int total = 0, n11 = 0, n22 = 0, n33 = 0;
    while (geo_next(p, &ref, iter)) {
        total++;
        if (ref == 11) n11++;
        else if (ref == 22) n22++;
        else if (ref == 33) n33++;
        else ASSERT(0);
    }
    ASSERT_EQ(total, 3);
    ASSERT_EQ(n11, 1);
    ASSERT_EQ(n22, 1);
    ASSERT_EQ(n33, 1);

    qmap_close(db);
    unlink(MV_F2);
}

int main(void) {
    test_suite_begin("Geo MV Persistence Integration Tests");
    RUN_TEST(mv_compat_unique_keys_reopen);
    RUN_TEST(mv_roundtrip_duplicates);
    return test_suite_end();
}
