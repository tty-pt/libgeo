/**
 * @file test_axis.c
 * @brief Unit tests for the rec_query "islet" axis registration:
 *        registered fill dispatch through the axis table must equal the
 *        direct rec_axis_fill_bbox_N() calls, and the decode fn must
 *        round-trip "dim=N s=x,... l=x,..." param strings.
 */

#include "../test_common.h"
#include "../../include/ttypt/islet.h"
#include <ttypt/rec.h>
#include <string.h>

static void setup_once(void) {
    static int initialized = 0;
    if (!initialized) {
        islet_init();
        initialized = 1;
    }
}

static int find_islet_slot(void)
{
    int i;
    const rec_axis_t *axis;

    for (i = 0; i < rec_axis_count(); i++) {
        axis = rec_axis_get(i);
        if (axis && !strcmp(axis->name, "islet"))
            return i;
    }
    return -1;
}

TEST(axis_registered) {
    setup_once();
    int slot = find_islet_slot();
    ASSERT(slot >= 0);

    const rec_axis_t *axis = rec_axis_get(slot);
    ASSERT_NOT_NULL(axis);
    ASSERT_NOT_NULL((void *)axis->fill);
    ASSERT_NULL((void *)axis->rank);
    ASSERT_NOT_NULL((void *)axis->decode);
}

TEST(axis_fill_equiv_3d) {
    setup_once();
    int slot = find_islet_slot();
    ASSERT(slot >= 0);

    uint32_t db = islet_open(NULL, "test_axis_fill_3d", 1023);
    ASSERT_EQ(rec_axis_set_ctx(slot, (void *)(uintptr_t)db), 0);

    int16_t at[3] = {5, 5, 5};
    islet_put_3(db, at, 11);
    islet_put_3(db, at, 22);
    int16_t at2[3] = {50, 50, 50};
    islet_put_3(db, at2, 99);

    int16_t s[3] = {0, 0, 0};
    uint16_t l[3] = {10, 10, 10};

    rec_set_t *direct = rec_set_new();
    ASSERT_EQ(rec_axis_fill_bbox_3(db, s, l, direct), 0);

    const rec_axis_t *axis = rec_axis_get(slot);
    void *params = rec_axis_decode(slot, "dim=3 s=0,0,0 l=10,10,10");
    ASSERT_NOT_NULL(params);

    rec_set_t *via = rec_set_new();
    ASSERT_EQ(axis->fill(axis->ctx, params, via), 0);

    ASSERT_EQ(rec_set_count(direct), rec_set_count(via));
    ASSERT_EQ(rec_set_count(direct), (size_t)2);
    ASSERT(memcmp(rec_set_at(direct), rec_set_at(via),
                  rec_set_count(direct) * sizeof(rec_ref_t)) == 0);

    rec_set_free(direct);
    rec_set_free(via);
}

TEST(axis_fill_equiv_2d) {
    setup_once();
    int slot = find_islet_slot();
    ASSERT(slot >= 0);

    uint32_t db = islet_open(NULL, "test_axis_fill_2d", 1023);
    ASSERT_EQ(rec_axis_set_ctx(slot, (void *)(uintptr_t)db), 0);

    int16_t at[2] = {3, 4};
    islet_put_2(db, at, 7);
    int16_t at2[2] = {100, 100};
    islet_put_2(db, at2, 8);

    int16_t s[2] = {0, 0};
    uint16_t l[2] = {16, 16};

    rec_set_t *direct = rec_set_new();
    ASSERT_EQ(rec_axis_fill_bbox_2(db, s, l, direct), 0);

    const rec_axis_t *axis = rec_axis_get(slot);
    void *params = rec_axis_decode(slot, "dim=2 s=0,0 l=16,16");
    ASSERT_NOT_NULL(params);

    rec_set_t *via = rec_set_new();
    ASSERT_EQ(axis->fill(axis->ctx, params, via), 0);

    ASSERT_EQ(rec_set_count(direct), rec_set_count(via));
    ASSERT_EQ(rec_set_count(direct), (size_t)1);
    ASSERT(memcmp(rec_set_at(direct), rec_set_at(via),
                  rec_set_count(direct) * sizeof(rec_ref_t)) == 0);

    rec_set_free(direct);
    rec_set_free(via);
}

TEST(axis_decode_bad_dim) {
    setup_once();
    int slot = find_islet_slot();
    ASSERT(slot >= 0);

    ASSERT_NULL(rec_axis_decode(slot, "dim=0 s=0 l=1"));
    ASSERT_NULL(rec_axis_decode(slot, "dim=5 s=0,0,0,0,0 l=1,1,1,1,1"));
    ASSERT_NULL(rec_axis_decode(slot, "s=0,0 l=1,1"));
    ASSERT_NULL(rec_axis_decode(slot, "dim=2 l=1,1"));
    ASSERT_NULL(rec_axis_decode(slot, "dim=2 s=0,0"));
    ASSERT_NULL(rec_axis_decode(slot, "dim=3 s=0,0 l=1,1,1"));
    ASSERT_NULL(rec_axis_decode(-1, "dim=2 s=0,0 l=1,1"));
}

TEST(axis_fill_bad_dim_guard) {
    setup_once();
    int slot = find_islet_slot();
    ASSERT(slot >= 0);

    uint32_t db = islet_open(NULL, "test_axis_fill_guard", 1023);
    ASSERT_EQ(rec_axis_set_ctx(slot, (void *)(uintptr_t)db), 0);

    const rec_axis_t *axis = rec_axis_get(slot);
    struct { int16_t s[4]; uint16_t l[4]; int dim; } bad_params;
    memset(&bad_params, 0, sizeof(bad_params));
    bad_params.dim = 0;

    rec_set_t *out = rec_set_new();
    ASSERT_EQ(axis->fill(axis->ctx, &bad_params, out), -1);
    ASSERT_EQ(axis->fill(axis->ctx, NULL, out), -1);

    rec_set_free(out);
}

TEST(axis_open_matches_direct_open) {
    setup_once();
    int slot = find_islet_slot();
    ASSERT(slot >= 0);

    void *ctx = rec_axis_open(":test_axis_open_db:1023");
    ASSERT_EQ(rec_axis_set_ctx(slot, ctx), 0);
    uint32_t db = (uint32_t)(uintptr_t)ctx;

    int16_t at[3] = {1, 2, 3};
    islet_put_3(db, at, 55);

    int16_t s[3] = {0, 0, 0};
    uint16_t l[3] = {10, 10, 10};

    rec_set_t *direct = rec_set_new();
    ASSERT_EQ(rec_axis_fill_bbox_3(db, s, l, direct), 0);

    const rec_axis_t *axis = rec_axis_get(slot);
    void *params = rec_axis_decode(slot, "dim=3 s=0,0,0 l=10,10,10");
    ASSERT_NOT_NULL(params);

    rec_set_t *via = rec_set_new();
    ASSERT_EQ(axis->fill(axis->ctx, params, via), 0);

    ASSERT_EQ(rec_set_count(direct), rec_set_count(via));
    ASSERT_EQ(rec_set_count(direct), (size_t)1);

    rec_set_free(direct);
    rec_set_free(via);
}

TEST(axis_open_empty_spec_defaults) {
    setup_once();

    /* empty spec -> filename=NULL, database=NULL, mask=0 (in-memory,
     * qmap-default mask) -- must not crash and must return a usable
     * handle. */
    void *ctx = rec_axis_open("");
    uint32_t db = (uint32_t)(uintptr_t)ctx;
    int16_t at[1] = {9};
    islet_put_1(db, at, 1);

    int16_t s[1] = {0};
    uint16_t l[1] = {20};
    rec_set_t *out = rec_set_new();
    ASSERT_EQ(rec_axis_fill_bbox_1(db, s, l, out), 0);
    ASSERT_EQ(rec_set_count(out), (size_t)1);
    rec_set_free(out);
}

int main(void) {
    test_suite_begin("Islet rec_query Axis Registration Tests");
    RUN_TEST(axis_registered);
    RUN_TEST(axis_fill_equiv_3d);
    RUN_TEST(axis_fill_equiv_2d);
    RUN_TEST(axis_decode_bad_dim);
    RUN_TEST(axis_fill_bad_dim_guard);
    RUN_TEST(axis_open_matches_direct_open);
    RUN_TEST(axis_open_empty_spec_defaults);
    return test_suite_end();
}
