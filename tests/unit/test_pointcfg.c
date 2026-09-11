/**
 * @file test_pointcfg.c
 * @brief Unit tests for the point config objects (pointcfg.h):
 *        Point1_2, Point2_2, Point3_2, Point4_2, Point2_4.
 */

#include "../test_common.h"
#include "../../include/ttypt/pointcfg.h"
#include "../../include/ttypt/point.h"
#include "../../include/ttypt/morton.h"
#include <unistd.h>
#include <string.h>

static void setup_once(void) {
    static int initialized = 0;
    if (!initialized) {
        geo_init();
        initialized = 1;
    }
}

static uint32_t fresh_db(void) {
    return geo_open(NULL, NULL, 0xFF);
}

static void check_i16(const int16_t *a, const int16_t *b, int dim) {
    for (int i = 0; i < dim; i++)
        ASSERT_EQ(a[i], b[i]);
}

/* All five objects must have every member populated. */
TEST(members_present) {
    ASSERT_NOT_NULL(Point1_2.morton_set);
    ASSERT_NOT_NULL(Point1_2.morton_get);
    ASSERT_NOT_NULL(Point2_2.morton_set);
    ASSERT_NOT_NULL(Point3_2.morton_set);
    ASSERT_NOT_NULL(Point4_2.morton_set);
    ASSERT_NOT_NULL(Point2_4.morton_set);

    const geo_point2b_t *p2b[4] = {
        &Point1_2, &Point2_2, &Point3_2, &Point4_2
    };
    for (int c = 0; c < 4; c++) {
        const geo_point2b_t *cfg = p2b[c];
        ASSERT_NOT_NULL(cfg->add);
        ASSERT_NOT_NULL(cfg->sub);
        ASSERT_NOT_NULL(cfg->min);
        ASSERT_NOT_NULL(cfg->max);
        ASSERT_NOT_NULL(cfg->copy);
        ASSERT_NOT_NULL(cfg->vol);
        ASSERT_NOT_NULL(cfg->set);
        ASSERT_NOT_NULL(cfg->debug);
        ASSERT_NOT_NULL(cfg->idx);
        ASSERT_NOT_NULL(cfg->put);
        ASSERT_NOT_NULL(cfg->get);
        ASSERT_NOT_NULL(cfg->replace);
        ASSERT_NOT_NULL(cfg->del);
        ASSERT_NOT_NULL(cfg->del_all);
        ASSERT_NOT_NULL(cfg->cell_count);
        ASSERT_NOT_NULL(cfg->get_multi);
        ASSERT_NOT_NULL(cfg->iter);
        ASSERT_NOT_NULL(cfg->fill_bbox);
        ASSERT_NOT_NULL(cfg->next);
    }
    ASSERT_NOT_NULL(Point2_4.add);
    ASSERT_NOT_NULL(Point2_4.put);
    ASSERT_NOT_NULL(Point2_4.iter);
    ASSERT_NOT_NULL(Point2_4.next);
    ASSERT_NOT_NULL(Point2_4.fill_bbox);
}

/* Shared exported symbols must be wired through intact. */
TEST(symbol_aliasing) {
    ASSERT(Point1_2.next == geo_next);
    ASSERT(Point3_2.next == geo_next);
    ASSERT(Point4_2.next == geo_next);
    ASSERT(Point2_4.next == geo_next32);

    ASSERT(Point1_2.iter == geo_iter_1);
    ASSERT(Point2_2.iter == geo_iter_2);
    ASSERT(Point3_2.iter == geo_iter_3);
    ASSERT(Point4_2.iter == geo_iter_4);
    ASSERT(Point2_4.iter == geo_iter_2_32);

    ASSERT(Point3_2.fill_bbox == rec_axis_fill_bbox_3);
    ASSERT(Point2_4.fill_bbox == rec_axis_fill_bbox_2_32);

    ASSERT(Point3_2.get_multi == geo_get_multi_3);
    ASSERT(Point2_4.get_multi == geo_get_multi_2_32);
}

/* Capture one stderr debug print and check its content. */
static void check_debug_p2b(void (*dbg)(char *, int16_t *)) {
    char tmpl[] = "/tmp/libgeo_pcfgd_XXXXXX";
    int16_t p[4] = { 10, 20, 30, 40 };
    int fd = fileno(stderr);
    int saved = dup(fd);
    int tfd = mkstemp(tmpl);
    ASSERT(tfd != -1);
    fflush(stderr);
    dup2(tfd, fd);
    dbg("pt", p);
    fflush(stderr);
    dup2(saved, fd);
    close(saved);
    close(tfd);
    FILE *f = fopen(tmpl, "r");
    char buf[128];
    size_t n = fread(buf, 1, sizeof buf - 1, f);
    buf[n] = '\0';
    fclose(f);
    unlink(tmpl);
    ASSERT(strstr(buf, "pt(") != NULL);
    ASSERT(strstr(buf, "(10") != NULL);
}

static void check_debug_p4b(void (*dbg)(char *, int32_t *)) {
    char tmpl[] = "/tmp/libgeo_pcfgd_XXXXXX";
    int32_t p[2] = { 7, -3 };
    int fd = fileno(stderr);
    int saved = dup(fd);
    int tfd = mkstemp(tmpl);
    ASSERT(tfd != -1);
    fflush(stderr);
    dup2(tfd, fd);
    dbg("pt", p);
    fflush(stderr);
    dup2(saved, fd);
    close(saved);
    close(tfd);
    FILE *f = fopen(tmpl, "r");
    char buf[128];
    size_t n = fread(buf, 1, sizeof buf - 1, f);
    buf[n] = '\0';
    fclose(f);
    unlink(tmpl);
    ASSERT(strstr(buf, "pt(7, -3)") != NULL);
}

/* Exercise every member of one int16 config against its flat inlines. */
static void exercise_p2b(const geo_point2b_t *cfg, int dim) {
    uint32_t db = fresh_db();

    /* --- codec: parity with the flat inlines, then decode round-trip */
    int16_t p[4] = { 1, 2, 3, 4 };
    uint64_t code;
    switch (dim) {
    case 1: code = morton_set_1(p); ASSERT_EQ(cfg->morton_set(p), code); break;
    case 2: code = morton_set_2(p); ASSERT_EQ(cfg->morton_set(p), code); break;
    case 3: code = morton_set_3(p); ASSERT_EQ(cfg->morton_set(p), code); break;
    case 4: code = morton_set_4(p); ASSERT_EQ(cfg->morton_set(p), code); break;
    }
    int16_t back[4];
    cfg->morton_get(back, code);
    check_i16(back, p, dim);

    /* --- vector utilities */
    int16_t a[4] = { 100, 200, 300, 400 };
    int16_t b[4] = { 5, 6, 7, 8 };
    int16_t r1[4], r2[4];

    cfg->add(r1, a, b);
    switch (dim) {
    case 1: point_add_1(r2, a, b); break;
    case 2: point_add_2(r2, a, b); break;
    case 3: point_add_3(r2, a, b); break;
    case 4: point_add_4(r2, a, b); break;
    }
    check_i16(r1, r2, dim);

    cfg->sub(r1, a, b);
    switch (dim) {
    case 1: point_sub_1(r2, a, b); break;
    case 2: point_sub_2(r2, a, b); break;
    case 3: point_sub_3(r2, a, b); break;
    case 4: point_sub_4(r2, a, b); break;
    }
    check_i16(r1, r2, dim);

    cfg->min(r1, a, b);
    switch (dim) {
    case 1: point_min_1(r2, a, b); break;
    case 2: point_min_2(r2, a, b); break;
    case 3: point_min_3(r2, a, b); break;
    case 4: point_min_4(r2, a, b); break;
    }
    check_i16(r1, r2, dim);

    cfg->max(r1, a, b);
    switch (dim) {
    case 1: point_max_1(r2, a, b); break;
    case 2: point_max_2(r2, a, b); break;
    case 3: point_max_3(r2, a, b); break;
    case 4: point_max_4(r2, a, b); break;
    }
    check_i16(r1, r2, dim);

    cfg->copy(r1, a);
    check_i16(r1, a, dim);

    cfg->set(r2, 7);
    for (int i = 0; i < dim; i++)
        ASSERT_EQ(r2[i], 7);

    int16_t dims[4] = { 2, 3, 4, 5 };
    switch (dim) {
    case 1: ASSERT_EQ(cfg->vol(dims), point_vol_1(dims)); break;
    case 2: ASSERT_EQ(cfg->vol(dims), point_vol_2(dims)); break;
    case 3: ASSERT_EQ(cfg->vol(dims), point_vol_3(dims)); break;
    case 4: ASSERT_EQ(cfg->vol(dims), point_vol_4(dims)); break;
    }

    int16_t s[4] = { 0, 0, 0, 0 };
    int16_t e[4] = { 3, 3, 3, 3 };
    switch (dim) {
    case 1: ASSERT_EQ(cfg->idx(p, s, e), point_idx_1(p, s, e)); break;
    case 2: ASSERT_EQ(cfg->idx(p, s, e), point_idx_2(p, s, e)); break;
    case 3: ASSERT_EQ(cfg->idx(p, s, e), point_idx_3(p, s, e)); break;
    case 4: ASSERT_EQ(cfg->idx(p, s, e), point_idx_4(p, s, e)); break;
    }

    check_debug_p2b(cfg->debug);

    /* --- db ops: multi-value CRUD through the object */
    ASSERT_EQ(cfg->get(db, p), GEO_MISS);
    cfg->put(db, p, 10);
    cfg->put(db, p, 20);
    cfg->put(db, p, 30);
    ASSERT_EQ(cfg->cell_count(db, p), 3);
    ASSERT_EQ(cfg->get(db, p), 10);

    uint32_t cur = cfg->get_multi(db, p);
    uint32_t got[3] = { 0, 0, 0 };
    uint32_t ref;
    int n = 0;
    while (geo_cell_next(&ref, cur) && n < 3)
        got[n++] = ref;
    ASSERT_EQ(n, 3);
    ASSERT_EQ(got[0], 10);
    ASSERT_EQ(got[1], 20);
    ASSERT_EQ(got[2], 30);

    cfg->replace(db, p, 99);
    ASSERT_EQ(cfg->cell_count(db, p), 1);
    ASSERT_EQ(cfg->get(db, p), 99);
    cfg->put(db, p, 7);
    cfg->del(db, p);
    ASSERT_EQ(cfg->cell_count(db, p), 1);
    ASSERT_EQ(cfg->get(db, p), 7);
    ASSERT_EQ(cfg->del_all(db, p), 1);
    ASSERT_EQ(cfg->cell_count(db, p), 0);

    /* --- iter + fused next: 2^dim grid, refs are linear indices */
    int cells = 1 << dim;
    int16_t start[4] = { 0, 0, 0, 0 };
    uint16_t len[4] = { 2, 2, 2, 2 };
    int16_t q[4];
    for (uint32_t refv = 0; refv < (uint32_t)cells; refv++) {
        uint32_t acc = refv;
        for (int i = 0; i < dim; i++) {
            start[i] = (int16_t)(acc & 1);
            acc >>= 1;
        }
        cfg->put(db, start, refv);
    }
    for (int i = 0; i < dim; i++) {
        start[i] = 0;
        len[i] = 2;
    }
    uint32_t it = cfg->iter(db, start, len);
    uint32_t refv;
    int seen = 0;
    uint64_t sum = 0;
    while (cfg->next(q, &refv, it)) {
        seen++;
        sum += refv;
    }
    ASSERT_EQ(seen, cells);
    ASSERT_EQ(sum, (uint64_t)cells * (cells - 1) / 2);

    /* --- fill_bbox into a sealed recall set */
    rec_set_t *out = rec_set_new();
    ASSERT_EQ(cfg->fill_bbox(db, start, len, out), 0);
    ASSERT_EQ((uint64_t)rec_set_count(out), (uint64_t)cells);
    rec_set_free(out);
}

/* 2D x 32-bit config: same member coverage on int32 lanes. */
static void exercise_p4b(const geo_point4b_t *cfg) {
    uint32_t db = fresh_db();

    int32_t p[2] = { 7, -3 };
    uint64_t code = cfg->morton_set(p);
    ASSERT_EQ(code, morton_set_2_32(p));
    int32_t backp[2];
    cfg->morton_get(backp, code);
    ASSERT_EQ(backp[0], p[0]);
    ASSERT_EQ(backp[1], p[1]);

    int32_t a[2] = { 1000, -2000 };
    int32_t b[2] = { 5, 6 };
    int32_t r[2];

    cfg->add(r, a, b);
    ASSERT_EQ(r[0], a[0] + b[0]);
    ASSERT_EQ(r[1], a[1] + b[1]);
    cfg->sub(r, a, b);
    ASSERT_EQ(r[0], a[0] - b[0]);
    ASSERT_EQ(r[1], a[1] - b[1]);
    cfg->min(r, a, b);
    ASSERT_EQ(r[0], 5);
    ASSERT_EQ(r[1], -2000);
    cfg->max(r, a, b);
    ASSERT_EQ(r[0], 1000);
    ASSERT_EQ(r[1], 6);
    cfg->copy(r, a);
    ASSERT_EQ(r[0], a[0]);
    ASSERT_EQ(r[1], a[1]);
    cfg->set(r, 42);
    ASSERT_EQ(r[0], 42);
    ASSERT_EQ(r[1], 42);

    int32_t dims[2] = { 4, 250 };
    ASSERT_EQ(cfg->vol(dims), 1000ULL);

    int32_t s[2] = { 0, 0 };
    int32_t e[2] = { 8, 8 };
    int32_t q[2] = { 3, 5 };
    ASSERT_EQ(cfg->idx(q, s, e), 5ULL * 8 + 3);

    check_debug_p4b(cfg->debug);

    /* db CRUD */
    ASSERT_EQ(cfg->get(db, p), GEO_MISS);
    cfg->put(db, p, 10);
    cfg->put(db, p, 20);
    cfg->put(db, p, 30);
    ASSERT_EQ(cfg->cell_count(db, p), 3);
    ASSERT_EQ(cfg->get(db, p), 10);
    uint32_t cur = cfg->get_multi(db, p);
    uint32_t ref;
    int n = 0;
    while (geo_cell_next(&ref, cur))
        n++;
    ASSERT_EQ(n, 3);
    cfg->replace(db, p, 99);
    ASSERT_EQ(cfg->get(db, p), 99);
    cfg->del(db, p);
    ASSERT_EQ(cfg->get(db, p), GEO_MISS);
    cfg->put(db, p, 5);
    ASSERT_EQ(cfg->del_all(db, p), 1);
    ASSERT_EQ(cfg->cell_count(db, p), 0);

    /* iter + fused next */
    for (int32_t y = 0; y < 4; y++) {
        for (int32_t x = 0; x < 4; x++) {
            int32_t cell[2] = { x, y };
            cfg->put(db, cell, (uint32_t)(y * 4 + x));
        }
    }
    int32_t start[2] = { 0, 0 };
    int32_t len[2] = { 4, 4 };
    uint32_t it = cfg->iter(db, start, len);
    n = 0;
    while (cfg->next(p, &ref, it))
        n++;
    ASSERT_EQ(n, 16);

    rec_set_t *out = rec_set_new();
    ASSERT_EQ(cfg->fill_bbox(db, start, len, out), 0);
    ASSERT_EQ((uint64_t)rec_set_count(out), 16ULL);
    rec_set_free(out);
}

TEST(point1_2_config)   { exercise_p2b(&Point1_2, 1); }
TEST(point2_2_config)   { exercise_p2b(&Point2_2, 2); }
TEST(point3_2_config)   { exercise_p2b(&Point3_2, 3); }
TEST(point4_2_config)   { exercise_p2b(&Point4_2, 4); }
TEST(point2_4_config)   { exercise_p4b(&Point2_4); }

int main(void) {
    setup_once();
    test_suite_begin("Point Config Objects Unit Tests");

    RUN_TEST(members_present);
    RUN_TEST(symbol_aliasing);
    RUN_TEST(point1_2_config);
    RUN_TEST(point2_2_config);
    RUN_TEST(point3_2_config);
    RUN_TEST(point4_2_config);
    RUN_TEST(point2_4_config);

    return test_suite_end();
}