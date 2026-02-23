/**
 * @file test_common.h
 * @brief Common testing utilities and macros for libgeo test suite.
 *
 * Provides a simple, lightweight testing framework with assertions,
 * colored output, and test statistics tracking.
 */

#ifndef TEST_COMMON_H
#define TEST_COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

/* ANSI color codes for test output */
#define COLOR_RESET   "\033[0m"
#define COLOR_RED     "\033[31m"
#define COLOR_GREEN   "\033[32m"
#define COLOR_YELLOW  "\033[33m"
#define COLOR_BLUE    "\033[34m"
#define COLOR_MAGENTA "\033[35m"
#define COLOR_CYAN    "\033[36m"
#define COLOR_BOLD    "\033[1m"

/* Test statistics */
static int test_count = 0;
static int test_passed = 0;
static int test_failed = 0;
static const char *current_test_name = NULL;
static int current_test_assertions = 0;

/* Timing */
static struct timeval test_start_time;

/* Get current time in microseconds */
static inline uint64_t get_time_usec(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (uint64_t)tv.tv_sec * 1000000 + tv.tv_usec;
}

/* Test failure handler */
static inline void test_fail(const char *file, int line, const char *expr) {
    fprintf(stderr, "%s    FAIL:%s %s:%d: Assertion failed: %s\n",
            COLOR_RED, COLOR_RESET, file, line, expr);
    test_failed++;
    exit(1); /* Exit immediately on first failure */
}

/* Test success marker */
static inline void test_pass_assertion(void) {
    current_test_assertions++;
}

/* Assertion macros */
#define ASSERT(expr) \
    do { \
        if (!(expr)) { \
            test_fail(__FILE__, __LINE__, #expr); \
        } \
        test_pass_assertion(); \
    } while(0)

#define ASSERT_EQ(a, b) \
    do { \
        if ((a) != (b)) { \
            fprintf(stderr, "%s    Expected: %lld, Got: %lld%s\n", \
                    COLOR_YELLOW, (long long)(b), (long long)(a), COLOR_RESET); \
            test_fail(__FILE__, __LINE__, #a " == " #b); \
        } \
        test_pass_assertion(); \
    } while(0)

#define ASSERT_NEQ(a, b) \
    do { \
        if ((a) == (b)) { \
            fprintf(stderr, "%s    Expected different values, both: %lld%s\n", \
                    COLOR_YELLOW, (long long)(a), COLOR_RESET); \
            test_fail(__FILE__, __LINE__, #a " != " #b); \
        } \
        test_pass_assertion(); \
    } while(0)

#define ASSERT_LT(a, b) \
    do { \
        if ((a) >= (b)) { \
            fprintf(stderr, "%s    Expected %lld < %lld%s\n", \
                    COLOR_YELLOW, (long long)(a), (long long)(b), COLOR_RESET); \
            test_fail(__FILE__, __LINE__, #a " < " #b); \
        } \
        test_pass_assertion(); \
    } while(0)

#define ASSERT_GT(a, b) \
    do { \
        if ((a) <= (b)) { \
            fprintf(stderr, "%s    Expected %lld > %lld%s\n", \
                    COLOR_YELLOW, (long long)(a), (long long)(b), COLOR_RESET); \
            test_fail(__FILE__, __LINE__, #a " > " #b); \
        } \
        test_pass_assertion(); \
    } while(0)

#define ASSERT_GE(a, b) \
    do { \
        if ((a) < (b)) { \
            fprintf(stderr, "%s    Expected %lld >= %lld%s\n", \
                    COLOR_YELLOW, (long long)(a), (long long)(b), COLOR_RESET); \
            test_fail(__FILE__, __LINE__, #a " >= " #b); \
        } \
        test_pass_assertion(); \
    } while(0)

#define ASSERT_LE(a, b) \
    do { \
        if ((a) > (b)) { \
            fprintf(stderr, "%s    Expected %lld <= %lld%s\n", \
                    COLOR_YELLOW, (long long)(a), (long long)(b), COLOR_RESET); \
            test_fail(__FILE__, __LINE__, #a " <= " #b); \
        } \
        test_pass_assertion(); \
    } while(0)

#define ASSERT_NULL(ptr) ASSERT((ptr) == NULL)
#define ASSERT_NOT_NULL(ptr) ASSERT((ptr) != NULL)

/* Point comparison */
#define ASSERT_POINT_EQ(p1, p2, dim) \
    do { \
        int _match = 1; \
        for (uint8_t _i = 0; _i < (dim); _i++) { \
            if ((p1)[_i] != (p2)[_i]) { \
                _match = 0; \
                fprintf(stderr, "%s    Point mismatch at dimension %d: %d != %d%s\n", \
                        COLOR_YELLOW, _i, (p1)[_i], (p2)[_i], COLOR_RESET); \
                break; \
            } \
        } \
        if (!_match) { \
            test_fail(__FILE__, __LINE__, "points equal"); \
        } \
        test_pass_assertion(); \
    } while(0)

/* String comparison */
#define ASSERT_STREQ(s1, s2) \
    do { \
        if (strcmp((s1), (s2)) != 0) { \
            fprintf(stderr, "%s    Expected: \"%s\", Got: \"%s\"%s\n", \
                    COLOR_YELLOW, (s2), (s1), COLOR_RESET); \
            test_fail(__FILE__, __LINE__, #s1 " == " #s2); \
        } \
        test_pass_assertion(); \
    } while(0)

/* Test definition macro */
#define TEST(name) \
    static void test_##name(void); \
    static void run_test_##name(void) { \
        current_test_name = #name; \
        current_test_assertions = 0; \
        test_count++; \
        test_start_time.tv_sec = 0; \
        gettimeofday(&test_start_time, NULL); \
        printf("%s[TEST]%s %s", COLOR_CYAN, COLOR_RESET, #name); \
        fflush(stdout); \
        test_##name(); \
        struct timeval end_time; \
        gettimeofday(&end_time, NULL); \
        uint64_t elapsed = (end_time.tv_sec - test_start_time.tv_sec) * 1000000 + \
                          (end_time.tv_usec - test_start_time.tv_usec); \
        printf(" %s✓%s (%d assertions, %.2f ms)\n", \
               COLOR_GREEN, COLOR_RESET, current_test_assertions, elapsed / 1000.0); \
        test_passed++; \
    } \
    static void test_##name(void)

/* Test runner macro */
#define RUN_TEST(name) run_test_##name()

/* Test suite header */
static inline void test_suite_begin(const char *suite_name) {
    printf("\n%s%s=== Test Suite: %s ===%s\n\n", 
           COLOR_BOLD, COLOR_BLUE, suite_name, COLOR_RESET);
    test_count = 0;
    test_passed = 0;
    test_failed = 0;
}

/* Test suite summary */
static inline int test_suite_end(void) {
    printf("\n%s%s=== Test Summary ===%s\n", COLOR_BOLD, COLOR_BLUE, COLOR_RESET);
    printf("Total:  %d tests\n", test_count);
    printf("%sPassed: %d%s\n", COLOR_GREEN, test_passed, COLOR_RESET);
    
    if (test_failed > 0) {
        printf("%sFailed: %d%s\n", COLOR_RED, test_failed, COLOR_RESET);
        return 1;
    } else {
        printf("%sAll tests passed!%s\n", COLOR_GREEN, COLOR_RESET);
        return 0;
    }
}

/* Benchmark utilities */
typedef struct {
    const char *name;
    uint64_t start_time;
    uint64_t ops_count;
} benchmark_t;

static inline void bench_start(benchmark_t *bench, const char *name) {
    bench->name = name;
    bench->start_time = get_time_usec();
    bench->ops_count = 0;
}

static inline void bench_end(benchmark_t *bench, uint64_t ops) {
    uint64_t end_time = get_time_usec();
    uint64_t elapsed = end_time - bench->start_time;
    double seconds = elapsed / 1000000.0;
    double ops_per_sec = ops / seconds;
    
    printf("%s[BENCH]%s %s: ", COLOR_MAGENTA, COLOR_RESET, bench->name);
    printf("%.2f M ops/sec (%.2f ms for %llu ops)\n",
           ops_per_sec / 1000000.0, elapsed / 1000.0, (unsigned long long)ops);
}

/* Simple random number generator for tests */
static uint64_t test_rng_state = 0x123456789ABCDEF0ULL;

static inline void test_seed_rng(uint64_t seed) {
    test_rng_state = seed;
}

static inline uint64_t test_rand64(void) {
    /* xorshift64* */
    test_rng_state ^= test_rng_state >> 12;
    test_rng_state ^= test_rng_state << 25;
    test_rng_state ^= test_rng_state >> 27;
    return test_rng_state * 0x2545F4914F6CDD1DULL;
}

static inline int16_t test_rand_coord(void) {
    return (int16_t)(test_rand64() & 0xFFFF);
}

static inline int16_t test_rand_coord_range(int16_t min, int16_t max) {
    int32_t range = (int32_t)max - min;
    return min + (int16_t)(test_rand64() % range);
}

#endif /* TEST_COMMON_H */
