# Testing libgeo

This document describes how to build and run the comprehensive test suite for libgeo.

## Overview

The test suite includes:
- **Unit Tests**: Test individual functions and components
- **Property-Based Tests**: Verify mathematical invariants and properties
- **Integration Tests**: End-to-end scenarios
- **Stress Tests**: Capacity and performance limits
- **Benchmarks**: Performance measurements for critical operations
- **Fuzz Tests**: Robustness testing with libFuzzer

## Quick Start

```bash
# Run all tests (unit + integration + property + stress)
cd tests
make clean test

# Run performance benchmarks
make bench

# Build fuzz test harnesses (requires libFuzzer for full functionality)
make fuzz
make fuzz-standalone  # Run in standalone mode
```

## Test Categories

### Unit Tests (`tests/unit/`)

Tests for individual components:
- **test_morton.c**: Morton code encoding/decoding (16 tests)
- **test_point.c**: Point arithmetic and utilities (32 tests)
- **test_geo_core.c**: Core geo API functions (17 tests)

Run with: `make test-unit`

### Integration Tests (`tests/integration/`)

End-to-end scenarios:
- **test_iteration.c**: Iterator correctness (5 tests)
- **test_persistence.c**: Database persistence (5 tests)
- **test_spatial_queries.c**: Spatial query behavior (3 tests)

Run with: `make test-integration`

### Property-Based Tests (`tests/property/`)

Mathematical property verification with thousands of random inputs:
- **test_morton_props.c**: Morton code invariants (reversibility, uniqueness, determinism)

Run with: `make test-property`

### Stress Tests (`tests/stress/`)

Capacity and performance limits:
- **test_capacity.c**: Hash table capacity and collision handling (5 tests)
- **test_large_datasets.c**: Large dataset operations (6 tests)

Run with: `make test-stress`

### Benchmarks (`tests/benchmark/`)

Performance measurements:
- **bench_morton.c**: Morton encode/decode throughput
- **bench_insertion.c**: Database creation performance
- **bench_queries.c**: Point lookup performance
- **bench_iteration.c**: Point insert performance

Run with: `make bench`

Example output:
```
Morton Encode (3D): 17.49 M ops/sec
Morton Decode (3D): 20.33 M ops/sec
Morton Round-Trip (3D): 10.13 M ops/sec
```

### Fuzz Tests (`tests/fuzz/`)

Robustness testing with randomized inputs:
- **fuzz_morton.c**: Morton encoding edge cases
- **fuzz_geo_api.c**: Random geo API operations
- **fuzz_queries.c**: Spatial query edge cases

Supports both standalone mode and libFuzzer integration.

Run with:
```bash
make fuzz-standalone  # Standalone mode (works now)
# For libFuzzer mode (requires libFuzzer library):
# CFLAGS="-DLIBFUZZER -fsanitize=fuzzer" make fuzz
```

## Advanced Testing

### Memory Leak Detection (Valgrind)

```bash
make valgrind
```

Runs all tests under Valgrind to detect memory leaks and invalid memory access.

### Address Sanitizer

```bash
make asan
make test
```

Builds tests with AddressSanitizer for runtime memory error detection.

### Undefined Behavior Sanitizer

```bash
make ubsan
make test
```

Detects undefined behavior at runtime.

### Thread Sanitizer

```bash
make tsan
make test
```

Detects data races (note: libgeo is NOT thread-safe by design).

### Code Coverage

```bash
make coverage
```

Generates HTML coverage report in `coverage_html/`.

## Test Structure

```
tests/
├── Makefile              # Build system for all tests
├── test_common.h         # Testing framework and utilities
├── unit/                 # Unit tests (65 tests)
├── integration/          # Integration tests (13 tests)
├── stress/               # Stress tests (11 tests)
├── benchmark/            # Performance benchmarks (4)
├── fuzz/                 # Fuzz testing harnesses (3)
└── property/             # Property-based tests (12K+ assertions)
```

## Writing Tests

The test framework provides simple macros:

```c
#include "../test_common.h"
#include "../../include/ttypt/geo.h"

TEST(my_test_name) {
    int16_t pos[3] = {10, 20, 30};
    uint32_t value = 42;
    
    uint32_t db = geo_open(NULL, "test_db", 1023);
    geo_put(db, pos, value, 3);
    uint32_t result = geo_get(db, pos, 3);
    
    ASSERT_EQ(result, value);
}

int main(void) {
    RUN_TEST(my_test_name);
    return test_suite_end();
}
```

Available assertions:
- `ASSERT(expr)` - Assert expression is true
- `ASSERT_EQ(a, b)` - Assert equality
- `ASSERT_NEQ(a, b)` - Assert inequality
- `ASSERT_LT(a, b)` - Assert less than
- `ASSERT_GT(a, b)` - Assert greater than
- `ASSERT_GE(a, b)` - Assert greater or equal
- `ASSERT_LE(a, b)` - Assert less or equal
- `ASSERT_POINT_EQ(p1, p2, dim)` - Assert points are equal
- `ASSERT_NULL(ptr)` / `ASSERT_NOT_NULL(ptr)` - Pointer checks

## Continuous Integration

Tests run automatically on every push and pull request via GitHub Actions.
See `.github/workflows/test.yml` for CI configuration.

## Current Test Coverage

- **Morton encoding/decoding**: 100% (16 unit tests + 10K+ property tests)
- **Point utilities**: 100% (32 tests covering all functions)
- **Geo core API**: ~95% (17 tests)
- **Iteration/queries**: 13 tests
- **Stress testing**: 11 tests
- **Fuzzing**: 3 harnesses

**Total**: 89 tests + 12,000+ property-based assertions + 4 benchmarks + 3 fuzz harnesses

## Known Limitations

- **Thread safety**: Not tested as libgeo is explicitly single-threaded

## Troubleshooting

### Tests fail to build

Ensure libgeo is built first:
```bash
cd .. && make && cd tests
```

### Valgrind reports errors

Check if errors are in libgeo code or dependencies (qmap, qsys). Libgeo-specific leaks should be investigated.

### Performance regression

Benchmark results may vary by system. Compare relative performance, not absolute numbers.
