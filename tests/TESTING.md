# Testing libgeo

This document describes how to build and run the comprehensive test suite for libgeo.

## Overview

The test suite includes:
- **Unit Tests**: Test individual functions and components
- **Property-Based Tests**: Verify mathematical invariants and properties
- **Benchmarks**: Performance measurements for critical operations
- **Integration Tests** (TODO): End-to-end scenarios
- **Stress Tests** (TODO): Capacity and performance limits
- **Fuzz Tests** (TODO): Robustness testing

## Quick Start

```bash
# Run all core tests (unit + property)
cd tests
make test

# Run only unit tests
make test-unit

# Run property-based tests
make test-property

# Run performance benchmarks
make bench
```

## Test Categories

### Unit Tests (`tests/unit/`)

Tests for individual components:
- **test_morton.c**: Morton code encoding/decoding (16 tests)
- **test_point.c**: Point arithmetic and utilities (32 tests)
- **test_geo_core.c**: Core geo API functions (17 tests)

Run with: `make test-unit`

### Property-Based Tests (`tests/property/`)

Mathematical property verification with thousands of random inputs:
- **test_morton_props.c**: Morton code invariants (reversibility, uniqueness, determinism)

Run with: `make test-property`

### Benchmarks (`tests/benchmark/`)

Performance measurements:
- **bench_morton.c**: Morton encode/decode throughput

Run with: `make bench`

Example output:
```
Morton Encode (3D): 17.49 M ops/sec
Morton Decode (3D): 20.33 M ops/sec
Morton Round-Trip (3D): 10.13 M ops/sec
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
├── unit/                 # Unit tests
├── integration/          # Integration tests (TODO)
├── stress/               # Stress tests (TODO)
├── benchmark/            # Performance benchmarks
├── fuzz/                 # Fuzz testing harnesses (TODO)
└── property/             # Property-based tests
```

## Writing Tests

The test framework provides simple macros:

```c
#include "../test_common.h"
#include "../../include/ttypt/geo.h"

TEST(my_test_name) {
    int16_t pos[3] = {10, 20, 30};
    uint32_t value = 42;
    
    geo_put(db, pos, value, 3);
    uint32_t result = geo_get(db, pos, 3);
    
    ASSERT_EQ(result, value);
}

int main(void) {
    test_suite_begin("My Test Suite");
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
- **Geo core API**: ~95% (17 tests, missing some edge cases)

**Total**: 65+ unit tests, 12,000+ property-based assertions

## Known Limitations

- **2D/1D support**: Tests currently focus on 3D operations as the library is optimized for 3D with `FAST_MORTON=1`
- **Thread safety**: Not tested as libgeo is explicitly single-threaded
- **Fuzzing**: Harnesses not yet implemented
- **Integration tests**: Not yet implemented
- **Stress tests**: Not yet implemented

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

## Future Work

- [ ] Integration tests for complex spatial queries
- [ ] Stress tests for capacity limits
- [ ] Fuzz testing with AFL/libFuzzer
- [ ] Memory tests for leak detection
- [ ] CI performance tracking
- [ ] Automated regression detection
