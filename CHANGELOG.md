## [0.4.0] - 2026-02-23
- Add comprehensive test suite
  - Unit tests: 65 tests for core functions (morton, point, geo)
  - Integration tests: 13 tests for iteration, persistence, queries
  - Property-based tests: 12K+ assertions for morton invariants
  - Stress tests: 11 tests for capacity and large datasets
  - Benchmarks: 4 benchmarks for performance measurement
  - Fuzz tests: 3 harnesses with libFuzzer support
- Add 'make test' target to root Makefile
- Update README with testing documentation

## [0.3.0] - 2026-02-23
- Update to libqmap 0.6.0
  - Improved pointer stability via allocation reuse optimization
  - File loading no longer requires QM_MIRROR flag
  - Enhanced documentation and bug fixes

## [0.2.0] - 2025-10-24
- Update to libqmap 0.5.0 (BTREE support)
