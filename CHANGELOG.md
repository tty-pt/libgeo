## [Unreleased]
- BMI2 PDEP/PEXT codec fast path (`GEO_USE_PDEP`, default 1): activates
  only when the TU is compiled with `-mbmi2` (`__BMI2__`); the scalar
  spread/compact bit-twiddle kernels stay as the portable fallback and
  as the reference for bit-identity (enforced by
  `tests/unit/test_codec_parity.c` and a bench-time cross-check).
  Measured 2-4x on isolated encode/decode/round-trip once the paired
  bench pre-generates inputs outside the timed loop (RNG cost otherwise
  swamps the codec signal — see `docs/PERF.md`). Opt out with
  `-DGEO_USE_PDEP=0` (e.g. on AMD Zen ≤3, where PDEP/PEXT are
  microcoded and slow).
- Bulk decode API: `morton_get_bulk()` / `morton_get_bulk4()`, the
  decode-side counterpart to `morton_set_bulk`/`morton_set_bulk4` (same
  `GEO_SIMD_MORTON` gate, same AVX2/NEON-stub/scalar three-tier
  structure). Closes the encode-only asymmetry in the bulk API.
  Measured win is thin/noisy under `-mavx2` (same scatter-cost story as
  the existing bulk encode); kept for API symmetry and because it's
  never worse than the scalar tail. Parity-tested against
  `morton_get_3`/`morton_get_4` in `test_codec_parity.c`.
- New primary API surface: point config objects (`include/ttypt/
  pointcfg.h`). Each config is one exported const struct — "a class with
  all-static methods" — named `Point<D>_<B>` (D = dims, B = bytes/lane):
  `Point1_2`, `Point2_2`, `Point3_2`, `Point4_2` (int16 lanes, shared
  `geo_point2b_t` type) and `Point2_4` (int32 lanes, `geo_point4b_t`).
  Each has ~20 members (codec, `add/sub/min/max/copy/vol/set/idx/debug`,
  `put/get/replace/del/del_all/cell_count`, `get_multi` + `geo_cell_next`,
  `iter` + fused `.next`, `fill_bbox`), so a language server
  autocompletes the whole config off one symbol. Struct members are
  addressable functions with the same bodies as the flat inlines; the
  flat API stays exported and documented as the tight-loop fast path /
  runtime-dim (`geo_ops[]`) escape hatch. Nothing removed.
- New 2D x 32-bit dense config (`geo_*_2_32`, `morton_set_2_32` /
  `morton_get_2_32`, `point_*_2_32`): int32_t lanes
  (-2147483648..2147483647), dense stride-2 codec filling all 64 key
  bits (2 x 32, no reserved bits). Box walker, Z-interval skip, fill,
  iterate, and multi-value chains all work on 32-bit lanes with
  int32_t start AND lengths; one cursor pool/idm serves both lane
  widths via per-cursor copy ops. `geo_ops[]` stays int16-lane-only;
  the new config is reached through its suffixed monomorphs, advanced
  with `geo_next32()`, multi-values with the shared `geo_cell_next()`.
  One config per database, never reopen with another config's API
  (the file carries no config tag). GEO_FILL_MAX_VOL (1M cells) caps
  every config alike. Walker internals (`inrange`, gap jump, box
  walk, iter/get_multi/fill stamps) are now stamped from
  config-parametric macros; the int16 1..4 instantiations emit the
  same expressions as before (full suite + oracle/scan-count tests
  confirm identical walks).
- BREAKING: per-dimension public API. The runtime-`dim` functions are
  gone — `morton_set(p, dim)` is now `morton_set_1..4(p)`,
  `morton_get` is now `morton_get_1..4`, every `point_*` is now
  `point_*_1..4`, and every database op is now `geo_put_N`,
  `geo_set_N`, `geo_get_N`, `geo_del_N`, `geo_del_all_N`,
  `geo_cell_count_N`, `geo_get_multi_N`, `geo_iter_N`,
  `rec_axis_fill_bbox_N` (N = 1..4, no dim argument). Invalid dims are
  unrepresentable: dim 0/5 rejections are deleted (dim-0 table slot is
  empty instead). Morton CODE VALUES unchanged (1D/2D/3D still
  bit-identical to v0.5.0).
- New `geo_ops[1..4]` runtime-dim table (`geo.h`): per-dim function
  pointers with no dim argument (the index is the dim), covering
  codec, point, scatter, iterate, and fill ops. `geo_ops[0]` is NULL.
  `geo_ops[dim].morton_set(p)` is the migration path for callers with
  a genuinely runtime dim.
- Monomorphized box walker: `geo_box_walk_1..4` stamped from one macro
  with literal dims (no per-iteration dim dispatch anywhere);
  `geo_jump_over_gap` forced inline so the literal folds its k-loop
  chain; cursor/collection carry a per-dim `point_copy` pointer.
  Measured parity vs the generic walker on this VM (see docs/PERF.md);
  kept for the structural win at ~19KB extra .text.
- 4D support (dim=4): dense stride-4 Morton codec using the full 64-bit
  key space; box walker, Z-interval skip, fill, and SIMD bulk all
  generalized (skip-cube span is now 2^(D*k)); 1D/2D/3D codes unchanged
- New tunables (default 0, measure-first): GEO_SMALLDIM_UNROLL (1D/2D
  fast paths), GEO_4D_UNROLL (4D fast paths incl. exact 8-byte
  point_copy); new morton_set_bulk4() batch encoder
- Retired GEO_PACKED_CURI (measured neutral-negative; 4D needs 4 slots
  anyway, struct stays 12 bytes either way)
- Retired GEO_3D_POINT_COPY (measured consistently below 1.0x in
  interleaved testing; plain loop restored)
- Promoted all remaining tunables to unconditional: inline Morton codec
  (duplicate non-inline implementation deleted, ABI wrappers kept),
  hoisted box bounds (single gap-jump implementation), clz kmax,
  dimension-unrolled inrange/gap-jump/point_copy (generic fallbacks
  deleted). Only GEO_SIMD_MORTON (batch API gate) remains tunable.
- New docs/PERF.md: per-flag verdicts, numbers, and the VM-drift
  benchmarking caveat; single-dimension-per-database convention
  documented (all dims share the uint64 keyspace)

## [0.5.0] - 2026-09-10
- Kernel form (requires libqmap >= 0.8.0): maps open QM_SORTED|QM_MULTIVALUE
  - Multi-value cells: geo_put appends; new geo_set (replace), geo_get_multi /
    geo_cell_next (chain read in insertion order), geo_cell_count, geo_del_all
  - geo_get returns the first value; geo_del removes the first value
  - New rec_axis_fill_bbox(): stream a bounding box into a sealed rec_set_t
    (sorted + deduped), never materializing the box volume; GEO_FILL_MAX_VOL
    (1M cells) cap with -1 errors; plain int return
- Raw iterator restructured: flat (point, value) list in morton-discovery
  order (matches the long-standing doc claim), grown on demand — no more
  volume-sized array, sparse queries over huge boxes stay cheap
- Z-interval skip: box walks ratchet a skip floor past the largest
  box-disjoint aligned cube at each false positive (provably sound; replaces
  the dead COMPUTE_BMLM path, whose constants were wrong); geo_last_scan_count
  diagnostic proves engagement
- Docs: capacity claims fixed (mask = initial table, map auto-grows, no
  termination); geo.pc version synced
- Tests: 10 MV unit + 11 fill unit + 7 range unit (brute-force oracle,
  mutation-validated) + parity integration + MV persistence + randomized
  fill/range property oracles + fill/sparse benches; asan/ubsan clean

## [0.4.1] - 2026-02-23
- Fix geo_iter segfaults: add bounds checking in geo_search()
- Fix 2D/1D Morton encoding: use dimension-specific pack functions
- Fix GEO_MISS: now returns UINT32_MAX (was incorrectly defined as 64-bit)

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
