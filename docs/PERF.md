# libgeo performance notes

Measured findings for the `GEO_*` optimization tunables (see `geo.h`),
the per-dimension API, plus the benchmarking methodology that produced
them. Last updated for the per-dimension API change (unreleased; see
CHANGELOG).

## Tunables

One flag remains, a `0/1` macro defaulted in `include/ttypt/geo.h`.
Override with `-DGEO_SIMD_MORTON=0/1` on the compiler command line.

> **Footgun:** the defaults block uses `#ifndef`, so `-U GEO_SIMD_MORTON`
> does **not** disable it — the block re-defines it to the default.
> Use `-DGEO_SIMD_MORTON=0` to force off.

| Flag | Default | Verdict |
|---|---|---|
| `GEO_SIMD_MORTON` | 1 | Proven win: `morton_set_bulk` / `morton_set_bulk4` |

All former per-path tunables were **promoted to unconditional** after
measurement showed the generic fallbacks never win (the unrolled /
hoisted / inlined forms won or tied in every paired and interleaved
comparison; identical `examined/walk` counts throughout). The ~140-line
duplicate non-inline codec, the non-hoisted gap-jump copy, the `kmax`
linear scan, and every unrolled-vs-loop guard are deleted. Retired
along the way:

| Flag | Fate |
|---|---|
| ~~`GEO_PACKED_CURI`~~ | **Retired.** Measured neutral-negative, and 4D support
requires 4 coordinate slots anyway (`p[4]` + `ref` = 12 bytes either
way — the "packed" form saved zero bytes). `geo_curi_t` is now
unconditionally `int16_t p[MAX_DIM]`. |
| ~~`GEO_3D_POINT_COPY`~~ | **Retired.** Consistently below 1.0 in both interleaved
sessions (BOX 0.93, LOOKUP 0.86, PUT 0.93 medians of per-round
ratios): the manual 4+2 copy pessimizes versus the compiler's loop
handling. Plain loop restored. |
| ~~`GEO_INLINE_MORTON`~~ | **Promoted.** Proven win; inline codec is now the only
codec (extern ABI wrappers retained in the `.so`). |
| ~~`GEO_HOIST_BOX`~~ / ~~`GEO_CLZ_KMAX`~~ / ~~`GEO_3D_UNROLL`~~ / ~~`GEO_SMALLDIM_UNROLL`~~ / ~~`GEO_4D_UNROLL`~~ | **Promoted.** Hoisted bounds, `clz` kmax, and all
dimension-specialized fast paths are now unconditional; the generic
fallbacks they beat-or-tied are deleted. |

## Methodology warning: this VM drifts

Sequential A-then-B benchmarking on the shared hypervisor VM used here
is **meaningless**: a clean baseline re-run measured ~3x faster than the
first baseline (point lookup 2.0 → 4.2M ops/sec) with no code change.
No frequency control is available (`/sys/.../cpufreq` absent).

Trusted numbers below come only from drift-proof designs:

- **Same-process paired legs** run back-to-back (scalar loop vs
  `morton_set_bulk`; direct vs forced-extern morton calls). Drift
  cancels; order bias noted where present.
- **Interleaved A/B** with per-config `.so` builds in separate dirs,
  run in shuffled order pinned with `taskset -c 0`, scored as
  per-round ratios against the baseline in the *same* round.

## Proven wins

- **SIMD bulk encode** (`morton_set_bulk`, AVX2 4-wide): beats the
  inlined scalar loop **7/7 paired rounds, 1.3–3.7x (median ~2.3x)**
  on 1M random 3D points. Correctness verified against scalar on
  1027 random points + edge cases (all `SHRT_MIN`/`SHRT_MAX`/0/±1).
  `morton_set_bulk4` is the 4D analogue (same structure, `spread4`
  sequence); scalar fallback on non-AVX2.
- **Inline morton** (`morton_set_N`/`morton_get_N` as `static inline`):
  paired direct-vs-forced-extern test: inline leg faster in **11/14
  legs despite always running first on a cold cache**, median
  ~1.5–1.8x on tight encode/decode loops. Removes call overhead on
  every `geo_put_N`/`geo_get_N`/`geo_del_N`/box-walk decode. Zero
  correctness risk (full suite passes both ways); the `.so` still
  exports thin ABI wrappers (`morton_set_1..4`, `morton_get_1..4`) so
  linking is unaffected.

## Per-dimension API + monomorphized walker (unreleased)

The runtime-`dim` API is gone: `morton_set_N`/`morton_get_N`,
`point_*_N`, `geo_*_N` (dim in the name), plus the `geo_ops[1..4]`
table for genuinely runtime dims. Internally the box walk is stamped
out 4× (`geo_box_walk_1..4`, one macro source) with the dim as a
compile-time literal, and `geo_jump_over_gap` is `always_inline` so
the literal reaches its body (maxd loop unrolls, the `dim==3/4/else`
chain in the k-loop folds to one path). No per-iteration dim dispatch
remains anywhere in the walker.

Interleaved baseline-vs-perdim A/B (separate processes, `taskset -c
0`, both libs at `-O2`, per-round ratios, 11 rounds ×2 sessions with
larger rep counts in session 2):

| bench | median B/A (s1) | median B/A (s2) | verdict |
|---|---|---|---|
| CODEC3 (encode+decode loop) | 0.88 | 0.98 | parity |
| PUT3 (scatter store) | 0.90 | 1.00 | parity |
| GET3 (scatter lookup) | 0.89 | 1.06 | parity |
| FILL3 (dense 3D box fill) | 0.97 | 1.03 | parity |
| FILL4 (dense 4D box fill) | 1.07 | 1.00 | parity |

Session 1's ~10% scatter lean evaporated with longer runs — it was
noise. Honest verdict: **monomorphization measures parity on this VM**
(all medians within ±6%, round spreads ±25–50%). Expected: the removed
dim branches were perfectly predicted (same direction every
iteration), so deleting them saves only a few µops per key against
`qmap_next`/decode costs. Kept anyway: zero regression, no
per-iteration dispatch left (the structural goal), one macro source,
and ~19KB extra `.text` (`.so` text 8.6KB → 27.4KB). Re-measure on
quiet bare metal before claiming more.

## 2D x 32-bit dense config (unreleased)

`bench_2d32` (`-O3`, same shared VM — absolutes, not A/B claims):

| bench | result |
|---|---|
| Morton Encode+Decode (2D32, 1M random wide lanes) | ~250 M ops/sec |
| Scatter Put (2D32, 100K) | ~4.0 M ops/sec |
| Scatter Get (2D32, 100K) | ~25 M ops/sec |
| Box Fill 256x256 (2D32, 65K cells) | ~6.8 ms/fill |

Z-interval skip engages on 32-bit lanes: 64x64 dense grid, 16x16
sub-box query finds all 289 cells (inclusive end) while decoding 397
of the 768-code interval — the gap jump skips the Z-spill without
dropping an in-box key (oracle-tested in `test_geo_2d32`).

## Neutral-measured paths (now unconditional)

Box-walker micro-opts (hoisted bounds, `clz` kmax, dimension-unrolled
`inrange`/gap-jump/`point_copy`) all measured ≈1.0x, inside the
±10–15% noise floor across two independent 7-round interleaved
sessions (medians of per-round ratios: 0.85–1.12, rank order flipping
between sessions). Since the generic fallbacks never won and the fast
paths reduce instruction counts with identical
`geo_last_scan_count()` behavior, the fast paths are now
unconditional and the fallbacks deleted — less code, same-or-better
speed. Re-measure on quiet bare metal with `perf stat` if you want
absolute numbers.

## Bug lessons

- An early 8-byte `point_copy` fast path for 3D read/wrote 2 bytes
  past bare `int16_t[3]` stack arrays: `test_point` caught it via
  stack-smashing abort. The per-dim `point_copy_3` is three exact
  scalar stores; `point_copy_4` is a clean 8-byte copy (4×int16 =
  exactly 8B). Never use a wider variant than the arrays hold.
- A dlopen-based A/B harness (two `libgeo.so` in one process) is
  invalid: libgeo's global state (documented not-thread-safe,
  not-multi-instance-safe) segfaults/hangs. Use separate processes.
- `tests/Makefile` does **not** rebuild the top-level lib: after any
  header or `src/libgeo.c` change, run top-level `make` first, or the
  tests will link/run against a stale `lib/libgeo.so` (symptom here:
  segfault through `geo_ops` entries whose layout shifted — the stale
  7-member table under new 13-member headers).

## Codec stability

1D/2D/3D Morton codes are **bit-identical** to v0.5.0 (sparse
stride-3 packing, 48-bit key space, upper 16 bits reserved). Proven by
dumping 5015 codes (15 edge + 5000 random points, dims 1–3) from a
pristine `HEAD` worktree build vs the current build: **zero diff**.
4D codes are new (dense stride-4 packing, full 64 bits). All dims
share the `uint64` keyspace: **one dimension per database** (already
the de facto convention; now documented in `geo.h`).

## 4D measurements

- **bulk4 correctness**: AVX2 `spread4_4x` path verified bit-exact vs
  scalar on 1M random points (paired harness) plus 1031 points in
  `test_dim4_props`.
- **bulk4 speed**: roughly at parity with an `-O2` scalar loop
  (paired speedups 0.57–1.64, median ~1.1). The scalar `spread4` is
  shorter than `spread3`, so GCC's auto-vectorizer already does well;
  the manual gather (`_mm256_set_epi64x` × 16) eats the SIMD margin.
  The AVX2 path still matters for `-O0` builds. Requires `-mavx2`
  at lib compile time — the default build has no AVX2 path.
- **4D box walker** (interleaved baseline-vs-all-on, 5 rounds,
  ~45k-point 4D cloud): `LOOKUP4` median **1.25x, all 5 rounds > 1**
  (inlining pays more on the larger 4-axis decode body); `BOX4_QUERY`
  median 0.88x, mixed — neutral, same as the 3D walker verdict.
  Identical checksums and `examined/walk` counts across configs.

## Test-harness trap: which lib are you measuring?

`tests/Makefile` links `-L../lib` but sets **no rpath**, so test
binaries load `/usr/lib/libgeo.so` at runtime — a stale system copy
will silently shadow your fresh build (symptom: new symbols like
`morton_set_bulk4` fail to resolve, or new behavior is absent).
Always run with:

```sh
export LD_LIBRARY_PATH=/home/quirinpa/libgeo/lib  # repo lib first
```

## Build-flag note

The default library build is `-g` with **no `-O` flag** (see root
`Makefile`/mk), i.e. effectively `-O0`, while `tests/benchmark`
binaries build with `-O3`. Absolute benchmark numbers therefore
understate an optimized build, and lib-vs-bench comparisons (e.g.
`bench_bulk4`'s scalar leg vs the lib's scalar fallback) mix
optimization levels. Paired same-process comparisons remain valid
directionally; for absolute numbers rebuild the lib with `-O2`
(and `-mavx2` for the AVX2 bulk paths).

## Reproducing

```sh
# default build + full suite (remember LD_LIBRARY_PATH, see above)
export LD_LIBRARY_PATH=/home/quirinpa/libgeo/lib
make && make -C tests clean && make test
# SIMD bulk API off (the only remaining tunable)
make clean && make CFLAGS="-g -DGEO_SIMD_MORTON=0"
# benchmarks (drift-prone on shared VMs; prefer paired/interleaved harnesses)
make bench
# sanitizers
make -C tests asan && make -C tests test-unit
make -C tests ubsan && make -C tests test-unit
```
