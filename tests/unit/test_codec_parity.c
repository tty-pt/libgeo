/**
 * @file test_codec_parity.c
 * @brief Bit-identity tests: PDEP/PEXT (when built with -mbmi2) vs
 * scalar reference, plus morton_get_bulk/bulk4 parity.
 *
 * On a non-BMI2 build this still runs the bulk-decode and round-trip
 * correctness checks — only the PDEP-vs-scalar cross-check is skipped.
 */

#include "../test_common.h"
#include "../../include/ttypt/morton.h"
#include "../../include/ttypt/geo.h"

/* ---- scalar reference implementations (always spread3/compact path) ---- */

static uint64_t ref_set_1(int16_t *p) {
	uint16_t u = geo_unsign(p[0]);
	return geo_spread3(u);
}
static uint64_t ref_set_2(int16_t *p) {
	return geo_spread3(geo_unsign(p[0]))
		| (geo_spread3(geo_unsign(p[1])) << 1);
}
static uint64_t ref_set_3(int16_t *p) {
	return geo_spread3(geo_unsign(p[0]))
		| (geo_spread3(geo_unsign(p[1])) << 1)
		| (geo_spread3(geo_unsign(p[2])) << 2);
}
static uint64_t ref_set_4(int16_t *p) {
	return geo_spread4(geo_unsign(p[0]))
		| (geo_spread4(geo_unsign(p[1])) << 1)
		| (geo_spread4(geo_unsign(p[2])) << 2)
		| (geo_spread4(geo_unsign(p[3])) << 3);
}
static uint64_t ref_set_2_32(int32_t *p) {
	return geo_spread2(geo_unsign32(p[0]))
		| (geo_spread2(geo_unsign32(p[1])) << 1);
}

static void ref_get_1(int16_t *pos, uint64_t code) {
	pos[0] = geo_sign((uint16_t)geo_compact_axis(code, 0));
}
static void ref_get_2(int16_t *pos, uint64_t code) {
	pos[0] = geo_sign((uint16_t)geo_compact_axis(code, 0));
	pos[1] = geo_sign((uint16_t)geo_compact_axis(code, 1));
}
static void ref_get_3(int16_t *pos, uint64_t code) {
	pos[0] = geo_sign((uint16_t)geo_compact_axis(code, 0));
	pos[1] = geo_sign((uint16_t)geo_compact_axis(code, 1));
	pos[2] = geo_sign((uint16_t)geo_compact_axis(code, 2));
}
static void ref_get_4(int16_t *pos, uint64_t code) {
	pos[0] = geo_sign((uint16_t)geo_compact_axis4(code, 0));
	pos[1] = geo_sign((uint16_t)geo_compact_axis4(code, 1));
	pos[2] = geo_sign((uint16_t)geo_compact_axis4(code, 2));
	pos[3] = geo_sign((uint16_t)geo_compact_axis4(code, 3));
}
static void ref_get_2_32(int32_t *pos, uint64_t code) {
	pos[0] = geo_sign32((uint32_t)geo_compact_axis2(code, 0));
	pos[1] = geo_sign32((uint32_t)geo_compact_axis2(code, 1));
}

#define NPAR 4096

static int16_t pa_i16[NPAR][4];
static int32_t pa_i32[NPAR][2];
static uint64_t pa_codes[NPAR];
static int16_t pa_decoded[4], pa_ref_decoded[4];

TEST(pdep_encode_parity) {
	test_seed_rng(0x0DE);
	for (int i = 0; i < NPAR; i++) {
		for (int d = 0; d < 4; d++)
			pa_i16[i][d] = test_rand_coord();
		pa_i32[i][0] = (int32_t)test_rand64();
		pa_i32[i][1] = (int32_t)test_rand64();
	}

	/* 1D encode */
	for (int i = 0; i < NPAR; i++) {
		uint64_t fast = morton_set_1(pa_i16[i]);
		uint64_t slow = ref_set_1(pa_i16[i]);
		if (fast != slow) { ASSERT_EQ(fast, slow); }
	}
	/* 2D encode */
	for (int i = 0; i < NPAR; i++) {
		uint64_t fast = morton_set_2(pa_i16[i]);
		uint64_t slow = ref_set_2(pa_i16[i]);
		if (fast != slow) { ASSERT_EQ(fast, slow); }
	}
	/* 3D encode */
	for (int i = 0; i < NPAR; i++) {
		uint64_t fast = morton_set_3(pa_i16[i]);
		uint64_t slow = ref_set_3(pa_i16[i]);
		if (fast != slow) { ASSERT_EQ(fast, slow); }
	}
	/* 4D encode */
	for (int i = 0; i < NPAR; i++) {
		uint64_t fast = morton_set_4(pa_i16[i]);
		uint64_t slow = ref_set_4(pa_i16[i]);
		if (fast != slow) { ASSERT_EQ(fast, slow); }
	}
	/* 2_32 encode */
	for (int i = 0; i < NPAR; i++) {
		uint64_t fast = morton_set_2_32(pa_i32[i]);
		uint64_t slow = ref_set_2_32(pa_i32[i]);
		if (fast != slow) { ASSERT_EQ(fast, slow); }
	}
}

TEST(pdep_decode_parity) {
	test_seed_rng(0xDE);
	for (int i = 0; i < NPAR; i++) {
		pa_codes[i] = test_rand64();
	}

	/* 1D decode */
	for (int i = 0; i < NPAR; i++) {
		morton_get_1(pa_decoded, pa_codes[i]);
		ref_get_1(pa_ref_decoded, pa_codes[i]);
		ASSERT_POINT_EQ(pa_decoded, pa_ref_decoded, 1);
	}
	/* 2D decode */
	for (int i = 0; i < NPAR; i++) {
		morton_get_2(pa_decoded, pa_codes[i]);
		ref_get_2(pa_ref_decoded, pa_codes[i]);
		ASSERT_POINT_EQ(pa_decoded, pa_ref_decoded, 2);
	}
	/* 3D decode */
	for (int i = 0; i < NPAR; i++) {
		morton_get_3(pa_decoded, pa_codes[i]);
		ref_get_3(pa_ref_decoded, pa_codes[i]);
		ASSERT_POINT_EQ(pa_decoded, pa_ref_decoded, 3);
	}
	/* 4D decode */
	for (int i = 0; i < NPAR; i++) {
		morton_get_4(pa_decoded, pa_codes[i]);
		ref_get_4(pa_ref_decoded, pa_codes[i]);
		ASSERT_POINT_EQ(pa_decoded, pa_ref_decoded, 4);
	}
	/* 2_32 decode */
	int32_t d32[2], r32[2];
	for (int i = 0; i < NPAR; i++) {
		morton_get_2_32(d32, pa_codes[i]);
		ref_get_2_32(r32, pa_codes[i]);
		ASSERT(d32[0] == r32[0] && d32[1] == r32[1]);
	}
}

TEST(bulk3_decode_parity) {
#if GEO_SIMD_MORTON
	int16_t orig[NPAR][3], bulk_out[NPAR][3], scalar_out[NPAR][3];

	test_seed_rng(0xB3);
	for (int i = 0; i < NPAR; i++) {
		for (int d = 0; d < 3; d++)
			orig[i][d] = test_rand_coord();
		pa_codes[i] = morton_set_3(orig[i]);
	}

	morton_get_bulk(bulk_out, pa_codes, NPAR);
	for (int i = 0; i < NPAR; i++)
		morton_get_3(scalar_out[i], pa_codes[i]);
	ASSERT(memcmp(bulk_out, scalar_out, sizeof scalar_out) == 0);

	/* also verify round-trip recovery */
	for (int i = 0; i < NPAR; i++) {
		ASSERT_POINT_EQ(bulk_out[i], orig[i], 3);
	}
#else
	printf("(bulk API disabled — skipped) ");
#endif
}

TEST(bulk4_decode_parity) {
#if GEO_SIMD_MORTON
	int16_t orig[NPAR][4], bulk_out[NPAR][4], scalar_out[NPAR][4];

	test_seed_rng(0xB4);
	for (int i = 0; i < NPAR; i++) {
		for (int d = 0; d < 4; d++)
			orig[i][d] = test_rand_coord();
		pa_codes[i] = morton_set_4(orig[i]);
	}

	morton_get_bulk4(bulk_out, pa_codes, NPAR);
	for (int i = 0; i < NPAR; i++)
		morton_get_4(scalar_out[i], pa_codes[i]);
	ASSERT(memcmp(bulk_out, scalar_out, sizeof scalar_out) == 0);

	for (int i = 0; i < NPAR; i++)
		ASSERT_POINT_EQ(bulk_out[i], orig[i], 4);
#else
	printf("(bulk API disabled — skipped) ");
#endif
}

int main(void) {
	test_suite_begin("Codec Parity");
	RUN_TEST(pdep_encode_parity);
	RUN_TEST(pdep_decode_parity);
	RUN_TEST(bulk3_decode_parity);
	RUN_TEST(bulk4_decode_parity);
	return test_suite_end();
}
