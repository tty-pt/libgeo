#ifndef POINT_H
#define POINT_H

#include <stdint.h>
#include <stdio.h>

static inline void
point_add(int16_t *tar, int16_t *orig,
		int16_t *tr, uint8_t dim)
{
	for (uint8_t i = 0; i < dim; i++)
		tar[i] = orig[i] + tr[i];
}

static inline void
point_sub(int16_t *tar, int16_t *orig,
		int16_t *tr, uint8_t dim)
{
	for (uint8_t i = 0; i < dim; i++)
		tar[i] = orig[i] - tr[i];
}

static inline void
point_min(int16_t *tar, int16_t *a,
		int16_t *b, uint8_t dim)
{
	for (uint8_t i = 0; i < dim; i++)
		tar[i] = a[i] < b[i] ? a[i] : b[i];
}

static inline void
point_max(int16_t *tar, int16_t *a,
		int16_t *b, uint8_t dim)
{
	for (uint8_t i = 0; i < dim; i++)
		tar[i] = a[i] > b[i] ? a[i] : b[i];
}

static inline void
point_copy(int16_t *tar, int16_t *orig, uint8_t dim)
{
	for (uint8_t i = 0; i < dim; i++)
		tar[i] = orig[i];
}

static inline int32_t
point_vol(int16_t *p, uint8_t dim)
{
	int32_t res = 1;

	for (uint8_t i = 0; i < dim; i++)
		res *= p[i];

	return res;
}

static inline void
point_set(int16_t *tar, int16_t value, uint8_t dim)
{
	for (uint8_t i = 0; i < dim; i++)
		tar[i] = value;
}

static inline void
point_debug(char *label, int16_t *p, uint8_t dim)
{
	fprintf(stderr, "%s ( ", label);
	for (uint8_t i = 0; i < dim; i++)
		fprintf(stderr, "%d ", p[i]);
	fprintf(stderr, ")\n");
}

static inline uint64_t
point_idx(int16_t *p, int16_t *s, int16_t *e, uint8_t dim)
{
	int16_t l[dim], ps[dim];

	point_sub(l, e, s, dim);
	point_sub(ps, p, s, dim);

	uint64_t res = 0, mul = 1;

	for (uint8_t i = 0; i < dim; mul *= l[i], i++)
		res += ps[i] * mul;

	return res;
}

#endif
