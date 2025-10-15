#ifndef MORTON_H
#define MORTON_H

#include <stdint.h>

uint64_t morton_set(int16_t *pos, uint8_t dim);
void morton_get(int16_t *pos, uint64_t morton, uint8_t dim);

#endif
