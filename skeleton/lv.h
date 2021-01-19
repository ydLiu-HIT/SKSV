#ifndef LV_H_
#define LV_H_

#include "stdio.h"
#include "stdint.h"


int32_t LandauVishkin_char(char *p, int32_t p_l, char *q, int32_t q_l, int32_t max_e);
int32_t LandauVishkin_64(uint8_t *p, int32_t p_l, int32_t p_off, \
        uint64_t *q, int32_t q_l, int32_t q_off, int32_t max_e, int flag, int *_e, int *_i);
int32_t LandauVishkin_64_su(uint8_t *p, int32_t p_l, int32_t p_off, \
        uint64_t *q, int32_t q_l, int32_t q_off, int32_t *memList, int32_t *mem, int32_t max_e, int flag, int *_e, int *_i);
int32_t LandauVishkin_64_log(uint8_t *p, int32_t p_l, int32_t p_off, \
        uint64_t *q, int32_t q_l, int32_t q_off, int32_t max_e, int flag, int *_e, int *_i);

#endif
