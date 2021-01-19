#ifndef _BSEQ_H
#define _BSEQ_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "bit_operation.h"
#include "kseq1.h"
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

// for read a block of seq from kseq
typedef struct {
	char *read_seq;
	char *read_seq_r;
	char *name;
	uint32_t read_length;
} rd_seq_dt;


typedef struct {
    int chr_name_i;
    int dir;
    uint32_t rd_l, rd_r;
    uint32_t rf_l, rf_r;
} rst_ent_dt;


typedef struct {
    int32_t rec_n;
    rst_ent_dt *rec;
}rst_info_dt;

typedef struct {
    gzFile fp;
    kseq_t *ks;
}bseq_file_t;


bseq_file_t *bseq_open(const char *fn);

void bseq_close(bseq_file_t *fp);

uint32_t bseq_read(bseq_file_t *fp, uint32_t batch, rd_seq_dt* query_info);

void free_rd_seq(uint32_t read_in, rd_seq_dt *query_info);

//uint32_t bseq_read(bseq_t *bp, uint32_t batch, read_type *rp);

//uint32_t bseq_read_2pass(bseq_file_t *fp, uint32_t batch, seq_io* seqio);

#endif
