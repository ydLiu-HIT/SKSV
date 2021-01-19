#ifndef __SDP_H_
#define __SDP_H_
#include <stdio.h>
#include <stdint.h>

#define RHT_L_LEN 5
#define CANON_SEED_S 10

//each read should not exit 64K for R_table at most is read len
//L_table save R_table index
typedef struct {
    //save kmer string "ATCG" mt most kmer 13 mer 5 + 8 mer
    uint32_t *kmer_s;
    uint16_t *kmer_ord;
    uint16_t *kmer_pos;
    //kmer first part [4^5] "AT", at most 64k element
    uint16_t *L_table;
    //kmer remain part "CG", at most 8 mer  
    uint16_t *R_table;
    uint16_t *R_table_n;

    uint8_t kmer_l, kmer_ll;
    uint16_t rec_n_mx;
}rht_t;


typedef struct {
	uint32_t head, tail, len;
}kmerp_t;

//full hash table, each key have maximum _mx value
typedef struct {
    uint8_t kmer_l;
    uint32_t _mx;    
	kmerp_t *ht;
	uint32_t *pos_v, *tmp_pos_v;
}kpht_t;


void init_KPHT(kpht_t *kht, uint8_t kmer_l, uint32_t _mx);
//void build_KPHT(kpht_t *kht, char *seq, uint32_t seq_l);
void build_KPHT(kpht_t *kht, uint8_t *seq, uint32_t seq_l);
void set_kmerp_v(kmerp_t *kp, uint32_t *pv, uint8_t kl, uint32_t v);
void t_build_NODE(kpht_t *kht, char *read, uint32_t read_l, uint32_t read_off, uint64_t *ref, uint32_t ref_l, uint32_t offset, uint16_t set_range);
//int run_greed_kht(kpht_t *kht, uint8_t *read, int32_t rd_l, uint32_t rd_off, uint64_t *ref, int32_t ref_l, uint32_t offset, uint16_t block_s, uint16_t *rd_ex, uint16_t *rf_ex, int dir_fg);

//int run_greed_kht(kpht_t *kht, uint8_t *rd, uint32_t rd_l, uint32_t rd_off, uint64_t *ref, uint32_t ref_l, uint32_t offset, int32_t block_s, double block_s_off, uint16_t *rd_ex, uint16_t *rf_ex, int dir_fg);
int run_greed_kht(kpht_t *kht, uint8_t *rd, uint32_t rd_l, uint32_t rd_off, uint64_t *ref, uint32_t ref_l, uint32_t offset, int32_t block_s, double block_s_off, uint32_t *rd_ex, uint32_t *rf_ex, int dir_fg, int dt_t);
typedef struct {
    uint16_t read_off;
    uint16_t ref_off;
    uint16_t len;
}node_t;

void init_RHT(rht_t *rht, uint8_t kmer_l, uint16_t rec_n_mx);
void build_RHT(rht_t *rht, char *seq, uint32_t seq_l);

void set_NODE(node_t *nd, uint16_t rd_o, uint16_t rf_o, uint16_t len);
void init_NODE_sp(node_t *nd_a, uint32_t max_n);
uint16_t build_NODE(rht_t *rht, node_t *nd_a, char *read, uint32_t read_l, uint32_t read_off, uint64_t *ref, uint32_t ref_l, uint32_t offset, uint16_t set_range);

int run_sdp(node_t *nd_a, uint16_t nd_n, uint32_t waiting_len, int *bt_path, int *bt_li, uint32_t r1, uint32_t r2);

int run_greed(rht_t *rht, char *read, uint32_t rd_l, uint32_t rd_off, uint64_t *ref, uint32_t ref_l, uint32_t offset, uint16_t block_s, uint16_t *rd_ex, uint16_t *rf_ex, int dir_fg);
double get_MAPQ_score(double n_match, double n_error, double ref_len, double error_p);
#endif
