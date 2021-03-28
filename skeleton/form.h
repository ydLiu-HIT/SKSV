#ifndef _FORM_H
#define _FORM_H

#include "desc.h"
#include <stdint.h>

#define PATH_LEN 1024
#define WAITING_LEN 448

typedef struct{
    int argc;
    char **argv;
    
    char read_path[PATH_LEN];
    char index_path[PATH_LEN];
    char sv_path[PATH_LEN];
    
    uint8_t db_k; //de bruijn graph index kmer len
    uint8_t seed_k; //seeding kmer len
    uint8_t sdp_k; //local sdp kmer len

    int32_t block_s;

    int data_type; 
    /*
        data_type = 0, unavaialble
        1, pacbio ccs
        2, pacbio clr 
        3, ont 2d
        4, ont 1d 
     */
    double error_rate;
    double er_ins, er_del, er_mis;
    
    uint8_t top_seed_extend; 
    int seed_step;
    int batch_size;
    int max_path_N;
    int rd_mx_gap;
    int rf_mx_gap;
    float secondary_ratio;
    float overlap_ratio;
    int min_chain_score;

    int the_sv_lim;

    //seed match max uniq path number
	int read_kmer_match_mx;
    //seed match uniq path's max ref pos number
    int read_kmer_ref_mx;
    int ref_kmer_hit_mx;
    // max lv extend error
    int max_lv_e;
    int ref_lv_bnd_len;

    uint32_t rht_n_max;
    uint32_t seed_hit_n_max;
    uint16_t rst_seed_hit_n_max;
    uint16_t rht_hs_range;
    int seed_Num;
	int thread_n;
}opt_t;

void opt_init(opt_t* P);
int opt_parse(int argc, char *argv[], opt_t* p);
void opt_log(opt_t* p);
#endif
