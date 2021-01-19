#ifndef VAR_ALIGNER_H_
#define VAR_ALIGNER_H_

#include "form.h"
#include "kseq1.h"
#include "bseq.h"
#include "sdp.h"
#include <zlib.h>

typedef struct {
    int chr_name_i;
    int dir;
    uint32_t rd_off;   
    uint32_t rf_off;
    uint32_t mem_len;
    //lv is read len, lv_i+lv is ref len
    //int lv_l, lv_l_i, lv_r, lv_r_i;     
    uint32_t rd_l, rd_r;
    uint32_t rf_l, rf_r;
    uint32_t sd_len;

    //for path check
    int path_range;
    int path_id;
    int is_valid;
	int is_exd;
    int dp_max;

    //path merge range
    uint32_t rd_l_rg, rd_r_rg, rf_l_rg, rf_r_rg;
    int64_t rd_l_bdy, rd_r_bdy, rf_l_bdy, rf_r_bdy;
    int32_t idx_l, idx_r;

}sd_hit_t;

typedef struct {
    kpht_t *kht;
    sd_hit_t **sd_hit_a, **sd_hit_rst;
    sd_hit_t *sd_hit_b;
    uint8_t *read_seq, *read_seqr;
    int **sd_bt_path, *sd_bt_path_ord;
    int *sd_hit_rst_n;
    int *sd_bt_n;
} rd_handle_dt;


//global variant
extern uint8_t _db_k;  	// deBGA kmer len
extern uint8_t _seed_k;	// deVar seed kmer len
extern uint8_t _sdp_k; 	// deVar greed extension kmer len

//opt and ref seq are gobal

void set_sd_hit_nd(sd_hit_t *sd, int32_t hit_n, uint32_t rd_len, uint32_t rd_off, uint32_t rf_off, int chr_name_i, uint32_t mem_len, int lv_l, int lv_l_i, int tl_e, int lv_r, int lv_r_i, int tr_e, int z);

void log_index();

void init_global_value();

int run_aligner(int argc, char* argv[]);

uint32_t get_seed_vl_len(uint8_t *rd, uint32_t read_l, uint8_t _seed_k);
//int run_seed_sdp(sd_hit_t *nd_a, uint32_t nd_n, int *bt_path, int *bt_l, double er_ins, double er_del);
int run_seed_sdp(sd_hit_t *nd_a, uint32_t nd_n, int *bt_path, int *bt_l, double er_ins, double er_del);
int compare_sd_len(const void *p, const void *q, void *t);

int get_range_cover(uint32_t a1, uint32_t a2, uint32_t  b1, uint32_t b2, \
                    uint32_t c1, uint32_t c2, uint32_t  d1, uint32_t d2, \
                    uint32_t *x1, uint32_t *x2, uint32_t  *y1, uint32_t *y2, \
                    int sv_lim, double er_ins, double er_del, double er_mis);

void get_sd_hit_rg(sd_hit_t *_sdt, uint32_t *c1, uint32_t *c2, uint32_t  *d1, uint32_t *d2);

void update_sd_hit_rg(sd_hit_t *_sdt, uint32_t a1, uint32_t a2, uint32_t  b1, uint32_t b2); 

void update_hit_vld(uint32_t *a1, uint32_t *a2, uint32_t  *b1, uint32_t *b2, \
                   uint32_t x1, uint32_t x2, uint32_t  y1, uint32_t y2, int lchr);


int single_core_aln(uint32_t rd_i,  rd_handle_dt * t_rd_d);

void log_rst_rec(uint32_t rd_i, double MAX_MQ);

static void *single_thread_core_aln(void *thread_i_p);
#endif
