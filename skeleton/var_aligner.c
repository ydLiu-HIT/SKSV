#define _GNU_SOURCE
#include <math.h>
#include <string.h>
#include <inttypes.h>
#include <pthread.h>
#include "bit_operation.h"
#include "var_aligner.h"
#include "index_loader.h"
#include "lv.h"
#include "stdlib.h"
#include "stdio.h"
#include "ktime.h"

#define MAX_READ_LEN 1000000

pthread_rwlock_t rwlock;

opt_t *opt;
FILE *fp_sv = NULL;

/*
	read deBGA index											stage 0
	for each read
	find split read, find read's hit to Uniq path as seed		stage 1
	Using seeds to run SDP and sort the result path and seed 	stage 2
	Run greedy extension on each seed and merge coverd seed		stage 3
*/
uint8_t _db_k;  	// deBGA kmer len
uint8_t _seed_k;	// deVar seed kmer len,	stage 1
uint8_t _sdp_k; 	// deVar greed extension kmer len, stage 3
int _max_lv_e;		// Seed extension maximun Landau Viskin error, stage 1
int _read_kmer_match_mx;
int _read_kmer_ref_mx;

int _rd_gap_mx;
int _rf_gap_mx;
int max_path_N;
int _min_chain_score;
float _secondary_ratio;
float _overlap_ratio;

int _thread_n;
rd_handle_dt *_rd_hd;

double _er_ins;		//read error rate
double _er_del;
double _er_mis;	

uint32_t _batch_n;
uint32_t _rht_n_max;
uint32_t _seed_hit_n_max;
int _seed_step;
uint16_t _rst_seed_hit_n_max;

int _ref_lv_bnd_len; //default 200, could be larger or smaller depend on the speed************
int THREAD_READ_I; 
uint32_t _rd_seq_i;
int32_t _rd_block_s; 
double _rf_block_s_off; 

int _sv_lim;
rd_seq_dt *_rd_seq;
rst_info_dt *_rst_info;

int _dt_type; //CCS 1, CLR 2
void log_index()
{
    //change in deBGA, START_POS_REF->0, init first ref space
    fprintf(stderr, "chrosome N: [%u] ", chr_file_n-1);
    //printf("reference_len: [%zu] \n", reference_len);
    fprintf(stderr, "reference_len: [%zu] uniq path len: [%zu]\n", reference_len-START_POS_REF, result_seq*32);
    if (0)
    for (int i=1; i<chr_file_n; ++i)
    {
        fprintf(stderr, ">%s ", chr_names[i]);
        fprintf(stderr, "len:[%zu], ", (chr_end_n[i]-chr_end_n[i-1]));
        size_t offset = chr_end_n[i-1]-1;
   
        fprintf(stderr, " | end:[%zu]\n", chr_end_n[i]-START_POS_REF);
    }
	
    fprintf(stderr, "1[%zu] 2[%zu] uniq_seq len: [%zu], uniq_path: [%zu]\n", result_hash_g, result_kmer_g, result_seq, result_seqf-1);
}

void del_deBGAmemory()
{
	//if (buffer_ref_seq)	free(buffer_ref_seq);   //free after the program finish
	free(buffer_seqf);
	free(buffer_seq);
	free(buffer_pp);
	free(buffer_p);
	free(buffer_hash_g);
	free(buffer_kmer_g);
	free(buffer_off_g);
}


int run_aligner(int argc, char* argv[])
{
    opt_t opt_ent;   
    opt = &opt_ent; //global opt
    opt_init(opt);
    if (opt_parse(argc, argv, opt) != 0)  return 1;
	// ---------------------------  set global variant
    _db_k = opt->db_k;  //22
	_seed_k = opt->seed_k; //18
    _sdp_k = opt->sdp_k; //7
    _max_lv_e = opt->max_lv_e;  //3
	_er_ins = opt->er_ins;
	_er_del = opt->er_del;
	_er_mis = opt->er_mis;
   	_rht_n_max = opt->rht_n_max; //regional read length for hash table 
	_seed_hit_n_max = opt->seed_hit_n_max;
    _seed_step = opt->seed_step; //10
	_rst_seed_hit_n_max = opt->rst_seed_hit_n_max;
	_read_kmer_match_mx = opt->read_kmer_match_mx;
	_read_kmer_ref_mx = opt->read_kmer_ref_mx;
    _rd_gap_mx = opt->rd_mx_gap;
    _rf_gap_mx = opt->rf_mx_gap;
    _thread_n = opt->thread_n;
    _sv_lim = opt->the_sv_lim;
    _dt_type = opt->data_type; //CCS 1, CLR 2
    _ref_lv_bnd_len = opt->ref_lv_bnd_len;
    max_path_N = opt->max_path_N;
    _min_chain_score = opt->min_chain_score;
    _secondary_ratio = opt->secondary_ratio;
    _overlap_ratio = opt->overlap_ratio;
    fp_sv = fopen(opt->sv_path, "w");
    if(fp_sv == NULL)
    {
        fprintf(stderr, "[Wrong] Falied to open file %s!!!\n", opt->sv_path);
        exit(0);
    }
    

	_rd_block_s = opt->block_s;
    _rf_block_s_off = 1.0 * (_er_del - _er_ins);

    double rht_time = 0, sd_time = 0, sd_sdp_time = 0, sd_greed_time = 0, ex_time = 0, x_time = 0;
    
    // --------------------------- 	load deBGA index 
    double tmap = realtime(), atmap, xtmap;
    fprintf(stderr, "loading index files\n");
    load_index_file(opt->index_path);
    //log_index();
    fprintf(stderr, "time used: %3.3fs\n", ((double)(realtime()-tmap)));
    tmap = realtime();
    
    int _cnt_suc_rd=0, _cnt_rd=0, _tot_rd=0, _sd_hit_tmp = 0, _sd_hit_cnt = 0;
    int cnt_sv_s=0, cnt_sv_a=0, cnt_sv_c=0, cnt_sv_c_all=0;
    int read_mx_len=0, read_mn_len=0, read_sum = 0, read_av_n = 0, read_av_patch=1000;
    double read_av_len = 0;
    
    _rd_hd= (rd_handle_dt *)malloc(sizeof(rd_handle_dt) * _thread_n);
    //init all required space for read handler
    kpht_t _kht_ent[_thread_n]; 
    for (int i=0; i<_thread_n; ++i)
    {
        _rd_hd[i].read_seq = (uint8_t*)malloc(sizeof(uint8_t) * (MAX_READ_LEN));
        _rd_hd[i].read_seqr = (uint8_t*)malloc(sizeof(uint8_t) * (MAX_READ_LEN));

        _rd_hd[i].kht = &(_kht_ent[i]);
        init_KPHT(_rd_hd[i].kht, _sdp_k, _rht_n_max);

        _rd_hd[i].sd_hit_a = (sd_hit_t **)malloc(sizeof(sd_hit_t*) * 2);
        _rd_hd[i].sd_hit_a[0] = (sd_hit_t *)malloc(sizeof(sd_hit_t) * (_seed_hit_n_max));
        _rd_hd[i].sd_hit_a[1] = (sd_hit_t *)malloc(sizeof(sd_hit_t) * (_seed_hit_n_max));
        
        //stroe sv segment from both strand
        _rd_hd[i].sd_hit_rst = (sd_hit_t **)malloc(sizeof(sd_hit_t *) * 2);
        _rd_hd[i].sd_hit_rst[0] = (sd_hit_t *)malloc(sizeof(sd_hit_t) * (_rst_seed_hit_n_max));
        _rd_hd[i].sd_hit_rst[1] = (sd_hit_t *)malloc(sizeof(sd_hit_t) * (_rst_seed_hit_n_max));
        _rd_hd[i].sd_hit_rst_n= (int *)malloc(sizeof(int) * 2);

        _rd_hd[i].sd_bt_path = (int **)malloc(sizeof(int*) * 2);
        _rd_hd[i].sd_bt_path[0] = (int *)malloc(sizeof(int) * (_seed_hit_n_max));
        _rd_hd[i].sd_bt_path[1] = (int *)malloc(sizeof(int) * (_seed_hit_n_max));
        _rd_hd[i].sd_bt_n = (int *)malloc(sizeof(int) * 2);
        _rd_hd[i].sd_bt_path_ord = (int *)malloc(sizeof(int) * (_seed_hit_n_max));
    } 
 
    //_batch_n = 655350;
    _batch_n = opt->batch_size;

    _rd_seq = (rd_seq_dt *)malloc(sizeof(rd_seq_dt) * _batch_n);
    _rst_info = (rst_info_dt *)malloc(sizeof(rst_info_dt) * _batch_n);
    for (uint32_t i=0; i<_batch_n; ++i)
    {
        _rst_info[i].rec = (rst_ent_dt *)malloc(sizeof(rst_ent_dt) * _rst_seed_hit_n_max);
    }
    fprintf(stderr, "time used for memory applying: %3.3fs\n", ((double)(realtime()-tmap)));


	bseq_file_t *bf;
    bf = bseq_open(opt->read_path);
    if(bf == 0)
    {
        fprintf(stderr, "wrong input file route or name: %s \n", opt->read_path);
        exit(1);
    }

    pthread_rwlock_init(&rwlock, NULL);

    _rd_seq_i= _batch_n;
    int batch_id = 1;
    while (_rd_seq_i >= _batch_n)
    {
        //read a batch of read into _rd_seq
        atmap = realtime();
        _rd_seq_i  = bseq_read(bf, _batch_n, _rd_seq);
        //fprintf(stderr, "time used for load batch reads: %3.3fs\n", (double)(realtime()-atmap));

		//handle each read
		THREAD_READ_I = 0;
        if (_thread_n <= 1)
        {
		    for (int _ti=0; _ti<_rd_seq_i; ++_ti)
            {
                fprintf(stderr, "%s\n", _rd_seq[_ti].name);
                //if (_ti % 500 == 0)
                //{
                //    fprintf(stderr, "run #%d reads ", _ti);
                //    fprintf(stderr, ">%s                     \r", _rd_seq[_ti].name);
                //}
                uint32_t rd_i = _ti;
                int thread_i = 0;

                single_core_aln(rd_i, _rd_hd+thread_i);
            }
        }
        else
		{
			pthread_attr_t attr;
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

			pthread_t* tid = (pthread_t* )calloc(_thread_n, sizeof(pthread_t));
			for (int t_i = 0; t_i < _thread_n; ++t_i)
			{
				int res = pthread_create(&tid[t_i], &attr, single_thread_core_aln, _rd_hd+t_i);
				if(res != 0)
				{
					fprintf(stderr, "create pthread error");
					exit(1);
				}
            }
			for (int t_i = 0; t_i < _thread_n; ++t_i) pthread_join(tid[t_i], 0);
			free(tid);
		}
        fprintf(stderr, "finish batch %d in %3.3fs\n", batch_id, (double)(realtime()-atmap));
        batch_id += 1;

#ifndef GEN_DATA
        double MAX_MQ = pow(10, -6);
		for (int _ti=0; _ti<_rd_seq_i; ++_ti)
        {
            uint32_t rd_i = _ti;
            log_rst_rec(rd_i, MAX_MQ);
        }
#endif
		//free space
		free_rd_seq(_rd_seq_i, _rd_seq);
    }

    pthread_rwlock_destroy(&rwlock);

	bseq_close(bf);
    fclose(fp_sv);

    //free memory
    for(int i = 0; i < _thread_n; ++i)
    {
        if(_rd_hd[i].read_seq != NULL) free(_rd_hd[i].read_seq);
        if(_rd_hd[i].read_seqr != NULL) free(_rd_hd[i].read_seqr); //can be remove by read_seq
        for(int tt = 0; tt < 2; ++tt)
        {
            if(_rd_hd[i].sd_hit_a[tt] != NULL) free(_rd_hd[i].sd_hit_a[tt]);
        }
        if(_rd_hd[i].sd_hit_a != NULL) free(_rd_hd[i].sd_hit_a);
        for(int t = 0; t < 2; ++t)
            if (_rd_hd[i].sd_bt_path[t] != NULL) free(_rd_hd[i].sd_bt_path[t]);
        if(_rd_hd[i].sd_bt_path != NULL) free(_rd_hd[i].sd_bt_path);
        if(_rd_hd[i].sd_bt_path_ord != NULL) free(_rd_hd[i].sd_bt_path_ord); //can remove outside, def in single_aln_core directly
        if(_rd_hd[i].sd_hit_rst_n != NULL) free(_rd_hd[i].sd_hit_rst_n);
        if(_rd_hd[i].sd_bt_n != NULL) free(_rd_hd[i].sd_bt_n);
        
        for(int tti = 0; tti < 2; ++tti)
        {
            if (_rd_hd[i].sd_hit_rst[tti] != NULL) free(_rd_hd[i].sd_hit_rst[tti]);
        }
        if(_rd_hd[i].sd_hit_rst != NULL) free(_rd_hd[i].sd_hit_rst);
        
        //free HPHT
        //if(_rd_hd[i].kht->ht != NULL) free(_rd_hd[i].kht->ht);
        //if(_rd_hd[i].kht->pos_v != NULL) free(_rd_hd[i].kht->pos_v);
        //if(_rd_hd[i].kht->tmp_pos_v != NULL) free(_rd_hd[i].kht->tmp_pos_v);
        //if(_rd_hd[i].kht != NULL) free(_rd_hd[i].kht);
    } 
    if(_rd_hd != NULL) free(_rd_hd);
    
    del_deBGAmemory();
    
    for(int i = 0; i < _batch_n; ++i)
    {
        if(_rst_info[i].rec != NULL) free(_rst_info[i].rec);
    }
    if(_rst_info != NULL) free(_rst_info);

    double aln_tot_time = ((double)(realtime()-tmap));
	fprintf(stderr, "aln time used: %3.3f                                 \n", aln_tot_time);
    return 0;
        
} 

int compare_sd_len(const void *p, const void *q, void *t)
{
	int f = *(int *)p;
	int h = *(int *)q;
	sd_hit_t *s = (sd_hit_t *)t;
    if (s[f].path_range != s[h].path_range) return s[f].path_range > s[h].path_range;
 	if (s[f].sd_len == s[h].sd_len) return f < h;
 	return s[h].sd_len - s[f].sd_len;
}


uint32_t get_seed_vl_len(uint8_t *read, uint32_t read_l, uint8_t _seed_k)
{
    uint32_t vl_seq_l=0;
    for (uint32_t i=0; i<read_l; ++i) 
    {
        if (read[i] == 4) continue;
        uint32_t start = i;
        while (read[i++] != 4 && i<read_l);
        if (i-start < _seed_k) continue;
        vl_seq_l += (i-start - _seed_k);
    }
    return vl_seq_l;
}


void set_sd_hit_nd(sd_hit_t *sd, int32_t hit_n, uint32_t rd_len, uint32_t rd_off, uint32_t rf_off, int chr_name_i, uint32_t mem_len, int lv_l, int lv_l_i, int tl_e, int lv_r, int lv_r_i, int tr_e, int z)
{
#ifdef ext_hit_nd
    if (hit_n == 0)
    {
        printf("uint32_t __hit_n = %u, __rd_len=%u, __rd_off=%lu, __rf_off=%lu, __mem_len=%u;\n", hit_n, rd_len, rd_off, rf_off, mem_len);
        printf("int __chr_name_i=%d, __lv_l=%d, __lv_l_i=%d, __tl_e=%d, __lv_r=%d, __lv_r_i=%d, __tr_e=%d, __z=%d;\n", chr_name_i, lv_l, lv_l_i, tl_e, lv_r, lv_r_i, tr_e, z);
    }
#endif

	sd->chr_name_i = chr_name_i;
    sd->dir = z;
    sd->rf_off = rf_off;
    sd->mem_len = mem_len;
    sd->rd_off = rd_off;
    //int v_l = lv_l-tl_e, v_r = lv_r-tr_e;
        
    sd->rd_l = rd_off - lv_l;
    sd->rd_r = rd_off + mem_len + lv_r - 1;
    sd->rf_l = rf_off - (lv_l + lv_l_i);
    sd->rf_r = rf_off + mem_len + (lv_r + lv_r_i) - 1;

    sd->sd_len = MN((sd->rd_r)-(sd->rd_l), (sd->rf_r)-(sd->rf_l));
    sd->rd_l_rg = sd->rd_l;
    sd->rd_r_rg = sd->rd_r;
    sd->rf_l_rg = sd->rf_l;
    sd->rf_r_rg = sd->rf_r;

    sd->rd_l_bdy = -1;
    sd->rd_r_bdy = -1;
    sd->rf_l_bdy = -1;
    sd->rf_r_bdy = -1;

    sd->is_exd = 0;
    //1 base
    sd->idx_l = hit_n;
    sd->idx_r = hit_n+1;
}
#define INF 0xfffff

int can_be_removed(sd_hit_t *sd_hit1, sd_hit_t *sd_hit2, int THRESHOLD)
{
    if(sd_hit2->rd_l>sd_hit1->rd_l && sd_hit2->rf_l>sd_hit1->rf_l && sd_hit2->rd_r>sd_hit1->rd_r && sd_hit2->rf_r>sd_hit1->rf_r && (int64_t)(sd_hit2->rf_l-sd_hit1->rf_r) < THRESHOLD)   
        return 1;
    else
        return 0;
}


int overlap_ratio(uint32_t s1, uint32_t e1, uint32_t s2, uint32_t e2)
{
    uint32_t start, end;
    start = MX(s1, s2);
    end = MN(e1, e2);
    if (start > end)
        return 0; //no overlap
    else
    {
        int sub = end - start;
        if(sub/(float)(e2-s2) > _overlap_ratio || sub/(float)(e1-s1) > _overlap_ratio) //0.4 consider?
            return 1;
        else
            return 0;
    }
}

int run_seed_sdp(sd_hit_t *nd_a, uint32_t nd_n, int *bt_path, int *bt_l, double er_ins, double er_del)
{
    int32_t dp[nd_n], dp_mx=-INF, dp_mx_pre = -INF;
    int bt[nd_n], dp_mx_i=0, dp_mx_i_s = 0, tpi=0;   
    int cpath[nd_n], tmp_path_1[nd_n], tmp_path_2[nd_n], *lpath, *mpath, tmp_l = 0;
    int cpath_i = 0, lpath_i = 0, mpath_i = 0, real_o = 0;
    uint8_t is_inp[nd_n];
    uint8_t is_upt[nd_n];
    memset(is_inp, 0, sizeof(is_inp));
    memset(is_upt, 0, sizeof(is_upt));
    for(int o=0; o<nd_n; ++o)
    {
        bt[o] = -1;
        dp[o] = nd_a[o].sd_len + 1;
        if (dp[o] > dp_mx)
        {
            dp_mx_i = o;
            dp_mx = dp[o];
        }
        nd_a[o].dp_max = -1;
    }

    //if (nd_n>3)
    float c;
    uint32_t cur_rds, cur_rde, last_rds, last_rde;
    int64_t rdl1, rdr1, rfl1, rfr1, rdl2, rdr2, rfl2, rfr2;

    int if_all_end = 0;
    for (int o=0; o<max_path_N; ++o)
    {
        if_all_end = 0;
        if (o)  
        {
            if (is_inp[0] == 0)
            {
                dp_mx = dp[0];
                dp_mx_i = 0;
            }
            else
            {
                dp_mx = -INF;
            }
        }

        for (int  i=1; i<(nd_n); ++i)
        {
            int MAT = 0;
            //if nd in path skip
            if (is_inp[i] > 0) continue;
            //when find 2nd path
            if (o)
            {
                //if prv nd not update, no need to run dp
                if (bt[i] == -1 || is_upt[bt[i]]==0)
                {
                    if (bt[i] != -1) if_all_end = 1;
                    if (dp[i] > dp_mx)
                    {
                        dp_mx_i = i;
                        dp_mx = dp[i];
                    }
                    continue;
                }
                else
                {
                    is_upt[i] = 1;
                    bt[i] = -1;
                    dp[i] = nd_a[i].sd_len + 1;
                    if (dp[i] > dp_mx)
                    {
                        dp_mx_i = i;
                        dp_mx = dp[i];
                    }
                }
            }
            
            for (int j=i-1; j>=0; --j)
            {
                if (MAT > 50) break;
                if (is_inp[j] > 0) continue;
               
                // in same chrosome 
                if (nd_a[j].chr_name_i != nd_a[i].chr_name_i)
                    continue;
                
                rdl1 = nd_a[j].rd_l; rdl2 = nd_a[i].rd_l;
                rdr1 = nd_a[j].rd_r; rdr2 = nd_a[i].rd_r;
                rfl1 = nd_a[j].rf_l; rfl2 = nd_a[i].rf_l;
                rfr1 = nd_a[j].rf_r; rfr2 = nd_a[i].rf_r;

                if (rdl1 < rdl2 && rfl1 < rfl2 && rdr1 < rdr2 && rfr1 < rfr2 && (rdl2-rdr1 < _rd_gap_mx) && (rfl2-rfr1 < _rf_gap_mx))
                { 
                    int64_t rd_dis, rf_dis, rd_len, rf_len, span_len = 0, diff_len = 0, gain_len = 0, diff_rg = 0, dp_wgt;
                    rd_dis = (rdl2 - rdr1);
                    rf_dis = (rfl2 - rfr1);
                    if (rd_dis > 0 && rf_dis > 0)
                        span_len = MN(rd_dis, rf_dis) + 1;
                    
                    diff_len = G_ABS(rd_dis-rf_dis);
                    rd_len = (rdl2 > rdr1) ? (rdr2-rdl2) : (rdr2-rdr1);
                    rf_len = (rfl2 > rfr1) ? (rfr2-rfl2) : (rfr2-rfr1);
                    gain_len = MN(rd_len, rf_len) + 1;

                    double ins_off = span_len * er_ins * 2;
                    double del_off = span_len * er_del * 2;
                    if (rd_dis == rf_dis) 
                        diff_rg = MX(ins_off, del_off);
                    else
                        diff_rg = rd_dis > rf_dis ? ins_off : del_off;
                    //if (diff_len > diff_rg+_sv_lim)  continue;
                    MAT += 1; 
                    
                    //c = (diff_len > 0)? 0.5 * log10(diff_len) + 0.01 * 20 * diff_len : 0;
                    //dp_wgt = gain_len - (int)(diff_len * c);
                    dp_wgt = gain_len - diff_len;
                    
                    //if ( i== 328 || i == 329 || i == 330 ||i==675 || i==676)
                    //{
                    //    printf("dp[%d]=%d, dp[%d]=%d, dp_wgt=%ld, gain_len = %d, diff_len = %ld, c= %f\n", i, dp[i], j, dp[j], dp_wgt, gain_len, diff_len, c);
                    //}

                    if (dp[i] < dp[j] + dp_wgt)
                    {
                        dp[i] = dp[j] + dp_wgt;
                        if (dp[i] > dp_mx)
                        {
                            dp_mx_i = i;
                            dp_mx = dp[i];
                        }
                        bt[i] = j;
                        if_all_end = 1;
                    }
                }
            }
        }

        int path_i=0, cur_i=dp_mx_i;
        if(dp_mx < _min_chain_score) break;

        cur_rde = nd_a[cur_i].rd_r;
        cpath_i = 0;
#ifdef print_result
        printf("path: %d----->, dp_mx = %d\n", o, dp_mx);
#endif
        while (cur_i != -1)
        {
#ifdef print_result
            printf("%d->", cur_i);
#endif
            cpath[cpath_i++] = cur_i;
            is_inp[cur_i] = 1;
            is_upt[cur_i] = 1;

            cur_rds = nd_a[cur_i].rd_l;
            cur_i = bt[cur_i];
        }
#ifdef print_result
        printf("\n");
#endif
        
        int olp = -1;
        if(o > 0)
        {
            olp = overlap_ratio(last_rds, last_rde, cur_rds, cur_rde);
            //fprintf(stderr, "o = %d, olp = %d, last_rds = %u, last_rde = %u, cur_rds = %u, cur_rde = %u, dp_mx = %d, dp_mx_pre = %d\n", o, olp, last_rds, last_rde, cur_rds, cur_rde, dp_mx, dp_mx_pre);
            if(olp && dp_mx / (float)dp_mx_pre < _secondary_ratio)
            //if(dp_mx / (float)dp_mx_pre < _secondary_ratio)
            {
                cpath_i = 0;
                continue;
            }
        }

        nd_a[cur_i].dp_max = dp_mx;
        if(o == 0)
        {
            last_rds = cur_rds;
            last_rde = cur_rde;
            dp_mx_pre = dp_mx;   
        }

        //mark the overaped pos in in path
        if (0)
        for (int  i=0; i<(nd_n); ++i)
        {
            if (is_inp[i] == 0)
                continue;
            rdl1 = nd_a[i].rd_l; 
            rdr1 = nd_a[i].rd_r; 

            //left
            int _j = i-1;
            while (_j >= 0 && is_inp[_j] == 0 && nd_a[_j].rd_r > rdl1)
            {
                is_inp[_j] = 2;
                is_upt[_j] = 2;
                _j -= 1;
            }

            //right
            _j = i+1;
            while (_j < nd_n && is_inp[_j] == 0 && nd_a[_j].rd_l < rdr1)
            {
                is_inp[_j] = 2;
                is_upt[_j] = 2;
                _j += 1;
            }
        }

        if (real_o%2 == 0) {lpath = tmp_path_1; mpath = tmp_path_2;}
        else    {lpath = tmp_path_2; mpath = tmp_path_1;}

        mpath_i = 0;
        int li = 0;
        if (cpath_i > 0)
            for (int i=0, ld, lf; i<(cpath_i); ++i)
            {
                int ti = cpath[cpath_i-1-i];
                nd_a[ti].path_range = real_o;
                while (li < lpath_i && lpath[li] < ti)  mpath[mpath_i++] = lpath[li++];
                mpath[mpath_i++] = ti; 
            }
        while (li < lpath_i)  mpath[mpath_i++] = lpath[li++];
 
        lpath_i = mpath_i; 
        real_o += 1;
    }

    *bt_l = mpath_i;
    for (int i=0, ld, lst_ti, lf; i<(mpath_i); ++i)
    {
        int ti = mpath[i];
        bt_path[i] = ti;
        nd_a[ti].path_id = i;
        nd_a[ti].is_valid = 1;
        if (i>0)
        {
            nd_a[lst_ti].idx_r = ti+1;
            nd_a[ti].idx_l = lst_ti+1;
        }
        else
            nd_a[ti].idx_l = 0;
        if (i==mpath_i)
            nd_a[ti].idx_r = nd_n+1;
        lst_ti = ti;
    }
    
    return dp_mx_pre;
    //return (dp_mx_pre<_min_chain_score)?_min_chain_score:dp_mx_pre;
}

int is_in_(uint32_t t, uint32_t l, uint32_t r)
{
    return (t>=l && t <= r);
}

float overlap_l(uint32_t s1, uint32_t e1, uint32_t s2, uint32_t e2)
{
    uint32_t start, end;
    start = MX(s1, s2);
    end = MN(e1, e2);
    if (start > end)
        return 0.0;
    else
        return (end-start)/(float)(e2-s2);
}

//check two block whether is SV gapped
//0 no, 1 have
int check_sv_gap(int rd_dis, int rf_dis, int sv_lim, double er_ins, double er_del, double er_mis)
{
    int span_len = 0, diff_len = 0, diff_off = 0;
    double tmp_v = 0;

    if (rd_dis>0 && rf_dis>0)
        span_len = MN(rd_dis, rf_dis);
    diff_len = G_ABS(rd_dis-rf_dis);

    if (span_len > 0)
        tmp_v = (rd_dis > rf_dis) ? (span_len * er_del * 2) : (span_len * er_ins * 2);
    diff_off = (int)(tmp_v);
    
    diff_off = sv_lim < diff_off ? 0 : sv_lim - diff_off;
#ifdef print_mg_info
    printf("(%d %d %d %d %d %d) ", rd_dis, rf_dis, span_len, diff_len, diff_off, \
           diff_len <= diff_off && span_len <= sv_lim );
#endif
    //if (diff_len <= diff_off && (int)(span_len*er_mis) + diff_len <= diff_off)
    if (diff_len <= diff_off && span_len <= 2*sv_lim)
        return 0;
    else
        return 1;
}

//if two range both have overlap merge
int get_range_cover(uint32_t a1, uint32_t a2, uint32_t  b1, uint32_t b2, \
                    uint32_t c1, uint32_t c2, uint32_t  d1, uint32_t d2, \
                    uint32_t *x1, uint32_t *x2, uint32_t  *y1, uint32_t *y2, \
                    int sv_lim, double er_ins, double er_del, double er_mis)
{
    int rd_dis=0, rf_dis=0, span_len=0, diff_len = 0, diff_off = 0;
    double tmp_v=0;

    if (a1 <= c1)
    {
        if (c1<=a2 && c2<=a2) //the next node included in the previous one
        {
            //subset
            *x1 = a1;   *y1 = b1;
            *x2 = a2;   *y2 = b2;
            return 1;
        }
        rd_dis = c1 - a2; rf_dis = d1 - b2;
    }
    else
    {
        if (a1<=c2 && a2<=c2) //the previous node included in the next one
        {
            //fprintf(stderr, "impossible!\n");
            //subset
            *x1 = c1;   *y1 = d1;
            *x2 = c2;   *y2 = d2;
            return 1;
        }
        rd_dis = a1 - c2; rf_dis = b1 - d2;
    } 
    
    if (rd_dis>0 && rf_dis>0)
        span_len = MN(rd_dis, rf_dis);
    diff_len = G_ABS(rd_dis-rf_dis);

    if (span_len > 0)
        tmp_v = (rd_dis > rf_dis) ? (span_len * er_ins * 2) : (span_len * er_del * 2);
    diff_off = (int)(tmp_v);
    
    //printf("diff_off = %d, sv_lim = %d\n", diff_off, sv_lim);
    diff_off = sv_lim < diff_off ? 0 : sv_lim - diff_off;
    if (diff_len > diff_off)  
    {
        return 0; //have variant
    }
    *x1 = MN(a1, c1);   *y1 = MN(b1, d1);
    *x2 = MX(a2, c2);   *y2 = MX(b2, d2);
#ifdef refine_bnd
    if (span_len > sv_lim * 2)
    {
        return 2;
    }
    else
    {
        return 1;
    }
#endif
    return 2;
    
    if (is_in_(a1, c1, c2))
    {
        rd_dis = c2 - a1; rf_dis = d2 - b1;
    }
    else if (is_in_(a2, c1, c2))
    {
        rd_dis = a2 - c1; rf_dis = b2 - d1;
    }
    else if (is_in_(c1, a1, a2))
    {
        if (is_in_(c2, a1, a2))
        {
            *x1 = a1;   *y1 = b1;
            *x2 = a2;   *y2 = b2;
            printf("(%d %d %d %d a)", a1, a2, c1, c2);
            return 1;
        }
        rd_dis = c1 - a1; rf_dis = d1 - b1;
    }
    else    
    {
        if (a2 < c1 && c1-a2 < sv_lim)
        {
            rd_dis = c1 - a2; rf_dis = d1 - b2;
        }
        else if (a1 > c2 && a1-c2 < sv_lim)
        {
            rd_dis = a1 - c2; rf_dis = b1 - d2;
        }
        else  
        {
            printf("(%d %d %d %d A)", a1, a2, c1, c2);
            return 0;
        }

    }
    if (rd_dis > 0 && rf_dis > 0)
        span_len = MN(rd_dis, rf_dis);
    diff_len = G_ABS(rd_dis-rf_dis);
    //diff_off = rd_dis > rf_dis ? (span_len * er_ins * 2) : (span_len * er_del * 2);
    //diff_off += sv_lim;

    diff_off = rd_dis < rf_dis ? (span_len * er_ins * 2) : (span_len * er_del * 2);
    diff_off = sv_lim < diff_off ? 0 : sv_lim - diff_off;

    if (diff_len > diff_off)  
    {
        printf("(%d %d %d %d %d %d %d %d B)", diff_len, diff_off, rd_dis, rf_dis, a1, a2, c1, c2);
        return 0;
    }
    *x1 = MN(a1, c1);   *y1 = MN(b1, d1);
    *x2 = MX(a2, c2);   *y2 = MX(b2, d2);
        printf("(%d %d %d %d %d %d b)", rd_dis, rf_dis, a1, a2, c1, c2);
    return 1;
}

void get_sd_hit_rg(sd_hit_t *_sdt, uint32_t *c1, uint32_t *c2, uint32_t  *d1, uint32_t *d2) 
{
    *c1 = _sdt->rd_l_rg; *d1 = _sdt->rf_l_rg;
    *c2 = _sdt->rd_r_rg; *d2 = _sdt->rf_r_rg;
}

void update_sd_hit_rg(sd_hit_t *_sdt, uint32_t a1, uint32_t a2, uint32_t  b1, uint32_t b2) 
{
    _sdt->rd_l = a1; _sdt->rf_l  = b1;
    _sdt->rd_r = a2; _sdt->rf_r  = b2;
    _sdt->rd_l_rg = a1; _sdt->rf_l_rg = b1;
    _sdt->rd_r_rg = a2; _sdt->rf_r_rg = b2;
}

void update_hit_vld(uint32_t *a1, uint32_t *a2, uint32_t  *b1, uint32_t *b2, \
                   uint32_t x1, uint32_t x2, uint32_t  y1, uint32_t y2, int lchr)
{
    if(y1 < chr_end_n[lchr-1] || y2 > chr_end_n[lchr])
        return ;
    
    *a1 = x1; *b1 = y1;
    *a2 = x2; *b2 = y2;
}


void do_sv_detection(uint32_t rd_i, rd_handle_dt * t_rd_d, rst_ent_dt *t_rst_rec, int32_t *t_rst_rec_n, sd_hit_t **sd_hit_rst, int *sd_hit_rst_n, sd_hit_t **sd_hit_a, int **sd_bt_path, int *sd_bt_n)
{
#ifdef print_seed
    printf("do sv detection\n");
#endif

    int *sd_bt_path_ord = t_rd_d->sd_bt_path_ord;
 
    uint32_t seq_l = _rd_seq[rd_i].read_length;
    kpht_t *kht = t_rd_d->kht;
    uint8_t *rd_seq;
    int sv_lim = _sv_lim;
    for(int z = 0; z < 2; ++z)
    {
#ifdef print_seed
        if(z==1)    printf("\n\nReverse--------------------------------------\n");
#endif
        rd_seq = (z == 0)? t_rd_d->read_seq : t_rd_d->read_seqr;
        build_KPHT(kht, rd_seq, seq_l);
        for (int o=0; o<sd_bt_n[z]; ++o)
        {
            sd_bt_path_ord[o] = sd_bt_path[z][o];
        }
#ifdef print_hit_exd
        printf("\nsd_bt_n = %d\n", sd_bt_n[z]);
        printf("before sort: ");
        for (int o=0; o<sd_bt_n[z]; ++o)
        {
            printf("%d(%d)--", sd_bt_path[z][o], sd_hit_a[z][sd_bt_path[z][o]].path_range);
        }
        printf("\n");
#endif
        //sort all path node
        qsort_r(sd_bt_path_ord, sd_bt_n[z], sizeof(int), compare_sd_len, sd_hit_a[z]);

#ifdef print_hit_exd
        printf("after sort: ");
        for(int o=0; o<sd_bt_n[z];++o)
            printf("%d--", sd_bt_path_ord[o]);
        printf("\n");

#endif

#ifdef merge_colinear
        
        //merge when co-linear
        int _cur_path_range, _cur_path_id, _cur_dir;
        int leftmost_id = -1, rightmost_id = sd_bt_n[z];
//#ifdef print_sv
//        for (int i = 0; i < rightmost_id; ++i )
//        {
//            int ti = sd_bt_path[z][i];
//            sd_hit_t *sdt = sd_hit_a[z]+ti;
//            printf("path_id = %d, node_id = %d, path_range = %d, is_valid = %d, rdl = %u, rdr = %u, rfl = %u, rfr = %u\n", i, ti, sdt->path_range, sdt->is_valid, sdt->rd_l, sdt->rd_r, sdt->rf_l, sdt->rf_r);
//        }
//#endif

        for (int o = 0; o < sd_bt_n[z]; ++o)
        {
            int nid = sd_bt_path_ord[o];
            sd_hit_t *sdt = sd_hit_a[z]+nid;
            if(sdt->is_valid == 0) continue;
            _cur_path_id = sdt->path_id;
            _cur_path_range = sdt->path_range;
            _cur_dir = sdt->dir;

            int64_t a, b, c, d;
            int merge_thre = 25;
            //right merge
            int _ti = _cur_path_id + 1;
            while (_ti < rightmost_id )
            {
                sd_hit_t *sdt1 = sd_hit_a[z] + sd_bt_path[z][_ti];
                if(sdt1->path_range != _cur_path_range || sdt1->is_valid == 0 || sdt1->dir != _cur_dir)
                {
                    _ti++; continue;
                }
                a = (int64_t)sdt->rd_r; b = (int64_t)sdt->rf_r; c = (int64_t)sdt1->rd_l; d = (int64_t)sdt1->rf_l;
                if (abs(c+b-a-d) < merge_thre)
                {
                    sdt->rd_r = sdt1->rd_r;
                    sdt->rf_r = sdt1->rf_r;
                    sdt1->is_valid = 0;

                    _ti ++;
                }
                else
                {
                    sdt->rd_r_bdy = c;
                    sdt->rf_r_bdy = d;
                    break;
                }
            }

            _ti = _cur_path_id - 1;
            while (_ti > leftmost_id)
            {
                sd_hit_t *sdt1 = sd_hit_a[z] + sd_bt_path[z][_ti];
                if(sdt1->path_range != _cur_path_range || sdt1->is_valid == 0 || sdt1->dir != _cur_dir)
                {
                    _ti--; continue;
                }
                a = (int64_t)sdt->rd_l; b = (int64_t)sdt->rf_l; c = (int64_t)sdt1->rd_r; d = (int64_t)sdt1->rf_r;
                if (abs(a+d-c-b) < merge_thre)
                {
                    sdt->rd_l = sdt1->rd_l;
                    sdt->rf_l = sdt1->rf_l;
                    sdt1->is_valid = 0;

                    _ti --;
                }
                else
                {
                    sdt->rd_l_bdy = c;
                    sdt->rf_l_bdy = d;
                    break;
                }
            }
        }
        
#ifdef print_sv
        for (int i = 0; i < rightmost_id; ++i )
        {
            int ti = sd_bt_path[z][i];
            sd_hit_t *sdt = sd_hit_a[z]+ti;
            if (sdt->is_valid == 1)
                printf("path_id = %d, node_id = %d, path_range = %d, is_valid = %d, rdl = %u, rdr = %u, rfl = %u, rfr = %u\n", i, ti, sdt->path_range, sdt->is_valid, sdt->rd_l, sdt->rd_r, sdt->rf_l, sdt->rf_r);
        }
#endif

#endif
        
        for (int o=0, ti; o<sd_bt_n[z]; ++o)
        {
            ti = sd_bt_path_ord[o];
            sd_hit_t *sdt = sd_hit_a[z]+ti;
            

            if (sdt->is_valid == 0)   continue;
            //printf("here, path_id = %d, node_id = %d, path_range = %d, is_valid = %d, rdl = %u, rdr = %u, rfl = %u, rfr = %u\n", o, ti, sdt->path_range, sdt->is_valid, sdt->rd_l, sdt->rd_r, sdt->rf_l, sdt->rf_r);

            uint32_t t_rd_l = sdt->rd_l, t_rd_r = sdt->rd_r;
            uint32_t t_rf_l = sdt->rf_l, t_rf_r = sdt->rf_r;
            if(t_rf_l < chr_end_n[sdt->chr_name_i-1])
            {
                t_rd_l += (chr_end_n[sdt->chr_name_i-1]-t_rf_l);
                t_rf_l = chr_end_n[sdt->chr_name_i-1];
            }
            if(t_rf_r > chr_end_n[sdt->chr_name_i])
            {
                t_rd_r -= (t_rf_r - chr_end_n[sdt->chr_name_i]);
                t_rf_r = chr_end_n[sdt->chr_name_i];
            }

            int64_t rb1 = sdt->rd_r_bdy, rfb1 = sdt->rf_r_bdy;
            uint32_t bi=0, rd_ex_r=0, rf_ex_r=0, rd_ex_l=0, rf_ex_l=0;
            uint32_t a1, a2, b1, b2, c1, c2, d1, d2, x1, x2, y1, y2;
            a1 = t_rd_l;   b1 = t_rf_l;  
            a2 = t_rd_r;   b2 = t_rf_r;
            int _cur_path_range, _cur_path_id, _cur_dir;
            _cur_path_id = sdt->path_id;
            _cur_path_range = sdt->path_range;
            _cur_dir = sdt->dir;
            int last_vi = -1, have_B = 0;

#ifdef print_hit_exd
            //printf("right extension with greedy strategy**********************************************\n");
#endif
            //right extend
            while (1)
            {
                if (rfb1 == -1)
                {
                    rfb1 = chr_end_n[ (sdt)->chr_name_i ]-START_POS_REF;
                    rb1 = seq_l-1;
                }
                /*
                input t_rd_r, rb1, a2, b2
                       _path_id, path_range, dir
                output rd_ex_r */
                bi = rd_ex_r = rf_ex_r = 0;
                
                double av_len = 0;
                if (rb1 > t_rd_r && rfb1 > t_rf_r)
                    av_len = MN((rb1-t_rd_r), (rfb1-t_rf_r));

                if ((av_len) > 1.0*_rd_block_s)
                {
                    if (_dt_type == 1)
#ifdef merge_colinear
                        bi = run_greed_kht(kht, rd_seq, rb1-t_rd_r, t_rd_r+1, buffer_ref_seq, MN(rfb1-t_rf_r, chr_end_n[sdt->chr_name_i]-t_rf_r), t_rf_r+1, _rd_block_s, _rf_block_s_off, &rd_ex_r, &rf_ex_r, 1, _dt_type);
#else
                        bi = run_greed_kht(kht, rd_seq, seq_l-1-t_rd_r, t_rd_r+1, buffer_ref_seq, MN(seq_l*2, chr_end_n[sdt->chr_name_i]-t_rf_r), t_rf_r+1, _rd_block_s, _rf_block_s_off, &rd_ex_r, &rf_ex_r, 1, _dt_type);
#endif
                    else
                        bi = run_greed_kht(kht, rd_seq, rb1-t_rd_r, t_rd_r+1, buffer_ref_seq, MN(rfb1-t_rf_r, chr_end_n[sdt->chr_name_i]-t_rf_r), t_rf_r+1, _rd_block_s, _rf_block_s_off, &rd_ex_r, &rf_ex_r, 1, _dt_type);
                }
                
                a2 = t_rd_r + rd_ex_r;   b2 = t_rf_r + rf_ex_r;
#ifdef print_hit_exd
                printf("rd_ex_r = %d, a2 = %d, rf_ex_r = %u, b2 = %u\n", rd_ex_r, a2, rf_ex_r, b2);
#endif
                
                last_vi = -1;
                int last_t = 1;
                if (_cur_path_id < sd_bt_n[z]-1)
                {
                    for (int _i=_cur_path_id+1, new_ti, _t; _i<sd_bt_n[z]; ++_i)
                    {
                        new_ti = sd_bt_path[z][_i];
                        sd_hit_t *_sdt = sd_hit_a[z]+new_ti;
                        if (_cur_dir != _sdt->dir) break;
                        if (_sdt->is_valid == 0) continue;

                        get_sd_hit_rg(_sdt, &c1, &c2, &d1, &d2);
                        //a1[rd_l], a2[rd_r+rd_ex_r], b1[rf_l], b2[rf_r+rf_ex_r]
                        //c1[rd_l_rg], c2[rd_r_rg], d1[rf_l_rg], d2[rf_r_rg]
                        _t = get_range_cover(a1, a2, b1, b2, c1, c2, d1, d2, \
                                            &x1, &x2, &y1, &y2, sv_lim, _er_ins, _er_del, _er_mis);

                        if (_t)
                        {
                            if (_t == 2) have_B = 1;
                            else _sdt->is_valid = 0;

                            if (_t==1 && _cur_path_range >= _sdt->path_range && sdt->chr_name_i == _sdt->chr_name_i) 
                            {
                                update_hit_vld(&a1, &a2, &b1, &b2, x1, x2, y1, y2, sdt->chr_name_i);
                                last_vi = _i;
                                last_t = _t;
                            }
                        }
                        else if(_cur_path_range >= _sdt->path_range) 
                        {
                            break;
                        }
                    }
                }
                    
                if (last_vi == -1) last_vi = _cur_path_id;
                
                if (_cur_path_id+1 < last_vi)
                    for (int _i=_cur_path_id+1, new_ti; _i<last_vi; ++_i)
                    {
                        new_ti = sd_bt_path[z][_i];
                        sd_hit_t *_sdt = sd_hit_a[z]+new_ti;
                        _sdt->is_valid = 0;
                    }

                for (int _i = last_vi; _i < sd_bt_n[z]; ++_i)
                {
                    sd_hit_t *_sdt = sd_hit_a[z]+sd_bt_path[z][_i];
                    
                    if (_i == _cur_path_id) continue;
                    if (_sdt->is_valid == 0 || sdt->chr_name_i != _sdt->chr_name_i) continue;
                    if (_sdt->rd_l_bdy == -1)
                    {
                        _sdt->rd_l_bdy = a2;    _sdt->rf_l_bdy = b2;
                        continue;
                    }
                    if (_sdt->rd_l_rg > a2)
                    {
                        _sdt->rd_l_bdy = MX(_sdt->rd_l_bdy, a2);
                        _sdt->rf_l_bdy = MX(_sdt->rf_l_bdy, b2); 
                    }
                    else    break;
                }
                if (last_vi != _cur_path_id && last_t > 0)
                {
                    //break;
                    sd_hit_t *_sdt = sd_hit_a[z]+sd_bt_path[z][last_vi];
                    if (_sdt->is_exd == 1 || _sdt->rf_r > chr_end_n[_sdt->chr_name_i])
                        break; //break while
                    uint32_t tmp_rd_r = _sdt->rd_r, tmp_rf_r = _sdt->rf_r;
                    int _rd_dis=(a2 - tmp_rd_r), _rf_dis=(b2 - tmp_rf_r);
                    if (_dt_type == 1)
                    {
                        if (check_sv_gap(_rd_dis, _rf_dis, sv_lim, _er_ins, _er_del, _er_mis) == 1)
                            break;
                    }  
                    t_rd_r = _sdt->rd_r; t_rf_r = _sdt->rf_r;
                    rb1 = _sdt->rd_r_bdy; rfb1 = _sdt->rf_r_bdy;
                    _cur_path_id = _sdt->path_id;
                } 
#ifdef gd_edge_cnt
                else
                    break;
#else
                break;
#endif
            }
            int64_t rb2 = sdt->rd_l_bdy, rfb2 = sdt->rf_l_bdy;
            _cur_path_id = sdt->path_id;
            _cur_path_range = sdt->path_range;
            _cur_dir = sdt->dir;
            
#ifdef print_hit_exd
            //printf("left extension with greedy strategy*********************************************************\n");
#endif
            while (1)
            {
                if (rfb2 == -1)
                {
                    rfb2 = chr_end_n[ (sdt)->chr_name_i - 1]-START_POS_REF;
                    rb2 = 0;
                }
                bi = rd_ex_l = rf_ex_l = 0;
                
                double av_len = 0;
                if ((rb2 < t_rd_l) && (rfb2 < t_rf_l))
                av_len = MN((1.0*t_rd_l - rb2), (1.0*t_rf_l - rfb2));

                if ((av_len) > 1.0*_rd_block_s)
                {
                    if (_dt_type == 1)
#ifdef merge_colinear
                        bi = run_greed_kht(kht, rd_seq, t_rd_l-rb2, t_rd_l-1, buffer_ref_seq, MN(t_rf_l-rfb2, t_rf_l - chr_end_n[sdt->chr_name_i-1]), t_rf_l-1, _rd_block_s, _rf_block_s_off, &rd_ex_l, &rf_ex_l, -1, _dt_type);
#else
                        bi = run_greed_kht(kht, rd_seq, t_rd_l, t_rd_l-1, buffer_ref_seq, MN(seq_l*2, t_rf_l - chr_end_n[sdt->chr_name_i-1]), t_rf_l-1, _rd_block_s, _rf_block_s_off, &rd_ex_l, &rf_ex_l, -1, _dt_type);
#endif
                    else
                        bi = run_greed_kht(kht, rd_seq, t_rd_l-rb2, t_rd_l-1, buffer_ref_seq, MN(t_rf_l-rfb2, t_rf_l - chr_end_n[sdt->chr_name_i-1]), t_rf_l-1, _rd_block_s, _rf_block_s_off, &rd_ex_l, &rf_ex_l, -1, _dt_type);
                }

                a1 = t_rd_l - rd_ex_l;   b1 = t_rf_l - rf_ex_l;  
#ifdef print_hit_exd
                printf("rd_ex_l = %d, a1 = %d, rf_ex_l = %u, b1 = %u\n", rd_ex_l, a1, rf_ex_l, b1);
#endif

                last_vi = -1;
                int last_t = 1;
                if (_cur_path_id > 0)
                {
                    for (int _i=_cur_path_id-1, _t, new_ti; _i>=0; --_i)
                    {
                        new_ti = sd_bt_path[z][_i];
                        sd_hit_t *_sdt = sd_hit_a[z]+new_ti;
                        if (_cur_dir != _sdt->dir) break;
                        if (_sdt->is_valid == 0) continue;
                        get_sd_hit_rg(_sdt, &c1, &c2, &d1, &d2);
                        //mask hit that covered all read region
                        _t = get_range_cover(a1, a2, b1, b2, c1, c2, d1, d2, \
                                            &x1, &x2, &y1, &y2, sv_lim, _er_ins, _er_del, _er_mis);

                        if (_t)
                        {
                            if (_t == 2) have_B = 1;
                            else _sdt->is_valid = 0;

                            if (_t==1 && _cur_path_range >= _sdt->path_range && sdt->chr_name_i == _sdt->chr_name_i) 
                            {
                                update_hit_vld(&a1, &a2, &b1, &b2, x1, x2, y1, y2, sdt->chr_name_i);
                                last_vi = _i;
                                last_t = _t;
                            }
                        }
                        else if(_cur_path_range >= _sdt->path_range)
                        {
                            break;
                        }
                    }
                }
//#ifdef print_hit_exd
//                printf("a1 = %u, a2 = %u\n", a1, a2);
//                printf("cur_path_id = %d, last_vi = %d\n", _cur_path_id, last_vi);
//                printf("order = %d\n", sd_bt_path[z][last_vi]);
//#endif
                if (last_vi == -1) last_vi = _cur_path_id;
                if (_cur_path_id-1 > last_vi)
                    for (int _i=_cur_path_id-1, new_ti; _i>last_vi; --_i)
                    {
                        new_ti = sd_bt_path[z][_i];
                        sd_hit_t *_sdt = sd_hit_a[z]+new_ti;
                        _sdt->is_valid = 0;
                    }
                for (int _i = last_vi; _i >= 0; --_i)
                {
                    sd_hit_t *_sdt = sd_hit_a[z]+sd_bt_path[z][_i];
                    if (_i == _cur_path_id) continue;
                    if (_sdt->is_valid == 0 || sdt->chr_name_i != _sdt->chr_name_i) continue;
                    if (_sdt->rd_r_bdy == -1) 
                    {
                        _sdt->rd_r_bdy = a1;    _sdt->rf_r_bdy = b1;
                        continue;
                    }
                    if (_sdt->rd_r_rg < a1)
                    {
                        _sdt->rd_r_bdy = MN(_sdt->rd_r_bdy, a1);
                        _sdt->rf_r_bdy = MN(_sdt->rf_r_bdy, b1); 
                    }
                    else break;
                }
                if (last_vi != _cur_path_id && last_t > 0)
                {
                    //break;
                    sd_hit_t *_sdt = sd_hit_a[z]+sd_bt_path[z][last_vi];
                    if (_sdt->is_exd == 1 || _sdt->rf_l < chr_end_n[_sdt->chr_name_i-1])
                        break;
                    uint32_t tmp_rd_l = _sdt->rd_l, tmp_rf_l = _sdt->rf_l;
                    int _rd_dis=(tmp_rd_l - a1), _rf_dis=(tmp_rf_l - b1);
                    if (_dt_type == 1)
                    if (check_sv_gap(_rd_dis, _rf_dis, sv_lim, _er_ins, _er_del, _er_mis) == 1)
                        break;

                    
                    t_rd_l = _sdt->rd_l; t_rf_l = _sdt->rf_l;
                    rb2 = _sdt->rd_l_bdy; rfb2 = _sdt->rf_l_bdy;
                    _cur_path_id = _sdt->path_id;
                }
#ifdef gd_edge_cnt
                else
                    break;
#else
                break;
#endif
            }

            sdt->is_exd = 1;
            update_sd_hit_rg(sdt, a1, a2, b1, b2);
        
#ifdef print_sv
            printf("Segment: [%c, %2lu,%2lu,%2lu,%2lu,%2lu]\n", \
                        sdt->dir==0?'+':'-', \
                        (unsigned long)sdt->rd_l_rg, (unsigned long)sdt->rd_r_rg, \
                        (unsigned long)sdt->rf_l_rg - chr_end_n[(sdt)->chr_name_i-1] + 1, (unsigned long)sdt->rf_r_rg - chr_end_n[(sdt)->chr_name_i-1] + 1, \
                        (unsigned long)(sdt->rd_r_rg-(sdt->rd_l_rg)));
#endif

                //if (have_B == 1)    continue;
            sdt->is_valid = 2;
        }

        int last_rst_n = sd_hit_rst_n[z], vk;
        for (int o=0, ti; o<sd_bt_n[z]; ++o)
        {
            ti = sd_bt_path[z][o];
            sd_hit_t *sdt = sd_hit_a[z]+ti;
            if (sdt->is_valid == 0)   continue;
            uint32_t a1, a2, b1, b2;
            get_sd_hit_rg(sdt, &a1, &a2, &b1, &b2);

            //filter bad block
            //if (MN((a2-a1), (b2-b1)) < sv_lim * 2)  continue;

            (sd_hit_rst[z]+sd_hit_rst_n[z])->dir = z;
            (sd_hit_rst[z]+sd_hit_rst_n[z])->chr_name_i = sdt->chr_name_i;
            (sd_hit_rst[z]+sd_hit_rst_n[z])->is_valid = 1;
            update_sd_hit_rg(sd_hit_rst[z]+(sd_hit_rst_n[z]++), a1, a2, b1, b2);
        }

        vk = sd_hit_rst_n[z] - last_rst_n;

#ifndef both_strand
        if (0)
#endif
        if (z == 1 && vk > 0)
        {
            for (int o=0; o<vk; ++o)
            {
                sd_hit_t *sdt1 = sd_hit_rst[z] + last_rst_n + o;
                sd_hit_t *sdt2 = sd_hit_rst[z] + last_rst_n + (vk-1-o);
                int chr_n1, chr_n2;
                uint32_t a1, a2, b1, b2, c1, c2, d1, d2;
                chr_n1 = sdt1->chr_name_i;
                chr_n2 = sdt2->chr_name_i;
                get_sd_hit_rg(sdt1, &a1, &a2, &b1, &b2);
                get_sd_hit_rg(sdt2, &c1, &c2, &d1, &d2);

                if (o > vk-1-o) break;
                update_sd_hit_rg(sdt1, seq_l-1-c2, seq_l-1-c1, d1, d2);
                update_sd_hit_rg(sdt2, seq_l-1-a2, seq_l-1-a1, b1, b2);
                sdt1->chr_name_i = chr_n2;
                sdt2->chr_name_i = chr_n1;
            }
        }
    } 
}


int single_core_aln(uint32_t rd_i, rd_handle_dt * t_rd_d)
{
    //load all space for aln
    uint8_t *read_seq, *read_seqr;
    read_seq = t_rd_d->read_seq;
    read_seqr = t_rd_d->read_seqr;
    sd_hit_t **sd_hit_a = t_rd_d->sd_hit_a;

    sd_hit_t **sd_hit_rst = t_rd_d->sd_hit_rst;
    int *sd_hit_rst_n = t_rd_d->sd_hit_rst_n;
    //---------
    //
    int **sd_bt_path = t_rd_d->sd_bt_path;
    int *sd_bt_n = t_rd_d->sd_bt_n;
    
    int32_t sd_hit_n = 0;

    //load taget read sequence
	char *t_seq = _rd_seq[rd_i].read_seq;  //read 
	char *t_seq_r = _rd_seq[rd_i].read_seq_r; //reverse read
	char *t_name = _rd_seq[rd_i].name;;

	uint32_t t_seq_l = _rd_seq[rd_i].read_length;
    
    //final storage
    rst_ent_dt *t_rst_rec = _rst_info[rd_i].rec;
    int32_t *t_rst_rec_n = &(_rst_info[rd_i].rec_n);
    *t_rst_rec_n = 0;

    uint32_t seq_l = t_seq_l, vl_seq_l = 0;
    trans_char_u8(t_seq, seq_l, read_seq, MAX_READ_LEN);
    trans_char_u8(t_seq_r, seq_l, read_seqr, MAX_READ_LEN);
    
	vl_seq_l = get_seed_vl_len(read_seq, seq_l, _seed_k);
    if (vl_seq_l <= _seed_k || vl_seq_l <= _sdp_k)
    {
        *t_rst_rec_n = 0;
        return 0;
    }
    if(seq_l < 100)
    {
        *t_rst_rec_n = 0;
        return 0;
    }

	// set variant for deBGA
    //default buffer_h_l is 14, buffer_k_l = 22 - 14 = 8
    //therefore _seed_k should > buffer_h_l which is 14
    uint8_t buffer_h_l = BUFFER_HASH_L, buffer_k_l = _db_k - buffer_h_l, buffer_k_r=_seed_k-buffer_h_l;
    uint8_t seed_h_off = (_seed_k-buffer_h_l) << 1;

    uint64_t _sd_mask = bits_right_1[_seed_k];
    uint64_t _sd_k_mask = bits_right_1[buffer_k_r];
    uint8_t *rd_seq;
    sd_bt_n[0] = sd_bt_n[1] = 0;
    sd_hit_rst_n[0] = sd_hit_rst_n[1] = 0;
    int dp_max[2] = {-1, -1};
    //int if_run_rev = 0;
    int if_run_rev = 1;

    for (int z=0; z<(1+if_run_rev); ++z)
    {
#ifdef print_seed
        if (z==1) printf("Reverse-------------------------------------------------------------------\n");
#endif
        rd_seq = (z==0) ? read_seq : read_seqr;
        sd_hit_n = 0;
 
        uint64_t kmer_s, kmer_h, kmer_k;

#ifdef debug_exd
        if (0)
#endif
        for (uint32_t i=0; i<seq_l; ++i)
        {
            if (rd_seq[i] == 4) continue;
            uint32_t start = i, rdb = 0, rde = 0, rd_VL = 0;
            while (i<seq_l && rd_seq[i++] != 4);
            if (i-start < _seed_k) continue;
            rdb = start; rde = i; rd_VL = i - start;
            //printf("rdb = %ld, rde = %ld, i = %ld, seq_l = %ld\n", rdb, rde, i, seq_l);

            //extract _seed_k and check if exist in deBGA
            //seeding strategy is divided into 3 part
            //[bg    100bp] [  ...   ]  [ed    100bp]   #0 base
            //rds_s_t == 1 : fix check time else fix jump len

            //rd_VL = 50200;
            uint32_t rds_r[3], rds_l[3], rds_t[3], rds_j[3], rds_b[3], rds_e[3], rds_n[3];
            //int32_t rds_o_n1 = 100, rds_o_c1 = 10, rds_o_K = 100, rds_o_J = 100;
            int32_t rds_o_n1 = 100, rds_o_c1 = 10, rds_o_K = 100, rds_o_J = _seed_step; //rds_o_c1 = 20 previous
            int rds_s_t = 2; 

            //if run
            rds_r[0] = rds_r[2] = (rd_VL >= 3 * rds_o_n1) ? 1 : 0;
            rds_r[1] = 1; 

            //part len
            rds_l[0] = rds_l[2] = (rds_r[0] == 1) ? rds_o_n1 : 0;
            rds_l[1] = rd_VL - rds_l[0] - rds_l[2];

            //check #
            rds_t[0] = rds_t[2] = rds_o_c1;
            rds_t[1] = (rds_s_t == 1) ? rds_o_K : (rds_l[1] / rds_o_J);

            //jump #
            rds_j[0] = rds_j[2] = MX(1, (rds_l[0] / rds_t[0]));
            rds_j[1] = (rds_s_t == 1) ? (rds_l[1] / rds_t[1]) : rds_o_J;
            rds_j[1] = MX(rds_j[1], 5);

            //maximum hit #
            rds_n[0] = rds_n[2] = 3;
            rds_n[1] = rds_t[1];

            //window begin & end
            rds_b[0] = 0;        rds_b[1] = rds_l[0];            rds_b[2] = rds_b[1] + rds_l[1];
            for (int _o=0; _o<3; ++_o)
                rds_b[_o] += (rds_j[_o] > _seed_k) ? (rds_j[_o] - _seed_k)/2 : 0;
            rds_e[0] = rds_l[0]; rds_e[1] = rds_e[0] + rds_l[1]; rds_e[2] = rds_e[1] + rds_l[2]; 
            
#ifdef  print_seed
            printf("rds_n: %u-%u-%u\n", rds_n[0], rds_n[1], rds_n[2]);
            printf("rds_j: %u-%u-%u\n", rds_j[0], rds_j[1], rds_j[2]);
#endif
            //generate seed and save
            uint32_t _ck_i = 20, r_b_v = 0;  //seeding from 20bp to filter polyA in the head
            for (int _r = 0; _r < 3; ++_r)
            {
                if (rds_r[_r] == 0) continue;
                int _r_cnt = 0;
                uint32_t j;
                //_ck_i = rds_b[_r];
#ifdef print_seed
                printf("_r = %d, rds_n = %u, rds_e = %u\n", _r, rds_n[_r], rds_e[_r]);
#endif
                //r_b_v = _ck_i;
                while (_r_cnt < rds_n[_r] && (_ck_i+_seed_k < rds_e[_r]))
                {
                    j = _ck_i + rdb;
                    _ck_i += rds_j[_r];
#ifdef print_seed
                    //printf("j = %d, _seed_k = %d, r_b_v = %d\n", j, _seed_k, r_b_v);
#endif
                    if (j + _seed_k < r_b_v)
                    {
                        continue;
                    }
                    if (sd_hit_n >= _seed_hit_n_max) break;
                    kmer_s = get_first_kmer_bits_u8(rd_seq + j, _seed_k);
                    kmer_h = kmer_s >> seed_h_off;
                    kmer_k = kmer_s & _sd_k_mask; 

                    //check if have hit
                    //get kmer offsets in uniq_seq 
                    int64_t kmer_uqs_off_bound[2] = {-1, -1};
                    uint64_t uqp_off, uqp_id, uqp_b, uqp_len;
                    uint64_t uqp_e, uqp_l_off, uqp_r_off, *uqp_seq;  
                    int flag_hit = 0;
		    	    //int tmp = multi_binsearch_offset64(kmer_k, buffer_kmer_g, buffer_hash_g[kmer_h + 1] - buffer_hash_g[kmer_h], buffer_hash_g[kmer_h], kmer_uqs_off_bound, (_db_k-_seed_k));
                    int tmp = binsearch_range(kmer_k, buffer_kmer_g+buffer_hash_g[kmer_h],\
                        buffer_hash_g[kmer_h+1]-buffer_hash_g[kmer_h],\
                        kmer_uqs_off_bound, (_db_k-_seed_k)<<1); 
                    if (tmp == -1) continue; 
                    kmer_uqs_off_bound[0] += buffer_hash_g[kmer_h];
		    	    kmer_uqs_off_bound[1] += buffer_hash_g[kmer_h];

                    if (kmer_uqs_off_bound[1] - kmer_uqs_off_bound[0] + 1> _read_kmer_match_mx)
                    {
                        continue;
                    }
                    
                    int uqp_ref_n = 0;
                    for (uint64_t o = kmer_uqs_off_bound[0]; o <= kmer_uqs_off_bound[1]; ++o)
                    {
                        uqp_off = buffer_off_g[o];

                        uqp_id = hash_get_uqp_id(uqp_off, map_uniq_id, map_uniq_id_n, buffer_seqf, result_seqf);

                        uqp_ref_n = buffer_pp[uqp_id+1] - buffer_pp[uqp_id];
                        if (uqp_ref_n  > _read_kmer_ref_mx || uqp_id + 1 > result_pp)
		    	    	{	
                            continue;
		    	    	}
                        uqp_b = buffer_seqf[uqp_id];
                        uqp_e = buffer_seqf[uqp_id+1]; 
                        uqp_len = uqp_e-uqp_b;

                        uqp_l_off = uqp_off;
                        uqp_r_off = uqp_off+_seed_k-1;
                        uqp_seq = buffer_seq;
                        
#ifdef print_seed
                        printf("before change to reference: uqp_l_off = %ld, uqp_b = %ld, uqp_r_off = %ld, uqp_e = %ld, uqp_ref_n|mx= %d|%d\n", uqp_l_off,uqp_b, uqp_r_off, uqp_e, uqp_ref_n, _read_kmer_ref_mx);
#endif

                        if (uqp_ref_n <= 1)
                        {
                            uqp_l_off = buffer_p[buffer_pp[uqp_id]+0]+(uqp_off - uqp_b) - 1;
                            uqp_r_off = uqp_l_off+_seed_k-1;
                            uqp_off = uqp_l_off;
                            uqp_b = (uqp_l_off > _ref_lv_bnd_len)? uqp_l_off - _ref_lv_bnd_len : 1;
                            uqp_e = uqp_r_off + _ref_lv_bnd_len + 1;
                            uqp_len = uqp_e - uqp_b;
                            uqp_seq = buffer_ref_seq;
#ifdef print_seed
                            printf("after change to reference: uqp_l_off = %ld, uqp_b = %ld, uqp_r_off = %ld, uqp_e = %ld\n", uqp_l_off,uqp_b, uqp_r_off, uqp_e); 
#endif
                        }

                        int read_b = r_b_v+1;
                        int read_e = rde,  read_l_off = j, read_r_off = j+_seed_k-1;
                        read_b = (read_b > 1)? read_b :0;
                        int tk_bnd_l=MN(uqp_l_off-uqp_b, read_l_off-read_b), tk_bnd_r=MN(uqp_e-uqp_r_off-1, read_e-read_r_off-1);
 
#ifdef print_seed
                        /*
                        printf("tk_bnd_l = %d, uqp_l_off = %ld, uqp_b = %ld, read_l_off = %ld, read_b = %ld\n", \
                                tk_bnd_l, uqp_l_off, uqp_b, read_l_off, read_b);
                        printf("tk_bnd_r = %d, uqp_e = %ld, uqp_r_off = %ld,  read_e = %ld, read_r_off = %ld\n", \
                                tk_bnd_r, uqp_e, uqp_r_off, read_e, read_r_off);
                        */
#endif
                        //handle each hit
                        int read_mem_len = 0, rd_block_off = 0, rf_block_off = 0; 
                        int lv_ex_l = 0, lv_ex_r = 0, tl_i=0, tl_e=0, tr_i=0, tr_e=0;

                        //mem left
#ifndef only_right_ext
                        if (uqp_off>uqp_b && j>read_b)
                        {
                            --uqp_l_off; --read_l_off;
                            while (uqp_l_off >= uqp_b && read_l_off >= read_b && (GET_BITS_A(uqp_seq, uqp_l_off) == rd_seq [read_l_off]))
                            {--uqp_l_off; --read_l_off;}
                            ++uqp_l_off; ++read_l_off; 
                        }
#endif
                        if (uqp_off+_seed_k < uqp_e && j+_seed_k < read_e)
                        {
                            ++uqp_r_off; ++read_r_off;
                            while (uqp_r_off < uqp_e && read_r_off < read_e && (GET_BITS_A(uqp_seq, uqp_r_off) == rd_seq [read_r_off]))
                            {++uqp_r_off; ++read_r_off;}
                            --uqp_r_off; --read_r_off;
                        }
                        int tk_mem_l=uqp_off-uqp_l_off, tk_mem_r=uqp_r_off-(uqp_off+_seed_k-1);

#ifdef print_seed
                        printf("after MEM                  :uqp_l_off = %ld, uqp_b = %ld, uqp_r_off = %ld, uqp_e = %ld\n", uqp_l_off,uqp_b, uqp_r_off, uqp_e); 
                        printf("First LV---------------\n");
#endif

                        read_mem_len = uqp_r_off - uqp_l_off + 1; 
                        rd_block_off = j - read_l_off;
                        rf_block_off = uqp_off - uqp_l_off;
                        
                        if (uqp_r_off < (uqp_e-1) && read_r_off < (read_e-1))
                        {
                            lv_ex_r = LandauVishkin_64_log(rd_seq , read_e-read_r_off-1, read_r_off+1, uqp_seq, uqp_e-uqp_r_off-1, uqp_r_off+1, _max_lv_e, 1, &tr_e, &tr_i); 
                            ++lv_ex_r;
                        }
#ifndef only_right_ext              
                        if (uqp_l_off > uqp_b && read_l_off > read_b)
                        {
                            lv_ex_l = LandauVishkin_64_log(rd_seq , read_l_off-read_b, read_l_off-1, uqp_seq, uqp_l_off-uqp_b, uqp_l_off-1, _max_lv_e, -1, &tl_e, &tl_i); 
                            ++lv_ex_l;
                        }
#endif
                        int sig = (uqp_ref_n > 1)? 1 : 0;
                        uint64_t ref_b, uqp_l_off_t=0, uqp_r_off_t, now_pos_l, now_begin_l, now_pos_r, now_end_r; 
                        int chr_name_i, tmp_rbv = 0;
#ifdef print_seed
                        printf("uqp_l_off = %ld, uqp_b = %ld, refpb = %ld, uqp_len = %ld, uqp_r_off = %ld uqp_e = %ld, refpe = %ld\n", uqp_l_off, uqp_b, uqp_l_off-(lv_ex_l + tl_i), uqp_len, uqp_r_off, uqp_e-1, uqp_r_off+lv_ex_r+tr_i);
                        printf("second LV-------------------\n");
#endif
                        for(int t = 0; t < uqp_ref_n; ++t)
                        {
                            if(uqp_ref_n > 1)
                                ref_b = buffer_p[buffer_pp[uqp_id]+t]+(uqp_off - uqp_b) - 1;
                            else
                                ref_b = uqp_off;
                            chr_name_i = chromosome_judge(ref_b);
                            
                            if(sig && (uqp_b == uqp_l_off-lv_ex_l-tl_i || uqp_e-1 == uqp_r_off+lv_ex_r+tr_i))
                            {
                                uqp_seq = buffer_ref_seq;
                                uqp_l_off_t = ref_b - tk_mem_l;
                                uqp_r_off_t = ref_b + _seed_k-1+tk_mem_r;
                                now_pos_l = uqp_l_off_t - lv_ex_l - tl_i;
                                now_begin_l = (now_pos_l > _ref_lv_bnd_len)? now_pos_l - _ref_lv_bnd_len : 1;
                                now_pos_r = uqp_r_off_t + lv_ex_r + tr_i;
                                now_end_r = now_pos_r + _ref_lv_bnd_len;
                            }

                            int tl_e_t = tl_e, tl_i_t = tl_i, lv_ex_l_t = lv_ex_l;
                            if(sig && uqp_b == uqp_l_off-lv_ex_l-tl_i && read_b < read_l_off-lv_ex_l) // left ext in ref for each copy
                            {
                                int l_e = 0, l_i = 0;
                                lv_ex_l_t = LandauVishkin_64_log(rd_seq, read_l_off-lv_ex_l-read_b, read_l_off-lv_ex_l-1, uqp_seq, now_pos_l - now_begin_l, now_pos_l-1,_max_lv_e - tl_e, -1, &l_e, &l_i);
                                lv_ex_l_t++;
                                lv_ex_l_t += lv_ex_l; tl_e_t += l_e; tl_i_t += l_i;
                            }
                            int tr_e_t = tr_e, tr_i_t = tr_i, lv_ex_r_t = lv_ex_r;
                            if(sig && uqp_e-1 == uqp_r_off+lv_ex_r+tr_i && read_e > read_r_off+lv_ex_r) //right ext in ref for each copy
                            {
                                int r_e = 0, r_i = 0;
                                lv_ex_r_t = LandauVishkin_64_log(rd_seq, read_e-read_r_off-lv_ex_r-1, read_r_off+lv_ex_r+1,uqp_seq, now_end_r - now_pos_r-1, now_pos_r+1, _max_lv_e - tr_e, 1, &r_e, &r_i);
                                lv_ex_r_t++; 
                                lv_ex_r_t += lv_ex_r; tr_e_t += r_e; tr_i_t += r_i;
                            }
                            
                            int tk_lv_l = MN(lv_ex_l_t, lv_ex_l_t+tl_i_t), tk_lv_r = MN(lv_ex_r_t, lv_ex_r_t+tr_i_t);
                            int tk_lv_li = G_ABS(tl_i_t), tk_lv_ri = G_ABS(tr_i_t);

                            double x_ref_l = 3000000000, x_tot_len, x_err_len, x_mat_len, x_Score;
                            double x_rnd_p = 0.25, x_hit_p = _er_ins +  _er_del + _er_mis;
                            x_ref_l = reference_len-START_POS_REF;
                            x_tot_len = read_mem_len + tk_lv_l + tk_lv_r;
                            x_err_len = tl_e_t + tr_e_t;
                            x_mat_len = x_tot_len - x_err_len + (tk_lv_li + tk_lv_ri);
                            x_Score = x_mat_len * -10 * (log10(x_rnd_p/(1-x_hit_p))) + \
                                      x_err_len * -10 * (log10((1-x_rnd_p)/x_hit_p)) + \
                                      -10 * log10(x_ref_l);
#ifndef all_seed 
#ifdef print_seed
                            printf("x_tot_len = %f, x_score = %f\n", x_tot_len, x_Score);
#endif
                            //if (x_tot_len < 30) continue;
                            //if (x_Score < 30) continue;
                            if (x_tot_len < _seed_k+2*_max_lv_e+7) continue;
                            if (x_Score < _seed_k+2*_max_lv_e+7) continue;
#endif
                            flag_hit = 1;
                            set_sd_hit_nd(sd_hit_a[z]+sd_hit_n, sd_hit_n, seq_l, read_l_off, ref_b-rf_block_off, chr_name_i, read_mem_len, lv_ex_l_t, tl_i_t, tl_e_t, lv_ex_r_t, tr_i_t, tr_e_t, z);
#ifdef  print_seed
                            printf("nd_id%2d, chr%s, %c, rd_p%4d| uid%2zu, %6zu, %4zu, %4zu, %2zu| %d, %4zu| ", \
                            (int)sd_hit_n, chr_names[chr_name_i], z==0?'+':'-', j, \
                            o - kmer_uqs_off_bound[0], uqp_id, uqp_l_off - uqp_b, uqp_r_off - uqp_b, uqp_r_off - uqp_l_off + 1, \
                            t, uqp_l_off_t - chr_end_n[chr_name_i-1] + 1);
                            printf("%3u, %4u, %4u, %6u, %6u| ", \
                            (sd_hit_a[z]+sd_hit_n)->sd_len, \
                            (sd_hit_a[z]+sd_hit_n)->rd_l, \
                            (sd_hit_a[z]+sd_hit_n)->rd_r, \
                            (sd_hit_a[z]+sd_hit_n)->rf_l, \
                            (sd_hit_a[z]+sd_hit_n)->rf_r);
                            printf("L:%d, %d, %d, %d R:%d, %d, %d, %d MAPQ:%.1f\n", tk_mem_l, tk_lv_l, tk_bnd_l, tl_e_t, \
                                            tk_mem_r, tk_lv_r, tk_bnd_r, tr_e_t, x_Score);
#endif

                            if (tmp_rbv < (sd_hit_a[z]+sd_hit_n)->rd_r) 
                            {
                                tmp_rbv = (sd_hit_a[z]+sd_hit_n)->rd_r;
                            }
                            sd_hit_n += 1;
                        }
                        if (tmp_rbv > r_b_v) r_b_v = tmp_rbv;
                    }
                    if (flag_hit != 0)
                    {
                        _r_cnt += 1;
                    }
                }
#ifdef  print_seed
                printf("_r = %d, _r_cnt = %d\n", _r, _r_cnt);
#endif
            }
        }

        if (sd_hit_n == 0) continue;
#ifdef print_result
        printf("before run sdp, sd_hit_n = %d\n", sd_hit_n);
        for (int o=0; o<sd_hit_n; ++o)
        {
            printf("%d, %d, %3u, %4u, %4u, %6u, %6u \n", o, (sd_hit_a[z]+o)->chr_name_i,\
                            (sd_hit_a[z]+o)->sd_len, \
                            (sd_hit_a[z]+o)->rd_l, \
                            (sd_hit_a[z]+o)->rd_r, \
                            (sd_hit_a[z]+o)->rf_l, \
                            (sd_hit_a[z]+o)->rf_r);
        }

#endif
       
        sd_bt_n[z] = 0;
        dp_max[z] = run_seed_sdp(sd_hit_a[z], sd_hit_n, sd_bt_path[z], &sd_bt_n[z], _er_ins, _er_del); 
        
#ifdef print_result
        printf("\nsd_bt_n = %d\n", sd_bt_n[z]);
        for(int o = 0; o <sd_bt_n[z]; ++o)
        {
            printf("%d(%d)--", sd_bt_path[z][o], sd_hit_a[z][sd_bt_path[z][o]].path_range);
        }
        printf("\n");
#endif
    }
     
    int which = -1, z = -1; 
    
    //if(dp_max[0] > dp_max[1] && dp_max[0] > _min_chain_score && dp_max[1] / (float)dp_max[0] < _secondary_ratio)
    //    which = 0, z = 0;
    //else if(dp_max[1] > dp_max[0] && dp_max[1] > _min_chain_score && dp_max[0] / (float)dp_max[1] < _secondary_ratio)
    //    which = 1, z = 1;
    //else if(dp_max[0] > _min_chain_score && dp_max[1] > _min_chain_score)
    //    which = 2, z = (dp_max[0] > dp_max[1])? 0 : 1;

    if(dp_max[0] > dp_max[1] && dp_max[1] < _min_chain_score)
        which = 0, z = 0;
    else if(dp_max[1] > dp_max[0] && dp_max[0] < _min_chain_score)
        which = 1, z = 1;
    else if(dp_max[0] > _min_chain_score && dp_max[1] > _min_chain_score)
        which = 2, z = (dp_max[0] > dp_max[1])? 0 : 1;

#ifdef print_result
    fprintf(stderr, "dp_max: %d-%d, which = %d, z = %d\n", dp_max[0], dp_max[1], which, z);
#endif
    
    
    if(which == -1)
    {
        *t_rst_rec_n = 0;
        return 0;
    }

    int Sig = 0, t_rec_n = 0; 
    if(Sig == 1)
    {
        return 0;
    }
    else
    {  
        *t_rst_rec_n = 0;
        do_sv_detection(rd_i, t_rd_d, t_rst_rec, t_rst_rec_n, sd_hit_rst, sd_hit_rst_n, sd_hit_a, sd_bt_path, sd_bt_n);
        
#ifndef best_strand 
        //best strand - > alt strand
        for (int ki=0; ki<sd_hit_rst_n[1-z]; ++ki)
        {
            sd_hit_t *sdt1 = sd_hit_rst[1-z] + (ki);
            uint32_t a1, a2, b1, b2, c1, c2, d1, d2;
            get_sd_hit_rg(sdt1, &a1, &a2, &b1, &b2);

            //fprintf(stderr, "A. chr=%d, a1 = %u, a2 = %u, b1 = %u, b2 = %u\n", sdt1->chr_name_i, a1, a2, b1, b2);

            int check_f = 1;
            for (int kj=0; kj<sd_hit_rst_n[z]; ++kj)
            {
                sd_hit_t *sdt2 = sd_hit_rst[z] + (kj);
                get_sd_hit_rg(sdt2, &c1, &c2, &d1, &d2);

                if (sdt2->is_valid == 0)
                    continue;
                //fprintf(stderr, "B. chr=%d, c1 = %u, c2 = %u, d1 = %u, d2 = %u\n", sdt2->chr_name_i, c1, c2, d1, d2);
                //fprintf(stderr, "ratio = %f\n", overlap_l(c1,c2,a1,a2));
                if(overlap_l(c1, c2, a1, a2) > 0.6)
                {
                    check_f = 0;
                    break;
                }
            }
            sdt1->is_valid = check_f;
        }
        // alt strand - > best strand
        for (int ki=0; ki<sd_hit_rst_n[z]; ++ki)
        {
            sd_hit_t *sdt1 = sd_hit_rst[z] + (ki);
            uint32_t a1, a2, b1, b2, c1, c2, d1, d2;
            get_sd_hit_rg(sdt1, &a1, &a2, &b1, &b2);

            //fprintf(stderr, "A. chr=%d, a1 = %u, a2 = %u, b1 = %u, b2 = %u\n", sdt1->chr_name_i, a1, a2, b1, b2);

            int check_f = 1;
            for (int kj=0; kj<sd_hit_rst_n[1-z]; ++kj)
            {
                sd_hit_t *sdt2 = sd_hit_rst[1-z] + (kj);
                get_sd_hit_rg(sdt2, &c1, &c2, &d1, &d2);

                if (sdt2->is_valid == 0)
                    continue;

                //fprintf(stderr, "B. chr=%d, c1 = %u, c2 = %u, d1 = %u, d2 = %u\n", sdt2->chr_name_i, c1, c2, d1, d2);
                //fprintf(stderr, "ratio = %f\n", overlap_l(c1,c2,a1,a2));
                if(overlap_l(c1, c2, a1, a2) > 0.8)
                {
                    check_f = 0;
                    break;
                }
            }
            sdt1->is_valid = check_f;
        }
#else 
        // only output skeletons in one contig
        for (int ki=0; ki<sd_hit_rst_n[1-z]; ++ki)
        {
            sd_hit_t *sdt1 = sd_hit_rst[1-z] + (ki);
            
            sdt1->is_valid = 0;
        }
#endif
        // remove overlaped segment in the same strand
        for (int ki=0; ki<sd_hit_rst_n[z]; ++ki)
        {
            sd_hit_t *sdt1 = sd_hit_rst[z] + (ki);
            if (sdt1->is_valid == 0) continue;

            uint32_t a1, a2, b1, b2, c1, c2, d1, d2;
            get_sd_hit_rg(sdt1, &a1, &a2, &b1, &b2);

            //fprintf(stderr, "C. chr=%d, a1 = %u, a2 = %u, b1 = %u, b2 = %u\n", sdt1->chr_name_i, a1, a2, b1, b2);

            int check_f = 1;
            for (int kj=ki+1; kj<sd_hit_rst_n[z]; ++kj)
            {
                sd_hit_t *sdt2 = sd_hit_rst[z] + (kj);
                get_sd_hit_rg(sdt2, &c1, &c2, &d1, &d2);

                if (sdt2->is_valid == 0)
                    continue;

                //fprintf(stderr, "D. chr=%d, c1 = %u, c2 = %u, d1 = %u, d2 = %u\n", sdt2->chr_name_i, c1, c2, d1, d2);
                //fprintf(stderr, "ratio = %f\n", overlap_l(c1,c2,a1,a2));
                if(overlap_l(c1, c2, a1, a2) > 0.8)
                {
                    check_f = 0;
                    break;
                }
            }
            sdt1->is_valid = check_f;
        }
        //
        for (int ki=0; ki<sd_hit_rst_n[1-z]; ++ki)
        {
            sd_hit_t *sdt1 = sd_hit_rst[1-z] + (ki);
            if (sdt1->is_valid == 0) continue;

            uint32_t a1, a2, b1, b2, c1, c2, d1, d2;
            get_sd_hit_rg(sdt1, &a1, &a2, &b1, &b2);

            //fprintf(stderr, "E. chr=%d, a1 = %u, a2 = %u, b1 = %u, b2 = %u\n", sdt1->chr_name_i, a1, a2, b1, b2);

            int check_f = 1;
            for (int kj=ki+1; kj<sd_hit_rst_n[1-z]; ++kj)
            {
                sd_hit_t *sdt2 = sd_hit_rst[1-z] + (kj);
                get_sd_hit_rg(sdt2, &c1, &c2, &d1, &d2);

                if (sdt2->is_valid == 0)
                    continue;

                //fprintf(stderr, "F. chr=%d, c1 = %u, c2 = %u, d1 = %u, d2 = %u\n", sdt2->chr_name_i, c1, c2, d1, d2);
                if(overlap_l(c1, c2, a1, a2) > 0.8)
                {
                    check_f = 0;
                    break;
                }
            }
            sdt1->is_valid = check_f;
        }
        
        // output segments
        int ki = 0, kj = 0;
        while (ki < sd_hit_rst_n[0] && kj < sd_hit_rst_n[1])
        {
            uint32_t a1, a2, b1, b2, c1, c2, d1, d2;
            sd_hit_t *sdt1 = sd_hit_rst[0] + (ki);
            sd_hit_t *sdt2 = sd_hit_rst[1] + (kj);
            if (sdt1->is_valid == 0)
            {
                ki++;
                continue;
            }
            if (sdt2->is_valid == 0)
            {
                kj++;
                continue;
            }
            get_sd_hit_rg(sdt1, &a1, &a2, &b1, &b2);
            get_sd_hit_rg(sdt2, &c1, &c2, &d1, &d2);
            
            sd_hit_t *sdt; 
            if (a1 > c1)
            {    
                sdt = sdt2;
                kj++;
            }
            else
            {
                sdt = sdt1;
                ki++;
            }
            //if(sdt->rd_r_rg - sdt->rd_l_rg < 500) continue;
            t_rst_rec[t_rec_n].chr_name_i = sdt->chr_name_i;
            t_rst_rec[t_rec_n].dir = sdt->dir;
            t_rst_rec[t_rec_n].rd_l = sdt->rd_l_rg;
            t_rst_rec[t_rec_n].rf_l = sdt->rf_l_rg;
            t_rst_rec[t_rec_n].rd_r = sdt->rd_r_rg;
            t_rst_rec[t_rec_n].rf_r = sdt->rf_r_rg;
            t_rec_n++;
        }
        while (ki < sd_hit_rst_n[0])
        {
            uint32_t a1, a2, b1, b2;
            sd_hit_t *sdt = sd_hit_rst[0] + (ki++);
            get_sd_hit_rg(sdt, &a1, &a2, &b1, &b2);

            //if (sdt->is_valid == 0 || (sdt->rd_r_rg - sdt->rd_l_rg < 500))
            if (sdt->is_valid == 0)
                continue;
            t_rst_rec[t_rec_n].chr_name_i = sdt->chr_name_i;
            t_rst_rec[t_rec_n].dir = sdt->dir;
            t_rst_rec[t_rec_n].rd_l = sdt->rd_l_rg;
            t_rst_rec[t_rec_n].rf_l = sdt->rf_l_rg;
            t_rst_rec[t_rec_n].rd_r = sdt->rd_r_rg;
            t_rst_rec[t_rec_n].rf_r = sdt->rf_r_rg;
            t_rec_n++;
        }
        while (kj < sd_hit_rst_n[1])
        {
            uint32_t a1, a2, b1, b2;
            sd_hit_t *sdt = sd_hit_rst[1] + (kj++);
            get_sd_hit_rg(sdt, &a1, &a2, &b1, &b2);

            //if (sdt->is_valid == 0 || (sdt->rd_r_rg - sdt->rd_l_rg < 500))
            if (sdt->is_valid == 0 )
                continue;

            t_rst_rec[t_rec_n].chr_name_i = sdt->chr_name_i;
            t_rst_rec[t_rec_n].dir = sdt->dir;
            t_rst_rec[t_rec_n].rd_l = sdt->rd_l_rg;
            t_rst_rec[t_rec_n].rf_l = sdt->rf_l_rg;
            t_rst_rec[t_rec_n].rd_r = sdt->rd_r_rg;
            t_rst_rec[t_rec_n].rf_r = sdt->rf_r_rg;
            t_rec_n++;
        }
        *t_rst_rec_n = t_rec_n;
    }

    return 0;
}


void log_rst_rec(uint32_t rd_i, double MAX_MQ)
{
	char *t_name = _rd_seq[rd_i].name;
	uint32_t t_seq_l = _rd_seq[rd_i].read_length;

    rst_ent_dt *t_rst_rec = _rst_info[rd_i].rec;
    int32_t t_rec_n = _rst_info[rd_i].rec_n;

    
    int sum_l = 0;
    uint32_t ld, lr;
    int tflag = 0;

    if(t_rec_n == 0) return ;
    fprintf(fp_sv, ">%s [%d]\n", t_name, t_seq_l);
    for (int i=0, lg, lchr; i<t_rec_n; ++i)
    {
        rst_ent_dt *sdt = t_rst_rec + i;

        lchr = sdt->chr_name_i;
        lg = sdt->dir;
        ld = sdt->rd_r;
        lr = sdt->rf_r;
        sum_l += (sdt->rd_r-(sdt->rd_l));

        if (tflag == 1)
            fprintf(fp_sv,"|");
        tflag = 0;
        fprintf(fp_sv, "%s,%c,%lu,%lu,%lu,%lu,%lu", \
                chr_names[lchr], \
                sdt->dir==0?'+':'-', \
                (unsigned long)sdt->rd_l, (unsigned long)sdt->rd_r, \
                (unsigned long)sdt->rf_l+2-chr_end_n[lchr-1], (unsigned long)sdt->rf_r+2-chr_end_n[lchr-1], \
                (unsigned long)(sdt->rd_r+1-(sdt->rd_l)));
        tflag=1;

        if (i==t_rec_n-1) fprintf(fp_sv, "\n");

    }
}
static void *single_thread_core_aln(void *d_p)
{
	rd_handle_dt *d = (rd_handle_dt* )d_p;
    int _read_lines;

	while(1)
	{
		pthread_rwlock_wrlock(&rwlock);
		_read_lines = THREAD_READ_I++;
		if (_read_lines < _rd_seq_i && _read_lines % 500 == 0)
		{
            fprintf(stderr, "[STC] run rd #%d reads >%.5s  \r", _read_lines, _rd_seq[_read_lines].name);
		}
		pthread_rwlock_unlock(&rwlock);

		if (_read_lines < _rd_seq_i)
		{
			single_core_aln(_read_lines, d);
		}
		else break;
	}
	return 0;

}
