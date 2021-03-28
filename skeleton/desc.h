#ifndef _DESC_H 
#define _DESC_H

#define PACKAGE_NAME    "SKSV-skeleton"
#define PACKAGE_VERSION "1.0.2"
#define CONTACT "<ydliu@hit.edu.cn>"

#define INDEX_KMER 22
#define SEEDING_KMER 22
#define HASH_KMER 7
#define SV_LIM 30
#define BATCH_SIZE 100000
#define TOP_N 4 
#define RD_GAP 1000
#define RF_GAP 1000
#define SEC_RATIO 0.8
#define OVERLAP_RATIO 0.4
#define GAP_FILL_EXT 50
#define MAPQ 0 
#define MIN_CHAIN_SCORE 100
#define ERROR_MODEL 1
#define SOFT 500
#define BLOCK 30
#define EDIT_DIS 3
#define SEED_STEP 20
#define EXTEND_MX 5000 
#define GAP_OPEN 6 
#define GAP_OPEN2 44
#define GAP_EXT 2
#define GAP_EXT2 1
#define MATCH_SCORE 1 
#define MISMATCH_SCORE 4
#define BLOCKS 30
#define OUTPUT_PREFIX "sk.svseg"


//" PRId64 "
//#define GEN_DATA
//#define log_mem_lv
//#define only_right_ext

//#define print_path_o

//#define debug_wf

//#define print_rst_path

//#define print_seed
//#define lv_log_s
//#define print_hit_exd
//#define print_result
//#define print_sv
//#define print_greed
//#define merge_colinear


//#define print_gd_detail
//#define print_exd_detail

//#define show_inter

//#define sdp_filter_chr
//#define print_mg_info

//#define run_only_one

//#define gd_use_bnd
//#define ext_hit_nd
//#define debug_exd

#define refine_bnd
//#define best_strand
#define both_strand 


#define filter_strand
//#define filter_strand_dp_score
#define filter_same_second_hit

//#define print_mg_info
#define gd_edge_cnt

#define run_mem_lv
#define run_hit_sdp 
#define run_greed_exd 

//#define all_seed
//#define run_o_rht
//#define only_seed

//#define no_build_RHT

#endif
// 1 <= chr_name_i < chr_file_n
// chr_name_i = chromosome_judge(ref_p);
// chr name: chr_names[chr_name_i]
// chr len: chr_end_n[chr_name_i] - chr_end_n[chr_name_i-1]
