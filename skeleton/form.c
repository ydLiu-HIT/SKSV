#include "form.h"

#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>

void opt_init(opt_t *opt)
{
    strcpy(opt->sv_path, "sk.svseg");
    opt->db_k = 22;
    opt->seed_k = 22;
    opt->sdp_k = 7;
    opt->rht_k = 7;
    
    opt->block_s = 30;

    opt->batch_size = 100000;
    opt->the_sv_lim = 50;
    opt->ref_lv_bnd_len = 100; //ccs, could be large or small depend on the speed

    opt->max_path_N = 3;
    opt->rd_mx_gap = 2000;
    opt->rf_mx_gap = 2000;
    opt->secondary_ratio = 0.9;
    opt->overlap_ratio = 0.4;
    opt->min_chain_score = 100;

    opt->thread_n = 4;
    opt->seed_step = 20;

    opt->max_lv_e = 0;
	opt->read_kmer_match_mx = 30;
	opt->read_kmer_ref_mx = 5;
	opt->ref_kmer_hit_mx = 20;

    //KHT Rtable max len, max read length
    opt->rht_n_max = 1000000;
    //MAX hit storing table sizeh
    opt->seed_hit_n_max = 25000; //could be smaller, 5000000 previous
    opt->rst_seed_hit_n_max = 2500; //could be smaller, 10000 previous

    opt->data_type = 1;
    //pacbio clr
    opt->error_rate = 0.15;
    opt->er_ins = 0.0592;
    opt->er_del = 0.0301;
    opt->er_mis = 0.013;

}
/*
void opt_log(opt_t *opt)
{
    printf("index path: %s\n", opt->index_path);
    printf("read path: %s\n", opt->read_path);
    printf("sv_path : %s\n", opt->sv_path);
    printf("db_k: %d\n", opt->db_k);
    printf("seed_k: %d\n", opt->seed_k);
    printf("sdp_k: %d\n", opt->sdp_k);
    printf("seeding step: %d\n", opt->seed_step);
    printf("batch_size = %d, top_n = %d, rd_gap = %d, rf_gap = %d, sec_ratio = %f\n", opt->batch_size, opt->max_path_N, opt->rd_mx_gap, opt->rf_mx_gap, opt->secondary_ratio);
}
*/


int help_usage()
{
    fprintf(stderr, "\n"); 
    fprintf(stderr, "Program:   %s, skeleton-based read alignment for variant calling\n", PACKAGE_NAME); 
    fprintf(stderr, "Version:   %s\n", PACKAGE_VERSION); 
    fprintf(stderr, "Contact:   %s\n", CONTACT); 
    fprintf(stderr, "Usage:     %s <command> [Options]\n\n", PACKAGE_NAME); 
    fprintf(stderr, "Command: \n");
	fprintf(stderr, "    index      index reference sequence\n");
	fprintf(stderr, "    aln        align long raw sequence to reference by index\n");
	fprintf(stderr, "\n");

    fprintf(stderr, "Usage:	%s index <ref.fa> <index_route>\n", PACKAGE_NAME);
	fprintf(stderr, "\n\n");
	fprintf(stderr, "Usage:	%s aln [options] <index_route> <read.fa/fq>\n\n", PACKAGE_NAME);

    fprintf(stderr, "    <index_route>                 The path of RdBG index.\n");
    fprintf(stderr, "    <read.fq/fa>                  The input reads in fasta or fastq format.\n\n");
    
    fprintf(stderr, "Algorithm options:\n\n");
	fprintf(stderr, "    -t --thread           [INT]    Number of threads. [4]\n");
    fprintf(stderr, "    -k --index_kmer       [INT]    RdBG-index kmer length. [%u]\n", INDEX_KMER);
	fprintf(stderr, "    -s --seeding_kmer     [INT]    K-mer length of seeding process (no long than RdBG-index kmer). [%u]\n", SEEDING_KMER);
	fprintf(stderr, "    -d --hash_kmer        [INT]    K-mer length of local hash process. [%u]\n", HASH_KMER);
    fprintf(stderr, "    -B --batch_size       [INT]    The number of reads to be processed in one loop. [%u]\n", BATCH_SIZE);
    fprintf(stderr, "    -l --sv_lim           [INT]    The minimum length to consider a structure variant. [%u]\n", SV_LIM);
    fprintf(stderr, "    -m --error_model      [INT]    Data type for different error model. [%u]\n", ERROR_MODEL);
    fprintf(stderr, "    -w --block            [INT]    The window length when extend node in the skeleton. [%d]\n", BLOCKS);
    fprintf(stderr, "    -x --edit_dis         [INT]    The maximum edit distance for Landua-Vishkin. [%u]\n", EDIT_DIS);
    fprintf(stderr, "    -e --seed_step        [INT]    The interval of seeding. [%u]\n", SEED_STEP);
    fprintf(stderr, "    -M --extend_mx        [INT]    The maximum length for an anchor. [%u]\n", EXTEND_MX);
    fprintf(stderr, "    -N --top_N            [INT]    Max allowed skeleton paths. [%u]\n", TOP_N);
    fprintf(stderr, "    -g --rd_mx_gap        [INT]    Maximum allowed distance in read that two anchor can be connected. [%u]\n", RD_GAP);
    fprintf(stderr, "    -G --rf_mx_gap        [INT]    Maximum allowed distance in reference that two anchor can be connected. [%u]\n", RF_GAP);
    fprintf(stderr, "    -Y --min_chain_score  [INT]    The minimum skeleton chaining score to be processed. [%u]\n", MIN_CHAIN_SCORE);
    fprintf(stderr, "    -p --secondary_ratio  [FLAOT]  Minimum secondary-to-primary alignment score ratio. [%.1f]\n", SEC_RATIO);
    fprintf(stderr, "    -P --overlap_ratio    [FLAOT]  Minimum overlap ratio to consider tow skeletons one as primary one as secondary. [%.1f]\n", OVERLAP_RATIO);
    
    //fprintf(stderr, "    -u --rd_match_mx      [INT]    ****\n");
    //fprintf(stderr, "    -f --rf_match_mx      [INT]    ****\n");

    fprintf(stderr, "\nOutput options:\n\n");
    fprintf(stderr, "    -o --output    [STR]    Output file for SV signatures in svseg format. [%s]\n", OUTPUT_PREFIX);
    fprintf(stderr, "    -h --help                      Show detailed usage.\n");

    return 1; 
}

int aln_usage(void)
{
    fprintf(stderr, "\n"); 
    fprintf(stderr, "Program:   %s: skeleton-based read alignment for variant calling\n", PACKAGE_NAME); 
    fprintf(stderr, "Usage:     %s aln [options] <index_route> <read.fa/fq>\n\n", PACKAGE_NAME);

    fprintf(stderr, "    <index_route>                 The path of RdBG index.\n");
    fprintf(stderr, "    <read.fq/fa>                  The input reads in fasta or fastq format.\n\n");

    fprintf(stderr, "Algorithm options:\n\n");
	fprintf(stderr, "    -t --thread           [INT]    Number of threads. [4]\n");
    fprintf(stderr, "    -k --index_kmer       [INT]    RdBG-index kmer length. [%u]\n", INDEX_KMER);
	fprintf(stderr, "    -s --seeding_kmer     [INT]    K-mer length of seeding process (no long than RdBG-index kmer). [%u]\n", SEEDING_KMER);
	fprintf(stderr, "    -d --hash_kmer        [INT]    K-mer length of local hash process. [%u]\n", HASH_KMER);
    fprintf(stderr, "    -B --batch_size       [INT]    The number of reads to be processed in one loop. [%u]\n", BATCH_SIZE);
    fprintf(stderr, "    -l --sv_lim           [INT]    The minimum length to consider a structure variant. [%u]\n", SV_LIM);
    fprintf(stderr, "    -m --error_model      [INT]    Data type for different error model. [%u]\n", ERROR_MODEL);
    fprintf(stderr, "    -w --block            [INT]    The window length when extend node in the skeleton. [%d]\n", BLOCKS);
    fprintf(stderr, "    -x --edit_dis         [INT]    The maximum edit distance for Landua-Vishkin. [%u]\n", EDIT_DIS);
    fprintf(stderr, "    -e --seed_step        [INT]    The interval of seeding. [%u]\n", SEED_STEP);
    fprintf(stderr, "    -M --extend_mx        [INT]    The maximum length for an anchor. [%u]\n", EXTEND_MX);
    fprintf(stderr, "    -N --top_N            [INT]    Max allowed skeleton paths. [%u]\n", TOP_N);
    fprintf(stderr, "    -g --rd_mx_gap        [INT]    Maximum allowed distance in read that two anchor can be connected. [%u]\n", RD_GAP);
    fprintf(stderr, "    -G --rf_mx_gap        [INT]    Maximum allowed distance in reference that two anchor can be connected. [%u]\n", RF_GAP);
    fprintf(stderr, "    -Y --min_chain_score  [INT]    The minimum skeleton chaining score to be processed. [%u]\n", MIN_CHAIN_SCORE);
    fprintf(stderr, "    -p --secondary_ratio  [FLAOT]  Minimum secondary-to-primary alignment score ratio. [%.1f]\n", SEC_RATIO);
    fprintf(stderr, "    -P --overlap_ratio    [FLAOT]  Minimum overlap ratio to consider tow skeletons one as primary one as secondary. [%.1f]\n", OVERLAP_RATIO);
    
    
    //fprintf(stderr, "    -u --rd_match_mx      [INT]    ****\n");
    //fprintf(stderr, "    -f --rf_match_mx      [INT]    ****\n");

    fprintf(stderr, "\nOutput options:\n\n");
    fprintf(stderr, "    -o --output    [STR]    Output file for SV signatures in svseg format. [%s]\n", OUTPUT_PREFIX);
    fprintf(stderr, "    -h --help                      Show detailed usage.\n");

    
    return 1;
}

char *const short_options = "u:f:x:vw:k:s:d:ht:e:o:O:m:l:M:B:N:g:G:p:P:t:Y:";
struct option long_options[] = {
    { "index_kmer", 1, NULL, 'k'},
    { "error_model", 1, NULL, 'm'},
    { "seeding_kmer", 1, NULL, 's'},
    { "hash_kmer", 1, NULL, 'd'},
    {"thread", 1, NULL, 't'},
    {"seed_step", 1, NULL, 'e'},
    {"sv_lim", 1, NULL, 'l'},
    {"batch_size", 1, NULL, 'B'},
    {"top_N", 1, NULL, 'N'},
    {"rd_mx_gap", 1, NULL, 'g'},
    {"rf_mx_gap", 1, NULL, 'G'},
    {"min_chain_score", 1, NULL, 'Y'},
    {"secondary_ratio", 1, NULL, 'p'},
    {"overlap_ratio", 1, NULL, 'P'},
    {"block", 1, NULL, 'w'},
    {"edit_dis", 1, NULL, 'x'},
    {"extend_mx", 1, NULL, 'M'},
    {"rd_match_mx", 1, NULL, 'u'},
    {"rf_match_mx", 1, NULL, 'f'},
    {"output", 1, NULL, 'o'},
    {"help", 0, NULL, 'h'},
    { 0, 0, 0, 0}
};


int opt_parse(int argc, char *argv[], opt_t *opt)
{
	int c; 
    char *p;
    int option_index=0;
    if (argc == 1) return help_usage();
    while((c = getopt_long(argc, argv, short_options, long_options, &option_index))>=0){
        switch(c){
            case 'u': opt->read_kmer_match_mx = atoi(optarg); break;
            case 'f': opt->read_kmer_ref_mx = atoi(optarg); break;
            case 'k': opt->db_k = atoi(optarg); break;
            case 'l': opt->the_sv_lim= atoi(optarg); break;
            case 'B': opt->batch_size = atoi(optarg); break;
            case 'N': opt->max_path_N = atoi(optarg); break;
            case 'g': opt->rd_mx_gap = atoi(optarg); break;
            case 'G': opt->rf_mx_gap = atoi(optarg); break;
            case 'p': opt->secondary_ratio = atof(optarg); break;
            case 'P': opt->overlap_ratio = atof(optarg); break;
            case 'Y': opt->min_chain_score = atoi(optarg); break;
            case 'm': opt->data_type= atoi(optarg); break;
            case 's': opt->seed_k = atoi(optarg); break;
            case 'd': opt->sdp_k = atoi(optarg); break;
            case 'w': opt->block_s = atoi(optarg); break;
            case 'h': return aln_usage(); break;
            case 'x': opt->max_lv_e = atoi(optarg); break;
            case 't': opt->thread_n = atoi(optarg); break;
            case 'e': opt->seed_step = atoi(optarg); break;
            case 'o': strcpy(opt->sv_path, optarg); break;
            case 'M': opt->ref_lv_bnd_len = atoi(optarg); break;
            default:
                fprintf(stderr,"Wrong parameters\n");
                return aln_usage();
        }
    }
    
    if(optind + 3 != argc){
        return aln_usage();
    }

    if (opt->seed_k > 22)  opt->seed_k = 22;
    else  if (opt->seed_k < 11) opt->seed_k = 11;

    if (opt->sdp_k > 11)  opt->sdp_k = 11;
    else  if (opt->sdp_k < 6) opt->sdp_k = 6;

    switch (opt->data_type)
    {
        //PacBio CLR
        case 2:
            opt->error_rate = 0.142;
            opt->er_ins = 0.0592;
            opt->er_del = 0.0301;
            opt->er_mis = 0.0527;
            break;
        //PacBio CCS
        case 1:
            opt->error_rate = 0.0172;
            opt->er_ins = 0.00087;
            opt->er_del = 0.0034;
            opt->er_mis = 0.0130;
            break;
        default:
            break;
    }
 
    opt->argv = argv;
    opt->argc = argc;
    optind++;
 	strncpy(opt->index_path, argv[optind++],sizeof(opt->index_path));
 	strncpy(opt->read_path,argv[optind++],sizeof(opt->read_path));

    //strcat(opt->sv_path, ".svseg") ;

    return 0;  
}
