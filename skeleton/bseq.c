#include "bseq.h"


bseq_file_t *bseq_open(const char *fn)
{
    bseq_file_t *fp;
    gzFile f;
    f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
    if (f == 0) return 0;
    fp = (bseq_file_t* )calloc(1, sizeof(bseq_file_t));
    fp->fp = f;
    fp->ks = kseq_init(fp->fp);
    return fp;
}

void bseq_close(bseq_file_t *fp)
{
    kseq_destroy(fp->ks);
    gzclose(fp->fp);
    free(fp);
}

uint32_t bseq_read(bseq_file_t *fp, uint32_t batch, rd_seq_dt* query_info)
{
    uint32_t seqii = 0;
    int64_t kr1 = 1;
    kseq_t *seq1 = fp->ks;
    uint32_t _max_readlen = 1000000;

    for(seqii = 0; (seqii < batch) && ((kr1 = kseq_read(seq1)) > 0); seqii++)
    {
        //strdup, later process like before and compare the memory
        query_info[seqii].name = strdup(seq1->name.s);
        query_info[seqii].read_seq = strdup(seq1->seq.s);
        if (seq1->seq.l > _max_readlen) seq1->seq.l = _max_readlen;
        query_info[seqii].read_length = seq1->seq.l;
        for (uint32_t i = 0; i < seq1->seq.l; ++i) //convert U to T 
            if (query_info[seqii].read_seq[i] == 'u' || query_info[seqii].read_seq[i] == 'U')
                --query_info[seqii].read_seq[i];
        query_info[seqii].read_seq_r = strdup(query_info[seqii].read_seq);
        revComRead(query_info[seqii].read_seq, query_info[seqii].read_seq_r, seq1->seq.l); 
    }

    return seqii;
}

void free_rd_seq(uint32_t read_in, rd_seq_dt *query_info)
{
	for (uint32_t r_i = 0; r_i < read_in; ++r_i)
	{
		if (query_info[r_i].name != NULL)
		{
			free(query_info[r_i].name);
			query_info[r_i].name = NULL;
		}

		if (query_info[r_i].read_seq != NULL)
		{
			free(query_info[r_i].read_seq);
			query_info[r_i].read_seq = NULL;
			free(query_info[r_i].read_seq_r);
			query_info[r_i].read_seq_r = NULL;
		}
	}
}


//uint32_t bseq_read(bseq_t *bp, uint32_t batch, read_type *rp)
//{
//    uint32_t seq_i = 0;
//    uint32_t i;
//    kseq_t *seqs = bp->seq;
//
//    for( ; (seq_i<batch) && (kseq_read(seqs)>=0); ++seq_i)
//    {
//        //strdup, later process like before and compare the memory
//        rp[seq_i].name = strdup(seq1->name.s);
//        query_info[seqii].read_seq = strdup(seq1->seq.s);
//        query_info[seqii].read_length = seq1->seq.l;
//        for (i = 0; i < seq1->seq.l; ++i) //convert U to T
//            if (query_info[seqii].read_seq[i] == 'u' || query_info[seqii].read_seq[i] == 'U')
//                --query_info[seqii].read_seq[i];
//    }
//
//    return seqii;
//}
//
//uint32_t bseq_read_2pass(bseq_file_t *fp, uint32_t batch, seq_io* seqio)
//{
//    uint32_t seqii = 0;
//    int64_t kr1 = 1;
//    uint32_t i;
//    kseq_t *seq1 = fp->ks;
//
//    for(seqii = 0; (seqii < batch) && ((kr1 = kseq_read(seq1)) > 0); seqii++)
//    {
//        seqio[seqii].name = strdup(seq1->name.s);
//        seqio[seqii].read_seq = strdup(seq1->seq.s);
//        seqio[seqii].qual = seq1->qual.l? strdup(seq1->qual.s):0;
//        seqio[seqii].read_length = seq1->seq.l;
//        for (i = 0; i < seq1->seq.l; ++i) //convert U to T
//            if (seqio[seqii].read_seq[i] == 'u' || seqio[seqii].read_seq[i] == 'U')
//                --seqio[seqii].read_seq[i];
//    }
//
//    return seqii;
//}
