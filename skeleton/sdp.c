#define _GNU_SOURCE
#include <stdlib.h>
#include <string.h>
#include "sdp.h"
#include "lv.h"
#include "desc.h"
#include <math.h>
#include "bit_operation.h"
#include <assert.h>


void init_KPHT(kpht_t *kht, uint8_t kmer_l, uint32_t _mx)
{
    kht->kmer_l = kmer_l;
    kht->_mx = _mx;

    //4^7 space
    kht->ht = (kmerp_t *)malloc((1<<(kmer_l<<1)) * sizeof(kmerp_t));
    kht->pos_v = (uint32_t *)malloc(_mx * sizeof(uint32_t));
    kht->tmp_pos_v = (uint32_t *)malloc(_mx * sizeof(uint32_t));
}

//todo, if kmer match pos with in kmer_l, only save one record
void set_kmerp_v(kmerp_t *kp, uint32_t *pv, uint8_t kl, uint32_t v)
{
    if (kp->len > 0) pv[(kp->tail)-1] = v;
    else kp->head = v; 
    ++(kp->len); 
    kp->tail = v;
    pv[v-1] = v;

    //v = 10, kp->head = 10, kp->len = 1, kp->tail = 10, pv[10-1]=10
    //v = 20, pv[10-1] = 20, kp->len = 2, kp->tail = 20, pv[20-1] = 20
    //v = 33, pv[20-1] = 33, kp->len = 3, kp->tail = 33, pv[33-1] = 33
    //pv[9] = 10, pv[19]=20, pv[32] = 33, kp->head = 10, kp->tail = 33
}

//void build_KPHT(kpht_t *kht, char *seq, uint32_t seq_l)
void build_KPHT(kpht_t *kht, uint8_t *seq, uint32_t seq_l)
{
    uint32_t kmer_tmp;
    uint8_t kl = kht->kmer_l; //_sdp_k  7
    uint32_t kmask = bits_right_1[kl];
    kmerp_t *_ht = kht->ht;
    uint32_t *_pos_v = kht->pos_v, *_tmp_pos_v = kht->tmp_pos_v;

    if (seq_l < kl || seq_l-kl+1 >= kht->_mx)
    {
        fprintf(stderr, "read length [%lu] exceed KHT spaces limited[%lu]\n",\
                (unsigned long)seq_l, (unsigned long)kht->_mx);
        assert(seq_l < kl || seq_l-kl+1 >= kht->_mx); 
        return;
    }
    
    uint8_t sc;
    uint32_t kmer_n = seq_l-kl+1;
    memset(_ht, 0, sizeof(kmerp_t) * (1<<(kl<<1))); 
    kmer_tmp = get_first_kmer_bits_u8(seq, kl);
    
    set_kmerp_v(_ht+kmer_tmp, _tmp_pos_v, kl, 1);
    for (uint32_t i=1; i<kmer_n; ++i)
    {
        sc = seq[i+kl-1];
        if(sc >= 4)
            kmer_tmp = ((kmer_tmp << 2) & kmask) | (rand()%4);
        else
            kmer_tmp = ((kmer_tmp << 2) & kmask) | sc;
        set_kmerp_v(_ht+kmer_tmp, _tmp_pos_v, kl, i+1);
    }
    // adjust;
    uint32_t tv = 0, ttv, R_index = 0;
    kmer_tmp = get_first_kmer_bits_u8(seq, kl);
    for (uint32_t i=0; i<kmer_n; ++i)
    {
        if (i>0)
        {
            sc = seq[i+kl-1];
            if(sc >= 4)
                kmer_tmp = ((kmer_tmp << 2) & kmask) | (rand()%4);
            else
                kmer_tmp = ((kmer_tmp << 2) & kmask) | sc;
        }
        if (_tmp_pos_v[i] == 0) continue;
        tv = i+1;
        (_ht+kmer_tmp)->head = R_index+1;
        while (_tmp_pos_v[tv-1] != 0)
        {
            _pos_v[R_index++] = tv;
            ttv = tv;
            tv = _tmp_pos_v[tv-1];
            _tmp_pos_v[ttv-1] = 0;
        }
    }
    /* if no adjust search by:
       int hit_n = (_ht+kmer_tmp)->len;
       uint32_t hit_i = (_ht+kmer_tmp)->head, tv = hit_i;
       for (int j=0; j<hit_n; ++j)
       {
            uint32_t k_i = tv-1; 
            tv = _tmp_pos_v[tv-1];
       }
       if adjust seach by:
            uint32_t k_i = _pos_v[tv-1]-1; 
            ++tv;
    */

}

uint32_t get_new_kmer_tmp(int i, int dir_fg, uint8_t kl, uint64_t *ref, uint32_t offset)
{
    uint32_t kmer_tmp = 0;
    if (dir_fg > 0)
    {
        for (int t=i; t<(i+kl); ++t) 
            kmer_tmp = (kmer_tmp<<2) | GET_BITS_A(ref, offset+t);
    }
    else
    {
        for (int t=(i+kl-1); t>=i; --t) 
            kmer_tmp = (kmer_tmp<<2) | GET_BITS_A(ref, offset-t);
    }
    return kmer_tmp;
}


int run_greed_kht(kpht_t *kht, uint8_t *rd, uint32_t rd_l, uint32_t rd_off, uint64_t *ref, uint32_t ref_l, uint32_t offset, int32_t block_s, double block_s_off, uint32_t *rd_ex, uint32_t *rf_ex, int dir_fg, int dt_t)
{
    //rd_block_s = 30, rf_block_s_off = (err_del-err_ins) * 100
    //dir_fg: 1 right ext  -1 left ext    
#ifdef print_greed
    printf("Inner greedy extension~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~, kmer = %d, block_s = %d, block_s_off = %f\n", kht->kmer_l, block_s, block_s_off);
    printf("rd_l = %d, rd_off = %d, ref_l = %d, ref_off = %d\n", rd_l, rd_off, ref_l, offset);
#endif

    uint32_t kmer_tmp;
    uint8_t kl = kht->kmer_l;
    uint32_t kmask = bits_right_1[kl];
    kmerp_t *_ht = kht->ht;
    uint32_t *_pos_v = kht->pos_v, *_tmp_pos_v = kht->tmp_pos_v;

    uint32_t hit_idx, hit_n, hit_i, last_rd_pos=0;
    int flag_n = 1,fd=0, global_fg_f = 0;
    uint32_t last_i=0, last_ki=0;
    int32_t bi = 0, bst_ki, blk_l, blk_r, ln_s, rf_block_b = 0;
    double rf_block_acc = 0;
    int bst_l_s;
    int32_t rf_block_o = 0, rd_block_o = 0;

    for (int i=0; i<(ref_l-kl+1); ++i)
    {
        fd = 0;
        bst_ki = block_s + 1;
        //minimum match offset
        bst_l_s = block_s + 1;
        blk_l = block_s * bi;
        blk_r = block_s * (bi+1);

        if (flag_n == 0)
        {
            if (dir_fg > 0)
                kmer_tmp =((kmer_tmp<<2 & kmask) | GET_BITS_A(ref, offset+(i+kl-1))); 
            else
                kmer_tmp =((kmer_tmp>>2) | (GET_BITS_A(ref, offset-(i+kl-1)) << ((kl-1)<<1))); 
        }
        else
        {
            kmer_tmp = get_new_kmer_tmp(i, dir_fg, kl, ref, offset);
            flag_n = 0;
        }
        //search 
        if ((_ht+kmer_tmp)->len > 0)
        {
            int hit_n = (_ht+kmer_tmp)->len;
#ifdef print_greed
            //printf("hit_n = %d\n", hit_n);
#endif
            uint32_t hit_i = (_ht+kmer_tmp)->head, tv = hit_i;
            uint16_t mx_e = 0;               
            for (int j=0; j<hit_n; ++j)
            {
                uint32_t k_i = _pos_v[tv-1]-1;
                ++tv;
                if (dir_fg > 0) //right ext
                {
                    if (k_i < rd_off)    continue;
                    k_i -= rd_off;
                    if (k_i < blk_l)    continue;
                    if (k_i >= blk_r)   break;
                    if (k_i+kl-1 >= rd_l) break; 
                }
                else
                {
                    if (k_i > (rd_off - kl + 1)) break;
                    k_i = (rd_off-kl+1-k_i);
                    if (k_i >= blk_r)   continue;
                    if (k_i+kl-1 >= rd_l) continue; 
                    if (k_i < blk_l)    break;
                }
                //if not fall in range continue
                int ln_s = G_ABS((k_i - blk_l)*1.0 - (i - rf_block_b));
#ifdef print_greed
                printf("k_i = %u, blk_l = %u, i = %u, rf_block_b = %u, ln_s = %d, comp = %f\n", k_i, blk_l, i, rf_block_b, ln_s, G_ABS(2*block_s*block_s_off*2*3));
#endif

                if (dt_t == 1)
                {
                    if (ln_s > block_s/2) continue;
                }
                else if (dt_t == 2)
                {
                    //within two block and , indel len is 2
                    if (ln_s > G_ABS(2 * block_s * block_s_off * 2 * 3))
                       continue;
                }
                if (ln_s > bst_l_s)  continue;
                bst_l_s = ln_s;
                bst_ki = k_i; 
                fd = 1;
            }
        }

        if (fd == 0)
        {
            if ((i - rf_block_b) >= block_s)   
            {
                break;
            }
        }
        else
        {
            rf_block_o = rf_block_b;
            rd_block_o = block_s * (bi);
            global_fg_f = 1;
            last_i = i;
            last_ki = bst_ki;
            double indel_s = 2.0;
            double tmp_rf_b;

            if (dt_t == 1)
            {
                //greed offset jump
                tmp_rf_b = i + blk_r - bst_ki;
            
                if (bst_l_s < 3)
                    tmp_rf_b = (i - rf_block_acc) - (bst_ki - blk_l) + rf_block_acc + block_s;
                else
                    tmp_rf_b = (rf_block_acc + block_s + indel_s * ((1) * block_s * block_s_off));
            }
            else
            {
                tmp_rf_b = i + blk_r - bst_ki;
            }

            rf_block_acc = tmp_rf_b;
            //update
            //if (tmp_rf_b+1 < (ref_l-kl+1))
            {
                rf_block_b = (int32_t)rf_block_acc;
                ++bi;
                i = rf_block_b - 1;
            }
            flag_n = 1;
        }
    } 
    if (global_fg_f <= 0)
    {
        *rd_ex = 0;
        *rf_ex = 0;
        return 0;
    }
    *rd_ex = last_ki + kl-1;
    *rf_ex = last_i + kl-1;
    

    //LV
    int rd_z_off = last_ki+kl, rf_z_off = last_i+kl, lv_s = 0, lv_e = 0, lv_i = 0;

    if (dir_fg > 0) {rd_z_off += rd_off; rf_z_off += offset;}
    else {rd_z_off = rd_off-rd_z_off; rf_z_off = offset - rf_z_off;}
    
    int av_e = 0;
    if (dt_t == 1)
        av_e = 0;
    else if (dt_t == 2)
        av_e = 3;

    //printf("inner greedy function\n");
    //lv_s = LandauVishkin_64_su(rd, rd_l-(last_ki+kl-1), rd_z_off, ref, ref_l-(last_i+kl-1), rf_z_off, av_e, dir_fg, &lv_e, &lv_i);
#ifdef print_gd_detail
    printf("lst ki%d i%d\n", last_ki, last_i);
    lv_s = LandauVishkin_64_log(rd, rd_l-(last_ki+kl-1), rd_z_off, ref, ref_l-(last_i+kl-1), rf_z_off, av_e, dir_fg, &lv_e, &lv_i);
    printf("%d %d\n", lv_s, lv_s+lv_i);
#endif 
    
    *rd_ex += (lv_s);
    *rf_ex += (lv_s + lv_i);

    return 1;
}
