#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include "bit_operation.h"

//ATCG -> 0123
int8_t Char2Bit[] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

uint64_t bits_right_1[] = {
						0,3,0Xf,
						0X3f,0Xff,
						0X3ff,0Xfff,
						0X3fff,0Xffff,//uint16_t:65536
						0X3ffff,0Xfffff,
						0X3fffff,0Xffffff,
						0X3ffffff,0Xfffffff,
						0X3fffffff,0Xffffffff, //uint32_t:4294967295s
						0X3ffffffff,0Xfffffffff,
						0X3fffffffff,0Xffffffffff,
						0X3ffffffffff,0Xfffffffffff,
						0X3fffffffffff,0Xffffffffffff,
						0X3ffffffffffff,0Xfffffffffffff,
						0X3fffffffffffff,0Xffffffffffffff,
						0X3ffffffffffffff,0Xfffffffffffffff,
						0X3fffffffffffffff,0Xffffffffffffffff
						};
/*
uint64_t bits_right_0[33] = {0Xffffffffffffffff,0X3fffffffffffffff,
                            0Xfffffffffffffff,0X3ffffffffffffff,
                            0Xffffffffffffff,0X3fffffffffffff,
                            0Xfffffffffffff,0X3ffffffffffff,
                            0Xffffffffffff,0X3fffffffffff,
                            0Xfffffffffff,0X3ffffffffff,
                            0Xffffffffff,0X3fffffffff,
                            0Xfffffffff,0X3ffffffff,
                            0Xffffffff,0X3fffffff,
                            0Xfffffff,0X3ffffff,
                            0Xffffff,0X3fffff,
                            0Xfffff,0X3ffff,
                            0Xffff,0X3fff,
                            0Xfff,0X3ff,
                            0Xff,0X3f,
                            0Xf,3,0
                            };
*/
uint8_t rev[128]={
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,84,4,71,4,4,4,67,4,4,4,4,4,4,78,4,
	4,4,4,4,65,4,4,4,4,4,4,4,4,4,4,4,4,
	84,4,71,4,4,4,67,4,4,4,4,4,4,78,4,4,
	4,4,4,65,4,4,4,4,4,4,4,4,4,4,4
};

//shoud avoid using like (GBA(A, p++)) it will become A[p++1>>5]>>(p++&0xf...
//get the p-th 2 bits in A array
//#define GET_BITS_A(A,p) ((((A)[(p)>>5]) >> ((31-((p)&0x1f))<<1))&0X3)
//get kmer K sub bits
//e.g. if kmer 22 mer k uint64_t, the K=k, O=(22-10)<<1, M=bits_right_1[10] will get first 10 mer
//if K=k, O=0, M=bits_right_1[10] will get the last 10 mer
#define GET_SUB_KMER(K,O) ((K)>>(O))
#define GET_SUB_KMER_M(K,O,M) (((K)>>(O))&(M))
#define MN(a,b) ((a<b)?(a):(b))
#define MX(a,b) ((a>b)?(a):(b))

uint64_t GET_BITS_A(uint64_t *A, uint32_t p)
{
    return ((((A)[(p)>>5]) >> ((31-((p)&0x1f))<<1))&0X3);
}

/*
uint64_t get_first_kmer_bits(char *seq, uint8_t klen)
{
	uint64_t value = 0;
	for (uint8_t i=0; i<klen; ++i) value = (value<<2) | Char2Bit[(int)seq[i]];
    return value;
}
*/


uint64_t get_first_kmer_bits_u8(uint8_t *seq, uint8_t klen)
{
	uint64_t value = 0;
	for (uint8_t i=0; i<klen; ++i) 
    {
        if(seq[i] >= 4)
            value = (value<<2) | (rand()%4);
        else
            value = (value<<2) | seq[i];
    }
    return value;
}

void log_kmer(uint64_t kmer, uint8_t len, int new_line)
{
	//for (uint8_t i=0; i<len; ++i)	printf(i==len-1?"%c\n":"%c", "ACGT"[kmer>>((len-1-i)<<1) & 0X3]);
	for (uint8_t i=0; i<len; ++i)	printf("%c", "ACGT"[kmer>>((len-1-i)<<1) & 0X3]);
    if (new_line) printf("\n");
    else printf(" ");
}
int multi_binsearch_offset64(uint32_t x, uint32_t v[], int64_t n, uint64_t offset, int64_t seed_binary[], int8_t k_r)
{
    int64_t low, high, mid;
    uint32_t temp = 0;
    uint64_t current_offset = 0;

    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = ((int64_t )(low + high)) >> 1;
        temp = v[mid + offset] >> (k_r << 1);
        if(x < temp)
        {
            high = mid - 1;
        }
        else if(x > temp)
        {
            low = mid + 1;
        }
        else  /*found match*/
        {
            //找到当前match的kmer，然后向两侧search，找到match的上下界
            if (mid == 0)
            {
                seed_binary[0] = offset;

                current_offset = mid + offset + 1;
                while(x == (v[current_offset] >> (k_r << 1)))
                {
                    current_offset++;
                }
                seed_binary[1] = current_offset - 1;
            }
            else if (mid == (n - 1))
            {
                seed_binary[1] = n + offset- 1;

                current_offset = mid + offset - 1;
                while(x == (v[current_offset] >> (k_r << 1)))
                {
                    current_offset--;
                }
                seed_binary[0] = current_offset + 1;
            }
            else
            {
                //down
                current_offset = mid + offset - 1;
                // while(((int)current_offset >= offset) && (x == (v[current_offset] >> (k_r << 1))))
                // {
                //     current_offset--;
                // }
                // seed_binary[0] = current_offset + 1;

                if (offset == 0)
                {
                    while((current_offset != 0) && (x == (v[current_offset] >> (k_r << 1))))
                    {
                        current_offset--;
                    }
                    seed_binary[0] = current_offset + 1;
                }
                else
                {
                    while(x == (v[current_offset] >> (k_r << 1)))
                    {
                        current_offset--;
                    }
                    seed_binary[0] = current_offset + 1;    
                }

                current_offset = mid + offset + 1;
                while((current_offset < offset + n) && (x == (v[current_offset] >> (k_r << 1))))
                {
                    current_offset++;
                }
                seed_binary[1] = current_offset - 1;
            }

            return 1;
        }
    }
    return -1;
}


int binsearch_range(uint64_t key, uint32_t *v, int64_t n,  int64_t *range, uint8_t k_off)
{
    int64_t l=0, r=n-1, m;
    uint32_t tmp = 0;
    range[0] = range[1] = -1;

    if (n < 5 && 0)
    {
        range[0] = range[1] = -1;
        for (m=0; m<=n-1; ++m)
        {
        tmp = v[m] >> k_off;
            if (tmp > key) break;
            if (tmp == key)
            {
                range[0] = range[1] = m;
                break;
            }
            if (tmp > key) break;
        }
        if (range[0] == -1) return -1;
        if (m+1 <= n-1)
        for (++m; m<=n-1; ++m)
        {
        tmp = v[m] >> k_off;
            if (tmp > key) break;
            if (tmp == key)
            {
                range[1] = m;
                continue;
            }
        }
        return 1;
        //printf("%zu %zu\n", range[0], range[1]);
    }

    while (l <= r)
    {
        m = (l+r)/2;
        tmp = v[m] >> k_off;
        if (tmp == key)
        {
            range[0] = range[1] = m;
            
            //run low bound
            int64_t sl=l, sr=m-1, sm;
            while (sl <= sr)
            {
                sm = (sl+sr)/2;
                tmp = v[sm] >> k_off;
                if (tmp == key)
                {
                    range[0] = sm;
                    sr = sm-1;
                }
                else if (tmp > key) sr = sm - 1;
                else    sl = sm + 1;
            }

            //run upper bound
            sl = m+1; sr = r;
            while (sl <= sr)
            {
                sm = (sl+sr)/2;
                tmp = v[sm] >> k_off;
                if (tmp == key)
                {
                    range[1] = sm;
                    sl = sm+1;
                }
                else if (tmp > key) sr = sm - 1;
                else    sl = sm + 1;
            }
            return 1;
        }
        else if (tmp > key) r = m - 1;
        else l = m + 1;
    }
    return -1;
}


int64_t hash_get_uqp_id(uint64_t x, uint32_t v[], uint64_t n, uint64_t  *K, uint64_t m)
{
    if (x >= K[m-1]) return n-1;
    if (x/16 < n-1)
    {
        uint32_t tmp_p = v[x/16];
        if ((x/16+1 < n-1) && (v[x/16+1] <= tmp_p)) return tmp_p;
        if (K[tmp_p+1] > x+16)
            return tmp_p;
        else
            return tmp_p+1;
    }
    else return v[n-1];
}

int64_t binsearch_low_bound_uqp(uint64_t x, uint64_t v[], uint64_t n)
{
    int64_t low=0, high=n-1, mid;

    while ( low <= high )
    {
        mid = ((int64_t )(low + high)) >> 1;
        if(x < v[mid])      high = mid - 1;
        else if(x > v[mid])     low = mid + 1;
        else  return mid;
    }
    return high;
}

uint16_t binsearch_rht(uint16_t x, uint16_t *A, uint16_t low, uint16_t high)
{
    //A = A + 1; ++low; ++high; 
    if (x < A[low] || x > A[high])
        return 0xffff;
    uint16_t mid;
    while (low <= high)
    {
        mid = (low+high) >> 1;
        if (x < A[mid]) high = mid - 1;
        else if (x > A[mid]) low = mid + 1;
        else return mid;
    }
    return 0xffff;
}

int trans_char_A_to_64_B(char *A, uint64_t a_len, uint64_t *B, uint64_t max_b)
{
    //memset(B, 0, max_b * sizeof(B));
    uint64_t tmp = 0, B_i = 0;
    if (a_len >= max_b) return -1;
    for (uint64_t i=0; i<a_len; ++i)
    {
        tmp = tmp << 2 |Char2Bit[(int)A[i]];
        if (i==31) {B[B_i++] = tmp; tmp = 0;}    
    }
    if (a_len & 0X1f)    B[B_i++] = (tmp << ((32 - (a_len & 0X1f))<<1));
    return 0;
}

int trans_char_u8(char *A, uint64_t a_len, uint8_t *B, uint64_t max_b)
{
    //memset(B, 0, max_b * sizeof(B));
    uint64_t tmp = 0, B_i = 0;
    //assert(a_len < max_b); 
    if (a_len >= max_b) return -1;
    for (uint64_t i=0; i<a_len; ++i)
        B[i] = Char2Bit[(int)A[i]];
    return 0;
}

int revComRead(char *read, char *rcRead, int len_read)
{
	for(int i=0;i<len_read;i++)
		rcRead[i] = rev[(int)read[len_read - 1 - i]];
	rcRead[len_read] = '\0';
	return 0;
}


void log_seq_char(char *seq, uint32_t len)
{
    for (uint32_t i=0; i<len; ++i) printf("%c", seq[i]);
    printf("\n");
}
void log_seq_64(uint64_t *seq, uint32_t len, uint32_t off)
{
    for (uint32_t i=0; i<len; ++i) printf("%c", "ACGT"[GET_BITS_A(seq, off+i)]);
    printf("\n");
}
