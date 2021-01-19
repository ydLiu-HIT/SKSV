#ifndef BIT_OPERATION_H_
#define BIT_OPERATION_H_
#include <stdint.h>

extern int8_t Char2Bit[];
extern uint64_t bits_right_1[];

#define GET_SUB_KMER(K,O) ((K)>>(O))
#define GET_SUB_KMER_M(K,O,M) (((K)>>(O))&(M))
#define MN(a,b) ((a<b)?(a):(b))
#define MX(a,b) ((a>b)?(a):(b))
#define G_ABS(x) (((x)>0)?(x):(-(x)))

uint64_t GET_BITS_A(uint64_t *A, uint32_t p);

uint64_t get_first_kmer_bits(char *seq, uint8_t klen);
uint64_t get_first_kmer_bits_u8(uint8_t *seq, uint8_t klen);

void log_kmer(uint64_t kmer, uint8_t len, int new_lne);

void log_seq_64(uint64_t *seq, uint32_t len, uint32_t off);
void log_seq_char(char *seq, uint32_t len);


int binsearch_range(uint64_t key, uint32_t *v, int64_t n,  int64_t *range, uint8_t k_off);


int64_t hash_get_uqp_id(uint64_t x, uint32_t v[], uint64_t n, uint64_t  *K, uint64_t m);
int64_t binsearch_low_bound_uqp(uint64_t x, uint64_t v[], uint64_t n);
uint16_t binsearch_rht(uint16_t x, uint16_t *A, uint16_t low, uint16_t high);


int multi_binsearch_offset64(uint32_t x, uint32_t v[], int64_t n, uint64_t offset, int64_t seed_binary[], int8_t k_r);
int trans_char_A_to_64_B(char *A, uint64_t a_len, uint64_t *B, uint64_t max_b);
int trans_char_u8(char *A, uint64_t a_len, uint8_t *B, uint64_t max_b);

int revComRead(char *read, char *rcRead, int len_read);
#endif
