#include "stdio.h"
#include "stdint.h"
#include "lv.h"
#include "desc.h"
#include "bit_operation.h"

//only for extension of small diverge region
//p is read, q is ref  
int32_t LandauVishkin_char(char *p, int32_t p_l, char *q, int32_t q_l, int32_t max_e)
{
    int32_t prv, cur, nxt;
    int32_t LV_E[99], *LV = LV_E + max_e, lv, mx_lv=-2;
    for (int i = -max_e; i < (3 + max_e); ++i) LV[i] = -2;
    for (int e = 0; e <= max_e; ++e)
    {
        prv = -2;
        cur = -2 + e;
        nxt = LV[-e+1];
        for (int i = -e; i <= e; ++i)
        {
            lv = MX(prv, cur+1);
            lv = MX(lv, nxt+1);
            lv = MN(lv, q_l-i-1);
            if ((lv+1) < p_l && (lv+i+1)< q_l)
            {
                //run LCE
                ++lv;
                while (lv<p_l && lv+i<q_l && p[lv] == q[lv+i])
                    ++lv;
                --lv;
            }
            mx_lv = MX(mx_lv, lv);
            if (lv >= p_l-1) return mx_lv;
            LV[i] = lv;
            prv = cur;
            cur = nxt;
            nxt = LV[i+2];
        } 
    }
    return mx_lv;
}

//flag = 1, right extension
//flag = -1, left extension
//x_off if the begin offset
//q saved as a huge uint64_t block, q_l and p+l is the len can extend
//p: read   q: ref

int32_t LandauVishkin_64_su(uint8_t *p, int32_t p_l, int32_t p_off, \
        uint64_t *q, int32_t q_l, int32_t q_off, int32_t *memList, int32_t *memi, int32_t max_e, int flag, int *_e, int *_i)
{
    if(flag == 1) printf("right LV:    ");
    else printf("left LV:     ");

    //printf("q_l = %d, q_off = %d, p_l = %d, p_off = %d\n", q_l, q_off, p_l, p_off);
#ifdef lv_log_s
    printf("ref:  ");
    for (int32_t k=0; k<MN(q_l, 50); ++k)
    {
        printf("%c", "ACGT"[GET_BITS_A(q, q_off+(flag*k))]);
    }
    printf("\n");
    //read
    printf("read: ");
    for (int32_t k=0; k<MN(p_l, 50); ++k)
    {
        printf("%c", "ACGT"[p[(p_off+(flag*k))]]);
    }
    printf("\n");
#endif
    int32_t prv, cur, nxt, e;
    int32_t LV_E[(1+max_e)<<1], *LV = LV_E + max_e, lv, mx_e, mx_i = 0, mx_lv=-2, fd=0, last_mx_i = 0, change_mx = 0, indel_l = 1;
    //int32_t memList[(max_e+2)<<1], mem_i = 0;
    int32_t mem_i = *memi;  
    for (int i = -max_e; i < (3 + max_e); ++i) LV[i] = -2;
    for (int e = 0; e <= max_e && fd==0; ++e)
    {
        prv = -2;
        cur = -2 + e;
        nxt = LV[-e+1];
        change_mx = 0;
        for (int i = -e; i <= e && fd==0; ++i)
        {   
            lv = MX(prv, cur+1);
            lv = MX(lv, nxt+1);
            lv = MN(lv, q_l-i-1);
            
            if ((lv+1) < p_l && (lv+i+1)< q_l)
            {
                //run LCE
                ++lv;
                while (lv<p_l && lv+i<q_l && p[(p_off+flag*lv)] == GET_BITS_A(q,(q_off+flag*(lv+i))))
                    ++lv, change_mx = 1;
                --lv;
            }

            if (lv > mx_lv || (lv == mx_lv && lv+i >= mx_lv+mx_i)) 
            { 
                mx_e = e;
                mx_i = i;
            }
            mx_lv = MX(mx_lv, lv);
            LV[i] = lv;
            prv = cur;
            cur = nxt;
            nxt = LV[i+2];
            if (lv >= p_l-1 || lv+i >= q_l-1) fd = 1;
        } 
        if(mem_i > *memi)
        {
            if (change_mx == 0)
            {
                indel_l += 1; 
            }
            else
            {
                int t = mx_i - last_mx_i;
                if (t < 0) memList[mem_i] = indel_l<<4 | 1; //I
                else if (t > 0) memList[mem_i] = indel_l<<4 | 2; //D
                else memList[mem_i] = indel_l<<4 | 8; //X
                mem_i += 1; last_mx_i = mx_i; indel_l = 1;
            } 
        }
        if(mem_i == *memi || change_mx == 1)
            memList[mem_i++] = mx_lv;
    }
#ifdef lv_log_s
    int tt = 0;
    while(tt < *memi-1)
    { 
        printf("memList[%d] = %d, ", tt, memList[tt]);
        tt++;
        printf("%d%c\n", memList[tt]>>4, "MIDNSHP=X"[memList[tt]&0xf]);
        tt++;
    }
    tt = *memi;
    while(tt < mem_i-1)
    { 
        printf("memList[%d] = %d, ", tt, memList[tt]);
        tt++;
        printf("%d%c\n", memList[tt]>>4, "MIDNSHP=X"[memList[tt]&0xf]);
        tt++;
    }
#endif

    *memi = mem_i; *_e = mx_e; *_i = mx_i;
    printf("max_lv = %d, mx_e = %d, mx_i = %d\n", mx_lv, mx_e, mx_i);
    return mx_lv;
}



//

int LCE(uint8_t *p, int32_t rd_off, uint64_t *q, int32_t rf_off, int n, int flag)
{
    int length = 0;
    while ((p[rd_off+flag*length] == GET_BITS_A(q,(rf_off+flag*length))) && (length < n))
        length++;

    return length;
}

void PRINT(uint8_t *p, int loff, int len, int flag)
{
    for(int i = 0; i < len; ++i)
    {
        printf("%c", "ACGT"[p[loff+flag*i]]);
    }
    printf("\n");
}

/*
int32_t LandauVishkin_64(uint8_t *p, int32_t p_l, int32_t p_off, \
        uint64_t *q, int32_t q_l, int32_t q_off, int32_t max_e, int flag, int *_e, int *_i)
{
    if(flag == 1) printf("right LV\n");
    else printf("left LV\n");
    printf("my, max_e = %d, q_l = %d, q_off = %d, p_l = %d, p_off = %d\n", max_e, q_l, q_off, p_l, p_off);

    printf("ref:  ");
    for (int32_t k=0; k<MN(q_l, 50); ++k)
    {
        printf("%c", "ACGT"[GET_BITS_A(q, q_off+(flag*k))]);
    }
    printf("\n");
    //read
    printf("read: ");
    for (int32_t k=0; k<MN(p_l, 50); ++k)
    {
        printf("%c", "ACGT"[p[(p_off+(flag*k))]]);
    }
    printf("\n");

    //printf("ref : "); PRINT(q, 0, q_l);
    //printf("read: "); PRINT(p, 0, p_l);
    //printf("\n");
    int32_t i = 0, j = 0, e = 0;
    int mx_e = 0, mx_i = 0, mem_i = 0;
    int32_t lrd, lrf;
    int s, s_max = 0;
    int32_t target_lrd, target_lrf, max_lv;
    //int *memList= (int *)calloc(max_e<<1, sizeof(int));
    
    while (i < p_l && j < q_l && e <= max_e)
    {
        if(p[p_off + flag*i] == GET_BITS_A(q,(q_off + flag*j)))
        {
            //printf("%c", "ACGT"[p[p_off + flag*i]]);
            ++i; ++j;
            continue;
        }
        //printf("\n");
        max_lv = i-1; *_e = 0; *_i = 0;
        s_max = 0;
        for(lrd = i; lrd < i+(max_e+1-e); ++lrd)
        {
            for(lrf = j; lrf < j+(max_e+1-e); ++lrf)
            {
                if(lrd == i && lrf == j) continue;
                if(lrd + lrf - i - j > max_e - e) continue;
                s = LCE(p, flag*lrd+p_off, q, flag*lrf+q_off, MN(p_l - lrd, q_l - lrf), flag);
                
                //printf("lrd = %d, lrf = %d, s = %d\n", lrd, lrf, s);
                if(s > s_max)
                {
                    s_max = s;
                    target_lrd = lrd;
                    target_lrf = lrf;
                }
            }
        }
        if(s_max == 0) break;
        //printf("\ni = %d, j = %d, target_lrd = %d, target_lrf = %d, s = %d\n", i, j, target_lrd, target_lrf, s_max);
        //PRINT(p, flag*target_lrd+p_off, s_max, flag);
        //judge the type: ins/del/snp
        int a = target_lrd - i, b = target_lrf - j;
        if ((a<=1) && (b<=1))//snp, single base indel
        {
            e += 1;
        }
        else
        {
            e += (a+b);
        }
        i = target_lrd + s_max;
        j = target_lrf + s_max;
        //memList[mem_i++] = target_lrf; memList[mem_i++] = j;
        max_lv = i-1; mx_e = e; mx_i = target_lrf - target_lrd;
    }
    
    printf("max_lv = %d, mx_e = %d, mx_i = %d\n", max_lv, mx_e, mx_i);
    *_e = mx_e; *_i = mx_i;
    
    return max_lv; 
}
*/



// log info for debug
int32_t LandauVishkin_64_log(uint8_t *p, int32_t p_l, int32_t p_off, \
        uint64_t *q, int32_t q_l, int32_t q_off, int32_t max_e, int flag, int *_e, int *_i)
{
#ifdef lv_log_s
    if(flag == 1) printf("right LV\n");
    else printf("left LV\n");
    printf("max_e = %d, q_l = %d, q_off = %d, p_l = %d, p_off = %d\n", max_e, q_l, q_off, p_l, p_off);
#endif

#ifdef lv_log_s
    printf("ref:  ");
    for (int32_t k=0; k<MN(q_l, 50); ++k)
    {
        printf("%c", "ACGT"[GET_BITS_A(q, q_off+(flag*k))]);
    }
    printf("\n");
    //read
    printf("read: ");
    for (int32_t k=0; k<MN(p_l, 50); ++k)
    {
        printf("%c", "ACGT"[p[(p_off+(flag*k))]]);
    }
    printf("\n");
#endif

    int32_t prv, cur, nxt;
    int32_t LV_E[(1+max_e)<<1], *LV = LV_E + max_e, lv, mx_e = 0, mx_i = 0, mx_lv=-2, fd=0;
    //for back trace
    for (int i = -max_e; i < (3 + max_e); ++i) LV[i] = -2;

    for (int e = 0; e <= max_e && fd==0; ++e)
    {
        prv = -2;
        cur = -2 + e;
        nxt = LV[-e+1];
        for (int i = -e; i <= e && fd==0; ++i)
        {   
            lv = MX(prv, cur+1);
            lv = MX(lv, nxt+1);
            lv = MN(lv, q_l-i-1);
            
            if ((lv+1) < p_l && (lv+i+1)< q_l)
            {
                ++lv;
                while (lv<p_l && lv+i<q_l)
                {
                    if ( p[(p_off+flag*lv)] == GET_BITS_A(q,(q_off+flag*(lv+i))))
                        ++lv;
                    else
                        break;
                }
                --lv;
            }
            if (lv > mx_lv || (lv == mx_lv && lv+i >= mx_lv+mx_i))
            { 
                mx_e = e;
                mx_i = i;
            }
            mx_lv = MX(mx_lv, lv);
            prv = cur;
            cur = nxt;
            nxt = LV[i+2];
            if (lv >= p_l-1 || lv+i >= q_l-1) fd = 1;
            LV[i] = lv;
        } 
    }

    *_e = mx_e;
    *_i = mx_i;
    return mx_lv;
}

