#!/usr/bin/env python

''' 
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title: cuteSV 
 * @author: tjiang
 * @data: Apr 26th 2020
 * @version V1.0.6
'''

from cuteSV_Description import parseArgs
from multiprocessing import Pool
from CommandRunner import *
from cuteSV_resolveINV import run_inv
from cuteSV_resolveTRA import run_tra
from cuteSV_resolveINDEL import run_ins, run_del
from cuteSV_resolveDUP import run_dup
from cuteSV_genotype import generate_output, load_valuable_chr
from collections import defaultdict
import shlex
import os
import re
import argparse
import logging
import sys
import time

DEBUG = False


def load_rst(path, candidate, mi, md):
    read_id = str()
    read_length = 0
    file = open(path, 'r')
    for line in file:
        if line[0] == '>':
            wd = [s for s in re.split(' |,|\[|\]', line.strip()) if s != '']
            read_id = wd[0][1:]
            read_length = int(wd[1])
        else:
            wd = [s.strip() for s in re.split('\||\<|\>', line) if s.strip() != '']
            if len(wd) < 2:
                continue
            #print(read_id, read_length, wd)
            split_read = list()
            for tp in wd:
                ele = [int(s) if s.isdigit() else s for s in re.split(',|\[|\]| ', tp) if s != '']
                if ele[6] < 300:
                    continue
                split_read.append([ele[2], ele[3], ele[4], ele[5], str(ele[0]), ele[1]])

            analysis_split_read(split_read, 30, read_length, read_id, candidate, 100000, mi, md)

    file.close()

def analysis_inv(ele_1, ele_2, read_name, candidate, SV_size):
    if ele_1[5] == '+':
        # +-
        if ele_1[3] - ele_2[3] >= SV_size:
            if ele_2[0] + 0.5 * (ele_1[3] - ele_2[3]) >= ele_1[1]:
                candidate["INV"][ele_1[4]].append(["++", ele_2[3], ele_1[3], read_name])
                #candidate["INV"][ele_1[4]].append(["--", ele_2[3], ele_1[3], read_name])
                # head-to-head
                # 5'->5'
        if ele_2[3] - ele_1[3] >= SV_size:
            if ele_2[0] + 0.5 * (ele_2[3] - ele_1[3]) >= ele_1[1]:
                candidate["INV"][ele_1[4]].append(["++", ele_1[3], ele_2[3], read_name])
                #candidate["INV"][ele_1[4]].append(["--", ele_1[3], ele_2[3], read_name])
                # head-to-head
                # 5'->5'
    else:
        # -+
        if ele_2[2] - ele_1[2] >= SV_size:
            if ele_2[0] + 0.5 * (ele_2[2] - ele_1[2]) >= ele_1[1]:
                candidate["INV"][ele_1[4]].append(["--", ele_1[2], ele_2[2], read_name])
                #candidate["INV"][ele_1[4]].append(["++", ele_1[2], ele_2[2], read_name])
                # tail-to-tail
                # 3'->3'
        if ele_1[2] - ele_2[2] >= SV_size:
            if ele_2[0] + 0.5 * (ele_1[2] - ele_2[2]) >= ele_1[1]:
                candidate["INV"][ele_1[4]].append(["--", ele_2[2], ele_1[2], read_name])
                #candidate["INV"][ele_1[4]].append(["++", ele_2[2], ele_1[2], read_name])
                # tail-to-tail
                # 3'->3'

def generate_combine(candidate, tmp_sigs, svtype, CHR, read_name, merge_dis):
    if len(tmp_sigs) == 0:
        pass
    elif len(tmp_sigs) == 1:
        if svtype == "INS":
            candidate[svtype][CHR].append([tmp_sigs[0][0], tmp_sigs[0][1], tmp_sigs[0][2], read_name, tmp_sigs[0][-1]])
        else:
            candidate[svtype][CHR].append([tmp_sigs[0][0], tmp_sigs[0][1], read_name])
    else:
        tmp_sigs = sorted(tmp_sigs, key=lambda x: x[0])
        tsig = tmp_sigs[0]
        if svtype == "INS":
            tsig += [tmp_sigs[0][0]]
            for i in tmp_sigs[1:]:
                if i[0] - tsig[4] <= merge_dis and i[3] == tsig[3]:
                    tsig[1] += i[1]
                    tsig[4] = i[0]
                else:
                    candidate[svtype][CHR].append([(tsig[0]+tsig[4])/2, tsig[1], tsig[2], read_name, tsig[3]])

                    tsig = i
                    tsig.append(i[0]) 
            candidate[svtype][CHR].append([(tsig[0]+tsig[4])/2, tsig[1], tsig[2], read_name, tsig[3]])
        else:
            tsig += [sum(tmp_sigs[0])]
            for i in tmp_sigs[1:]:
                if i[0] - tsig[2] <= merge_dis:
                    tsig[1] += i[1]
                    tsig[2] = sum(i)
                else:
                    candidate[svtype][CHR].append([tsig[0], tsig[1], read_name])

                    tsig = i
                    tsig.append(sum(i)) # tsig.append(i[0])
            candidate[svtype][CHR].append([tsig[0], tsig[1], read_name])

    
def analysis_split_read(split_read, SV_size, RLength, read_name, candidate, MaxSize, merge_dis_ins, merge_dis_del):
    '''
    read_start	read_end	ref_start	ref_end	chr	strand
    #0		#1		#2		#3	#4	#5
    '''
    SP_list = sorted(split_read, key=lambda x: x[0])

    if DEBUG:
        print(read_name, len(SP_list))
        for i in SP_list:
            print(i)
        print("==============")

    # detect INS involoved in a translocation
    trigger_INS_TRA = 0

    # Store Strands of INV

    if len(SP_list) == 2:
        ele_1 = SP_list[0]
        ele_2 = SP_list[1]
        if ele_1[4] == ele_2[4] and ele_1[5] != ele_2[5]:
            analysis_inv(ele_1, ele_2, read_name, candidate, SV_size)
    else:
        for a in range(len(SP_list[1:-1])):
            ele_1 = SP_list[a]
            ele_2 = SP_list[a + 1]
            ele_3 = SP_list[a + 2]

            if ele_1[4] == ele_2[4] == ele_3[4]:
                if ele_1[5] == ele_3[5] and ele_1[5] != ele_2[5]:
                    if ele_2[5] == '-':
                        # +-+
                        #print(ele_2[0], ele_3[2], ele_1[3], ele_1[1], 0.5*(ele_3[2]-ele_1[3]))
                        #print(ele_3[0], ele_3[2], ele_1[3], ele_2[1])
                        if ele_2[0] + 0.5 * (ele_3[2] - ele_1[3]) >= ele_1[1] and ele_3[0] + 0.5 * (ele_3[2] - ele_1[3]) >= ele_2[1]:
                            # No overlaps in split reads
                            candidate["INV"][ele_1[4]].append(["++", ele_1[3], ele_2[3], read_name])
                            # head-to-head
                            # 5'->5'
                            candidate["INV"][ele_1[4]].append(["--", ele_2[2], ele_3[2], read_name])
                            # tail-to-tail
                            # 3'->3'
                    else:
                        # -+-
                        #print(ele_1[1], ele_2[0], ele_1[2], ele_3[3], 0.5*(ele_1[2]-ele_3[3]))
                        #print(ele_3[0], ele_2[1], ele_1[2], ele_3[3])
                        if ele_1[1] <= ele_2[0] + 0.5 * (ele_1[2] - ele_3[3]) and ele_3[0] + 0.5 * (ele_1[2] - ele_3[3]) >= ele_2[1]:
                            # No overlaps in split reads
                            candidate["INV"][ele_1[4]].append(["++", ele_3[3], ele_2[3], read_name])
                            # head-to-head
                            # 5'->5'
                            candidate["INV"][ele_1[4]].append(["--", ele_2[2], ele_1[2], read_name])
                            # tail-to-tail
                            # 3'->3'

                if ele_1[5] != ele_3[5]:
                    if ele_2[5] == ele_1[5]:
                        # ++-/--+
                        analysis_inv(ele_2, ele_3, read_name, candidate, SV_size)
                    else:
                        # +--/-++
                        analysis_inv(ele_1, ele_2, read_name, candidate, SV_size)

    tmp_sigs_ins = defaultdict(list)  
    tmp_sigs_del = defaultdict(list) 
    CHR = str()
    for a in range(len(SP_list[:-1])):
        ele_1 = SP_list[a]
        ele_2 = SP_list[a + 1]
        if ele_1[4] == ele_2[4]:
            if ele_1[5] == ele_2[5]:
                CHR = ele_2[4]
                reverse = False
                # dup & ins & del
                if ele_1[5] == '-':
                    ele_1 = [RLength - SP_list[a + 1][1], RLength - SP_list[a + 1][0]] + SP_list[a + 1][2:]
                    ele_2 = [RLength - SP_list[a][1], RLength - SP_list[a][0]] + SP_list[a][2:]
                    reverse = True
                    if DEBUG:
                        print(ele_1)
                        print(ele_2)
                        print("****************")

                #if ele_1[3] - ele_2[2] >= SV_size:
                #    candidate["DUP"][ele_2[4]].append([ele_2[2], ele_1[3], read_name])

                if ele_1[3] - ele_2[2] >= SV_size and ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1] <= MaxSize:
                    #candidate["DUP"][ele_2[4]].append([ele_2[2], ele_1[3], read_name])
                    if ele_1[1] < ele_2[0]:
                        tmp_sigs_ins[CHR].append([ele_2[2], ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1], ele_1[1] - int(ele_1[3]-ele_2[2]), reverse]) # add read start pos
                    else:
                        if ele_2[0]-ele_1[1] > ele_2[2]-ele_1[3]:
                            tmp_sigs_ins[CHR].append([(ele_2[2] + ele_1[3]) / 2, ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1], ele_1[1] + int((ele_2[2]-ele_1[3])/2), reverse])
                        else:
                            tmp_sigs_del[CHR].append([(ele_2[2] + ele_1[3]) / 2, ele_2[2] + ele_1[1] - ele_2[0] - ele_1[3]])

                if ele_1[3] - ele_2[2] < SV_size:
                    if ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1] >= SV_size:
                        if ele_2[2] - ele_1[3] <= 200 and ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1] <= MaxSize:
                            #tmp_sigs_ins[CHR].append([(ele_1[3] + ele_2[2])/2, ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1]])
                            if ele_1[3] < ele_2[2]:
                                tmp_sigs_ins[CHR].append([ele_1[3], ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1], ele_1[1], reverse])
                            else:
                                tmp_sigs_ins[CHR].append([ele_2[2], ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1], ele_1[1] - int(ele_1[3]-ele_2[2]), reverse])
                    if ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3] >= SV_size:
                        if ele_2[0] - ele_1[1] <= 200 and ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3] <= MaxSize:
                            #tmp_sigs_del[CHR].append([ele_1[3], ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3]]) 
                            if ele_1[1] < ele_2[0]:
                                tmp_sigs_del[CHR].append([ele_1[3], ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3]]) 
                            else:
                                tmp_sigs_del[CHR].append([ele_1[3] + ele_2[0] - ele_1[1] - 1, ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3]]) 
        else:
            trigger_INS_TRA = 1
            '''
			*********Description*********
			*	TYPE A:		N[chr:pos[	*
			*	TYPE B:		N]chr:pos]	*
			*	TYPE C:		[chr:pos[N	*
			*	TYPE D:		]chr:pos]N	*
			*****************************
			'''
            if ele_2[0] - ele_1[1] <= 100:
                if ele_1[5] == '+':
                    if ele_2[5] == '+':
                        # +&+
                        if ele_1[4] < ele_2[4]:
                            candidate["TRA"][ele_1[4]].append(['A', ele_1[3], ele_2[4], ele_2[2], read_name])
                            # N[chr:pos[
                        else:
                            candidate["TRA"][ele_2[4]].append(['D', ele_2[2], ele_1[4], ele_1[3], read_name])
                            # ]chr:pos]N
                    else:
                        # +&-
                        if ele_1[4] < ele_2[4]:
                            candidate["TRA"][ele_1[4]].append(['B', ele_1[3], ele_2[4], ele_2[3], read_name])
                            # N]chr:pos]
                        else:
                            candidate["TRA"][ele_2[4]].append(['B', ele_2[3], ele_1[4], ele_1[3], read_name])
                            # N]chr:pos]
                else:
                    if ele_2[5] == '+':
                        # -&+
                        if ele_1[4] < ele_2[4]:
                            candidate["TRA"][ele_1[4]].append(['C', ele_1[2], ele_2[4], ele_2[2], read_name])
                            # [chr:pos[N
                        else:
                            candidate["TRA"][ele_2[4]].append(['C', ele_2[2], ele_1[4], ele_1[2], read_name])
                            # [chr:pos[N
                    else:
                        # -&-
                        if ele_1[4] < ele_2[4]:
                            candidate["TRA"][ele_1[4]].append(['D', ele_1[2], ele_2[4], ele_2[3], read_name])
                            # ]chr:pos]N
                        else:
                            candidate["TRA"][ele_2[4]].append(['A', ele_2[3], ele_1[4], ele_1[2], read_name])
                            # N[chr:pos[ 

    
    if len(SP_list) >= 3 and trigger_INS_TRA == 1:
        for i in range(len(SP_list)-2):
            for j in range(i+2, len(SP_list)):
                A = SP_list[i]; B = SP_list[i+1]; C = SP_list[j]
                reverse = False
                if A[4] != B[4] and A[4] == C[4] and A[5] == C[5]:
                    if A[5] == '+':
                        ele_1 = A
                        ele_2 = C
                    else:
                        ele_1 = [RLength - C[1], RLength - C[0]] + C[2:]
                        ele_2 = [RLength - A[1], RLength - A[0]] + A[2:]
                        reverse = True
                    
                    dis_ref = ele_2[2] - ele_1[3]
                    dis_read = ele_2[0] - ele_1[1]
                    if DEBUG:
                        print(ele_1)
                        print(ele_2)
                        print(dis_ref, dis_read)
                    if dis_ref < 100 and dis_read - dis_ref >= SV_size and dis_read - dis_ref <= MaxSize:
                        tmp_sigs_ins[ele_2[4]].append([min(ele_2[2], ele_1[3]), dis_read - dis_ref, ele_1[1]+int(dis_ref/2), reverse])

                    if dis_ref <= -SV_size:
                        candidate["DUP"][ele_2[4]].append([ele_2[2], ele_1[3], read_name])
                        
                    break
    if DEBUG:
        print("tmp_sig_ins:", tmp_sigs_ins)
        print("tmp_sig_del:", tmp_sigs_del)
    for chr in tmp_sigs_ins.keys():
        generate_combine(candidate, tmp_sigs_ins[chr], "INS", chr, read_name, merge_dis_ins)
    for chr in tmp_sigs_del.keys():
        generate_combine(candidate, tmp_sigs_del[chr], "DEL", chr, read_name, merge_dis_del)



def CheckFileExist(fn, sfx=""):
    if not os.path.isfile(fn+sfx):
        sys.exit("Error: %s not found" %(fn+sfx))
    return os.path.abspath(fn+sfx)


#def getTotalReadCnt(sig_path):
#    fin = open(sig_path, 'r')
#    svsegINFO = defaultdict(list)
#    for line in fin:
#        if line.startswith('>'):
#            seq_name = line.strip('>').split(' ')[0]
#            continue
#        
#        line = line.strip().split('|')
#        if len(line) == 1:
#            record = line[0].split(',')
#            svsegINFO[record[0]].append([int(record[4]), int(record[5]), seq_name])
#        else:
#            record = [l.split(',') for l in line]
#            record = sorted(record, key=lambda x: (x[0], int(x[4]), int(x[5])))
#            #record = sorted(record, key=lambda x: x[0])
#            CHR = record[0][0]; start = int(record[0][4]); end = int(record[0][5])
#            tmpc = CHR; tmps = start; tmpe = end; tmpl = 0
#            i = 1
#            while i < len(record):
#                if int(record[i][6]) < 300: 
#                    i += 1
#                    continue
#                if record[i][0] == CHR:
#                    if int(record[i][4]) - end < 50000:
#                        end = int(record[i][5])
#                    else:
#                        if end - start > tmpl:
#                            tmps = start; tmpe = end; tmpc = record[i][0]
#                            tmpl = end - start
#                        CHR = record[i][0]; start = int(record[i][4]); end = int(record[i][5])
#                else:
#                    # may be tra
#                    #svsegINFO[tmpc].append([tmps, tmpe, seq_name])
#
#                    if end - start > tmpl:
#                        tmps = start; tmpe = end; tmpc = record[i][0]
#                        tmpl = end - start
#
#                    CHR = record[i][0]; start = int(record[i][4]); end = int(record[i][5])
#                i += 1
#            else:
#                if end - start > tmpl:
#                    tmps = start; tmpe = end
#                svsegINFO[tmpc].append([tmps, tmpe, seq_name])
#
#    fin.close()
#
#    for chr in svsegINFO:
#        svsegINFO[chr] = sorted(svsegINFO[chr], key=lambda x:(x[0], x[1]))
#
#    return svsegINFO
#

def getTotalReadCnt(sig_path):
    fin = open(sig_path, 'r')
    svsegINFO = defaultdict(list)

    for line in fin:
        if line.startswith('>'):
            seq_name = line.strip('>').split(' ')[0]
            continue
        
        line = line.strip().split('|')
        if len(line) == 1:
            record = line[0].split(',')
            svsegINFO[record[0]].append([int(record[4]), int(record[5]), seq_name])
        else:
            record = [l.split(',') for l in line]
            record = sorted(record, key=lambda x: (x[0], int(x[2])))
            #record = sorted(record, key=lambda x: x[0])
            CHR = record[0][0]; start = int(record[0][4]); end = int(record[0][5])
            tmpc = CHR; tmps = start; tmpe = end; tmpl = 0
            i = 1
            while i < len(record):
                if int(record[i][6]) < 300: 
                    i += 1
                    continue
                if record[i][0] == CHR:
                    if record[i][1] == "+":
                        if int(record[i][4]) - end < 50000:
                            end = int(record[i][5])
                        else:
                            if end - start > tmpl:
                                tmps = start; tmpe = end; tmpc = record[i][0]
                                tmpl = end - start
                            CHR = record[i][0]; start = int(record[i][4]); end = int(record[i][5])
                    else:
                        if start - int(record[i][5]) < 50000:
                            start = int(record[i][4])
                        else:
                            if end - start > tmpl:
                                tmps = start; tmpe = end; tmpc = record[i][0]
                                tmpl = end - start
                            CHR = record[i][0]; start = int(record[i][4]); end = int(record[i][5])
                else:
                    # may be tra
                    #svsegINFO[tmpc].append([tmps, tmpe, seq_name])

                    if end - start > tmpl:
                        tmps = start; tmpe = end; tmpc = record[i][0]
                        tmpl = end - start

                    CHR = record[i][0]; start = int(record[i][4]); end = int(record[i][5])
                i += 1
            else:
                if end - start > tmpl:
                    tmps = start; tmpe = end
                svsegINFO[tmpc].append([tmps, tmpe, seq_name])

    fin.close()

    for chr in svsegINFO:
        svsegINFO[chr] = sorted(svsegINFO[chr], key=lambda x:(x[0], x[1]))

    #for i in svsegINFO:
    #    for j in svsegINFO[i]:
    #        print(i, j)

    return svsegINFO


def rst_ctrl(args, argv):
    starttime = time.time()
    candidate = dict()
    candidate["DEL"] = defaultdict(list) 
    candidate["INS"] = defaultdict(list)
    candidate["INV"] = defaultdict(list)
    candidate["DUP"] = defaultdict(list)
    candidate["TRA"] = defaultdict(list)
    if args.print_allele_seq:
        if not os.path.exists(args.read):
            logging.info("[ERROR] --read is need to print allele sequence!")
            exit()
    if args.work_dir[-1] == '/':
        temporary_dir = args.work_dir
    else:
        temporary_dir = args.work_dir+'/'
    if not os.path.exists(temporary_dir):
        os.makedirs(temporary_dir)
    
    load_rst(args.input, candidate, args.merge_ins_threshold, args.merge_del_threshold)
    file = open("%ssignals.txt"%temporary_dir, 'w')
    for chr in candidate["DEL"]:
        for ele in candidate["DEL"][chr]:
            file.write("%s\t%s\t%d\t%d\t%s\n" %("DEL", chr, ele[0], ele[1], ele[2]))

    for chr in candidate["INS"]:
        for ele in candidate["INS"][chr]:
            file.write("%s\t%s\t%d\t%d\t%d\t%s\t%s\n" % ("INS", chr, ele[0], ele[1], ele[2], ele[3], '-' if ele[4] else '+'))

    for chr in candidate["INV"]:
        for ele in candidate["INV"][chr]:
            file.write("%s\t%s\t%s\t%d\t%d\t%s\n" % ("INV", chr, ele[0], ele[1], ele[2], ele[3]))

    for chr in candidate["TRA"]:
        for ele in candidate["TRA"][chr]:
            file.write("%s\t%s\t%s\t%d\t%s\t%d\t%s\n" % ("TRA", chr, ele[0], ele[1], ele[2], ele[3], ele[4]))

    for chr in candidate["DUP"]:
        for ele in candidate["DUP"][chr]:
            file.write("%s\t%s\t%d\t%d\t%s\n" %("DUP", chr, ele[0], ele[1], ele[2]))
    
    file.close()

    analysis_pools = Pool(processes=5)
    cmd_del = ("cat %ssignals.txt | grep DEL | sort -u | sort -k 2,2 -k 3,3n > %sDEL.sigs"%(temporary_dir, temporary_dir))
    cmd_ins = ("cat %ssignals.txt | grep INS | sort -u | sort -k 2,2 -k 3,3n > %sINS.sigs"%(temporary_dir, temporary_dir))
    cmd_inv = ("cat %ssignals.txt | grep INV | sort -u | sort -k 2,2 -k 3,3 -k 4,4n > %sINV.sigs"%(temporary_dir, temporary_dir))
    cmd_tra = ("cat %ssignals.txt | grep TRA | sort -u | sort -k 2,2 -k 5,5 -k 3,3 -k 4,4n > %sTRA.sigs"%(temporary_dir, temporary_dir))
    cmd_dup = ("cat %ssignals.txt | grep DUP | sort -u | sort -k 1,1r -k 2,2 -k 3,4n > %sDUP.sigs"%(temporary_dir, temporary_dir))
    #cmd_del = ("cat %ssignals.txt | grep DEL | sort -T /home/ydliu/ -u | sort -T /home/ydliu/ -k 2,2 -k 3,3n > %sDEL.sigs"%(temporary_dir, temporary_dir))
    #cmd_ins = ("cat %ssignals.txt | grep INS | sort -T /home/ydliu/ -u | sort -T /home/ydliu/ -k 2,2 -k 3,3n > %sINS.sigs"%(temporary_dir, temporary_dir))
    #cmd_inv = ("cat %ssignals.txt | grep INV | sort -T /home/ydliu/ -u | sort -T /home/ydliu/ -k 2,2 -k 3,3 -k 4,4n > %sINV.sigs"%(temporary_dir, temporary_dir))
    #cmd_tra = ("cat %ssignals.txt | grep TRA | sort -T /home/ydliu/ -u | sort -T /home/ydliu/ -k 2,2 -k 5,5 -k 3,3 -k 4,4n > %sTRA.sigs"%(temporary_dir, temporary_dir))
    #cmd_dup = ("cat %ssignals.txt | grep DUP | sort -T /home/ydliu/ -u | sort -T /home/ydliu/ -k 1,1r -k 2,2 -k 3,4n > %sDUP.sigs"%(temporary_dir, temporary_dir))


    for i in [cmd_del, cmd_ins, cmd_inv, cmd_tra, cmd_dup]:
        analysis_pools.map_async(exe, (i,))
    analysis_pools.close()
    analysis_pools.join()
    logging.info("Loading svseg information in %0.2f seconds." % (time.time() - starttime))

def main_ctrl(args, argv):
    fin_ref_fai = CheckFileExist(args.ref+".fai")
    contigINFO = list()
    if args.work_dir[-1] == '/':
        temporary_dir = args.work_dir
    else:
        temporary_dir = args.work_dir+'/'
    with open(fin_ref_fai, 'r') as fin:
        for row in fin:
            column = row.strip().split("\t")
            contig_name = column[0]
            local_ref_len = int(column[1])
            contigINFO.append([contig_name, local_ref_len])
            

    # get read count for genotyping
    svsegINFO = defaultdict(list)
    if args.genotype:
        logging.info("Collecting alignment skeletons for genotyping.")
        svsegINFO = getTotalReadCnt(args.input)

    valuable_chr = load_valuable_chr(temporary_dir)

    logging.info("Clustering structural variants.")
    analysis_pools = Pool(processes=int(args.threads))
    result = list()

    # +++++DEL+++++
    for chr in valuable_chr["DEL"]:
        para = [("%s%s.sigs" % (temporary_dir, "DEL"),
                 chr,
                 "DEL",
                 args.min_support,
                 args.diff_ratio_merging_DEL,
                 args.max_cluster_bias_DEL,
                 # args.diff_ratio_filtering_DEL,
                 min(args.min_support, 5),
                 svsegINFO[chr],
                 args.genotype,
                 args.gt_round)]
        result.append(analysis_pools.map_async(run_del, para))
        logging.info("Finished %s:%s." % (chr, "DEL"))

    # +++++INS+++++
    for chr in valuable_chr["INS"]:
        para = [("%s%s.sigs" % (temporary_dir, "INS"),
                 chr,
                 "INS",
                 args.min_support,
                 args.diff_ratio_merging_INS,
                 args.max_cluster_bias_INS,
                 # args.diff_ratio_filtering_INS,
                 min(args.min_support, 5),
                 svsegINFO[chr],
                 args.genotype,
                 args.gt_round,
                 args.read,
                 False,
                 args.print_allele_seq)]
        result.append(analysis_pools.map_async(run_ins, para))
        logging.info("Finished %s:%s." % (chr, "INS"))

    # +++++INV+++++
    for chr in valuable_chr["INV"]:
        para = [("%s%s.sigs" % (temporary_dir, "INV"),
                 chr,
                 "INV",
                 args.min_support,
                 args.max_cluster_bias_INV,
                 args.min_size,
                 svsegINFO[chr],
                 args.genotype,
                 args.max_size,
                 args.gt_round)]
        result.append(analysis_pools.map_async(run_inv, para))
        logging.info("Finished %s:%s." % (chr, "INV"))

    ## +++++DUP+++++
    for chr in valuable_chr["DUP"]:
        para = [("%s%s.sigs" % (temporary_dir, "DUP"),
                 chr,
                 args.min_support,
                 args.max_cluster_bias_DUP,
                 args.min_size,
                 svsegINFO[chr],
                 args.genotype,
                 args.max_size,
                 args.gt_round)]
        result.append(analysis_pools.map_async(run_dup, para))
        logging.info("Finished %s:%s." % (chr, "DUP"))

    # +++++TRA+++++
    for chr in valuable_chr["TRA"]:
        for chr2 in valuable_chr["TRA"][chr]:
            para = [("%s%s.sigs" % (temporary_dir, "TRA"),
                     chr,
                     chr2,
                     args.min_support,
                     args.diff_ratio_filtering_TRA,
                     args.max_cluster_bias_TRA,
                     svsegINFO[chr],
                     args.genotype,
                     args.gt_round)]
            result.append(analysis_pools.map_async(run_tra, para))
    logging.info("Finished %s." % ("TRA/BND"))

    analysis_pools.close()
    analysis_pools.join()

    semi_result = list()
    for res in result:
        try:
            semi_result += res.get()[0]
        except:
            pass

    logging.info("Writing to your output file.")

    # sort SVs by [chr] and [pos]
    semi_result = sorted(semi_result, key=lambda x: (x[0], int(x[2])))
    
    generate_output(args, semi_result, contigINFO, argv, args.print_allele_seq)

    if args.retain_work_dir:
        pass
    else:
        logging.info("Cleaning temporary files.")
        cmd_remove_tempfile = ("rm -r %ssignatures" % (temporary_dir))
        exe(cmd_remove_tempfile)

def setupLogging(debug=False):
    logLevel = logging.DEBUG if debug else logging.INFO
    logFormat = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(stream=sys.stderr, level=logLevel, format=logFormat)
    logging.info("Running %s" % " ".join(sys.argv))


def run(argv):
    args = parseArgs(argv)
    setupLogging(False)

    rst_ctrl(args, argv)

    starttime = time.time()
    main_ctrl(args, argv)
    logging.info("Finished in %0.2f seconds." % (time.time() - starttime))


if __name__ == '__main__':
    run(sys.argv[1:])
