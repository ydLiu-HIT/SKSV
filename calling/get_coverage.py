#!/usr/bin/env python
# coding=utf-8
import sys
from collections import defaultdict
from intervaltree import IntervalTree

def count_coverage(chr, s, e, f, read_count, up_bound, itround):
    status = 0
    iteration = 0
    primary_num = 0

    for i in f:
        if i[0] > e: break
        if i[0] < s and i[1] > e:
            read_count.add(i[2])
            if len(read_count) >= up_bound:
                status = 1
                break

    return status


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
            record = sorted(record, key=lambda x: x[0])
            CHR = record[0][0]; start = int(record[0][4]); end = int(record[0][5])
            tmpc = None; tmps = 0; tmpe = 0; tmpl = 0
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

    for i in svsegINFO:
        print(i, svsegINFO[i])
    return svsegINFO


sig_path = sys.argv[1]
getTotalReadCnt(sig_path)
