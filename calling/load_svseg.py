import sys
from collections import defaultdict
import re
import time

DEBUG = False 

def load_rst(path, candidate):
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

            analysis_split_read(split_read, 30, read_length, read_id, candidate, 100000)

    file.close()


def analysis_inv(ele_1, ele_2, read_name, candidate, SV_size):
    if ele_1[5] == '+':
        # +-
        if ele_1[3] - ele_2[3] >= SV_size:
            if ele_2[0] + 0.5 * (ele_1[3] - ele_2[3]) >= ele_1[1]:
                candidate["INV"][ele_1[4]].append(["++", ele_2[3], ele_1[3], read_name])
                # head-to-head
                # 5'->5'
        if ele_2[3] - ele_1[3] >= SV_size:
            if ele_2[0] + 0.5 * (ele_2[3] - ele_1[3]) >= ele_1[1]:
                candidate["INV"][ele_1[4]].append(["++", ele_1[3], ele_2[3], read_name])
                # head-to-head
                # 5'->5'
    else:
        # -+
        if ele_2[2] - ele_1[2] >= SV_size:
            if ele_2[0] + 0.5 * (ele_2[2] - ele_1[2]) >= ele_1[1]:
                candidate["INV"][ele_1[4]].append(["--", ele_1[2], ele_2[2], read_name])
                # tail-to-tail
                # 3'->3'
        if ele_1[2] - ele_2[2] >= SV_size:
            if ele_2[0] + 0.5 * (ele_1[2] - ele_2[2]) >= ele_1[1]:
                candidate["INV"][ele_1[4]].append(["--", ele_2[2], ele_1[2], read_name])
                # tail-to-tail
                # 3'->3'

def generate_combine(tmp_sigs, svtype, CHR, read_name, merge_dis):
    if len(tmp_sigs) == 0:
        pass
    elif len(tmp_sigs) == 1:
        candidate[svtype][CHR].append([tmp_sigs[0][0], tmp_sigs[0][1], read_name])
    else:
        tmp_sigs = sorted(tmp_sigs, key=lambda x: x[0])
        tsig = tmp_sigs[0]
        if svtype == "INS":
            tsig += [tmp_sigs[0][0]]
        else:
            tsig += [sum(tmp_sigs[0])]

        for i in tmp_sigs[1:]:
            if i[0] - tsig[2] <= merge_dis:
                tsig[1] += i[1]
                if svtype == "INS":
                    tsig[2] = i[0]
                else:
                    tsig[2] = sum(i)
            else:
                if svtype == "INS":
                    candidate[svtype][CHR].append([(tsig[0]+tsig[2])/2, tsig[1], read_name])
                else:
                    candidate[svtype][CHR].append([tsig[0], tsig[1], read_name])

                tsig = i
                tsig.append(i[0])
        if svtype == "INS":
            candidate[svtype][CHR].append([(tsig[0]+tsig[2])/2, tsig[1], read_name])
        else:
            candidate[svtype][CHR].append([tsig[0], tsig[1], read_name])

    
def analysis_split_read(split_read, SV_size, RLength, read_name, candidate, MaxSize):
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
    merge_dis_ins = 500 
    merge_dis_del = 500
    for a in range(len(SP_list[:-1])):
        ele_1 = SP_list[a]
        ele_2 = SP_list[a + 1]
        if ele_1[4] == ele_2[4]:
            if ele_1[5] == ele_2[5]:
                CHR = ele_2[4]
                # dup & ins & del
                if ele_1[5] == '-':
                    ele_1 = [RLength - SP_list[a + 1][1], RLength - SP_list[a + 1][0]] + SP_list[a + 1][2:]
                    ele_2 = [RLength - SP_list[a][1], RLength - SP_list[a][0]] + SP_list[a][2:]
                    if DEBUG:
                        print(ele_1)
                        print(ele_2)
                        #print(ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1])
                        print("****************")

                #if ele_1[3] - ele_2[2] >= SV_size:
                #    candidate["DUP"][ele_2[4]].append([ele_2[2], ele_1[3], read_name])

                if ele_1[3] - ele_2[2] >= SV_size and ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1] <= MaxSize:
                    #candidate["DUP"][ele_2[4]].append([ele_2[2], ele_1[3], read_name])
                    if ele_1[1] < ele_2[0]:
                        tmp_sigs_ins[CHR].append([ele_2[2], ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1]])
                    else:
                        if ele_2[0]-ele_1[1] > ele_2[2]-ele_1[3]:
                            tmp_sigs_ins[CHR].append([(ele_2[2] + ele_1[3]) / 2, ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1]])
                        else:
                            tmp_sigs_del[CHR].append([(ele_2[2] + ele_1[3]) / 2, ele_2[2] + ele_1[1] - ele_2[0] - ele_1[3]])

                if ele_1[3] - ele_2[2] < SV_size:
                    if ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1] >= SV_size:
                        if ele_2[2] - ele_1[3] <= 200 and ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1] <= MaxSize:
                            #tmp_sigs_ins[CHR].append([(ele_1[3] + ele_2[2])/2, ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1]])
                            if ele_1[3] < ele_2[2]:
                                tmp_sigs_ins[CHR].append([ele_1[3], ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1]])
                            else:
                                tmp_sigs_ins[CHR].append([ele_2[2], ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1]])
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

    #if DEBUG:
    #    print("tmp_sig_ins:", tmp_sigs_ins)
    #    print("tmp_sig_del:", tmp_sigs_del)
    #for chr in tmp_sigs_ins.keys():
    #    generate_combine(tmp_sigs_ins[chr], "INS", chr, read_name, merge_dis_ins)
    #for chr in tmp_sigs_del.keys():
    #    generate_combine(tmp_sigs_del[chr], "DEL", chr, read_name, merge_dis_del)


    #if len(SP_list) >= 3 and trigger_INS_TRA == 1:
    #    if SP_list[0][4] == SP_list[-1][4]:
    #        if SP_list[0][5] != SP_list[-1][5]:
    #            pass
    #        else:
    #            if SP_list[0][5] == '+':
    #                ele_1 = SP_list[0]
    #                ele_2 = SP_list[-1]
    #            else:
    #                ele_1 = [RLength - SP_list[-1][1], RLength -
    #                         SP_list[-1][0]] + SP_list[-1][2:]
    #                ele_2 = [RLength - SP_list[0][1], RLength -
    #                         SP_list[0][0]] + SP_list[0][2:]
    #            
    #            dis_ref = ele_2[2] - ele_1[3]
    #            dis_read = ele_2[0] - ele_1[1]
    #            if DEBUG:
    #                print(ele_1)
    #                print(ele_2)
    #                print(dis_ref, dis_read)
    #            if dis_ref < 100 and dis_read - dis_ref >= SV_size and dis_read - dis_ref <= MaxSize:
    #                # print(min(ele_2[2], ele_1[3]), dis_read - dis_ref, read_name)
    #                #tmp_sigs_ins.append([min(ele_2[2], ele_1[3]), dis_read - dis_ref])
    #                candidate["INS"][ele_2[4]].append([min(ele_2[2], ele_1[3]), dis_read - dis_ref, read_name])

    #            if dis_ref <= -SV_size:
    #                candidate["DUP"][ele_2[4]].append([ele_2[2], ele_1[3], read_name])

    if len(SP_list) >= 3 and trigger_INS_TRA == 1:
        for i in range(len(SP_list)-2):
            for j in range(i+2, len(SP_list)):
                A = SP_list[i]; B = SP_list[i+1]; C = SP_list[j]
                if A[4] != B[4] and A[4] == C[4] and A[5] == C[5]:
                    if A[5] == '+':
                        ele_1 = A
                        ele_2 = C
                    else:
                        ele_1 = [RLength - C[1], RLength - C[0]] + C[2:]
                        ele_2 = [RLength - A[1], RLength - A[0]] + A[2:]
                    
                    dis_ref = ele_2[2] - ele_1[3]
                    dis_read = ele_2[0] - ele_1[1]
                    if DEBUG:
                        print(ele_1)
                        print(ele_2)
                        print(dis_ref, dis_read)
                    if dis_ref < 100 and dis_read - dis_ref >= SV_size and dis_read - dis_ref <= MaxSize:
                        # print(min(ele_2[2], ele_1[3]), dis_read - dis_ref, read_name)
                        tmp_sigs_ins[ele_2[4]].append([min(ele_2[2], ele_1[3]), dis_read - dis_ref])
                        #candidate["INS"][ele_2[4]].append([min(ele_2[2], ele_1[3]), dis_read - dis_ref, read_name])

                    if dis_ref <= -SV_size:
                        candidate["DUP"][ele_2[4]].append([ele_2[2], ele_1[3], read_name])
                        
                    break
    if DEBUG:
        print("tmp_sig_ins:", tmp_sigs_ins)
        print("tmp_sig_del:", tmp_sigs_del)
    for chr in tmp_sigs_ins.keys():
        generate_combine(tmp_sigs_ins[chr], "INS", chr, read_name, merge_dis_ins)
    for chr in tmp_sigs_del.keys():
        generate_combine(tmp_sigs_del[chr], "DEL", chr, read_name, merge_dis_del)


if __name__ == '__main__':
    starttime = time.time()
    candidate = dict()
    candidate["DEL"] = defaultdict(list) 
    candidate["INS"] = defaultdict(list)
    candidate["INV"] = defaultdict(list)
    candidate["DUP"] = defaultdict(list)
    candidate["TRA"] = defaultdict(list)
    load_rst(sys.argv[1], candidate)
    file = open(sys.argv[2], 'w')
    for sv_type in ["DEL", "INS", "INV", "DUP", 'TRA']:
        try:
            for chr in candidate[sv_type]:
                for ele in candidate[sv_type][chr]:
                    if len(ele) == 3:
                        file.write("%s\t%s\t%d\t%d\t%s\n" %
                                   (sv_type, chr, ele[0], ele[1], ele[2]))
                    elif len(ele) == 5:
                        file.write("%s\t%s\t%s\t%d\t%s\t%d\t%s\n" % (sv_type, chr, ele[0],
                                                                     ele[1], ele[2], ele[3], ele[4]))
                    elif len(ele) == 4:
                        file.write("%s\t%s\t%s\t%d\t%d\t%s\n" % (
                            sv_type, chr, ele[0], ele[1], ele[2], ele[3]))
                        # INV chr strand pos1 pos2 read_ID
        except:
            pass
    file.close()
    print("[INFO]: Finished in %0.2f seconds." % (time.time() - starttime))
