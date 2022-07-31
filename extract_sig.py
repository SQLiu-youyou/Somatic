from cgi import print_environ
from telnetlib import WILL
from turtle import pos
from description import load_sigs_chr
from cluster_INDEL import run_ins,resolution_INS,run_del
import statistics
import sys
from numpy import append, min_scalar_type, piecewise, promote_types
import pysam
import cigar
from python_utils import listify
from multiprocessing import Pool
from parameter import parseArgs
import os
import argparse
import time
import gc
from Bio import SeqIO
# from Bio import SeqIO


def generate_candidate_in_a_read(chr_name,sig_list,candidate,read_name, combine_min_size, svtype):
    # combine_min_size = 200
    #cuteSV是500
    
    #对同一条read中，很接近的两个ins进行合并
    if len(sig_list) == 0:
        pass
    elif len(sig_list) == 1:
        if chr_name not in candidate[svtype]:
            candidate[svtype][chr_name] = list()
        #后续还需添加序列信息，目前存储：变异位点，变异长度，read_name,read来源
        if svtype == 'INS':
            candidate['INS'][chr_name].append([sig_list[0][0],sig_list[0][1],read_name,sig_list[0][2]])
        if svtype == 'DEL':
            candidate['DEL'][chr_name].append([sig_list[0][0],sig_list[0][1],read_name,sig_list[0][2]])
    else:
        #用第一列保存ins开始的位置，第二列保存INS的长度，用最后一列待比较的ins位点
        if svtype == 'INS':
            temp_sig = sig_list[0]
            temp_sig += [sig_list[0][0]]
            for sig in sig_list[1:]:
                # temp_sig[2] = sig[0]
                if (sig[1] + temp_sig[1])/2 > combine_min_size:
                    if sig[0] - temp_sig[-1] <= combine_min_size:
                        temp_sig[1] = temp_sig[1] + sig[1] 
                        temp_sig[-1] = sig[0]
                    else:
                        if chr_name not in candidate['INS']:
                            candidate['INS'][chr_name] = list()
                        candidate['INS'][chr_name].append([temp_sig[0],temp_sig[1],read_name,temp_sig[2]])
                        temp_sig = sig
                        temp_sig.append(sig[0])
                else:
                    if sig[0] - temp_sig[-1] <= (sig[1] + temp_sig[1])/2:
                        temp_sig[1] = temp_sig[1] + sig[1] 
                        temp_sig[-1] = sig[0]
                    else:
                        if chr_name not in candidate['INS']:
                            candidate['INS'][chr_name] = list()
                        candidate['INS'][chr_name].append([temp_sig[0],temp_sig[1],read_name,temp_sig[2]])
                        temp_sig = sig
                        temp_sig.append(sig[0])
            
            #处理最后一个merge的ins
            if chr_name not in candidate['INS']:
                    candidate['INS'][chr_name] = list()
            candidate['INS'][chr_name].append([temp_sig[0],temp_sig[1],read_name,temp_sig[2]])
        else:
            #第一列保存del开始的位置，第二列保存del的长度，最后一列使用del开始+del长度的和作为待比较的位点
            temp_sig = sig_list[0]
            temp_sig += [sum(sig_list[0])]
            for sig in sig_list[1:]:
                
                if (sig[1] + temp_sig[1])/2 > combine_min_size:
                    if sig[0] - temp_sig[-1] <= combine_min_size:
                        temp_sig[1] = temp_sig[1] + sig[1] 
                        temp_sig[-1] = sum(sig)
                    else:
                        if chr_name not in candidate['DEL']:
                            candidate['DEL'][chr_name] = list()
                        candidate['DEL'][chr_name].append([temp_sig[0],temp_sig[1],read_name,temp_sig[2]])
                        temp_sig = sig
                        # temp_sig.append(sig[0])
                        temp_sig.append(sum(sig))
                else:
                    if sig[0] - temp_sig[-1] <= (sig[1] + temp_sig[1])/2:
                        temp_sig[1] = temp_sig[1] + sig[1] 
                        temp_sig[-1] = sum(sig)
                    else:
                        if chr_name not in candidate['DEL']:
                            candidate['DEL'][chr_name] = list()
                        candidate['DEL'][chr_name].append([temp_sig[0],temp_sig[1],read_name,temp_sig[2]])
                        temp_sig = sig
                        # temp_sig.append(sig[0])
                        temp_sig.append(sum(sig))           
                
                
                # temp_sig[2] = sig[0]
                #注释开始，上一句不是
                # if sig[0] - temp_sig[2] <= combine_min_size:
                #     temp_sig[1] = temp_sig[1] + sig[1] 
                #     temp_sig[2] = sum(sig)
                # else:
                #     if chr_name not in candidate['DEL']:
                #         candidate['DEL'][chr_name] = list()
                #     candidate['DEL'][chr_name].append([temp_sig[0],temp_sig[1],read_name])
                #     temp_sig = sig
                #     temp_sig.append(sum(sig))
                    # print(temp_sig)
            #处理最后一个merge的ins
            if chr_name not in candidate['DEL']:
                    candidate['DEL'][chr_name] = list()
            candidate['DEL'][chr_name].append([temp_sig[0],temp_sig[1],read_name,temp_sig[2]])
            # print(candidate)

signal = {1 << 2: 0, \
			1 >> 1: 1, \
			1 << 4: 2, \
			1 << 11: 3, \
			1 << 4 | 1 << 11: 4}
'''
	1 << 2 means unmapped read
    1 >> 1 means normal_foward read (正链)
	1 << 4 means reverse_complement read
	1 << 11 means supplementary alignment read
	1 << 4 | 1 << 11 means supplementary alignment with reverse_complement read
'''
def detect_flag(Flag):
	back_sig = signal[Flag] if Flag in signal else 0
	return back_sig
dic_starnd = {1: '+', 2: '-'}

def generate_candidate_intra_read(split_read,read_length,read_name,candidate,chase_ins_min_size, \
    chase_ins_max_size):
    #cuteSV 30
    # chase_ins_min_size = 20 
    # chase_ins_max_size = 100000
    chase_del_min_size = 20
    chase_del_max_size = 3000000
    '''
	read_start	read_end	ref_start	ref_end     chr	    strand
	#0			#1			#2			#3		    #4	    #5
	'''
    primary = split_read[0]
    sp_list = sorted(split_read,key = lambda x:x[0])
    # print(len(sp_list)-1)
    # print(len(sp_list[:-1]))
    for i in range(0,len(sp_list)-1):
        ele1 = sp_list[i]
        ele2 = sp_list[i+1]
        if ele1[4] == ele2[4] and ele1[5] == ele2[5]:
            #负链转换成正链
            if ele1[5] == '-':
                ele1 = [read_length-sp_list[i+1][1],read_length-sp_list[i+1][0]]+sp_list[i+1][2:]
                ele2 = [read_length-sp_list[i][1],read_length-sp_list[i][0]]+sp_list[i][2:]
            #与dup相反
            if ele1[3] - ele2[2] < chase_ins_min_size:
                #INS
                #sv长度（read差距-ref差距）大于最小设定长度并小于最大sv长度限制, ref上差距 < 100, 
                if ele2[0] - ele1[1] - (ele2[2] - ele1[3]) >= chase_ins_min_size:
                    if ele2[2] - ele1[3] <= 100 and \
                            (ele2[0] - ele1[1] - (ele2[2] - ele1[3]) <= chase_ins_max_size or chase_ins_max_size == -1):
                        if ele1[4] not in candidate['INS']:
                            candidate['INS'][ele1[4]] = list()
                        if ele1 == primary or ele2 == primary:
                            candidate['INS'][ele1[4]].append([(ele1[1]+ele2[3])/2, \
                            ele2[0] - ele1[1] - (ele2[2] - ele1[3]),read_name,1])
                        else:
                            candidate['INS'][ele1[4]].append([(ele1[1]+ele2[3])/2, \
                            ele2[0] - ele1[1] - (ele2[2] - ele1[3]),read_name,0])
                #DEL
                if ele2[2] - ele1[3] - (ele2[0] - ele1[1]) >= chase_del_min_size:
                    if ele2[0] - ele1[1] <= 100 and \
                            (ele2[2] - ele1[3] - (ele2[0] - ele1[1]) <= chase_del_max_size or chase_ins_max_size == -1):
                        if ele1[4] not in candidate['DEL']:
                            candidate['DEL'][ele1[4]] = list()
                        if ele1 == primary or ele2 == primary:
                            candidate['DEL'][ele1[4]].append([ele1[3],\
                            ele2[2] - ele1[3] - (ele2[0] - ele1[1]), read_name,1])
                        else:
                            candidate['DEL'][ele1[4]].append([ele1[3],\
                            ele2[2] - ele1[3] - (ele2[0] - ele1[1]), read_name,0])


def organize_supple_info(primary_info,supple_info,read_length,read_name,candidate,\
    max_split_parts, chase_ins_min_size, chase_ins_max_size):
    # max_split_parts = 7

    split_read = list()
    if len(primary_info) > 0:
        split_read.append(primary_info)
    for supple in supple_info:
        info = supple.split(',')
        local_cigar = info[3]
        local_mapq = int(info[4]) 
        ref_start = int(info[1])
        chr = info[0]
        if local_mapq > 0:
            seq = list(cigar.Cigar(local_cigar).items())
            #supple里没有hardclip？
            if seq[0][1] == 'S':
                supple_first = seq[0][0]
            else :
                supple_first = 0
            if seq[-1][1] == 'S':
                supple_end = seq[-1][0]
            else:
                supple_end = 0
            # 计算ref上的偏差，I对ref不产生偏差
            bias = 0
            for item in seq:
                if item[1] == 'M' or item[1] == 'D' or item[1] == '=' or item[1] == 'X':
                    bias += item[0]
            #这里不是负链转换成正链，只是根据clip以及正负链信息，确定alignment的位置区间
            if info[2] == '+':
                split_read.append([supple_first,read_length-supple_end,ref_start,ref_start+bias,\
                    chr,'+'])
            else:
                split_read.append([supple_end,read_length-supple_first,ref_start,ref_start+bias,\
                    chr,'-'])
    # or max_split_parts == -1:
    if len(split_read) > 0 and len(split_read) <= max_split_parts:
        generate_candidate_intra_read(split_read,read_length,read_name,candidate, chase_ins_min_size, \
            chase_ins_max_size)


def parse_read(chr_name, read, candidate, min_mapq, sig_min_cigar_size, max_split_parts, \
    chase_ins_min_size, chase_ins_max_size, combine_min_size):
    # min_mapq = 20
    # sig_min_cigar_size = 20    
    process_signal = detect_flag(read.flag)
    Combine_sig_in_same_read_ins = list()
    Combine_sig_in_same_read_del = list()
    #找到所有ins的变异（inter）(reference染色体，起始位置，片段长度)
    read_name = read.query_name
    '''
    INTER
    1.通过mapq过滤
    2.对所有read(P and S)都通过cigar找INS
    3.generate_candidate_in_a_read对一条read中,相邻500bp之内的INS进行合并
    '''
    if read.mapq >= min_mapq:
        read_align_start = read.reference_start
        read_align_end = read.reference_end
        #确定S和H clip的长度
        shift_ins = 0
        shift_del = 0
        soft_left = 0
        hard_left = 0
        soft_right = 0
        hard_right = 0
        if read.cigar[0][0] == 4:
            soft_left = read.cigar[0][1]
        if read.cigar[0][0] == 5:
            hard_left = read.cigar[0][1]
        for cigar_pair in read.cigar:
            # DEL
            if cigar_pair[0] in [0, 7, 8]:
                shift_del += cigar_pair[1]
            if cigar_pair[0] == 2 and cigar_pair[1] < sig_min_cigar_size:
                shift_del += cigar_pair[1]
            if cigar_pair[0] == 2 and cigar_pair[1] >= sig_min_cigar_size:
                if process_signal == 1 or process_signal == 2:
                    Combine_sig_in_same_read_del.append([read_align_start + shift_del, cigar_pair[1],1])
                    # print(Combine_sig_in_same_read_del)
                    shift_del += cigar_pair[1]
                else:
                    Combine_sig_in_same_read_del.append([read_align_start + shift_del, cigar_pair[1],0])
                    # print(Combine_sig_in_same_read_del)
                    shift_del += cigar_pair[1]


                
            #ins处理的时候，也需要加上read来源
            #INS
            if cigar_pair[0] in [0, 2, 7, 8]:
                shift_ins += cigar_pair[1]
            if cigar_pair[0] == 1 and cigar_pair[1] >=  sig_min_cigar_size:
                if process_signal == 1 or process_signal == 2:
                # print(Combine_sig_in_same_read_del)
                    Combine_sig_in_same_read_ins.append([read_align_start+shift_ins, cigar_pair[1],1])
                # print(Combine_sig_in_same_read_ins)
                else:
                    Combine_sig_in_same_read_ins.append([read_align_start+shift_ins, cigar_pair[1],0])

        if read.cigar[-1][0] == 4:
            soft_right= read.cigar[-1][1]
        if read.cigar[-1][0] == 5:
            hard_right = read.cigat[-1][1]
        #SH clip只会存在一种，因此全部使用S，而不是H形式的clip
        if hard_left != 0 :
            soft_left = hard_left
        if hard_right != 0:
            soft_right = hard_right   
    # for ele in Combine_sig_in_same_read_ins:
    #     if int(ele[0]) > 21263000 and int(ele[0] < 21264000):
    #         print(Combine_sig_in_same_read_ins)
    #         print(read_name)
    
    #inter
    generate_candidate_in_a_read(chr_name,Combine_sig_in_same_read_ins,candidate,read_name, combine_min_size, 'INS')
    
    # if Combine_sig_in_same_read_del != []:
    #     print(Combine_sig_in_same_read_del)
    generate_candidate_in_a_read(chr_name,Combine_sig_in_same_read_del,candidate,read_name, combine_min_size, 'DEL')
    '''
    INTRA
    1.找出所有的P_read(区分正负链)[read_start,read_end,ref_start,ref_end,chr,strand]
    2.根据tag信息,找S_read
    3.organize_supple_info处理S_read的信息[read_start,read_end,ref_start,ref_end,chr,strand]
    4.generate_candidate_intra_read,通过相邻read之间的相邻关系,确定candiadte
    '''
    if process_signal == 1 or process_signal == 2:
        if read.mapq >= min_mapq:
            if process_signal == 1: #正链
                primary_info=[soft_left,read.query_length-soft_right,read_align_start,read_align_end,\
            chr_name,dic_starnd[process_signal]]
            else: #负链
                primary_info=[soft_right,read.query_length-soft_left,read_align_start,read_align_end,\
            chr_name,dic_starnd[process_signal]]
        else:
            primary_info = []
        Tags = read.get_tags()
        for tag in Tags:
            if tag[0] == 'SA':
                supple_info = tag[1].split(';')[:-1]
                # print(supple_info)
                #intra
                organize_supple_info(primary_info,supple_info,read.query_length,read_name,\
                candidate, max_split_parts, chase_ins_min_size, chase_ins_max_size)
       

def single_pipe(sam_path,task, min_mapq, sig_min_cigar_size, max_split_parts, chase_ins_min_size,\
    chase_ins_max_size, combine_min_size):
    candidate = dict()
    candidate['INS'] = dict()
    candidate['DEL'] = dict()
    reads_list = list()
    samfile = pysam.AlignmentFile(sam_path)
    for read in samfile.fetch(task[0],task[1],task[2]):
        parse_read(task[0],read,candidate, min_mapq, sig_min_cigar_size, max_split_parts,\
        chase_ins_min_size, chase_ins_max_size, combine_min_size)
        if read.reference_start < task[1]:
            continue
        is_primary = 0
        if read.flag in [0, 16]:
            is_primary = 1
        pos_start = read.reference_start # 0-based
        pos_end = read.reference_end
        reads_list.append([pos_start, pos_end, is_primary, read.query_name])

            # print(read_coverage[int(pos_start/100000)])
    # if len(candidate['INS']) != 0 :
    #     print(candidate)
    samfile.close()
    # print(read_coverage)

    '''
    1.生成进程区间的bed文件
    2.将candidate中的信息存到bed文件中 [type,chr,ins_pos,ins_len,read_id]
    '''
    output = "signatures/_%s_%d_%d.bed"%(task[0], task[1], task[2])
    reads_output = "signatures/_%s_%d_%d.reads.bed"%(task[0], task[1], task[2])
    file = open(output, 'w')
    reads_file = open(reads_output, 'w')
    #INS后续是4
    for ele in reads_list:
        reads_file.write("%s\t%d\t%d\t%d\t%s\n"%(task[0], ele[0], ele[1], ele[2], ele[3]))
        # print(ele)
    for svtype in ['INS','DEL']:
        for chr in candidate[svtype]:
            for ele in candidate[svtype][chr]:
                if len(ele) == 4:
                    file.write("%s\t%s\t%d\t%d\t%s\t%d\n"%(svtype, chr, ele[0], ele[1], ele[2], ele[3]))
                    
    file.close()
    reads_file.close()
    print("Finished %s:%d-%d."%(task[0], task[1], task[2]))	


def multi_run_wrapper(args):
	return single_pipe(*args)

    
def solve_bam(batch, max_cluster_bias_INS, min_support, min_size, bam_path,tumor_or_normal,\
    lenth_ratio, min_mapq, sig_min_cigar_size, max_split_parts, chase_ins_min_size, chase_ins_max_size,\
    combine_min_size, threshold_gloab, detailed_length_ratio):
    # 使用pysam将read读进程序，并保存
    normal_sam_file = pysam.AlignmentFile(bam_path,"rb")
    # tumor_sam_file = pysam.AlignmentFile(sys.argv[2],"rb")
    chr_number = len(normal_sam_file.get_index_statistics())
    #区块划分基因组
    task_list = list()
    # read_coverage = dict()
    ref_name_list = list()
    reads_info_dict = dict() # reads_dict_pool["chr1"] = [start, end, is_primary, read_name]

    for chr in normal_sam_file.get_index_statistics():
        ref_name_list.append(chr[0])
        #区域划分基因组，存储到task_list中e
        interval = int(normal_sam_file.get_reference_length(chr[0])/batch)
        for i in range(0,interval+1):
            if (i+1) * batch > normal_sam_file.get_reference_length(chr[0]):
                task_list.append([chr[0],i*batch,normal_sam_file.get_reference_length(chr[0])])
            else:
                task_list.append([chr[0],i*batch,(i+1)*batch])
    
    analysis_pools = Pool(processes=24)
    os.mkdir("%ssignatures"%"./")
    for task in task_list:
        para = [(bam_path, task, min_mapq, sig_min_cigar_size, max_split_parts, chase_ins_min_size,\
            chase_ins_max_size, combine_min_size)]
        
        # 使用多进程保存read_list
        # if  not in reads_info_dict:
        #     reads_info_dict[task[0]] = list()
        # reads_dict_pool[task[0]].append(analysis_pools.map_async(multi_run_wrapper, para))
        
        analysis_pools.map_async(multi_run_wrapper, para)
    analysis_pools.close()
    analysis_pools.join()


    # reads_info_dict = dict()
    # 把进程里的所有candidate放到最终处理的列表中
    # for chr in reads_dict_pool:
    #     reads_info_dict[chr] = list()
    #     for item in reads_dict_pool[chr]:
    #         try:
    #             reads_info_dict[chr].extend(item.get()[0])
    #         except:
    #             pass
    # for chr in reads_info_dict:
    #     print(chr + '\t' + str(len(reads_info_dict[chr])))

    print("Rebuilding signatures of structural variants.")
    analysis_pools = Pool(processes=24)
    cmd_ins = ("cat %ssignatures/*.bed | grep -w INS | sort -u -T %s | sort -k 2,2 -k 3,3n -T %s > %s%sINS.sigs"%("./", "./", "./", "./", tumor_or_normal))
    cmd_del = ("cat %ssignatures/*.bed | grep -w DEL | sort -u -T %s | sort -k 2,2 -k 3,3n -T %s > %s%sDEL.sigs"%("./", "./", "./", "./", tumor_or_normal))
    cmd_reads = ("cat %ssignatures/*.reads.bed > %s%sreads.sigs"%("./", "./", tumor_or_normal))
    # for i in [cmd_ins,cmd_del]:
    for i in [cmd_ins,cmd_del,cmd_reads]:
        analysis_pools.map_async(os.system, (i,))
    analysis_pools.close()
    analysis_pools.join()


    # result_sv = list()
    # # '''
    # # 聚类
    # # '''
    # value_sv = load_sigs_chr('./',tumor_or_normal)

    # reads_info_dict = dict()
    # for chr in value_sv['DEL']:
    #     reads_info_dict[chr] = list()
    # readsfile = open("%s%sreads.sigs"%("./", tumor_or_normal), 'r')
    
    # for line in readsfile:
    #     seq = line.strip().split('\t')
    #     chr = seq[0]
    #     if chr not in reads_info_dict:
    #         reads_info_dict[chr] = list()
    #     reads_info_dict[chr].append([int(seq[1]), int(seq[2]), int(seq[3]), seq[4]])
    
    # readsfile.close()

    # analysis_pools = Pool(processes=8)
    # # # for chr in value_sv['INS']:
    # # #     para = [("%s%s%s.sigs"%('./', tumor_or_normal,'INS'),chr,max_cluster_bias_INS,min_support,min_size,\
    # # #         lenth_ratio, threshold_gloab, detailed_length_ratio)]
    # # #     result_sv.append(analysis_pools.map_async(run_ins,para))
    # # #     print("Finish %s:%s"%(chr,'INS'))
    # for chr in value_sv['DEL']:
    #     para = [("%s%s%s.sigs"%('./', tumor_or_normal,'DEL'), chr, min_support,lenth_ratio, reads_info_dict[chr],bam_path)]
    #     result_sv.append(analysis_pools.map_async(run_del,para))
    #     print("Finish %s:%s"%(chr,'DEL'))

        
    # analysis_pools.close()
    # analysis_pools.join()
     
    # # print("Writing to your output file.")

    # #把进程里的所有candidate放到最终处理的列表中
    # semi_result = list()
    # for res in result_sv:
    #     try:
    #         semi_result += res.get()[0]
    #     except:
    #         pass
    # #按照染色体和位置排序
    # semi_result = sorted(semi_result, key = lambda x:(x[0], int(x[2])))
    
    # # print(semi_result)

    # # print(semi_result)
    # # print("Loading reference genome...")
    # # ref_g = SeqIO.to_dict(SeqIO.parse(sys.argv[3], "fasta"))
    os.system("rm -r %ssignatures "%("./"))
    
    # return semi_result

def main_ctrl(args, argv):
    # batch = 10000000
    # max_cluster_bias_INS = 50
    # min_support = 3（49x）
    # min_size_INS = 20
    # lenth_ratio = 0
    print("resolve normal bam")
    # normal_sv_list = solve_bam(10000000, 50, 3, 20, sys.argv[1], "normal",0)
    normal_sv_list = solve_bam(args.batches, args.max_cluster_bias_INS, args.normal_min_support, \
        args.normal_min_size, sys.argv[1], "normal", args.normal_length_ratio_flag, args.min_mapq, \
        args.sig_min_cigar_size, args.max_split_parts, args.chase_ins_min_size, args.chase_ins_max_size,\
        args.combine_min_size, args.threshold_gloab, args.detailed_length_ratio)
    
    file1 = open("normal_sv","w")
    for ele in normal_sv_list:
        for i in range(0,11):
            file1.write(str(ele[i])+'\t')
        file1.write(str(ele[11])+'\n')
    file1.close()
    

    # batch = 10000000
    # max_cluster_bias_INS = 50
    # min_support = 20（45x）
    # min_size = 50
    # len_ratio
    print("resolve tumor bam")
    # tumor_sv_list = solve_bam(10000000, 50, 20, 50, sys.argv[2], "tumor",1)
    tumor_sv_list = solve_bam(args.batches, args.max_cluster_bias_INS, args.tumor_min_support, \
        args.tumor_min_size, sys.argv[2], "tumor", args.tumor_length_ratio_flag, args.min_mapq, \
        args.sig_min_cigar_size, args.max_split_parts, args.chase_ins_min_size, args.chase_ins_max_size,\
        args.combine_min_size, args.threshold_gloab, args.detailed_length_ratio)

    # file1 = open("normal_sv","w")
    # for ele in normal_sv_list:
    #     for i in range(0,11):
    #         file1.write(str(ele[i])+'\t')
    #     file1.write(str(ele[11])+'\n')
    # file1.close()

    file2 = open("tumor_sv","w")
    for ele in tumor_sv_list:
        for i in range(0,11):
            file2.write(str(ele[i])+'\t')
        file2.write(str(ele[11])+'\n')
    file2.close()

def run(argv):
    args = parseArgs(argv)
    main_ctrl(args, argv)

if __name__=='__main__':
    run(sys.argv[1:])