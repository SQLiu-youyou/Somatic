from audioop import lin2ulaw
import sqlite3
import sys

from numpy import size
from genotype import *


def new_count_coverage(reads_list, search_start, search_end, up_bound, itround,read_id_num):
    status = 0
    iteration = 0
    primary_num = 0
    read_count = set()
    iteration_set = set()


    for read in reads_list:
        if (read[0] <= search_start and read[1] > search_start) or (search_start < read[0] < search_end):
            iteration += 1
            iteration_set.add(read[3])
            if read[2] == 0:
                continue
            primary_num += 1
            if read[0] < search_start and read[1] > search_end:
                read_count.add(read[3])
                if len(read_count) >= up_bound:
                    status = 1
                    break
            if iteration >= itround:
                if float(primary_num / iteration) <= 0.2:
                    status = 1
                else:
                    status = -1
                break
        
    if len(iteration_set) > 10*read_id_num:

        return status,iteration_set,1
    else:
        return status,iteration_set,0


def call_new_gt(reads_list, search_threshold, read_id_num, max_cluster_bias, gt_round):
    querydata = set()
    search_start = max(int(search_threshold) - max_cluster_bias, 0)
    # search_end = min(int(search_threshold) + max_cluster_bias, bamfile.get_reference_length(chr))
    search_end = int(search_threshold) + max_cluster_bias
    # if search_threshold == 118516185:
    #     print(search_start,search_end)
    up_bound = threshold_ref_count(read_id_num)
    
    status,querydata,flag = new_count_coverage(reads_list, 
                            search_start, 
                            search_end,  
                            up_bound, 
                            gt_round, read_id_num)


    if status == -1:
        # if search_threshold == 14962957:
        #     print("hi")
        DR = '.'
        GT = "./."
        GL = ".,.,."
        GQ = "."
        QUAL = "."

    
    else:
        # print(search_threshold)
        # print(querydata)
        if len(querydata) == 0 or len(querydata) < read_id_num:
            GT = '0/0' 
        else:
            DR = len(querydata) - read_id_num

            GT, GL, GQ, QUAL = cal_GL(DR, read_id_num)

    return GT,flag    

if __name__=='__main__':
    
    reads_info_dict = dict()
    
    
    readsfile = open("tumorreads.sigs",'r')
    
    for line in readsfile:
        seq = line.strip().split('\t')
        chr = seq[0]
        if chr not in reads_info_dict:
            reads_info_dict[chr] = list()
        reads_info_dict[chr].append([int(seq[1]), int(seq[2]), int(seq[3]), seq[4]])
    
    readsfile.close()

    tumor_ans = open(sys.argv[1],'r')
    for line in tumor_ans:
        seq = line.strip().split('\t')
        chr = line.strip().split('\t')[0]
        support = int(line.strip().split('\t')[4])
        pos = int(line.strip().split('\t')[2])
        old_gt = line.strip().split('\t')[8]
        if support <= 5:
            new_gt,flag = call_new_gt(reads_info_dict[chr],pos,support,1000,500)
            # print("new_gt")
            # print(new_gt)
            if new_gt != old_gt:
                if old_gt == '1/1':
                    continue
                elif old_gt == '0/1':
                    if flag:
                        continue
                    else:
                        print('\t'.join(seq[0:11])+'\t'+'low')
            else:
                print('\t'.join(seq[0:11])+'\t'+'low')
        else:
            print(line.strip())


            
            # def call_gt(reads_list, search_threshold, chr, read_id_list, max_cluster_bias, gt_round, gt_primary):

