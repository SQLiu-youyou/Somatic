from re import S, X
from tkinter import E
from turtle import right
import numpy as np
from sortedcontainers import SortedDict, SortedList
from genotype import *

def call_gt(reads_list, search_threshold, chr, read_id_list, max_cluster_bias, gt_round, gt_primary):
    # print("in")
    # return '1/1'
    # status = 2 
    # gt_round = 500
    # max_cluster_bias = 1000
    querydata = set()
    search_start = max(int(search_threshold) - max_cluster_bias, 0)
    search_end = int(search_threshold) + max_cluster_bias

    # if search_threshold == 89700301:
    #     print(search_start,search_end)
    up_bound = threshold_ref_count(len(read_id_list))

    status,querydata = count_coverage(reads_list, 
                            search_start, 
                            search_end,  
                            up_bound, 
                            gt_round, read_id_list, gt_primary)

    if search_threshold == 89700301:
        print(status)
    if status == -1:
        DR = '.'
        GT = "./."
        GL = ".,.,."
        GQ = "."
        QUAL = "."

    
    else:
        DR = 0
        for query in querydata:
            if query not in read_id_list:
                DR += 1
        if len(querydata) == 0:
            GT = './.'
        else:
            # if search_start == 42372687:
            #     print(querydata)
            #     print(search_threshold)
            #     print(DR,len(read_id_list))
            GT, GL, GQ, QUAL = cal_GL(DR, len(read_id_list))

    # GT = '1/1'
    return GT
# def call_gt(reads_list, search_threshold, chr, read_id_list, max_cluster_bias, gt_round):
#     # print("in")
#     GT = '1/1'
#     return GT

def generate_final_ins(cluster_list,candidate_single_SV,chr,min_size):
    # max_support = 500
    ins_start = int(np.mean([i[0] for i in cluster_list]))
    ins_len = int(np.mean([i[1] for i in cluster_list]))
    if ins_len >= min_size :
        CIPOS = cal_CIPOS(np.std([i[0] for i in cluster_list]), len([i[0] for i in cluster_list]))
        CILEN = cal_CIPOS(np.std([i[1] for i in cluster_list]), len([i[1] for i in cluster_list]))
        DR = '.'
        GT = './.'
        GL = '.,.,.'
        GQ = "."
        QUAL = "."
        candidate_single_SV.append([chr, "INS", str(int(ins_start)), str(int(ins_len)), str(len(cluster_list)), \
            str(CIPOS), str(CILEN), str(DR), str(GT), str(GL), str(GQ), str(QUAL)])
    # print(candidate_single_SV)

def generate_final_del(cluster_list,candidate_single_SV,chr, reads_info_dict):
    new_max_cluster_bias = 1000
    gt_round = 500
    gt_thres = 0.7
    del_start = int(np.mean([i[0] for i in cluster_list]))
    del_len = int(np.mean([i[1] for i in cluster_list]))
    # print("cluster_list")
    # print(cluster_list)
    if del_len >= 20 :
        # print("yes")
        CIPOS = cal_CIPOS(np.std([i[0] for i in cluster_list]), len([i[0] for i in cluster_list]))
        CILEN = cal_CIPOS(np.std([i[1] for i in cluster_list]), len([i[1] for i in cluster_list]))
        DR = '.'
        # GT = './.'
        GL = '.,.,.'
        GQ = "."
        QUAL = "."
        read_id_list = list()
        for ele in cluster_list:
            read_id_list.append(ele[2])
        if np.sum(i[3] for i in cluster_list)/len(cluster_list) < 0.7:
            gt_primary = 0
        else:
            gt_primary = 1
        GT = call_gt(reads_info_dict, del_start, chr, read_id_list, new_max_cluster_bias, gt_round, gt_primary)
        # DR, GT, GL, GQ, QUAL = call_gt(reads_info_dict, del_start, read_id_list, new_max_cluster_bias, gt_round)
            
        candidate_single_SV.append([chr, "DEL", str(int(del_start)), str(int(del_len)), str(len(cluster_list)), \
            str(CIPOS), str(CILEN), str(DR), str(GT), str(GL), str(GQ), str(QUAL)])
        
        # flag = 0
        # for ele in candidate_single_SV:
        #     if ele[2] == '2893750' and flag == 0:
        #         print(cluster_list)
        #         print("length")
        #         print(len(cluster_list))
        #         flag = 1
                

def cluster_through_len_ins(cluster_list,chr,min_support,candidate_single_SV,min_size,lenth_ratio, threshold_gloab,\
    detailed_length_ratio):
    # threshold_gloab = 0.2
    # detailed_length_ratio = 0.8
    '''
    第二轮按INS_len聚类
    1.长度相邻的两条read,len差值小于30%的平均长度,聚成一类
    2.类中元素个数>min_support,才输出最终ins结果
    3.generate_final_ins 输出最终的ins信息
    '''    
    # 对于同一read 里的变异，只保留更长的那一个
    reserved_read = dict()
    for ele in cluster_list:
        if ele[2] not in reserved_read:
            reserved_read[ele[2]] = ele
        else:
            if ele[1] > reserved_read[ele[2]][1]:
                reserved_read[ele[2]] = ele
    if len(reserved_read) < min_support:
        return
    
    
    '''
    ***********************************************************************
        #SortedList 		#standard_len		    #global_len		
    -------------------------------------------------------------
        按ins_pos排序的list   目前待比较的标准ins_len 	cluster中所有len的值
    ***********************************************************************
    '''
    # ins_pos ins_len read_id
    SortedList = sorted(list(reserved_read.values()), key = lambda x:x[1])
    # SortedList = sorted(cluster_list, key=lambda x:x[1])
    # print(SortedList)
    standard_len = SortedList[0][1]
    global_len = [i[1] for i in SortedList]
    # 按ins_len聚类是，允许的len差值的最大值
    DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP = threshold_gloab * np.mean(global_len)
    
    # for ele in cluster_list:
    #     if ele[0]>= 13468446 and ele[0] <= 13468483:
    #         print(SortedList)
    #         print(len(SortedList))
    #         print(DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP)
    to_SV_list = list()
    to_SV_list.append(SortedList[0])
    for ele in SortedList[1:]:
        #根据len判断，是一个cluster里的，放进list里
        if ele[1] - standard_len < DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP:
            to_SV_list.append(ele)
            standard_len = ele[1]
        #不是一个堆了，判断当前是否足以支持一个变异，并重新申请cluster_list
        else:
            if len(to_SV_list) >= 0.9 * min_support:
                # if ele[0] == 129771776:
                #     print(to_SV_list)
                if lenth_ratio:
                    if len(to_SV_list) / len(cluster_list) >= detailed_length_ratio:
                        generate_final_ins(to_SV_list,candidate_single_SV,chr,min_size)
                else:
                    generate_final_ins(to_SV_list,candidate_single_SV,chr,min_size)
            to_SV_list = list()
            standard_len = ele[1]
            to_SV_list.append(ele)
    #对最后一个堆的判定
    if len(to_SV_list) >= 0.9 * min_support:
        if lenth_ratio:
            if len(to_SV_list) / len(cluster_list) >= detailed_length_ratio:
                generate_final_ins(to_SV_list,candidate_single_SV,chr,min_size)
        else:
            generate_final_ins(to_SV_list,candidate_single_SV,chr,min_size)

def cluster_through_len_del(cluster_list,chr,candidate_single_SV, min_support, length_ratio, \
    reads_info_dict):
    # print("in")
    #clr:0.3 ont 0.1
    threshold_gloab = 0.1
    # min_support = 14
    detailed_length_ratio = 0.8
    
    # for ele in cluster_list:
    #     if ele[0]>= 72446676 and ele[0] <= 72446856:
    #         print(cluster_list)
    '''
    第二轮按DEL_len聚类
    1.长度相邻的两条read,len差值小于50%的平均长度,聚成一类
    2.类中元素个数>min_support,才输出最终ins结果
    3.generate_final_ins 输出最终的ins信息
    '''    
    #对于同一read 里的变异，只保留更长的那一个
    reserved_read = dict()
    for ele in cluster_list:
        if ele[2] not in reserved_read:
            reserved_read[ele[2]] = ele
        else:
            if ele[1] > reserved_read[ele[2]][1]:
                reserved_read[ele[2]] = ele
    if len(reserved_read) < min_support:
        return
    
    '''
    ***********************************************************************
        #SortedList 		#standard_len		    #global_len		
    -------------------------------------------------------------
        按ins_pos排序的list   目前待比较的标准ins_len 	cluster中所有len的值
    ***********************************************************************
    '''
    # ins_pos ins_len read_id
    # SortedList = sorted(list(reserved_read.values()), key = lambda x:x[1])
    SortedList = sorted(cluster_list, key = lambda x:x[1])
    # for ele in cluster_list:
    #     if ele[0]>= 13468446 and ele[0] <= 13468483:
    #         print(SortedList)
    #         print(len(SortedList))
    # SortedList = sorted(cluster_list, key = lambda x:x[1])
    
    # print(SortedList)
    standard_len = SortedList[0][1]
    # global_len = [i[1] for i in SortedList]
    # 按ins_len聚类是，允许的len差值的最大值
    # DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP = threshold_gloab * np.mean(global_len)
    DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP = threshold_gloab * standard_len
    # print(DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP)

    to_SV_list = list()
    to_SV_list.append(SortedList[0])
    # for ele in cluster_list:
    #     if ele[0]>= 134948801 and ele[0] <= 134952894:
    #         print(SortedList)
    #         print(len(SortedList))
    #         print(DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP)

    for ele in SortedList[1:]:
        #根据len判断，是一个cluster里的，放进list里
        if ele[1] - standard_len < DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP:
            
            to_SV_list.append(ele)
            # if ele[1] == 2426:
            #     print(to_SV_list)
            standard_len = ele[1]
            DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP = np.mean([i[1] for i in to_SV_list]) * threshold_gloab
        #不是一个堆了，判断当前是否足以支持一个变异，并重新申请cluster_list
        else:
            if len(to_SV_list) >= min_support:
                # if to_SV_list[0][0] == 134951781:
                # print("to_sv_list")
                # print(to_SV_list)
                #     print(DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP)
                #     print(ele[1])
                # if length_ratio:
                #     if len(to_SV_list) / len(cluster_list) >= detailed_length_ratio:
                #         generate_final_del(to_SV_list,candidate_single_SV,chr)
                # else:
                # print("tp_sv_list")
                # print(to_SV_list)
                generate_final_del(to_SV_list,candidate_single_SV,chr, reads_info_dict)
                
                # generate_final_del(to_SV_list,candidate_single_SV,chr)
            to_SV_list = list()
            standard_len = ele[1]
            
            DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP = threshold_gloab * standard_len
            to_SV_list.append(ele)
            # print(to_SV_list)
    #对最后一个堆的判定
    if len(to_SV_list) >= min_support:
        # if length_ratio:
        #     if len(to_SV_list) / len(cluster_list) >= detailed_length_ratio:
        #         generate_final_del(to_SV_list,candidate_single_SV,chr)
        # else:
        #     generate_final_del(to_SV_list,candidate_single_SV,chr)
        # if to_SV_list[0][0] == 134951781:
        #     print("to_sv_list")
        #     print(to_SV_list)
        # print("sv_list")
        # print(to_SV_list)
        generate_final_del(to_SV_list,candidate_single_SV,chr, reads_info_dict)
        # generate_final_del(to_SV_list,candidate_single_SV,chr)
    
    
    # for ele in SortedList[1:]:
    #     #根据len判断，是一个cluster里的，放进list里
    #     if ele[1] - standard_len < DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP:
    #         to_SV_list.append(ele)
    #         standard_len = ele[1]
    #     #不是一个堆了，判断当前是否足以支持一个变异，并重新申请cluster_list
    #     else:
    #         if len(to_SV_list) >= min_support:
    #             # if to_SV_list[0][0] == 1702496:
    #             #     print("to_sv_list")
    #             #     print(to_SV_list)
    #             if length_ratio:
    #                 if len(to_SV_list) / len(cluster_list) >= detailed_length_ratio:
    #                     generate_final_del(to_SV_list,candidate_single_SV,chr)
    #             else:
    #                 generate_final_del(to_SV_list,candidate_single_SV,chr)
    #             # print(to_SV_list)
    #             # generate_final_del(to_SV_list,candidate_single_SV,chr)
    #         to_SV_list = list()
    #         standard_len = ele[1]
    #         to_SV_list.append(ele)
    # #对最后一个堆的判定
    # if len(to_SV_list) >= min_support:
    #     if length_ratio:
    #         if len(to_SV_list) / len(cluster_list) >= detailed_length_ratio:
    #             generate_final_del(to_SV_list,candidate_single_SV,chr)
    #     else:
    #         generate_final_del(to_SV_list,candidate_single_SV,chr)
    #     # print(to_SV_list)
    #     # generate_final_del(to_SV_list,candidate_single_SV,chr)

def resolution_INS(sigs_path,chr,max_cluster_bias,min_support,min_size, lenth_ratio, threshold_gloab, detailed_length_ratio):
    #首先按照起始位点坐标，sigs聚成一类
    '''
    第一轮按INS_pos聚类
    1.将相近的ins sigs聚成一个cluster(max_cluster_bias = 100)
    2.cluster中的元素个数大于min_support,才进入第二轮聚类
    3.cluster_through_len进行二轮聚类
    '''
    
    candidate_single_SV = list()
    cluster_list = list()
    file = open(sigs_path,'r')
    for line in file:
        seq = line.strip('\n').split('\t')
        if seq[1] != chr:
            continue
        ins_pos = int(seq[2])
        ins_len = int(seq[3])
        read_id = seq[4]
        #初始的第一个sig放进cluster中进行处理
        if len(cluster_list) == 0:
            cluster_list.append([ins_pos, ins_len, read_id])
            # print(cluster_list)
            continue
        #在一定位置阈值内的均放在cluster当中
        if ins_pos - cluster_list[-1][0] <= max_cluster_bias:
            cluster_list.append([ins_pos, ins_len, read_id])
        #否则 判断cluster中的条数，是否大于-s
        else:
            if len(cluster_list) >= min_support:
                # if ins_pos == 2584331:
                #     print(cluster_list)
                #     print(len(cluster_list))
                cluster_through_len_ins(cluster_list,chr,min_support,candidate_single_SV,min_size,lenth_ratio,\
                    threshold_gloab, detailed_length_ratio)
            #无论继续处理与否，都需要将当前sig放入新一轮的处理中
            cluster_list = list()
            cluster_list.append([ins_pos, ins_len, read_id])
            # print(cluster_list)
    #对最后一个cluster的处理
    if len(cluster_list) >= min_support:
        # print(cluster_list)
        cluster_through_len_ins(cluster_list,chr,min_support,candidate_single_SV,min_size,lenth_ratio,\
            threshold_gloab, detailed_length_ratio)
    file.close()
    return candidate_single_SV

def resolution_DEL(sigs_path,chr,min_support,length_ratio, reads_info_dict):
    #clr 200 ont 100
    max_cluster_bias = 100
    # min_support = 14
    # print(reads_info_dict)
    # print("ok")
    #首先按照起始位点坐标，sigs聚成一类
    '''
    第一轮按DEL_pos聚类
    1.将相近的ins sigs聚成一个cluster(max_cluster_bias = 200)
    2.cluster中的元素个数大于min_support,才进入第二轮聚类
    3.cluster_through_len进行二轮聚类
    '''
    
    candidate_single_SV = list()
    cluster_list = list()
    file = open(sigs_path,'r')
    for line in file:
        seq = line.strip('\n').split('\t')
        if seq[1] != chr:
            continue
        ins_pos = int(seq[2])
        ins_len = int(seq[3])
        read_id = seq[4]
        sig_origin = int(seq[5])
        # print(ins_pos,chr)
        #初始的第一个sig放进cluster中进行处理
        if len(cluster_list) == 0:
            cluster_list.append([ins_pos, ins_len, read_id, sig_origin])
            # print(cluster_list)
            continue
        #在一定位置阈值内的均放在cluster当中
        if ins_pos - cluster_list[-1][0] <= max_cluster_bias:
            cluster_list.append([ins_pos, ins_len, read_id, sig_origin])
        #否则 判断cluster中的条数，是否大于-s
        else:
            if len(cluster_list) >= min_support:
                # if ins_pos == 134948801:
                #     print(cluster_list)
                #     print(len(cluster_list))
                # print("cluster")
                # print(cluster_list)
                # if cluster_list[0][0] == 134948801:
                #     print(cluster_list)
                #     print(len(cluster_list))
                # print(cluster_list)
                cluster_through_len_del(cluster_list,chr,candidate_single_SV,min_support,\
                    length_ratio, reads_info_dict)
            #无论继续处理与否，都需要将当前sig放入新一轮的处理中
            cluster_list = list()
            cluster_list.append([ins_pos, ins_len, read_id, sig_origin])
            # print(cluster_list)
    #对最后一个cluster的处理
    if len(cluster_list) >= min_support:
        # print(cluster_list)
        cluster_through_len_del(cluster_list,chr,candidate_single_SV, min_support,\
            length_ratio, reads_info_dict)
    file.close()
    return candidate_single_SV

def run_ins(args):
   '''
    对sigs进行聚类(两轮)
    1.resolution_INS 对ins_pos聚类
    2.cluster_through_len对ins_len聚类
    '''
   return resolution_INS(*args)

def run_del(args):
    # print("in")
    return resolution_DEL(*args)

