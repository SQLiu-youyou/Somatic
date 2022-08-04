from hashlib import new
from re import S, X
from tkinter import E
from turtle import right
import numpy as np
from sortedcontainers import SortedDict, SortedList
from genotype import *


# def call_gt(bam_path, search_threshold, chr, read_id_list, max_cluster_bias, gt_round):
# 	import pysam
# 	# print(search_threshold)
# 	querydata = set()
# 	bamfile = pysam.AlignmentFile(bam_path)
# 	search_start = max(int(search_threshold) - max_cluster_bias, 0)
# 	search_end = min(int(search_threshold) + max_cluster_bias, bamfile.get_reference_length(chr))

# 	up_bound = threshold_ref_count(len(read_id_list))

# 	status = count_coverage(chr, 
# 							search_start, 
# 							search_end, 
# 							bamfile, 
# 							querydata, 
# 							up_bound, 
# 							gt_round)
# 	bamfile.close()

# 	if status == -1:
# 		DR = '.'
# 		GT = "./."
# 		GL = ".,.,."
# 		GQ = "."
# 		QUAL = "."

# 	# elif status == 1:
# 	# 	pass
# 	else:
# 		DR = 0
# 		for query in querydata:
# 			if query not in read_id_list:
# 				DR += 1
# 		GT, GL, GQ, QUAL = cal_GL(DR, len(read_id_list))
	
# 	# if search_threshold == 30353767 :
# 	# 	print("here")
# 	# 	print(len(read_id_list))
# 	# 	print(DR)

# 	return GT
#     # return len(read_id_list), DR, GT, GL, GQ, QUAL


# def call_gt(reads_list, search_threshold, chr, read_id_list, max_cluster_bias, gt_round, gt_primary):
#     # print("in")
#     # return '1/1'
#     # status = 2 
#     # gt_round = 500
#     # max_cluster_bias = 1000
#     # import pysam
#     # bamfile = pysam.AlignmentFile(bam_path)
#     querydata = set()
#     search_start = max(int(search_threshold) - max_cluster_bias, 0)
#     # search_end = min(int(search_threshold) + max_cluster_bias, bamfile.get_reference_length(chr))
#     search_end = int(search_threshold) + max_cluster_bias
#     # if search_threshold == 118516185:
#     #     print(search_start,search_end)
#     up_bound = threshold_ref_count(len(read_id_list))
    
#     status,querydata = count_coverage(reads_list, 
#                             search_start, 
#                             search_end,  
#                             up_bound, 
#                             gt_round, read_id_list, gt_primary)

#     # if search_threshold == 89700301:
#     #     print(status)
#     # if search_threshold == 118516185:
#     #     print(status)
#     #     print(len(querydata))
#     #     print(len(read_id_list))
#     # if len(read_id_list)==32:
#     #     print(search_threshold)
#     #     print("status")
#     #     print(status)

#     if status == -1:
#         # if search_threshold == 14962957:
#         #     print("hi")
#         DR = '.'
#         GT = "./."
#         GL = ".,.,."
#         GQ = "."
#         QUAL = "."

    
#     else:
#         DR = 0
#         for query in querydata:
#             if query not in read_id_list:
#                 DR += 1
#         # if search_threshold == 33759248:
#         #     print("DR")
#         #     print(DR)
#         #     print(len(read_id_list))
        
#         if len(querydata) == 0:
#             GT = './.'
#         # else:
#             # if search_start == 42372687:
#             #     print(querydata)
#             #     print(search_threshold)
#             #     print(DR,len(read_id_list))
#         GT, GL, GQ, QUAL = cal_GL(DR, len(read_id_list))
    
#     # if search_threshold == 33759248:
#     #         print("DR")
#     #         print(DR)
#     #         print(len(read_id_list))
    
#     # if len(read_id_list) <= 5:
#     #     new_status,new_querydata = count_coverage(reads_list, 
#     #                         search_start, 
#     #                         search_end,  
#     #                         up_bound, 
#     #                         gt_round, read_id_list, gt_primary)

#     # GT = '1/1'
#     return DR, GT, GL, GQ, QUAL
#     # return GT
# # def call_gt(reads_list, search_threshold, chr, read_id_list, max_cluster_bias, gt_round):
# #     # print("in")
# #     GT = '1/1'
# #     return GT

def generate_final_ins(cluster_list,candidate_single_SV,chr, reads_info_dict):
    # max_support = 500
    new_max_cluster_bias = 1000
    gt_round = 500
    gt_thres = 0.7
    ins_start = int(np.mean([i[0] for i in cluster_list]))
    ins_len = int(np.mean([i[1] for i in cluster_list]))
    if ins_len >= 20 :
        CIPOS = cal_CIPOS(np.std([i[0] for i in cluster_list]), len([i[0] for i in cluster_list]))
        CILEN = cal_CIPOS(np.std([i[1] for i in cluster_list]), len([i[1] for i in cluster_list]))
        # DR = '.'
        # GT = './.'
        # GL = '.,.,.'
        # GQ = "."
        # QUAL = "."
        read_id_list = list()
        for ele in cluster_list:
            read_id_list.append(ele[2])
        if np.sum(i[3] for i in cluster_list)/len(cluster_list) < 0.3:
            gt_primary = 0
        else:
            gt_primary = 1
        candidate_single_SV.append([chr, "INS", int(ins_start), int(ins_len), len(cluster_list), \
            str(CIPOS), str(CILEN),read_id_list])
        # candidate_single_SV.append([chr, "INS", str(int(ins_start)), str(int(ins_len)), str(len(cluster_list)), \
        #     str(CIPOS), str(CILEN), str(DR), str(GT), str(GL), str(GQ), str(QUAL)])
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
        # DR = '.'
        # GT = './.'
        # GL = '.,.,.'
        # GQ = "."
        # QUAL = "."
        read_id_list = list()
        for ele in cluster_list:
            read_id_list.append(ele[2])
        if np.sum(i[3] for i in cluster_list)/len(cluster_list) < 0.3:
            gt_primary = 0
        else:
            gt_primary = 1
        # if len(cluster_list) == 17:
        #     print(del_start)
        # DR, GT, GL, GQ, QUAL = call_gt(reads_info_dict, del_start, chr, read_id_list, new_max_cluster_bias, gt_round, gt_primary)
        # GT= call_gt(reads_info_dict, del_start, chr, read_id_list, new_max_cluster_bias, gt_round, gt_primary)

        # candidate_single_SV.append([chr, "DEL", str(int(del_start)), str(int(del_len)), str(len(cluster_list)), \
        #     str(CIPOS), str(CILEN), str(DR), str(GT), str(GL), str(GQ), str(QUAL)])

        candidate_single_SV.append([chr, "DEL", int(del_start), int(del_len), len(cluster_list), \
            str(CIPOS), str(CILEN),read_id_list])
        
        
        # flag = 0
        # for ele in candidate_single_SV:
        #     if ele[2] == '2893750' and flag == 0:
        #         print(cluster_list)
        #         print("length")
        #         print(len(cluster_list))
        #         flag = 1
                

def cluster_through_len_ins(cluster_list,chr,candidate_single_SV, min_support, length_ratio, \
    reads_info_dict,bam_path):
    threshold_gloab = 0.1
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
    # SortedList = sorted(list(reserved_read.values()), key = lambda x:x[1])
    SortedList = sorted(cluster_list, key=lambda x:x[1])
    # print(SortedList)
    standard_len = SortedList[0][1]
    # global_len = [i[1] for i in SortedList]
    # 按ins_len聚类是，允许的len差值的最大值
    # DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP = threshold_gloab * np.mean(global_len)
    DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP = threshold_gloab * standard_len
    
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
            DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP = np.mean([i[1] for i in to_SV_list]) * threshold_gloab
        #不是一个堆了，判断当前是否足以支持一个变异，并重新申请cluster_list
        else:
            if len(to_SV_list) >= min_support:
                # if ele[0] == 129771776:
                #     print(to_SV_list)
                # if lenth_ratio:
                #     if len(to_SV_list) / len(cluster_list) >= detailed_length_ratio:
                #         generate_final_ins(to_SV_list,candidate_single_SV,chr,min_size)
                # else:
                #     generate_final_ins(to_SV_list,candidate_single_SV,chr,min_size)
                generate_final_ins(to_SV_list,candidate_single_SV,chr, reads_info_dict)
            to_SV_list = list()
            standard_len = ele[1]
            
            DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP = threshold_gloab * standard_len
            to_SV_list.append(ele)
    #对最后一个堆的判定
    if len(to_SV_list) >= min_support:
        # if lenth_ratio:
        #     if len(to_SV_list) / len(cluster_list) >= detailed_length_ratio:
        #         generate_final_ins(to_SV_list,candidate_single_SV,chr,min_size)
        # else:
        #     generate_final_ins(to_SV_list,candidate_single_SV,chr,min_size)
        generate_final_ins(to_SV_list,candidate_single_SV,chr, reads_info_dict)

def cluster_through_len_del(cluster_list,chr,candidate_single_SV, min_support, length_ratio, \
    reads_info_dict,bam_path):
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

def resolution_INS(sigs_path,chr,min_support,length_ratio, reads_info_dict,bam_path):
    #首先按照起始位点坐标，sigs聚成一类
    max_cluster_bias = 100
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
        sig_origin = int(seq[5])
        #初始的第一个sig放进cluster中进行处理
        if len(cluster_list) == 0:
            cluster_list.append([ins_pos, ins_len, read_id,sig_origin])
            # print(cluster_list)
            continue
        #在一定位置阈值内的均放在cluster当中
        if ins_pos - cluster_list[-1][0] <= max_cluster_bias:
            cluster_list.append([ins_pos, ins_len, read_id,sig_origin])
        #否则 判断cluster中的条数，是否大于-s
        else:
            if len(cluster_list) >= min_support:
                # if ins_pos == 2584331:
                #     print(cluster_list)
                #     print(len(cluster_list))
                cluster_through_len_ins(cluster_list,chr,candidate_single_SV,min_support,\
                    length_ratio, reads_info_dict,bam_path)
            #无论继续处理与否，都需要将当前sig放入新一轮的处理中
            cluster_list = list()
            cluster_list.append([ins_pos, ins_len, read_id, sig_origin])
            # print(cluster_list)
    #对最后一个cluster的处理
    if len(cluster_list) >= min_support:
        # print(cluster_list)
        cluster_through_len_ins(cluster_list,chr,candidate_single_SV,min_support,\
                    length_ratio, reads_info_dict,bam_path)
    file.close()

    iteration_dict, primary_num_dict, cover_dict = overlap_cover(candidate_single_SV, reads_info_dict) # both key(sv idx), value(set(read id))
    read_id_dict = list()
    for i in candidate_single_SV:
        read_id_dict.append(i[7])
    assign_list = assign_gt(iteration_dict, primary_num_dict, cover_dict,read_id_dict)
    # print(assign_list)
# assign_list.append([len(read_id_dict[idx]), DR, GT, GL, GQ, QUAL])
    # candidate_single_SV.append([chr, "DEL", str(int(del_start)), str(int(del_len)), str(len(cluster_list)), \
        #     str(CIPOS), str(CILEN), str(DR), str(GT), str(GL), str(GQ), str(QUAL)])
    
    final_sv_list = list()
    for i in range(len(assign_list)):

        final_sv_list.append([chr,"INS",str(candidate_single_SV[i][2]),str(candidate_single_SV[i][3]),\
            str(assign_list[i][0]),str(candidate_single_SV[i][5]),str(candidate_single_SV[i][6]),\
            str(assign_list[i][1]),str(assign_list[i][2]),str(assign_list[i][3]),str(assign_list[i][4]),str(assign_list[i][5])])
    return final_sv_list
    
    # return candidate_single_SV

def resolution_DEL(sigs_path,chr,min_support,length_ratio, reads_info_dict,bam_path):
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
                    length_ratio, reads_info_dict,bam_path)
            #无论继续处理与否，都需要将当前sig放入新一轮的处理中
            cluster_list = list()
            cluster_list.append([ins_pos, ins_len, read_id, sig_origin])
            # print(cluster_list)
    #对最后一个cluster的处理
    if len(cluster_list) >= min_support:
        # print(cluster_list)
        cluster_through_len_del(cluster_list,chr,candidate_single_SV, min_support,\
            length_ratio, reads_info_dict,bam_path)
    file.close()
    # print(candidate_single_SV)
    # print(len(candidate_single_SV))
    iteration_dict, primary_num_dict, cover_dict = overlap_cover(candidate_single_SV, reads_info_dict) # both key(sv idx), value(set(read id))
    read_id_dict = list()
    for i in candidate_single_SV:
        read_id_dict.append(i[7])
    assign_list = assign_gt(iteration_dict, primary_num_dict, cover_dict,read_id_dict)
    # print(assign_list)
# assign_list.append([len(read_id_dict[idx]), DR, GT, GL, GQ, QUAL])
    # candidate_single_SV.append([chr, "DEL", str(int(del_start)), str(int(del_len)), str(len(cluster_list)), \
        #     str(CIPOS), str(CILEN), str(DR), str(GT), str(GL), str(GQ), str(QUAL)])
    
    final_sv_list = list()
    for i in range(len(assign_list)):
        final_sv_list.append([chr,"DEL",str(candidate_single_SV[i][2]),str(candidate_single_SV[i][3]),\
            str(assign_list[i][0]),str(candidate_single_SV[i][5]),str(candidate_single_SV[i][6]),\
            str(assign_list[i][1]),str(assign_list[i][2]),str(assign_list[i][3]),str(assign_list[i][4]),str(assign_list[i][5])])

            # print(final_sv_list)
    return final_sv_list
    # return candidate_single_SV

def overlap_cover(svs_list, reads_list):
    # [(10024, 12024), (89258, 91258), ...]
    # [[10000, 10468, 0, 'm54238_180901_011437/52298335/ccs'], [10000, 17490, 1, 'm54238_180901_011437/44762027/ccs'], ...]
    sort_list = list()
    idx = 0
    for i in reads_list:
        sort_list.append([i[0], 1, idx, i[2], i[3]])
        sort_list.append([i[1], 2, idx, i[2], i[3]])
        idx += 1
    idx = 0
    for i in svs_list:
        sort_list.append([i[2]-1000, 3, idx])
        sort_list.append([i[2]+1000, 0, idx])
        idx += 1
    sort_list = sorted(sort_list, key = lambda x:(x[0], x[1]))
    svs_set = set()
    read_set = set()
    overlap_dict = dict()
    cover_dict = dict()
    for node in sort_list:
        if node[1] == 1: # set2(read) left
            read_set.add(node[2])
            for x in svs_set:
                if svs_list[x][1] == node[0]:
                    continue
                if x not in overlap_dict:
                    overlap_dict[x] = set()
                overlap_dict[x].add(node[2])
        elif node[1] == 2: # set2(read) right
            read_set.remove(node[2])
        elif node[1] == 3: # set1(sv) left
            svs_set.add(node[2])
            overlap_dict[node[2]] = set()
            for x in read_set:
                overlap_dict[node[2]].add(x)
            cover_dict[node[2]] = set()
            for x in read_set:
                cover_dict[node[2]].add(x)
        elif node[1] == 0: # set1(sv) right
            svs_set.remove(node[2])
            temp_set = set()
            for x in read_set:
                temp_set.add(x)
            cover_dict[node[2]] = cover_dict[node[2]] & temp_set
    cover2_dict = dict()
    iteration_dict = dict()
    primary_num_dict = dict()
    for idx in cover_dict:
        iteration_dict[idx] = len(overlap_dict[idx])
        primary_num_dict[idx] = 0
        for x in overlap_dict[idx]:
            if reads_list[x][2] == 1:
                primary_num_dict[idx] += 1
        cover2_dict[idx] = set()
        for x in cover_dict[idx]:
            if reads_list[x][2] == 1:
                cover2_dict[idx].add(reads_list[x][3])
    # duipai(svs_list, reads_list, iteration_dict, primary_num_dict, cover2_dict)
    return iteration_dict, primary_num_dict, cover2_dict

def assign_gt(iteration_dict, primary_num_dict, cover_dict, read_id_dict):
    assign_list = list()
    for idx in range(len(read_id_dict)):
        iteration = iteration_dict[idx]
        primary_num = primary_num_dict[idx]
        read_count = cover_dict[idx]
        DR = 0
        for query in read_count:
            if query not in read_id_dict[idx]:
                DR += 1
        GT, GL, GQ, QUAL = cal_GL(DR, len(read_id_dict[idx]))
        
        assign_list.append([len(read_id_dict[idx]), DR, GT, GL, GQ, QUAL])
    return assign_list


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

