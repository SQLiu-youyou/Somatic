import numpy as np
from genotype import cal_CIPOS

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

def generate_final_del(cluster_list,candidate_single_SV,chr):
    del_start = int(np.mean([i[0] for i in cluster_list]))
    del_len = int(np.mean([i[1] for i in cluster_list]))
    if del_len >= 20 :
        CIPOS = cal_CIPOS(np.std([i[0] for i in cluster_list]), len([i[0] for i in cluster_list]))
        CILEN = cal_CIPOS(np.std([i[1] for i in cluster_list]), len([i[1] for i in cluster_list]))
        DR = '.'
        GT = './.'
        GL = '.,.,.'
        GQ = "."
        QUAL = "."
        candidate_single_SV.append([chr, "DEL", str(int(del_start)), str(int(del_len)), str(len(cluster_list)), \
            str(CIPOS), str(CILEN), str(DR), str(GT), str(GL), str(GQ), str(QUAL)])

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
    SortedList = sorted(list(reserved_read.values()), key = lambda x:x[1])
    # print(SortedList)
    standard_len = SortedList[0][1]
    global_len = [i[1] for i in SortedList]
    # 按ins_len聚类是，允许的len差值的最大值
    DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP = threshold_gloab * np.mean(global_len)
    
    # for ele in cluster_list:
    #     if ele[0]>= 2583425 and ele[0] <= 2584224:
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
            if len(to_SV_list) > 0.9 * min_support:
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
    if len(to_SV_list) > 0.9 * min_support:
        if lenth_ratio:
            if len(to_SV_list) / len(cluster_list) >= detailed_length_ratio:
                generate_final_ins(to_SV_list,candidate_single_SV,chr,min_size)
        else:
            generate_final_ins(to_SV_list,candidate_single_SV,chr,min_size)

def cluster_through_len_del(cluster_list,chr,candidate_single_SV):
    threshold_gloab = 0.5
    min_support = 14
    # detailed_length_ratio = 0.8
    
    
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
    SortedList = sorted(list(reserved_read.values()), key = lambda x:x[1])
    # print(SortedList)
    standard_len = SortedList[0][1]
    global_len = [i[1] for i in SortedList]
    # 按ins_len聚类是，允许的len差值的最大值
    DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP = threshold_gloab * np.mean(global_len)
    
    # for ele in cluster_list:
    #     if ele[0]>= 110393325 and ele[0] <= 110393334:
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
            if len(to_SV_list) > 0.9 * min_support:
                # if ele[0] == 129771776:
                #     print(to_SV_list)
                # if lenth_ratio:
                #     if len(to_SV_list) / len(cluster_list) >= detailed_length_ratio:
                #         generate_final_ins(to_SV_list,candidate_single_SV,chr,min_size)
                # else:
                #     generate_final_ins(to_SV_list,candidate_single_SV,chr,min_size)
                # print(to_SV_list)
                generate_final_del(to_SV_list,candidate_single_SV,chr)
            to_SV_list = list()
            standard_len = ele[1]
            to_SV_list.append(ele)
    #对最后一个堆的判定
    if len(to_SV_list) > 0.9 * min_support:
        # if lenth_ratio:
        #     if len(to_SV_list) / len(cluster_list) >= detailed_length_ratio:
        #         generate_final_ins(to_SV_list,candidate_single_SV,chr,min_size)
        # else:
        #     generate_final_ins(to_SV_list,candidate_single_SV,chr,min_size)
        # print(to_SV_list)
        generate_final_del(to_SV_list,candidate_single_SV,chr)
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

def resolution_DEL(sigs_path,chr):
    max_cluster_bias = 200
    min_support = 14
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
                # print(cluster_list)
                cluster_through_len_del(cluster_list,chr,candidate_single_SV)
            #无论继续处理与否，都需要将当前sig放入新一轮的处理中
            cluster_list = list()
            cluster_list.append([ins_pos, ins_len, read_id])
            # print(cluster_list)
    #对最后一个cluster的处理
    if len(cluster_list) >= min_support:
        # print(cluster_list)
        cluster_through_len_del(cluster_list,chr,candidate_single_SV)
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
    return resolution_DEL(*args)