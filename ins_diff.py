import sys
from multiprocessing import Pool

def query_sigs(normal_list,window_number,tumor_start,tumor_len,bias,len_similatity):
    # if tumor_start == 114293036:
    #     print("ok")
    for ele in normal_list[window_number]:
        if 1 >= ele[1]/tumor_len > len_similatity or 1>= tumor_len/ele[1] > len_similatity :
            if  tumor_start - bias < ele[0] < tumor_start + bias:
                return True
        if ele[1] > (1+(1-len_similatity)+0.1)*tumor_len:
            break
    try:
        for ele in normal_list[window_number-1]:
            if 1 >= ele[1]/tumor_len > len_similatity or 1>= tumor_len/ele[1] > len_similatity :
                if  tumor_start - bias < ele[0] < tumor_start + bias:
                    return True
            if ele[1] > (1+(1-len_similatity)+0.1)*tumor_len:
                break
    except:
        pass
    try:
        for ele in normal_list[window_number+1]:
            if 1 >= ele[1]/tumor_len > len_similatity or 1>= tumor_len/ele[1] > len_similatity :
                if  tumor_start - bias < ele[0] < tumor_start + bias:
                    return True
            if ele[1] > (1+(1-len_similatity)+0.1)*tumor_len:
                break
    except:
        pass
    return False


def query_sigs_pro(normal_list,window_number,tumor_start,tumor_len,bias,len_similatity):
    # if tumor_start == 114293036:
    #     print("ok")
    for ele in normal_list[window_number]:
        if 1 >= ele[1]/tumor_len > len_similatity or 1>= tumor_len/ele[1] > len_similatity :
            if  tumor_start - bias < ele[0] < tumor_start + bias:
                # if tumor_start == 114293036:
                #     print(normal_list[window_number])
                return True
        if ele[1] > (1+(1-len_similatity)+0.1)*tumor_len:
            break
    
    return False

def load_sigs_chr(path):
    sigs_chr = dict()
    # sigs_chr['INS'] = list()
    sigs_chr['INS'] = list()
    # for svtype in ['INS']:
    #     file = open("%s%s%s.sigs"%(path,tumor_or_normal,svtype),'r')
    #     for line in file:
    #         chr = line.strip('\n').split('\t')[1]
    #         if chr not in sigs_chr[svtype]:
    #             sigs_chr[svtype].append(chr)
    #     file.close()
    #     sigs_chr[svtype].sort()
    for svtype in ['INS']:
        file = open(path,'r')
        for line in file:
            chr = line.strip('\n').split('\t')[1]
            if chr not in sigs_chr[svtype]:
                sigs_chr[svtype].append(chr)
        file.close()
        sigs_chr[svtype].sort()
    return sigs_chr


def resolution_INS(path,chr):
    # print("resolution")
    candidate = list()
    normal_sv_file = open(sys.argv[1],'r')
    normal_sv = dict()
    normal_sv['INS'] = dict()
    normal_sv_200 = dict()
    normal_sv_200['INS'] = dict()
    normal_sv['INS'] = dict()
    j = 0  
    for line in normal_sv_file:
        svtype = line.strip().split('\t')[0]    
        if svtype == 'INS':
            
            chrom = line.strip().split('\t')[1]
            if chrom != chr :
                continue
            del_start = int(line.strip().split('\t')[2])
            del_len = int(line.strip().split('\t')[3])
            if j != 0:
                if last_start == del_start  and last_len == del_len:
                    continue
                else:
                    # del_end = int(del_start) + int(del_len)
                    if chr not in normal_sv['INS']:
                        normal_sv['INS'][chr] = dict()
                        # sorted_normal_sv['INS'][chr] = dict()
                    if del_start < 1000:
                        if 0 not in normal_sv['INS'][chr]:
                            normal_sv['INS'][chr][0] = list()
                        normal_sv['INS'][chr][0].append([del_start,del_len])
                    else:
                        if int(del_start/1000) not in normal_sv['INS'][chr]:
                            normal_sv['INS'][chr][int(del_start/1000)] = list()
                        normal_sv['INS'][chr][int(del_start/1000)].append([del_start,del_len])
                    
                    if chr not in normal_sv_200['INS']:
                        normal_sv_200['INS'][chr] = dict()
                        # sorted_normal_sv['INS'][chr] = dict()
                    if del_start < 200:
                        if 0 not in normal_sv_200['INS'][chr]:
                            normal_sv_200['INS'][chr][0] = list()
                        normal_sv_200['INS'][chr][0].append([del_start,del_len])
                    else:
                        if int(del_start/200) not in normal_sv_200['INS'][chr]:
                            normal_sv_200['INS'][chr][int(del_start/200)] = list()
                        normal_sv_200['INS'][chr][int(del_start/200)].append([del_start,del_len])
            else:
                    if chr not in normal_sv['INS']:
                        normal_sv['INS'][chr] = dict()
                        # sorted_normal_sv['INS'][chr] = dict()
                    if del_start < 1000:
                        if 0 not in normal_sv['INS'][chr]:
                            normal_sv['INS'][chr][0] = list()
                        normal_sv['INS'][chr][0].append([del_start,del_len])
                    else:
                        if int(del_start/1000) not in normal_sv['INS'][chr]:
                            normal_sv['INS'][chr][int(del_start/1000)] = list()
                        normal_sv['INS'][chr][int(del_start/1000)].append([del_start,del_len])
                    
                    if chr not in normal_sv_200['INS']:
                        normal_sv_200['INS'][chr] = dict()
                        # sorted_normal_sv['INS'][chr] = dict()
                    if del_start < 200:
                        if 0 not in normal_sv_200['INS'][chr]:
                            normal_sv_200['INS'][chr][0] = list()
                        normal_sv_200['INS'][chr][0].append([del_start,del_len])
                    else:
                        if int(del_start/200) not in normal_sv_200['INS'][chr]:
                            normal_sv_200['INS'][chr][int(del_start/200)] = list()
                        normal_sv_200['INS'][chr][int(del_start/200)].append([del_start,del_len])

            last_start = del_start
            last_len = del_len
            j = j+1
            # normal_sv['INS'][chr].append([del_start, del_len])
        # sorted_normal_sv = dict()
        # if svtype == 'INS':
        #     for ele in normal_sv['INS']:
        #         sorted_list = sorted(normal_sv['INS'][ele].items(), key=lambda x: x[1], reverse=True)
        #         # sorted_normal_sv[]
        # sorted_list = sorted(normal_sv.items(), key=lambda x: x[1], reverse=True)
        # print(normal_sv)
    for svtype in normal_sv:
        for chrome in normal_sv[svtype]:
            for position in normal_sv[svtype][chrome]:
                normal_sv[svtype][chrome][position].sort(key=lambda x: x[1])
    for svtype in normal_sv_200:
        for chrome in normal_sv_200[svtype]:
            for position in normal_sv_200[svtype][chrome]:
                normal_sv_200[svtype][chrome][position].sort(key=lambda x: x[1])

    print(normal_sv)
    print(normal_sv_200)
    
    i = 0
    flag = 0
    tumor_sv_file = open(path,'r')
    for line in tumor_sv_file:
        svtype = line.strip().split('\t')[0]

        if svtype == 'INS':
            if line.strip().split('\t')[1] == chr:
            # chr = line.strip().split('\t')[1]
                if chr in normal_sv['INS']:
                    del_start = int(line.strip().split('\t')[2])
                    # print(del_start)
                    del_len = int(line.strip().split('\t')[3])
                    if del_len < 100:
                        bias = 200
                        len_similarity = 0.8
                        if i != 0:
                            if del_start == last_pos and del_len == last_len:
                                if  flag == 1:
                                    candidate.append(line)
                                    # print(line.strip())
                            else:

                                if int(del_start/200) in normal_sv_200['INS'][chr]:
                                    if not query_sigs(normal_sv_200['INS'][chr],int(del_start/200),del_start,del_len,bias,len_similarity):
                                        flag = 1
                                        candidate.append(line)
                                        # print(line.strip())
                                    else:
                                        flag = 0
                                elif int(del_start/200)-1 in normal_sv_200['INS'][chr] or int(del_start/200)+1 in normal_sv_200['INS'][chr]:
                                    if int(del_start/200)-1 in normal_sv_200['INS'][chr] and int(del_start/200)+1 in normal_sv_200['INS'][chr]:
                                        if not query_sigs_pro(normal_sv_200['INS'][chr],int(del_start/200-1),del_start,del_len,bias,len_similarity):
                                            # print(line.strip())
                                        
                                            if not query_sigs_pro(normal_sv_200['INS'][chr],int(del_start/200+1),del_start,del_len,bias,len_similarity):
                                                flag = 1
                                                candidate.append(line)
                                                # print(line.strip())
                                            else:
                                                flag = 0
                                        else:
                                            flag = 0 
                                            
                                
                                    elif int(del_start/200)-1 in normal_sv_200['INS'][chr]:
                                        if not query_sigs_pro(normal_sv_200['INS'][chr],int(del_start/200-1),del_start,del_len,bias,len_similarity):
                                            flag = 1
                                            candidate.append(line)
                                            # print(line.strip())
                                        else:
                                            flag = 0
                                    else:
                                        if not query_sigs_pro(normal_sv_200['INS'][chr],int(del_start/200+1),del_start,del_len,bias,len_similarity):
                                            flag = 1
                                            candidate.append(line)
                                            # print(line.strip())
                                        else:
                                            flag = 0
                                else:
                                    candidate.append(line)
                                    # print(line.strip())
                                    flag = 1
                        else:

                            if int(del_start/200) in normal_sv_200['INS'][chr]:
                                if not query_sigs(normal_sv_200['INS'][chr],int(del_start/200),del_start,del_len,bias,len_similarity):
                                    flag = 1
                                    candidate.append(line)
                                    # print(line.strip())
                                else:
                                    flag = 0
                            elif int(del_start/200)-1 in normal_sv_200['INS'][chr] or int(del_start/200)+1 in normal_sv_200['INS'][chr]:
                                if int(del_start/200)-1 in normal_sv_200['INS'][chr] and int(del_start/200)+1 in normal_sv_200['INS'][chr]:
                                    if not query_sigs_pro(normal_sv_200['INS'][chr],int(del_start/200-1),del_start,del_len,bias,len_similarity):
                                        # print(line.strip())
                                    
                                        if not query_sigs_pro(normal_sv_200['INS'][chr],int(del_start/200+1),del_start,del_len,bias,len_similarity):
                                            flag = 1
                                            candidate.append(line)
                                            # print(line.strip())
                                        else:
                                            flag = 0
                                    else:
                                        flag = 0 
                                        
                            
                                elif int(del_start/200)-1 in normal_sv_200['INS'][chr]:
                                    if not query_sigs_pro(normal_sv_200['INS'][chr],int(del_start/200-1),del_start,del_len,bias,len_similarity):
                                        flag = 1
                                        candidate.append(line)
                                        # print(line.strip())
                                    else:
                                        flag = 0
                                else:
                                    if not query_sigs_pro(normal_sv_200['INS'][chr],int(del_start/200+1),del_start,del_len,bias,len_similarity):
                                        flag = 1
                                        candidate.append(line)
                                        # print(line.strip())
                                    else:
                                        flag = 0
                            else:
                                candidate.append(line)
                                # print(line.strip())
                                flag = 1
                    else:
                        bias = 1000
                        len_similarity = 0.5
                        if i != 0:
                            if del_start == last_pos and del_len == last_len:
                                if  flag == 1:
                                    candidate.append(line)
                                    # print(line.strip())
                            else:

                                if int(del_start/1000) in normal_sv['INS'][chr]:
                                    if not query_sigs(normal_sv['INS'][chr],int(del_start/1000),del_start,del_len,bias,len_similarity):
                                        flag = 1
                                        candidate.append(line)
                                        # print(line.strip())
                                    else:
                                        flag = 0
                                elif int(del_start/1000)-1 in normal_sv['INS'][chr] or int(del_start/1000)+1 in normal_sv['INS'][chr]:
                                    if int(del_start/1000)-1 in normal_sv['INS'][chr] and int(del_start/1000)+1 in normal_sv['INS'][chr]:
                                        if not query_sigs_pro(normal_sv['INS'][chr],int(del_start/1000-1),del_start,del_len,bias,len_similarity):
                                            # print(line.strip())
                                        
                                            if not query_sigs_pro(normal_sv['INS'][chr],int(del_start/1000+1),del_start,del_len,bias,len_similarity):
                                                flag = 1
                                                candidate.append(line)
                                                # print(line.strip())
                                            else:
                                                flag = 0
                                        else:
                                            flag = 0 
                                            
                                
                                    elif int(del_start/1000)-1 in normal_sv['INS'][chr]:
                                        if not query_sigs_pro(normal_sv['INS'][chr],int(del_start/1000-1),del_start,del_len,bias,len_similarity):
                                            flag = 1
                                            candidate.append(line)
                                            # print(line.strip())
                                        else:
                                            flag = 0
                                    else:
                                        if not query_sigs_pro(normal_sv['INS'][chr],int(del_start/1000+1),del_start,del_len,bias,len_similarity):
                                            flag = 1
                                            candidate.append(line)
                                            # print(line.strip())
                                        else:
                                            flag = 0
                                else:
                                    candidate.append(line)
                                    # print(line.strip())
                                    flag = 1
                        else:

                            if int(del_start/1000) in normal_sv['INS'][chr]:
                                if not query_sigs(normal_sv['INS'][chr],int(del_start/1000),del_start,del_len,bias,len_similarity):
                                    flag = 1
                                    candidate.append(line)
                                    # print(line.strip())
                                else:
                                    flag = 0
                            elif int(del_start/1000)-1 in normal_sv['INS'][chr] or int(del_start/1000)+1 in normal_sv['INS'][chr]:
                                if int(del_start/1000)-1 in normal_sv['INS'][chr] and int(del_start/1000)+1 in normal_sv['INS'][chr]:
                                    if not query_sigs_pro(normal_sv['INS'][chr],int(del_start/1000-1),del_start,del_len,bias,len_similarity):
                                        # print(line.strip())
                                    
                                        if not query_sigs_pro(normal_sv['INS'][chr],int(del_start/1000+1),del_start,del_len,bias,len_similarity):
                                            flag = 1
                                            candidate.append(line)
                                            # print(line.strip())
                                        else:
                                            flag = 0
                                    else:
                                        flag = 0 
                                        
                            
                                elif int(del_start/1000)-1 in normal_sv['INS'][chr]:
                                    if not query_sigs_pro(normal_sv['INS'][chr],int(del_start/1000-1),del_start,del_len,bias,len_similarity):
                                        flag = 1
                                        candidate.append(line)
                                        # print(line.strip())
                                    else:
                                        flag = 0
                                else:
                                    if not query_sigs_pro(normal_sv['INS'][chr],int(del_start/1000+1),del_start,del_len,bias,len_similarity):
                                        flag = 1
                                        candidate.append(line)
                                        # print(line.strip())
                                    else:
                                        flag = 0
                            else:
                                candidate.append(line)
                                # print(line.strip())
                                flag = 1

                
                        

            last_pos = del_start
            last_len = del_len
            i = i + 1
    return candidate

def run_del(args):
    # print("in")
    return resolution_INS(*args)

if __name__=='__main__':

    # sorted_normal_sv = dict()
    # sorted_normal_sv['INS'] = dict()

    
    # print(normal_sv)
    
    analysis_pools = Pool(processes=16)
    
    value_sv = load_sigs_chr(sys.argv[2])
    # print(value_sv)
    result_sv = list()

    for chr in value_sv['INS']:
        # print(chr)
        para = [(sys.argv[2], chr)]
        # run_del(para[0])
        # analysis_pools.map_async(run_del,para)
        result_sv.append(analysis_pools.map_async(run_del,para))
        # print("Finish %s:%s"%(chr,'INS'))

        
    analysis_pools.close()
    analysis_pools.join()
     
    # # print("Writing to your output file.")

    #把进程里的所有candidate放到最终处理的列表中
    semi_result = list()
    for res in result_sv:
        try:
            semi_result += res.get()[0]
        except:
            pass
    #按照染色体和位置排序
    semi_result = sorted(semi_result, key = lambda x:(x[0]))
    # for ele in semi_result:
    #     print(ele)
    ans = open(sys.argv[3],'w')
    for ele in semi_result:
        ans.write(ele)