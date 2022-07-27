import sys


# def query_sv_set(normal_sv_list,breakpoint,break_len):
#     bias = 1500
#     len_similarity = 0.7
    
#     left = 0
#     right = len(normal_sv_list) - 1
#     mid = 0
#     while left <= right:
#         mid = (left + right + 1) >> 1
#         if int(breakpoint) >= (int(normal_sv_list[mid][0]) - bias) and \
#             int(breakpoint) <= (int(normal_sv_list[mid][0]) + bias):
            
#             if  1 >= (int(normal_sv_list[mid][1])/int(break_len)) > len_similarity or \
#                 1 >= (int(break_len)/int(normal_sv_list[mid][1])) > len_similarity:
#                 return True
              
        
#         if int(normal_sv_list[mid][0]) < int(breakpoint):
#             if mid > 0 and int(breakpoint) >= (int(normal_sv_list[mid-1][0]) - bias) and \
#                 int(breakpoint) <= (int(normal_sv_list[mid-1][0]) + bias):
#                 if  1 >= (int(normal_sv_list[mid-1][1])/int(break_len)) > len_similarity or \
#                 1 >= (int(break_len)/int(normal_sv_list[mid-1][1])) > len_similarity:
#                     return True
#             left = mid + 1
            
#         else:
#             if mid < len(normal_sv_list)-1 and int(breakpoint) >= (int(normal_sv_list[mid+1][0]) - bias) \
#                 and int(breakpoint) <= (int(normal_sv_list[mid+1][0]) + bias):
#                     if  1 >= (int(normal_sv_list[mid+1][1])/int(break_len)) > len_similarity or \
#                     1 >= (int(break_len)/int(normal_sv_list[mid+1][1])) > len_similarity:
#                         return True
#             right = mid - 1
#     return False


# if __name__=='__main__':
#     normal_sv_file = open(sys.argv[1],'r')
#     normal_sv = dict()
#     normal_sv['INS'] = dict()
#     normal_sv['DEL'] = dict()
#     for line in normal_sv_file:
#         svtype = line.strip().split('\t')[1]
#         if svtype == 'INS':
#             chr = line.strip().split('\t')[0]
#             ins_pos = line.strip().split('\t')[2]
#             ins_len = line.strip().split('\t')[3]
#             if chr not in normal_sv['INS']:
#                 normal_sv['INS'][chr] = list()
#             normal_sv['INS'][chr].append([ins_pos,ins_len])
#         if svtype == 'DEL':
#             chr = line.strip().split('\t')[0]
#             del_start = line.strip().split('\t')[2]
#             del_len = line.strip().split('\t')[3]
#             del_end = int(del_start) + int(del_len)
#             if chr not in normal_sv['DEL']:
#                 normal_sv['DEL'][chr] = list()
#             normal_sv['DEL'][chr].append([del_start, del_len, del_end])

#     tumor_sv_file = open(sys.argv[2],'r')
#     for line in tumor_sv_file:
#         svtype = line.strip().split('\t')[1]
#         support = int(line.strip().split('\t')[4])
#         if svtype == 'INS':
#             chr = line.strip().split('\t')[0]
#             if chr in normal_sv['INS']:
        
#                 breakpoint = line.strip().split('\t')[2]
#                 break_len = line.strip().split('\t')[3]
#                 if not query_sv_set(normal_sv['INS'][chr],breakpoint,break_len):
#                     print(line.strip())
#         if svtype == 'DEL':
#             chr = line.strip().split('\t')[0]
#             if chr in normal_sv['DEL']:
#                 del_start = line.strip().split('\t')[2]
#                 del_len = line.strip().split('\t')[3]
#                 if not query_sv_set(normal_sv['DEL'][chr], del_start, del_len):
#                     if int(del_len) > 30 :
#                         print(line.strip())



# bias = 1000
# len_similatity = 0.5

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

if __name__=='__main__':
    normal_sv_file = open(sys.argv[1],'r')
    normal_sv = dict()
    normal_sv['INS'] = dict()
    normal_sv['DEL'] = dict()
    # sorted_normal_sv = dict()
    # sorted_normal_sv['DEL'] = dict()
    for line in normal_sv_file:
        svtype = line.strip().split('\t')[0]
        # if svtype == 'INS':
        #     chr = line.strip().split('\t')[1]
        #     ins_pos = line.strip().split('\t')[2]
        #     ins_len = line.strip().split('\t')[3]
        #     if chr not in normal_sv['INS']:
        #         normal_sv['INS'][chr] = list()
        #     normal_sv['INS'][chr].append([ins_pos,ins_len])
        
        if svtype == 'DEL':
            chr = line.strip().split('\t')[1]
            del_start = int(line.strip().split('\t')[2])
            del_len = int(line.strip().split('\t')[3])
            # del_end = int(del_start) + int(del_len)
            if chr not in normal_sv['DEL']:
                normal_sv['DEL'][chr] = dict()
                # sorted_normal_sv['DEL'][chr] = dict()
            if del_start < 1000:
                if 0 not in normal_sv['DEL'][chr]:
                    normal_sv['DEL'][chr][0] = list()
                normal_sv['DEL'][chr][0].append([del_start,del_len])
            else:
                if int(del_start/1000) not in normal_sv['DEL'][chr]:
                    normal_sv['DEL'][chr][int(del_start/1000)] = list()
                normal_sv['DEL'][chr][int(del_start/1000)].append([del_start,del_len])
            # normal_sv['DEL'][chr].append([del_start, del_len])
        # sorted_normal_sv = dict()
        # if svtype == 'DEL':
        #     for ele in normal_sv['DEL']:
        #         sorted_list = sorted(normal_sv['DEL'][ele].items(), key=lambda x: x[1], reverse=True)
        #         # sorted_normal_sv[]
        # sorted_list = sorted(normal_sv.items(), key=lambda x: x[1], reverse=True)
        # print(normal_sv)
    for svtype in normal_sv:
        for chrome in normal_sv[svtype]:
            for position in normal_sv[svtype][chrome]:
                normal_sv[svtype][chrome][position].sort(key=lambda x: x[1])
    # print(normal_sv)
    
    
    i = 0
    flag = 0
    tumor_sv_file = open(sys.argv[2],'r')
    for line in tumor_sv_file:
        svtype = line.strip().split('\t')[0]

        if svtype == 'DEL':
            chr = line.strip().split('\t')[1]
            if chr in normal_sv['DEL']:
                del_start = int(line.strip().split('\t')[2])
                del_len = int(line.strip().split('\t')[3])
                if del_len < 100:
                    bias = 200
                    len_similarity = 0.8
                else:
                    bias = 1000
                    len_similarity = 0.5


                # if del_start == 114293036:
                #     print("in")
                if chr in normal_sv['DEL']:
                    if i != 0:
                        if del_start == last_pos and del_len == last_len:
                            if  flag == 1:
                                 print(line.strip())
                        else:

                            if int(del_start/1000) in normal_sv['DEL'][chr]:
                                if not query_sigs(normal_sv['DEL'][chr],int(del_start/1000),del_start,del_len,bias,len_similarity):
                                    flag = 1
                                    print(line.strip())
                                else:
                                    flag = 0
                            elif int(del_start/1000)-1 in normal_sv['DEL'][chr] or int(del_start/1000)+1 in normal_sv['DEL'][chr]:
                                if int(del_start/1000)-1 in normal_sv['DEL'][chr] and int(del_start/1000)+1 in normal_sv['DEL'][chr]:
                                    if not query_sigs_pro(normal_sv['DEL'][chr],int(del_start/1000-1),del_start,del_len,bias,len_similarity):
                                        # print(line.strip())
                                    
                                        if not query_sigs_pro(normal_sv['DEL'][chr],int(del_start/1000+1),del_start,del_len,bias,len_similarity):
                                            flag = 1
                                            print(line.strip())
                                        else:
                                            flag = 0
                                    else:
                                        flag = 0 
                                        
                            
                                elif int(del_start/1000)-1 in normal_sv['DEL'][chr]:
                                    if not query_sigs_pro(normal_sv['DEL'][chr],int(del_start/1000-1),del_start,del_len,bias,len_similarity):
                                        flag = 1
                                        print(line.strip())
                                    else:
                                        flag = 0
                                else:
                                    if not query_sigs_pro(normal_sv['DEL'][chr],int(del_start/1000+1),del_start,del_len,bias,len_similarity):
                                        flag = 1
                                        print(line.strip())
                                    else:
                                        flag = 0
                            else:
                                print(line.strip())
                                flag = 1
                    else:

                        if int(del_start/1000) in normal_sv['DEL'][chr]:
                            if not query_sigs(normal_sv['DEL'][chr],int(del_start/1000),del_start,del_len,bias,len_similarity):
                                flag = 1
                                print(line.strip())
                            else:
                                flag = 0
                        elif int(del_start/1000)-1 in normal_sv['DEL'][chr] or int(del_start/1000)+1 in normal_sv['DEL'][chr]:
                            if int(del_start/1000)-1 in normal_sv['DEL'][chr] and int(del_start/1000)+1 in normal_sv['DEL'][chr]:
                                if not query_sigs_pro(normal_sv['DEL'][chr],int(del_start/1000-1),del_start,del_len,bias,len_similarity):
                                    # print(line.strip())
                                
                                    if not query_sigs_pro(normal_sv['DEL'][chr],int(del_start/1000+1),del_start,del_len,bias,len_similarity):
                                        flag = 1
                                        print(line.strip())
                                    else:
                                        flag = 0
                                else:
                                    flag = 0 
                                    
                        
                            elif int(del_start/1000)-1 in normal_sv['DEL'][chr]:
                                if not query_sigs_pro(normal_sv['DEL'][chr],int(del_start/1000-1),del_start,del_len,bias,len_similarity):
                                    flag = 1
                                    print(line.strip())
                                else:
                                    flag = 0
                            else:
                                if not query_sigs_pro(normal_sv['DEL'][chr],int(del_start/1000+1),del_start,del_len,bias,len_similarity):
                                    flag = 1
                                    print(line.strip())
                                else:
                                    flag = 0
                        else:
                            print(line.strip())
                            flag = 1
                        
                    

        last_pos = del_start
        last_len = del_len
        i = i + 1
        

                    
                
          