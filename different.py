from re import T
import sys


def query_sv_set(normal_sv_list,breakpoint,break_len):
    bias = 1500
    len_similarity = 0.7
    
    left = 0
    right = len(normal_sv_list) - 1
    # print(len(normal_sv_list))
    mid = 0
    while left <= right:
        mid = (left + right + 1) >> 1
        if int(breakpoint) >= (int(normal_sv_list[mid][0]) - bias) and int(breakpoint) <= (int(normal_sv_list[mid][0]) + bias):
            # if int(breakpoint) == 16208164:
            #     print("hree")
            # # print("yes")
            #     print(break_len)
            #     print(normal_sv_list[mid][0])
            #     print(int(normal_sv_list[mid][1])/int(break_len))
            #     print(int(break_len)/int(normal_sv_list[mid][1]))
            #     print(int(normal_sv_list[mid][0]) + bias)
            if int(normal_sv_list[mid][1]) > int(break_len) and int(break_len)/int(normal_sv_list[mid][1]) > len_similarity:
                # if int(breakpoint) == 16208164:
                #     print("1")
                return True
            elif int(normal_sv_list[mid][1]) < int(break_len) and int(normal_sv_list[mid][1])/int(break_len) > len_similarity:
                # print("2")
                return True
            elif int(normal_sv_list[mid][1]) == int(break_len):
                return True
            # if ((int(normal_sv_list[mid][1]) >= int(break_len) and (int(break_len)/int(normal_sv_list[mid][1])) > len_similarity)) or \
            #     ((int(normal_sv_list[mid][1]) < int(break_len) and (int(normal_sv_list[mid][1])/int(break_len)) > len_similarity)):
                # print(normal_sv_list[mid][1])
                # print(int)
                # print(int(normal_sv_list[mid][1])/int(break_len))
                # print(int(break_len)/int(normal_sv_list[mid][1]))
                # return True
            
            
            # if (int(normal_sv_list[mid][1])/int(break_len)) > len_similarity or \
            #     (int(break_len)/int(normal_sv_list[mid][1])) > len_similarity:
            #     print(int(normal_sv_list[mid][1])/int(break_len))
            #     # print("yes")
            #     return True
              
        
        if int(normal_sv_list[mid][0]) < int(breakpoint):
            left = mid + 1
        else:
            right = mid - 1
    return False


if __name__=='__main__':
    normal_sv_file = open(sys.argv[1],'r')
    normal_sv = dict()
    normal_sv['INS'] = dict()
    normal_sv['DEL'] = dict()
    for line in normal_sv_file:
        svtype = line.strip().split('\t')[1]
        if svtype == 'INS':
            chr = line.strip().split('\t')[0]
            ins_pos = line.strip().split('\t')[2]
            ins_len = line.strip().split('\t')[3]
            if chr not in normal_sv['INS']:
                normal_sv['INS'][chr] = list()
            normal_sv['INS'][chr].append([ins_pos,ins_len])
        if svtype == 'DEL':
            chr = line.strip().split('\t')[0]
            del_start = line.strip().split('\t')[2]
            del_len = line.strip().split('\t')[3]
            del_end = int(del_start) + int(del_len)
            if chr not in normal_sv['DEL']:
                normal_sv['DEL'][chr] = list()
            normal_sv['DEL'][chr].append([del_start, del_len, del_end])
    # print(normal_sv['DEL'])
            
    # print(normal_sv)
    tumor_sv_file = open(sys.argv[2],'r')
    for line in tumor_sv_file:
        svtype = line.strip().split('\t')[1]
        support = int(line.strip().split('\t')[4])
        if svtype == 'INS':
            chr = line.strip().split('\t')[0]
            if chr in normal_sv['INS']:
        
                breakpoint = line.strip().split('\t')[2]
                break_len = line.strip().split('\t')[3]
                # print("check")
                # print(normal_sv['INS'][chr])
                # print(breakpoint)
                if not query_sv_set(normal_sv['INS'][chr],breakpoint,break_len):
                    # if support <= max_support:
                    print(line.strip())
        if svtype == 'DEL':
            chr = line.strip().split('\t')[0]
            if chr in normal_sv['DEL']:
                del_start = line.strip().split('\t')[2]
                del_len = line.strip().split('\t')[3]
                # print(normal_sv['DEL'][chr],del_start, del_len)
                if not query_sv_set(normal_sv['DEL'][chr], del_start, del_len):
                    if int(del_len) > 30 :
                        print(line.strip())


        
