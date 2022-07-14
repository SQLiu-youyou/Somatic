from math import log10
import rlcompleter
import numpy as np

def cal_CIPOS(std, num):
    pos = int(1.96 * std / num ** 0.5)
    return "-%d,%d"%(pos,pos)

err = 0.1
prior = float(1/3)
Genotype = ["0/0", "0/1", "1/1"]

def log10sumexp(log10_probs):
    # Normalization of Genotype likelihoods
    m = max(log10_probs)
    return m + log10(sum(pow(10.0, x-m) for x in log10_probs))

def normalize_log10_probs(log10_probs):
    # Adjust the Genotype likelihoods
    log10_probs = np.array(log10_probs)
    lse = log10sumexp(log10_probs)
    return np.minimum(log10_probs - lse, 0.0)

def rescale_read_counts(c0, c1, max_allowed_reads=100):
    """Ensures that n_total <= max_allowed_reads, rescaling if necessary."""
    Total = c0 + c1
    if Total > max_allowed_reads:
        c0 = int(max_allowed_reads * float(c0/Total))
        c1 = max_allowed_reads - c0
    return c0, c1

def cal_GL(c0, c1):
    # Approximate adjustment of events with larger read depth
    c0, c1 = rescale_read_counts(c0, c1)
    # original genotype likelihood
    # ori_GL00 = np.float64(pow((1-err), c0)*pow(err, c1)*comb(c0+c1,c0)*(1-prior)/2)
    # ori_GL11 = np.float64(pow(err, c0)*pow((1-err), c1)*comb(c0+c1,c0)*(1-prior)/2)
    # ori_GL01 = np.float64(pow(0.5, c0+c1)*comb(c0+c1,c0)*prior)
    
    ori_GL00 = np.float64(pow((1-err), c0)*pow(err, c1)*(1-prior)/2)
    ori_GL11 = np.float64(pow(err, c0)*pow((1-err), c1)*(1-prior)/2)
    ori_GL01 = np.float64(pow(0.5, c0+c1)*prior)

    # normalized genotype likelihood
    prob = list(normalize_log10_probs([log10(ori_GL00), log10(ori_GL01), log10(ori_GL11)]))
    GL_P = [pow(10, i) for i in prob]
    PL = [int(np.around(-10*log10(i))) for i in GL_P]
    GQ = [int(-10*log10(GL_P[1] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[1]))]
    QUAL = abs(np.around(-10*log10(GL_P[0]), 1))

    return Genotype[prob.index(max(prob))], "%d,%d,%d"%(PL[0], PL[1], PL[2]), max(GQ), QUAL

def count_coverage(reads_list, search_start, search_end, up_bound, itround,\
    read_id_list, gt_primary):
    status = 0
    iteration = 0
    primary_num = 0
    read_count = set()
    iteration_set = set()

    # if search_start == 248554444:
    #     print(gt_primary)

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
    
    if status == 0 and len(read_count) == 0:
        read_count = iteration_set

    # if search_start == 248554444:
    #     print(read_count)
    #     print(len(read_count))
    #     print("iteration=%d"%(iteration))
    #     print("primary_num=%d"%(primary_num))

    # if search_start == 14961957:
    #     print(status)
    #     print(read_count)
    #     print(len(read_count))
    #     print(gt_primary)
    
    # if len(read_count)*100 < iteration:
    #     return status,iteration_set

    if gt_primary == 1:
        # if search_start == 248554444:
        #     print(status)
        #     print(read_count)
        if len(read_count)*100 < iteration:
            return status,iteration_set
        else:
            return status,read_count
    else:
        # if search_start == 248554444:
        #     print(status)
        #     print(len(iteration_set))
        return status,iteration_set

    # if len(read_count) < len(read_id_list):
    #     return status,iteration_set         
    
    # if search_start == 42372687:
    #     print(read_count)
    #     print(len(read_count))
    #     print("iteration=%d"%(iteration))
    #     print("primary_num=%d"%(primary_num))

    # return status,read_count

def threshold_ref_count(num):
    if num <= 2:
        return 10*num
    elif 3 <= num <= 5:
        return 5*num 
    elif 6 <= num <= 15:
        return 4*num
    else:
        return 3*num