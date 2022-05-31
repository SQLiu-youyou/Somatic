import argparse



def parseArgs(argv):
    parser = argparse.ArgumentParser(prog="cuteSV", \
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("normal_input", 
        metavar="[BAM]", 
        type = str, 
        help ="Sorted normal.bam file from NGMLR or Minimap2.")
    parser.add_argument("tumor_input", 
        metavar="[BAM]", 
        type = str, 
        help ="Sorted tumor.bam file from NGMLR or Minimap2.")
    parser.add_argument('-b', '--batches', 
        help = "Batch of genome segmentation interval.[%(default)s]", 
        default = 10000000, 
        type = int)

    GroupSignaturesCollect = parser.add_argument_group('Collection of SV signatures')  
    GroupSignaturesCollect.add_argument('-p', '--max_split_parts', 
        help = "Maximum number of split segments a read may be aligned before it is ignored. All split segments are considered when using -1. \
            (Recommand -1 when applying assembly-based alignment.)[%(default)s]", 
        default = 7, 
        type = int)
    GroupSignaturesCollect.add_argument('-q', '--min_mapq', 
        help = "Minimum mapping quality value of alignment to be taken into account.[%(default)s]", 
        default = 20, 
        type = int)
    GroupSignaturesCollect.add_argument('-sm', '--sig_min_cigar_size', 
        help = "min size of gain ins from cigar.[%(default)s]", 
        default = 20, 
        type = int)
    GroupSignaturesCollect.add_argument('-cmin', '--chase_ins_min_size', 
        help = "min size of combine ins intra read.[%(default)s]", 
        default = 20, 
        type = int)
    GroupSignaturesCollect.add_argument('-cmax', '--chase_ins_max_size', 
        help = "max size of combine ins intra read.[%(default)s]", 
        default = 100000, 
        type = int)
    GroupSignaturesCollect.add_argument('-cm', '--combine_min_size', 
        help = "min size of combine ins in same read.[%(default)s]", 
        default = 200, 
        type = int)


    GroupSVCluster = parser.add_argument_group('Generation of SV clusters')
    GroupSVCluster.add_argument('-ns', '--normal_min_support', 
        help = "Minimum number of reads that support a SV to be reported.[%(default)s]", 
        default = 2, 
        type = int)
    GroupSVCluster.add_argument('-ts', '--tumor_min_support', 
        help = "Minimum number of reads that support a SV to be reported.[%(default)s]", 
        default = 1, 
        type = int)
    GroupSVCluster.add_argument('-nl', '--normal_min_size', 
        help = "Minimum size of SV to be reported.[%(default)s]", 
        default = 20, 
        type = int)
    GroupSVCluster.add_argument('-tl', '--tumor_min_size', 
        help = "Minimum size of SV to be reported.[%(default)s]", 
        default = 50, 
        type = int)
    GroupSVCluster.add_argument('-nf', '--normal_length_ratio_flag', 
        help = "Minimum size of SV to be reported.[%(default)s]", 
        default = 0, 
        type = int)
    GroupSVCluster.add_argument('-tf', '--tumor_length_ratio_flag', 
        help = "Minimum size of SV to be reported.[%(default)s]", 
        default = 1, 
        type = int)
    GroupSVCluster.add_argument('-tg', '--threshold_gloab', 
        help = "Minimum size of SV to be reported.[%(default)s]", 
        default = 0.2, 
        type = int)
    GroupSVCluster.add_argument('-lr', '--detailed_length_ratio', 
        help = "Minimum size of SV to be reported.[%(default)s]", 
        default = 0.8, 
        type = int)
    GroupSVCluster.add_argument('--max_cluster_bias_INS', 
        help = "Maximum distance to cluster read together for insertion.[%(default)s]", 
        default = 50, 
        type = int)

    args = parser.parse_args(argv)
    return args