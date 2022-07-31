# import imp
import sys
import pysam
from Bio import SeqIO

def Generation_VCF_header(file, contiginfo):
    # General header
    file.write("##fileformat=VCFv4.2\n")
    # import time
    # file.write("##fileDate=%s\n"%(time.strftime('%Y-%m-%d %H:%M:%S %w-%Z',time.localtime())))
    for i in contiginfo:
        file.write("##contig=<ID=%s,length=%d>\n"%(i[0], i[1]))

    # Specific header
    # ALT
    file.write("##ALT=<ID=INS,Description=\"Insertion of novel sequence relative to the reference\">\n")
    file.write("##ALT=<ID=DEL,Description=\"Deletion relative to the reference\">\n")
    file.write("##ALT=<ID=DUP,Description=\"Region of elevated copy number relative to the reference\">\n")
    file.write("##ALT=<ID=INV,Description=\"Inversion of reference sequence\">\n")
    file.write("##ALT=<ID=BND,Description=\"Breakend of translocation\">\n")

    # INFO
    file.write("##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variant\">\n")
    file.write("##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variant\">\n")
    file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    file.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
    file.write("##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">\n")
    file.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n")
    file.write("##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n")
    file.write("##INFO=<ID=CILEN,Number=2,Type=Integer,Description=\"Confidence interval around inserted/deleted material between breakends\">\n")
    # file.write("##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">\n")
    file.write("##INFO=<ID=RE,Number=1,Type=Integer,Description=\"Number of read support this record\">\n")
    file.write("##INFO=<ID=STRAND,Number=A,Type=String,Description=\"Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)\">\n")
    file.write("##INFO=<ID=RNAMES,Number=.,Type=String,Description=\"Supporting read names of SVs (comma separated)\">\n")
    file.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency.\">\n")
    file.write("##FILTER=<ID=q5,Description=\"Quality below 5\">\n")
    # FORMAT
    # file.write("\n")
    file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    file.write("##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# High-quality reference reads\">\n")
    file.write("##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# High-quality variant reads\">\n")
    file.write("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"# Phred-scaled genotype likelihoods rounded to the closest integer\">\n")
    file.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"# Genotype quality\">\n")

    # file.write("##CommandLine=\"cuteSV %s\"\n"%(" ".join(argv)))


if __name__=='__main__':

    svid = dict()
    svid["INS"] = 0
    svid["DEL"] = 0

    file = open("tumor.vcf",'w')

    samfile = pysam.AlignmentFile(sys.argv[1])
    contigINFO = list()
    ref_ = samfile.get_index_statistics()
    for i in ref_:
        local_ref_len = samfile.get_reference_length(i[0])
        contigINFO.append([i[0], local_ref_len])
    ref_g = SeqIO.to_dict(SeqIO.parse(sys.argv[3], "fasta"))
    
    Generation_VCF_header(file, contigINFO)
    file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    
    inputfile = open(sys.argv[2],'r')
    for line in inputfile:
        i = line.strip().split('\t')
        if i[1] in ["DEL", "INS"]:
            if i[1] == "INS":
                cal_end = int(i[2])
            else:
                cal_end = int(i[2]) + abs(int(float(i[3])))
            info_list = "{PRECISION};SVTYPE={SVTYPE};SVLEN={SVLEN};END={END};CIPOS={CIPOS};CILEN={CILEN};RE={RE};".format(
                PRECISION = "IMPRECISE" if i[8] == "0/0" else "PRECISE", 
                SVTYPE = i[1], 
                SVLEN = i[3], 
                END = str(cal_end), 
                CIPOS = i[5], 
                CILEN = i[6], 
                RE = i[4])
            # if action:
            #     try:
            #         info_list += ";AF=" + str(round(int(i[4]) / (int(i[4]) + int(i[7])), 4))
            #     except:
            #         info_list += ";AF=."
            if i[1] =="DEL":
                info_list += ";STRAND=+-"
            if i[11] == "." or i[11] == None:
                filter_lable = "PASS"
            else:
                filter_lable = "PASS" if float(i[11]) >= 5.0 else "q5"
            file.write("{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}\t{GT}:{DR}:{RE}:{PL}:{GQ}\n".format(
                CHR = i[0], 
                POS = str(int(i[2])), 
                ID = "%s.%d"%(i[1], svid[i[1]]),
                # REF = str(ref_g[i[0]].seq[max(int(i[2])-1, 0)]) if i[1] == 'INS' else str(ref_g[i[0]].seq[max(int(i[2])-1, 0):int(i[2])+int(i[3])]),
                REF = str(ref_g[i[0]].seq[max(int(i[2])-1, 0)]) if i[1] == 'INS' else "<DEL>",
                ALT = "%s"%(str(ref_g[i[0]].seq[max(int(i[2])-1, 0)])+i[13] if i[1] == 'INS' else str(ref_g[i[0]].seq[max(int(i[2])-1, 0)])), 
                INFO = info_list, 
                FORMAT = "GT:DR:DV:PL:GQ", 
                GT = i[8],
                DR = i[7],
                RE = i[4],
                PL = i[9],
                GQ = i[10],
                QUAL = i[11],
                PASS = filter_lable))
            svid[i[1]] += 1