def load_sigs_chr(path,tumor_or_normal):
    sigs_chr = dict()
    sigs_chr['INS'] = list()
    for svtype in ['INS']:
        file = open("%s%s%s.sigs"%(path,tumor_or_normal,svtype),'r')
        for line in file:
            chr = line.strip('\n').split('\t')[1]
            if chr not in sigs_chr[svtype]:
                sigs_chr[svtype].append(chr)
        file.close()
        sigs_chr[svtype].sort()
    
    return sigs_chr
