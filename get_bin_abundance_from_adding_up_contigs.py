# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 11:49:07 2022

@author: acer
"""

import json
import sys




def invert_bin_dict(bin_dict):
    contig2bin={}
    for b in bin_dict:
        for contig in bin_dict[b]:
            contig2bin[contig]=b
    return contig2bin



if __name__=='__main__':
    dastool_abundance={}
    
    path2bindict=sys.argv[1]
    sample=path2bindict.split('/')[-1].split('.')[0]
    
    
    bin2contig=json.load(open(path2bindict))
    contig2bin=invert_bin_dict(bin2contig)
    contig_abundances=open('/net/phage/linuxhome/mgx/people/tina/wageningen/2018/{}/{}.RAT.contig.abundance.txt'.format(sample,sample)).read().rstrip().split('\n')[3:]
    for c in contig_abundances:
        c=c.split('\t')
        if c[0] in contig2bin:
            if contig2bin[c[0]] not in dastool_abundance:
                dastool_abundance[contig2bin[c[0]]]=0
            dastool_abundance[contig2bin[c[0]]]+=float(c[2])
    
    with open('/net/phage/linuxhome/mgx/people/tina/wageningen/2018/{}vs{}.DASTool.bin.reads.txt'.format(sample,sample), 'w') as outf:
        for mag in dastool_abundance:
            outf.write('{}\t{}\n'.format(mag, dastool_abundance[mag]))
    