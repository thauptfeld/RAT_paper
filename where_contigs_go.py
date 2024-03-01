# -*- coding: utf-8 -*-
"""
Created on Fri May 15 14:33:23 2020

@author: Tina
"""

import json

CAT='/hosts/linuxhome/chaperone/tina/Griftpark/CAT/'
bin_per_contig='/hosts/linuxhome/chaperone/tina/Griftpark/metagenomics/bin_per_contig.json'

def where_contigs_go(lineage, C2C_abundance_file, bin_per_contig_dict):
    '''
    Takes a lineage, in the format of a string: e.g. Betaproteobacteria_NA_NA_NA_NA
    Then, for each contig whose classification ends with that string, it checks in what bin it ends
    up, and the sum of read fractions that end up in that bin for that classification.
    '''
    OP_DICT={}
    errors=[]
    bin_info=read_in_dict_all(CAT+'BAT/bin_info_table.txt')
    C2C=read_in_dict_all(C2C_abundance_file)
    tmp=json.load(open(bin_per_contig_dict))
    C2B=tmp[C2C_abundance_file.split('.')[0]]
    del tmp
    for contig in C2C:
        try:
            if C2C[contig][0].endswith(lineage):
                if C2B[contig] not in OP_DICT:
                    OP_DICT[C2B[contig]]={'contigs': [], 'sum_of_fractions':0}
                OP_DICT[C2B[contig]]['contigs']+=[contig]
                OP_DICT[C2B[contig]]['sum_of_fractions']+=float(C2C[contig][1])
        except KeyError:
            errors += [contig]
   
    with open('{}.betaproteobacteria.to.bins.json'.format(C2C_abundance_file.split('.')[0]), 'w+') as op:
        op.write(json.dumps(OP_DICT, indent=4))
    with open('{}.betaproteobacteria.bins.txt'.format(C2C_abundance_file.split('.')[0]), 'w+') as op2:
        op2.write('{} total: {}\n'.format(C2C_abundance_file.split('.')[0], 
                  sum([float(OP_DICT[i]['sum_of_fractions']) for i in OP_DICT])))
        for b in OP_DICT:
            op2.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(b, OP_DICT[b]['sum_of_fractions'], 
                      bin_info[b][2], bin_info[b][3], bin_info[b][4], bin_info[b][5]))




def read_in_dict_all(dict_file, relevant_pathways=None, delimiter='\t', skip=[]):
    OP_DICT={}
    lines=open(dict_file).read().strip().split('\n')
    for line in lines:
        if not line=='':
            key=line.split(delimiter)[0]
            if relevant_pathways!=None:
                if key in relevant_pathways:
                    OP_DICT[key]=line.strip().split(delimiter)[1:]
                    if skip!=[]:
                        skip.sort(reverse=True)
                        for i in skip:
                            del OP_DICT[key][i-1]
            else:
                OP_DICT[key]=line.strip().split(delimiter)[1:]
                if skip!=[]:
                    skip.sort(reverse=True)
                    for i in skip:
                        del OP_DICT[key][i-1]
    return OP_DICT


if __name__=='__main__':
    samples=['2-100807', '2-100808', '2-100809', '2-100810', '2-100811', '2-100812',
             '2-100813', '2-100814', '2-100815', '2-100816', '2-100817', '2-100818',
             'S1','S2','S3','S4','S6','S7','S8','S9','S10','S11','S12', 'S13', 'S14']

    for sample in samples:
        print('Sample: '+sample)
        where_contigs_go('Betaproteobacteria_NA_NA_NA_NA', '{}.contig.abundance.txt'.format(sample),
                         bin_per_contig)