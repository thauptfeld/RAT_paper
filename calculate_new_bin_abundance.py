# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 16:20:50 2022

@author: acer
"""

import json
import os

samples=['W19-1','W19-2','W19-3','W19-4','W19-5','W19-6',
          'W22-1','W22-2','W22-3','W22-4','W22-5','W22-6',
          'W23-1','W23-2','W23-3','W23-4','W23-5','W23-6'
    ]
# path_to_reads='/home/tina/Documents/RAT/bio_example/all.reads.txt'
# path_to_json='/net/phage/linuxhome/mgx/people/tina/wageningen/2018/all_vs_all_treebins/'
# path_to_dRep='/home/tina/Documents/RAT/bio_example/dRep_clusters_useful.json'


# sample2reads={}
# with open(path_to_reads) as inf:
#     for line in inf:
#         line=line.split()
#         smp=line[1].split('/')[9]
#         print(smp)
#         sample2reads[smp]=int(line[0])-1

# mag2cluster=json.load(open(path_to_dRep))

# tree_abundance={}
# files=[path_to_json+f for f in os.listdir(path_to_json)]

# for f in files:
#     filename=os.path.split(f)[-1][:-5]
#     read=filename.split('vs')[0]
#     contig=filename.split('vs')[1]
    
#     print('{} vs {}'.format(read,contig))
    
#     if contig not in tree_abundance:
#         tree_abundance[contig]={}
#     if read not in tree_abundance[contig]:
#         tree_abundance[contig][read]={}
    
#     smp=json.load(open(f))
#     for mag in smp:
#         if mag in mag2cluster and mag2cluster[mag].startswith(contig):
#             if mag not in tree_abundance[contig][read]:
#                 tree_abundance[contig][read][mag2cluster[mag]]=0
#             tree_abundance[contig][read][mag2cluster[mag]]+=smp[mag]/sample2reads[read]
    
# with open('/home/tina/Documents/RAT/bio_example/tree_abundances_all_vs_all_mapped.json', 'w') as outf:
#     outf.write(json.dumps(tree_abundance, indent=4))
    
tree_genomes=open('/home/tina/Documents/RAT/bio_example/tree_genomes.txt').read().strip().split()
abundances=json.load(open('/home/tina/Documents/RAT/bio_example/tree_abundances_all_vs_all_mapped.json'))

abund_alm={}

for genome in tree_genomes:
    smp=genome.split('_')[0]
    if genome not in abund_alm:
        abund_alm[genome]=[]
    for s in samples:
        if genome in abundances[smp][s]:
            abund_alm[genome].append(abundances[smp][s][genome])
        else:
            abund_alm[genome].append(0)

abund_norm={}
for genome in abund_alm:
    abund_norm[genome]=[]
    divisor=max(abund_alm[genome])
    for value in abund_alm[genome]:
        abund_norm[genome].append(value/divisor)
        
        
with open('/home/tina/Documents/RAT/bio_example/abundances_normalized.json', 'w') as outf:
    outf.write(json.dumps(abund_norm, indent=4))
        
        