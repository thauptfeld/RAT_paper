# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 10:36:39 2021

@author: acer
"""

import json
import numpy as np

clusters=json.load(open('../bio_example/dRep_clusters.json'))
genomes=open('../bio_example/tree_genomes.txt').read().strip().split()
samples=open('../bio_example/samples.txt').read().strip().split()    
# clusters2={}
# for c in clusters:
#     for b in clusters[c]:
#         if b in genomes:
#             clusters2[b]=clusters[c]
# del b, c, clusters


# clusters_t={}
# for c in clusters2:
#     for ci in clusters2[c]:
#         clusters_t[ci]=c
# del c, ci

# abundance_cluster={}
# with open('../bio_example/all.DASTool.bin.reads.txt') as inp:
#     for line in inp:
#         mag= line.split('\t')[0]
#         frac=line.split('\t')[1].rstrip()
#         if mag in clusters_t:
#             abundance_cluster[mag]=float(frac)
            
# genomes_c_abundance={}
# for g in genomes:
#     genomes_c_abundance[g]={}
#     for s in samples:
#         genomes_c_abundance[g][s]=0

# del g,s,line,inp,mag,frac


# for c in abundance_cluster:
#     smp=c.split('_')[0]
#     genomes_c_abundance[clusters_t[c]][smp]=abundance_cluster[c]
    
    
# scores={}    
# with open('../bio_example/tree/dRep_cluster_abundances_not_normalized.txt', 'w') as outf:
#     outf.write('ID\t{}\n'.format('\t'.join(samples)))
#     for c in genomes_c_abundance:
#         values=[genomes_c_abundance[c][s] for s in samples]
#         scores[c]=sum(values)
#         values_norm=[]
#         for v in values:
#             values_norm.append(v)
#         outf.write('{}\t{}\n'.format(c, '\t'.join([str(v) for v in values_norm])))

# with open('../bio_example/tree/stack_score_clusters.json', 'w') as outf:
#     outf.write(json.dumps(scores, indent=4))
        
        
        
# ### Abundances per phylum

# abundance_per_phylum={}
# for smp in samples:
    
#     with open('/net/phage/linuxhome/mgx/people/tina/wageningen/2018/{}/{}.GTDB.complete.abundance.txt'.format(smp, smp)) as inp:
#         for line in inp:
#             if not line.startswith('#') and not line.startswith('un'):
#                 try:
#                     phylum=line.split('\t')[0].split(';')[2].replace('*','')
#                 except IndexError:
#                     phylum='NA'
#                 frac=float(line.split('\t')[2])
#                 if phylum not in ['NA', 'not classified']:
#                     if phylum not in abundance_per_phylum:
#                         abundance_per_phylum[phylum]={}
#                         for s in samples:
#                             abundance_per_phylum[phylum][s]=0
#                     abundance_per_phylum[phylum][smp]+=frac
#     json.dump(abundance_per_phylum, open('/home/tina/Documents/RAT/bio_example/phylum_abundances_GTDB.json', 'w'))
        
        
# abundance_per_phylum=json.load(open('../bio_example/phylum_abundances_GTDB.json'))
# abundance_for_itol={}
# b2c=open('../bio_example/tree/all.GTDB.bin2classification.txt').read().strip().split('\n')[1:]
# for b in b2c:
#     if b.split()[0] in genomes:
#         try:
#             phylum=b.split('\t')[3].split(';')[2]
#         except IndexError:
#             phylum='no support'
#         if phylum in abundance_per_phylum:
#             abundance_for_itol[b.split()[0]]=abundance_per_phylum[phylum]

    
# with open('../bio_example/tree/abundances_perphylum_GTDB_normalized.txt', 'w') as outf:
#     outf.write('ID\t{}\n'.format('\t'.join(samples)))
#     for c in abundance_for_itol:
#         values=[abundance_for_itol[c][s] for s in samples]
#         values_norm=[]
#         for v in values:
#             values_norm.append(v/max(values))
#         outf.write('{}\t{}\n'.format(c, '\t'.join([str(v) for v in values_norm])))
        
        
        
# abundance=open('../bio_example/tree/abundances_perphylum_GTDB_not_normalized.txt').read().strip().split('\n')[1:]
# scores={}
# for line in abundance:
#     mag=line.split()[0]
#     score=sum([float(s) for s in line.split()[1:]])
#     scores[mag]=score
# with open('../bio_example/tree/stack_score_phylum.json', 'w') as outf:
#     outf.write(json.dumps(scores, indent=4))
        
b2p={}
phylum_counts={}
b2c=open('../bio_example/tree/all.GTDB.bin2classification.txt').read().strip().split('\n')[1:]
for b in b2c:
    if b.split()[0] in genomes:
        try:
            phylum=b.split('\t')[3].split(';')[2]
        except IndexError:
            phylum='no support'
        if phylum not in phylum_counts:
            phylum_counts[phylum]=0
        phylum_counts[phylum]+=1
        b2p[b.split()[0]]=phylum