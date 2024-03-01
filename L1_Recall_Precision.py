#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 11:40:07 2020

@author: tina
"""
import json
ranks=["superkingdom" ,"phylum","class","order","family","genus","species"]



def evaluate_RAT(RAT_json, CAMI_profile):
    stats={}
    cami=make_cami_dict(CAMI_profile)
    rat=json.load(open(RAT_json))
    
    for rank in ranks:
        del rat[rank]['unclassified'],rat[rank]['unmapped'] # don't count unmapped and unclassified as taxa
        stats[rank]={'l1': None, 
                 'TP': None, 
                 'precision': None,
                 'FP': None,
                 'FN': None,
                 'sensitivity': None}
        
        n_taxa_rat=len([t for t in rat[rank] if not t.startswith('un')])    
        TP=len([t for t in cami[rank] if t in rat[rank]])
        FN=len(cami[rank])-TP
        FP=n_taxa_rat-TP
        
        stats[rank]['TP']=TP
        stats[rank]['FN']=FN
        stats[rank]['FP']=FP
        stats[rank]['precision']=TP/n_taxa_rat
        stats[rank]['sensitivity']=TP/(TP+FN)

        
        total_diff=0
        for s in rat[rank]:
            xr=rat[rank][s][1]
            xi=0
            if s in cami[rank]:
                xi=cami[rank][s]
            total_diff+=abs(xr-xi)
        for u in [v for v in cami[rank] if v not in rat[rank]]:
            xi=cami[rank][u]
            total_diff+=xi
        stats[rank]['l1']=total_diff
        
    return stats
    
  
        
    

    
def make_cami_dict(CAMI_profile):
    cami={}
    with open(CAMI_profile, 'r') as inp:
        for line in inp:
            if not line.startswith('@') and not line.startswith('\n'):
                line=line.rstrip().split('\t')
                rank, lineage, abundance = line[1], line[2], float(line[4])/100
                
                if rank not in cami and abundance>0:
                    cami[rank]={}
                if abundance>0:
                    cami[rank][lineage]=abundance
    return cami

def write_evaluation(stats, outfile):
    with open(outfile, 'w+') as outf:
        for metric in ['TP', 'FN', 'FP','sensitivity', 'precision', 'l1']:
            for rank in ranks:
                outf.write('{}, {}, '.format(metric, rank))
                for sample in stats:
                    outf.write('{},'.format(stats[sample][rank][metric]))
                outf.write('\n')



# if __name__=='__main__':
#     stats={}
#     # stats2={}
#     for i in range(0,64):
#         print(i)
#         stats['smp{}'.format(i)]=evaluate_RAT('/net/phage/linuxhome/tina/RAT/benchmark/CAMI_II_mousegut/profiles/smp{}.RAT.json'.format(i), 
#              '/net/phage/linuxhome/tina/CAMI/CAMI_II/mousegut/19122017_mousegut_scaffolds/taxonomic_profile_{}.txt'.format(i))
    

#     write_evaluation(stats, '/net/phage/linuxhome/tina/RAT/benchmark/CAMI_II_mousegut/ALL_samples_stats.csv')

