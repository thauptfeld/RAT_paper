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
    # cami=CAMI_profile
    # rat=RAT_json
    
    for rank in ranks:
        # del rat[rank]['unclassified'],rat[rank]['unmapped'] # don't count unmapped and unclassified as taxa
        stats[rank]={'l1': None, 
                 'TP': None, 
                 'precision': None,
                 'FP': None,
                 'FN': None,
                 'sensitivity': None}
        
        n_taxa_rat=len([t for t in rat[rank] if not t.startswith('un')])    
        TP=len([t for t in cami[rank] if t.split('|')[-1] in rat[rank]])
        FN=len(cami[rank])-TP
        FP=n_taxa_rat-TP
        
        stats[rank]['TP']=TP
        stats[rank]['FN']=FN
        stats[rank]['FP']=FP
        stats[rank]['precision']=TP/n_taxa_rat
        stats[rank]['sensitivity']=TP/(TP+FN)

        
        total_diff=0
        for s in rat[rank]:
            rat_lineage=rat[rank][s]['lineage']
            xr=rat[rank][s]['percent']
            xi=0
            if rat_lineage in cami[rank]:
                xi=cami[rank][rat_lineage]
            total_diff+=abs(xr-xi)
        for u in [v for v in cami[rank] if v.split('|')[-1] not in rat[rank]]:
            xi=cami[rank][u]
            total_diff+=xi
        stats[rank]['l1']=total_diff
        
    return stats
    
  
def evaluate_RAT_old(RAT_json, CAMI_profile):
    stats={}
    cami=make_cami_dict(CAMI_profile)
    rat=json.load(open(RAT_json))
    
    for rank in ranks:
        # del rat[rank]['unclassified'],rat[rank]['unmapped'] # don't count unmapped and unclassified as taxa
        stats[rank]={'l1': None, 
                 'TP': None, 
                 'precision': None,
                 'FP': None,
                 'FN': None,
                 'sensitivity': None}
        
        n_taxa_rat=len([t for t in rat[rank] if not t.startswith('un')])    
        TP=len([t for t in cami[rank] if t.split('|')[-1] in rat[rank]])
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
        for u in [v for v in cami[rank] if v.split('|')[-1] not in rat[rank]]:
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
                for sample in stats:
                    outf.write('{},{},{},{}\n'.format(sample,
                                                      rank,
                                                      metric,
                                                      stats[sample][rank][metric]))



if __name__=='__main__':
    # for m in ['wbin', 'wmetabat', 'nobin']:
    #     stats={}
    #     for s in [6,13,23,25,26,30,33,34,38,53]:
    #         stats['smp{}'.format(s)]=evaluate_RAT(
    #             '/net/phage/linuxhome/mgx/people/tina/RAT/benchmark/CAMI_II_mousegut_w_dm/'
    #             'profiles_all/smp{}.{}.json'.format(s,m), 
    #             '/net/phage/linuxhome/mgx/people/tina/CAMI/CAMI_II/mousegut/19122017_mousegut_scaffolds/'
    #             'taxonomic_profile_{}.txt'.format(s))
        
    
    #     write_evaluation(stats, '/net/phage/linuxhome/mgx/people/tina/RAT/benchmark/'
    #                      'CAMI_II_mousegut_w_dm/10_samples_{}_stats.csv'.format(m))
    
    
    # for m in ['kraken2', 'kaiju', 'centrifuge', 'bracken']:
    #     stats={}
    #     for s in [6,13,23,25,26,30,33,34,38,53]:
    #         stats['smp{}'.format(s)]=evaluate_RAT(
    #             '/net/phage/linuxhome/mgx/people/tina/RAT/new_read_classifier_results/'
    #             '20230321_profiles/20230321.smp{}.{}.0.001.json'.format(s,m), 
    #             '/net/phage/linuxhome/mgx/people/tina/CAMI/CAMI_II/mousegut/19122017_mousegut_scaffolds/'
    #             'taxonomic_profile_{}.txt'.format(s))
        
    
    #     write_evaluation(stats, '/net/phage/linuxhome/mgx/people/tina/RAT/'
    #                       'new_read_classifier_results/'
    #                       '20230321.10_samples_{}_stats_t.0.001.csv'.format(m))
        
    # for m in ['wbin', 'wmetabat', 'nobin','robust']:
    #     stats={}
    #     for s in [6,13,23,25,26,30,33,34,38,53]:
    #         stats['smp{}'.format(s)]=evaluate_RAT(
    #                 '/net/phage/linuxhome/mgx/people/tina/RAT/benchmark/CAMI_II_mousegut_w_dm/'
    #                 'profiles_all/only_classified/smp{}.{}_only_classified.json'.format(s,m), 
    #                 '/net/phage/linuxhome/mgx/people/tina/CAMI/CAMI_II/mousegut/19122017_mousegut_scaffolds/'
    #                 'taxonomic_profile_{}.txt'.format(s))
            
        
    #         write_evaluation(stats, '/net/phage/linuxhome/mgx/people/tina/RAT/benchmark/'
    #                           'CAMI_II_mousegut_w_dm/stats_all/10_samples_{}_oc_stats.csv'.format(m))

    
    ### Revisions
    
    ## Read classifiers
    
    for m in ['kraken2', 'kaiju', 'centrifuge', 'bracken']:
        stats={}
        for s in [2,3,5,8,10,12,14,15,18,19]:
            stats['plant{}'.format(s)]=evaluate_RAT(
                '../profiles/20230927.plant{}.{}.0.001.json'.format(s,m), 
                '../reference/gs_rhizosphere.noplasmids_{}.profile'.format(s))

        write_evaluation(stats, '../evaluations/'
                          '20231024.plant_{}_stats.0.001.csv'.format(m))
        for s in range(0,10):
            stats['marine{}'.format(s)]=evaluate_RAT(
                '../profiles/20230927.marine{}.{}.0.001.json'.format(s,m), 
                '../reference/gs_marine.noplasmids_{}.profile'.format(s))
        write_evaluation(stats, '../evaluations/'
                          '20231024.marine_{}_stats.0.001.csv'.format(m))          
        
    ## RAT
    for m in ['robust', 'sensitive', 'contig']:
        stats={}
        for s in [2,3,5,8,10,12,14,15,18,19]:
            stats['plant{}'.format(s)]=evaluate_RAT(
                '../profiles/20231020.plant{}.{}.1e-05.json'.format(s,m), 
                '../reference/gs_rhizosphere.noplasmids_{}.profile'.format(s))

        write_evaluation(stats, '../evaluations/'
                          '20231024.plant_{}_stats.0.001.csv'.format(m))
        for s in range(0,10):
            stats['marine{}'.format(s)]=evaluate_RAT(
                '../profiles/20231020.marine{}.{}.1e-05.json'.format(s,m), 
                '../reference/gs_marine.noplasmids_{}.profile'.format(s))
        write_evaluation(stats, '../evaluations/'
                          '20231024.marine_{}_stats.0.001.csv'.format(m))       
            
            