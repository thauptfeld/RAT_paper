# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 11:18:18 2021

@author: acer
"""

tools=['bracken', 'centrifuge', 'kaiju', 'kraken2',
       'robust', 'wbin', 'wmetabat', 'nobin']
samples=['6','13','23','25', '26','30','33','34','38','53']
ranks=['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

sens_prec_dict={}
for t in tools:
    sens_prec_dict[t]={}
    
for t in sens_prec_dict:
    for s in samples:
        sens_prec_dict[t][s]={}
        
for t in sens_prec_dict:
    for s in samples:
        for r in ranks:
            sens_prec_dict[t][s][r]={}

all_tools=open('../all_tools_stats_0.1.csv').read().strip().split('\n')[1:]           
for t in all_tools:
    t=t.split(',')
    tool, metric, sample, rank = t[0], t[1], t[2], t[3]
    if metric in ['precision', 'sensitivity']:
        sens_prec_dict[tool][sample][rank][metric]=t[4]
        
with open('../all_tools_CAMI2_sens_prec_0.1.csv', 'w') as outf:
    outf.write('tool,sample,rank,precision,sensitivity\n')
    for t in sens_prec_dict:
        for s in sens_prec_dict[t]:
            for r in sens_prec_dict[t][s]:
                outf.write('{0},{1},{2},{3},{4}\n'.format(t,s,r,
                                                          sens_prec_dict[t][s][r]['precision'],
                                                          sens_prec_dict[t][s][r]['sensitivity']))