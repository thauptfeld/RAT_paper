#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 14:38:28 2020

@author: tina
"""

ranks=['superkingdom','phylum','class','order','family','genus','species']
samples=['6','13','23','25','26','30','33','34','38','53']

import json
def summarize_stats(list_d2c, list_r2c):
    stats={'d': {},
           'r': {}}
    for f in list_d2c:
        stats['d'][f.split('.')[0].split('_')[-1]]=json.load(open(f))
    for f in list_r2c:
        stats['r'][f.split('.')[0].split('_')[-1]]=json.load(open(f))
    return stats

def compare_d2r(stats):
    comp={}
    for r in range(len(ranks)):
        comp[ranks[r]]={}
        avg_d_cor=0
        avg_r_cor=0
        avg_d_incor=0
        avg_r_incor=0
        avg_d_uncl=0
        avg_r_uncl=0
        for m in stats['d']:
            avg_d_cor+=stats['d'][m][r]

def write_stats(stats):
    with open('./20211117_mousegut_w_dm.csv', 'w+') as outf:
        outf.write('measure')
        for state in stats:
            for sample in samples:
                outf.write(',{}_{}'.format(state,sample))
        outf.write('\n')
        for measure in ['correct', 'corr_frac', 'incorrect', 'unclassified']:
            for r in range(len(ranks)):
                outf.write('{}_{}'.format(measure, ranks[r]))
                for state in stats:
                    for sample in samples:
                        outf.write(','+str(stats[state][sample][measure][r]))
                outf.write('\n')
                
def write_melted_stats(stats, outfile):
    with open(outfile, 'w') as outf:
        outf.write('measure,rank,value,group,sample\n')
        for group in stats:
            for sample in stats[group]:
                for r in range(len(stats[group][sample]['corr_frac'])):
                    outf.write('Correctly classified,{0},{1},{2},{3}\n'
                               ''.format(ranks[r], 
                                         stats[group][sample]['corr_frac'][r],
                                         group, sample))
                for r in range(len(stats[group][sample]['incorrect'])):
                    incorr_frac=stats[group][sample]['incorrect'][r]/stats[group][sample]['total_reads']
                    outf.write('Incorrectly classified,{0},{1},{2},{3}\n'
                               ''.format(ranks[r], 
                                         incorr_frac,
                                         group, sample))