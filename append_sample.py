# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 10:08:49 2021

@author: acer
"""

import json

ranks=['superkingdom','phylum','class','order','family','genus','species']


def add_sample_to_table(json_file, outfile, sample, group):
    sample_dict=json.load(open(json_file))
    with open(outfile, 'a') as outf:
        for r in range(len(ranks)):
            rank=ranks[r]

            total_reads=sample_dict['total_reads']
            measure='Correctly classified'
            value=sample_dict['correct'][r]/total_reads
            line='{0},{1},{2},{3},{4}\n'.format(measure,
                                                rank,
                                                value,
                                                group,
                                                sample)
            outf.write(line)
            measure='Incorrectly classified'
            value=sample_dict['incorrect'][r]/total_reads
            line='{0},{1},{2},{3},{4}\n'.format(measure,
                                                rank,
                                                value,
                                                group,
                                                sample)
            outf.write(line)
            measure='Unclassified'
            value=(total_reads-sample_dict['correct'][r]-sample_dict['incorrect'][r])/total_reads
            line='{0},{1},{2},{3},{4}\n'.format(measure,
                                                rank,
                                                value,
                                                group,
                                                sample)
            # print(group)
            # print (sample)
            # print(total_reads-sample_dict['correct'][r]-sample_dict['incorrect'][r])
            outf.write(line)

if __name__=='__main__':
    # for i in [13]:
    for i in [6,13,23,25,26,30,33,34,38,53]:
        # for t in ['diamond']:
        for t in ["kaiju", "centrifuge", "kraken", "wbin", "wmetabat", 'nobin', 'robust']:
            json_file='../benchmark/r2c_mousegut_{}_ncbi{}.json'.format(t,i)
            # json_file='../benchmark/d2c_stats_mousegut_{}.txt'.format(i)
    
            outfile='../stats_CAMI2/20230127_CAMI2_w_unclassified_CORRECT.csv'
            sample=i
            if t in ["wbin", "wmetabat", 'nobin', 'robust']:
                group='r_{}'.format(t)
            else:
                group=t
            # group='diamond'
            add_sample_to_table(json_file, outfile, sample, group)