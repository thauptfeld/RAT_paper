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
    for i in [19]:
    # for i in [0,1,2,3,4,5,6,7,8,9]:
    # for i in [2,3,5,8,10,12,14,15,18,19]:
        # for t in ['diamond']:
        for t in ["kaiju", "centrifuge", "kraken", "contig", "sensitive", "robust"]:
            # json_file='../correctIncorrect_2/r2c_marine{}_nr_{}.json'.format(i,t)
            json_file='../correctIncorrect_2/r2c_plant{}_nr_{}.json'.format(i,t)
    
            # outfile='../evaluations/20231020.CAMI2_marine.csv'
            outfile='../evaluations/20231020_CAMI2_plant.csv'
            sample=i
            if t in ["contig", "sensitive", 'robust']:
                group='r_{}'.format(t)
            else:
                group=t
            add_sample_to_table(json_file, outfile, sample, group)
        # json_file_diamond='../diamond_2/d2c_stats_marine{}_against_nr_lca.txt'.format(i)
        # group='diamond'
        json_file_diamond='../diamond_2/d2c_stats_plant{}_against_nr_lca.txt'.format(i)
        add_sample_to_table(json_file_diamond, outfile, sample, 'diamond')