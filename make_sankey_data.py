# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 12:04:47 2021

@author: acer
"""

import os
import json

ranks=['superkingdom','phylum','class','order','family','genus','species']


output_data={'d':{'correct':7*[0],
                  'incorrect':7*[0],
                  'unclassified':7*[0]}, 
             'r':{'correct':7*[0],
                  'incorrect':7*[0],
                  'unclassified':7*[0]}}
total={'d':0,
       'r':0}


for f in os.listdir('./benchmark/'):
    group=f[0]
    smp_dict=json.load((open('./benchmark/{}'.format(f))))
    total_reads=smp_dict['total_reads']
    total[group]+=total_reads
    for r in range(len(ranks)):
        output_data[group]['correct'][r]+=smp_dict['correct'][r]
        output_data[group]['incorrect'][r]+=smp_dict['incorrect'][r]
        unclassified=total_reads-smp_dict['correct'][r]-smp_dict['incorrect'][r]
        output_data[group]['unclassified'][r]+=unclassified

with open('./20211118_sankey.csv') as outf:
    outf.write('')
    for group in output_data:
        for state in output_data[group]:
            for r in range(len(ranks)):
                