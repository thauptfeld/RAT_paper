# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 17:51:37 2021

@author: acer
"""

import os

tool_tables=[f for f in os.listdir('../revision/evaluations/') if f.startswith('20231009.plant')]
# smps=[6,13,23,25,26,30,33,34,38,53]
smps=[0,1,2,3,4,5,6,7,8,9]
# smps=[2,3,5,8,10,12,14,15,18,19]
outfile='../revision/evaluations/20231009.all_tools.plant.0.001.csv'
with open(outfile, 'w') as o:
    o.write('tool,metric,sample,rank,value\n')


for f in tool_tables:
    tool=f.split('_')[-2]
    with open('../revision/evaluations/'+f) as inf, open(outfile,'a') as outf:
        for line in inf:
            line=line.strip().split(',')
            metric=line[2]
            rank=line[1]
            sample=line[0]
            value=line[3]
            outf.write('{0},{1},{2},{3},{4}\n'.format(tool,
                                                      metric,
                                                      sample,
                                                      rank,
                                                      value))

