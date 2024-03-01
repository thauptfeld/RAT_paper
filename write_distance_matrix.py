# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 18:25:20 2021

@author: acer
"""

import json

json_file='/home/tina/Documents/RAT/unifrac/wageningen_robust.json'
distance_matrix=json.load(open(json_file))
samples_cami=[t for t in distance_matrix if t.startswith('taxon')]
samples_rest=[t for t in distance_matrix if t not in samples_cami]
samples_cami.sort()
samples_rest.sort()
samples_all=samples_cami+samples_rest

with open('/home/tina/Documents/RAT/unifrac/wageningen_distance_matrix_robust.csv', 'w') as outf:
    outf.write(','+','.join(samples_all)+'\n')
    for s in samples_all:
        line=s
        for t in samples_all:
            line+=','+str(distance_matrix[s][t])
        line+='\n'
        outf.write(line)
        