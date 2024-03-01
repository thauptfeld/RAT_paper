# -*- coding: utf-8 -*-
"""
Created on Sun Sep 11 11:24:12 2022

@author: acer
"""

import json

ranks=['superkingdom','phylum','class','order','family','genus','species']

with open('../per_rank_classified_cami.csv', 'w') as outf:
    outf.write('sample,rank,classified\n')

    for i in [6,13,23,25,26,30,33,34,38,53]:

        sample='smp{}'.format(i)
    
        profile=json.load(open('../profiles/{}.wbin_only_classified.json'.format(sample)))
            
        for rank in ranks:
            classified_total=0
            for taxon in profile[rank]:
                classified_total+=profile[rank][taxon]['percent']
            outf.write('{},{},{}\n'.format(sample,rank,classified_total))
            
            
with open('../per_rank_classified_wageningen.csv', 'w') as outf:
    outf.write('sample,rank,classified\n')

    for i in ['W19-1','W19-2','W19-3','W19-4','W19-5','W19-6',
              'W22-1','W22-2','W22-3','W22-4','W22-5','W22-6',
              'W23-1','W23-2','W23-3','W23-4','W23-5','W23-6']:

        sample=i
    
        profile=json.load(open('../profiles_wageningen/{}.0.001.json'.format(sample)))
            
        for rank in ranks:
            classified_total=0
            for taxon in profile[rank]:
                classified_total+=profile[rank][taxon]['percent']
            outf.write('{},{},{}\n'.format(sample,rank,classified_total))