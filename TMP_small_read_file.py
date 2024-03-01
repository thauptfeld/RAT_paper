# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 11:42:15 2023

@author: haup0007
"""

import sys


potentials=[str(p) for p in range(0,10000)]


with open(sys.argv[1]) as inf, open('{}.10000.tsv'.format(sys.argv[:-4]), 'w') as outf:
    i=0
    for line in inf:
        if i%1000000==0:
            print(f'Processed {i} lines...')
        read_id=line.split('/')[0][3:]
        if read_id in potentials:
            outf.write(line)
        i+=1