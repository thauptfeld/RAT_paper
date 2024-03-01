# -*- coding: utf-8 -*-
"""
Created on Sat Apr  9 12:37:59 2022

@author: acer
"""
import os

folder='../rarefaction/CAMI/'

list_of_files=[f for f in os.listdir(folder) if f.endswith('.csv')]

with open('../rarefaction/all.cami.csv', 'w') as outf:
    outf.write('sample,iteration,percent,tool,taxa\n')
    for f in list_of_files:
        content=open(folder+f).read().strip().split('\n')[1:]
        for line in content:
            outf.write(line + '\n')
        