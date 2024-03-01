# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 12:32:18 2022

@author: acer
"""

with open('../time_ram/monitor_small.txt') as inf:
    # text=inf.read().strip().split('\n')
    for line in inf:
        line=line.strip()
        # print(line)
        if 'python' in line:
            print('\t'.join(line.split()[0:3]))
            
            
