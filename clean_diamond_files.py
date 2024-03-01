# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 15:40:29 2021

@author: acer
"""
import sys 

filename=sys.argv[1]


with open(filename) as inf, open('{}.clean'.format(filename), 'w') as outf:
    for line in inf:
        outf.write(line.replace('/1', '').replace('/2', ''))