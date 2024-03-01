# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 15:39:50 2021

@author: acer
"""

import json
import sys
import random

inf=sys.argv[1]
reads=json.load(open(inf))

count={}


for r in reads:
    if not r=='total_reads':
        pattern=''.join(list(reads[r]))
        if pattern not in count:
            count[pattern]=0
        count[pattern]+=1



with open('/net/phage/linuxhome/mgx/people/tina/RAT/benchmark/'
          'CAMI_II_mousegut_w_dm/evaluation/{}_patterns.json'.format(sys.argv[2]),
          'w') as outf:
    outf.write(json.dumps(count))