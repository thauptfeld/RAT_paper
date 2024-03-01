# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 11:02:59 2023

@author: Ernestina Hauptfeld
"""

import os
import json
from collections import Counter


marine=['../checkdata/'+f for f in os.listdir('../checkdata/') if 'marine' in f]
plant=['../checkdata/'+f for f in os.listdir('../checkdata/') if 'plant' in f]
mouse=['../checkdata/'+f for f in os.listdir('../checkdata/') if 'mousegut' in f]

m={}
p={}
mo={}
for f in marine:
    tmp=Counter(json.loads('{' + open(f).read()))
    m=dict(Counter(m)+tmp)
    



marine



with open('../evaluations/20231028.centrifuge_data_summarized.tsv') as outf:
    outf.write()
    


dict1 = {}
dict2 = {'b': 5, 'c': 15, 'd': 25}

counter1 = Counter(dict1)
counter2 = Counter(dict2)

merged_counter = counter1 + counter2

merged_dict = dict(merged_counter)

print(merged_dict)
