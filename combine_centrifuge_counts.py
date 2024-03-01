# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 11:18:44 2023

@author: haup0007
"""


import os


files=['../checkdata/' + f for f in os.listdir('../checkdata/') if f.endswith('txt')]

all_values={'marine': {},
            'plant': {},
            'mousegut': {},
            'groundwater': {}}
env=''

for f in files:
    if 'marine' in f:
        env='marine'
    elif 'plant' in f:
        env='plant'
    elif 'smp' in f:
        env='mousegut'
    elif '-' in f or 'W2' in f:
        env='groundwater'
            
    with open(f) as inf:        
        for line in inf:
            line=line.strip().split()
            
            if int(line[1]) not in all_values[env]:
                all_values[env][int(line[1])]=0
            all_values[env][int(line[1])]+=int(line[0])

total={'marine': 0,
       'plant': 0,
       'mousegut': 0,
       'groundwater': 0}

for env in all_values:
    for key in all_values[env]:
        total[env]+=all_values[env][key]


            
with open('../checkdata/centrifuge_values_normalized.csv','w') as outf:
    outf.write('environment,matchLength,counts\n')
    for env in all_values:
        for key in all_values[env]:
            value=all_values[env][key]/total[env]
            outf.write(f'{env},{key},{value}\n')