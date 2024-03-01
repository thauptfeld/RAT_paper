# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 12:42:51 2022

@author: acer
"""

# import subprocess
import os
place='../rarefaction/wageningen/raw/'
# place='../rarefaction/CAMI/raw/'

# all_files={}
# for i in [1,5,10,20,30,40,50,60,70,80,90,100]:
#     tmp_files=os.listdir(('/hosts/linuxhome/phage/mgx/people/tina/RAT/'
#                              'benchmark/wageningen/rarefaction_profiles/{}/'
#                              'W19-1.rat.rf1.it001.csv'.format(i)))
    
#     all_files[i]=[f for f in tmp_files if 'W19-1' in f and 'min' in f]

# for i in all_files:
#     for file in all_files[i]:
#         filename=file.split('/')[-1]
#         file.stats=filename.split('.')

# cmd=['wc','-l',file]
# proc=subprocess.Popen(cmd, stdout=subprocess.PIPE)

# result=proc.stdout
# for p in proc.stdout:
#     p=str(p)
#     p=p[2:].split()[0]
#     if p.isdigit():
#         print('This file has {} lines'.format(p))

info=dict()
files_kaiju=[f for f in os.listdir(place) if 'rat_r.' in f and not 'nomin' in f]
for f in files_kaiju:
    with open(place+f) as inf:
        sample=f.split('.')[0]
        for line in inf:
            if line.startswith('#'):
                percent=int(line[1:].strip()[:-1])
                if sample not in info:
                    info[sample]={}
                info[sample][percent]={}
            else:
                lines=int(line.split()[0])
                iteration=line.split()[1].split('.')[-2]
                info[sample][percent][iteration]=lines

for sample in info:
    with open('../rarefaction/wageningen/{}.rat_r.min0.001.csv'.format(sample), 'w') as outf:
        outf.write('sample,iteration,percent,tool,taxa\n')
        for percent in info[sample]:
            for iteration in info[sample][percent]:
                text2write='{0},{1},{2},{3},{4}\n'.format(sample,
                                                iteration,
                                                percent,
                                                'rat_r',
                                                info[sample][percent][iteration])
                outf.write(text2write)
            