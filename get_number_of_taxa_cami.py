# -*- coding: utf-8 -*-
"""
Created on Sat Apr  9 09:42:10 2022

@author: acer
"""

import os
import sys

folder=('/hosts/linuxhome/phage/mgx/people/tina/RAT/benchmark/'
        'CAMI_II_mousegut_w_dm/rarefaction/wmin/')

outdict={}
smp=sys.argv[1]

for i in [1,5,10,20,30,40,50,60,70,80,90]:
    outdict[i]={}
    files=[f for f in os.listdir(folder+str(i)) if 'cami' in f and '{}_reads'.format(smp) in f]
    for f in files:
        path=folder+str(i)+'/'+f
        outdict[i][f]=int(open(path).read().strip())

with open(folder + '{}.taxa.cami.txt'.format(smp), 'w') as outf:
    for i in outdict:
        outf.write('#{}%\n'.format(i))
        for f in outdict[i]:
            outf.write('{}\t{}\n'.format(outdict[i][f],f))
    