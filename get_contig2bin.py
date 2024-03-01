# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 10:45:17 2022

@author: acer
"""
import os
import json

samples=['W19-1','W19-2','W19-3','W19-4','W19-5','W19-6',
         'W22-1','W22-2','W22-3','W22-4','W22-5','W22-6',
         'W23-1','W23-2','W23-3','W23-4','W23-5','W23-6'
    ]


for s in samples:
    print(s)
    contig2bin={}
    path_to_dir=('/net/phage/linuxhome/mgx/people/tina/wageningen/2018/'
              '{}/bins/DASTool/'.format(s))
    list_of_bins=os.listdir(path_to_dir)
    for b in list_of_bins:
        with open(path_to_dir+b) as inf:
            for line in inf:
                if line.startswith('>'):
                    contig_id=line[1:].rstrip()
                    contig2bin[contig_id]=b
    json.dump(contig2bin, open('/net/phage/linuxhome/mgx/people/tina/wageningen/'
                               '2018/{}.contig2bin.json'.format(s), 'w'))
    