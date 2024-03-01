# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 09:59:51 2022

@author: acer
"""

import json
import subprocess

samples=['W19-1','W19-2','W19-3','W19-4','W19-5','W19-6',
          'W22-1','W22-2','W22-3','W22-4','W22-5','W22-6',
          'W23-1','W23-2','W23-3','W23-4','W23-5','W23-6'
    ]
comples=['W22-1','W22-2','W22-3','W22-4','W22-5','W22-6'
    ]

dRep_clusters=json.load(open('/home/tina/Documents/RAT/bio_example/dRep_clusters.json'))
dRep_genomes={}
# print(len(dRep_clusters))
for d in dRep_clusters:
    for e in dRep_clusters[d]:
        dRep_genomes[e]=d
        
# print(len(dRep_genomes))

for comp in comples:
    contig2bin=json.load(open('/net/phage/linuxhome/mgx/people/tina/wageningen/'
                              '2018/{}.contig2bin.json'.format(comp)))
    for smp in samples:
        bam_file=('/net/phage/linuxhome/mgx/people/tina/wageningen/'
                  '2018/bam/{}.fasta.{}_R1.fq.gz.bwamem.sorted'.format(comp,smp))
        smp_vs_comp={}
        print('{} vs {}'.format(smp,comp))
        # open alignment file with samtools, so that we don't need pysam dependency
        cmd=['samtools', 'view', bam_file]
        proc=subprocess.Popen(cmd, stdout=subprocess.PIPE)
        contig=''
        # check whether it is a paired bam file by checking the flags for "read paired"
        # flag and setting paired to true if it finds it
        while contig!='*':
            for read in proc.stdout:
                read=read.decode("utf-8").rstrip().split('\t')
                read_id, flag, contig, score=read[0], int(read[1]), read[2], int(read[4])
                if contig in contig2bin and score>=2 and len(bin(flag))<14:
                    if contig in contig2bin and contig2bin[contig] in dRep_genomes:
                        if contig2bin[contig] not in smp_vs_comp:
                            smp_vs_comp[contig2bin[contig]]=0
                        smp_vs_comp[contig2bin[contig]]+=1
        with open('/home/tina/tmp/test{}vs{}.json'.format(smp,comp), 'w') as outf:
            outf.write(json.dumps(smp_vs_comp, indent=4))
                        