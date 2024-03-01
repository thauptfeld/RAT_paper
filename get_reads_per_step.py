# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 21:19:23 2021

@author: acer
"""

import json

# samples=['W19-1','W19-2','W19-3','W19-4','W19-5','W19-6',
#           'W22-1','W22-2','W22-3','W22-4','W22-5','W22-6',
#           'W23-1','W23-2','W23-3','W23-4','W23-5','W23-6']

samples=['smp13','smp23','smp25','smp26','smp30',
         'smp33','smp34','smp38','smp53','smp6']

output={}

for sample in samples:
    print(sample)
    output[sample]={'total_reads':0,
        'bin_reads':0,
        'contig_reads':0,
        'contig_dm':0,
        'read_dm':0,
        'weird':0
        }
    path=('/net/phage/linuxhome/mgx/people/tina/RAT/benchmark/'
          'CAMI_II_mousegut_w_dm/wbin/')
    with open(path+'{}_clean.read2classification.txt'.format(sample)) as inf:
        
        for line in inf:
            if not line.startswith('#'):
                output[sample]['total_reads']+=1
                line=line.rstrip().split('\t')
                if len(line)>2:
                    if line[2]!='':
                        output[sample]['bin_reads']+=1
                        if len(line)>3 and line[3]=='':
                            output[sample]['weird']+=1
                    elif line[2]=='' and line[3]!='':
                        output[sample]['contig_reads']+=1
                    elif line[2]=='' and line[3]=='' and line[4]!='':
                        output[sample]['contig_dm']+=1
                    elif line[2]=='' and line[3]=='' and line[4]=='' and line[5]!='':
                        output[sample]['read_dm']+=1
                    
with open('/home/tina/Documents/RAT/RAT_wbin_reads_per_step.json', 'w') as outf:
    outf.write(json.dumps(output, indent=4))

# samples=json.load(open('../benchmark/wageningen_reads_per_step.json'))
# with open('../benchmark/wageningen_reads_per_step.csv', 'w') as outf:
#     outf.write('sample,bin,contig,contig_dm,read_dm\n')
#     for s in samples:
#         outf.write('{0},{1},{2},{3},{4}\n'.format(s,
#                                                   samples[s]['bin_reads'],
#                                                   samples[s]['contig_reads'],
#                                                   samples[s]['contig_dm'],
#                                                   samples[s]['read_dm']))