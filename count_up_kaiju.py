#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 13:34:07 2022

@author: tina
"""

### KAIJU

# import os
import sys

# path_to_rcs='/hosts/linuxhome/phage/mgx/people/tina/RAT/benchmark/wageningen/read_classifiers/'
path_to_rcs=('/hosts/linuxhome/phage/mgx/people/tina/RAT/benchmark/'
             'CAMI_II_mousegut_w_dm/kraken_kaiju_centrifuge/')

abundances={}

# dir_list=[d for d in os.listdir(path_to_rcs) if d.startswith('W')]
dir_list=[sys.argv[1]]
weird_ids=[]


print('Making taxon dict...')
for s in dir_list:
    abundances[s]={'u': {'num_reads':0,
                          'frac_reads':0}}
    rat=open('{}{}/kaiju/{}.kaiju.names.out'.format(path_to_rcs,s,s)).read().strip().split('\n')
    total=len(rat)
    for r in rat:
        if r.startswith('U'):
            abundances[s]['u']['num_reads']+=1
            abundances[s]['u']['frac_reads']+=1/total
        else:
            try:
                taxon=r.split('\t')[7]
                if taxon not in abundances[s]:
                    abundances[s][taxon]={'num_reads':0,
                                          'frac_reads':0}
                abundances[s][taxon]['num_reads']+=1
                abundances[s][taxon]['frac_reads']+=1/total
            except IndexError:
                weird_ids.append(r.split('\t')[2])


    print('Writing files...')
    with open('{}{}.kaiju.manually_counted.csv'.format(path_to_rcs,s), 'w') as outf1, open('{}{}.kaiju.manually_counted.0.001.csv'.format(path_to_rcs,s), 'w') as outf2:
        for taxon in abundances[s]:
            outf1.write('{},{},{}\n'.format(taxon,
                                            abundances[s][taxon]['num_reads'],
                                            abundances[s][taxon]['frac_reads']))
            if abundances[s][taxon]['frac_reads']>0.00001:
                outf2.write('{},{},{}\n'.format(taxon,
                                            abundances[s][taxon]['num_reads'],
                                            abundances[s][taxon]['frac_reads']))
                




# ### RAT
# #/hosts/linuxhome/phage/mgx/people/tina/wageningen/2018/W19-1/W19-1.RAT.complete.abundance.txt
# import os
# path_to_rcs='/hosts/linuxhome/phage/mgx/people/tina/wageningen/2018/'
# dir_list=[d for d in os.listdir(path_to_rcs) if d.startswith('W')]
# # dir_list=['W19-1']
# dir_list.sort()

# for s in dir_list:
#     rat=open('{}{}/{}.RAT.complete.abundance.txt'.format(path_to_rcs,s,s)).read().strip().split('\n')[2:]
#     with open('{}{}.RAT.complete.abundance.0.001.csv'.format(path_to_rcs,s), 'w') as outf1:
#         for r in rat:
#             if float(r.split('\t')[2])>=0.00001:
#                 outf1.write('{}\n'.format(','.join(r.split('\t')[0:3])))
        








# ### CENTRIFUGE

# import os

# path_to_rcs='/hosts/linuxhome/phage/mgx/people/tina/RAT/benchmark/wageningen/read_classifiers/'

# abundances={}

# dir_list=[d for d in os.listdir(path_to_rcs) if d.startswith('W')]
# # dir_list=['W23-6']
# weird_ids=[]

# print('Making taxon dict...')
# for s in dir_list:
#     abundances[s]={'taxa': 0,
#                    'all':0}
#     rat=open('{}{}/centrifuge/{}.centrifuge.kreport.txt'.format(path_to_rcs,s,s)).read().strip().split('\n')
#     total=0
#     for r in rat:
#         taxon=r.split()[4]
        
#         if taxon in ['0','1']:
#             total += int(r.split()[1])
#             abundances[s]['taxa']+=1
#             if taxon=='0':
#                 abundances[s]['U']=int(r.split()[2])
#         else:
#             uniq_reads=float(r.split()[2])
#             if uniq_reads>0:
#                 abundances[s]['all']+=1
#             if uniq_reads/total>0.00001:
#                 abundances[s]['taxa']+=1
#         abundances[s]['uncl']=abundances[s]['U']/total
            


# print('Writing file...')
# with open('{}all.centrifuge.manually_counted.csv'.format(path_to_rcs), 'w') as outf1:
#     outf1.write('sample,taxa_all,taxa_0.001,unclassified\n')
#     for s in abundances:
#         outf1.write('{},{},{},{}\n'.format(s,
#                                         abundances[s]['all'],
#                                         abundances[s]['taxa'],
#                                         abundances[s]['uncl']))




# ### Kraken

# import os

# path_to_rcs='/hosts/linuxhome/phage/mgx/people/tina/RAT/benchmark/wageningen/read_classifiers/'

# abundances={}

# dir_list=[d for d in os.listdir(path_to_rcs) if d.startswith('W')]
# # dir_list=['W23-6']
# weird_ids=[]

# print('Making taxon dict...')
# for s in dir_list:
#     abundances[s]={'taxa': 0,
#                    'all':0}
#     rat=open('{}{}/kraken2/{}.kraken2_report.out'.format(path_to_rcs,s,s)).read().strip().split('\n')
#     total=0
#     for r in rat:
#         taxon=r.split()[4]
        
#         if taxon in ['0','1']:
#             total += int(r.split()[1])
#             abundances[s]['taxa']+=1
#             if taxon=='0':
#                 abundances[s]['U']=int(r.split()[2])
#         else:
#             uniq_reads=float(r.split()[2])
#             if uniq_reads>0:
#                 abundances[s]['all']+=1
#             if uniq_reads/total>0.00001:
#                 abundances[s]['taxa']+=1
#         abundances[s]['uncl']=abundances[s]['U']/total
            


# print('Writing file...')
# with open('{}all.kraken2.manually_counted.csv'.format(path_to_rcs), 'w') as outf1:
#     outf1.write('sample,taxa_all,taxa_0.001,unclassified\n')
#     for s in abundances:
#         outf1.write('{},{},{},{}\n'.format(s,
#                                         abundances[s]['all'],
#                                         abundances[s]['taxa'],
#                                         abundances[s]['uncl']))
        
        
        
        
        
        


### bracken

# import os

# path_to_rcs='/hosts/linuxhome/phage/mgx/people/tina/RAT/benchmark/wageningen/read_classifiers/'

# abundances={}

# dir_list=[d for d in os.listdir(path_to_rcs) if d.startswith('W')]
# # dir_list=['W23-6']
# weird_ids=[]

# print('Making taxon dict...')
# for s in dir_list:
#     abundances[s]={'taxa': 0,
#                    'all':0}
    
#     rat=open('{}{}/kraken2/{}.kraken2_report.out'.format(path_to_rcs,s,s)).read().strip().split('\n')
#     total=0
#     for r in rat:
#         taxon=r.split()[4]
        
#         if taxon in ['0','1']:
#             total += int(r.split()[1])
#     # print(total)
#     rat=open('{}{}/kraken2/{}.kraken2_report_bracken_species.out'.format(path_to_rcs,s,s)).read().strip().split('\n')
    
#     for r in rat:
#         taxon=r.split()[4]
#         all_reads=float(r.split()[1])
        
#         uniq_reads=float(r.split()[2])
        
#         if taxon=='1':
#             # print(total)
#             uncl=(total-all_reads)/total
            
#         if uniq_reads>0:
#             abundances[s]['all']+=1
#         if uniq_reads/total>0.00001:
#             abundances[s]['taxa']+=1
#         abundances[s]['uncl']=uncl
            


# print('Writing file...')
# with open('{}all.bracken.manually_counted.csv'.format(path_to_rcs), 'w') as outf1:
#     outf1.write('sample,taxa_all,taxa_0.001,unclassified\n')
#     for s in abundances:
#         outf1.write('{},{},{},{}\n'.format(s,
#                                         abundances[s]['all'],
#                                         abundances[s]['taxa'],
#                                         abundances[s]['uncl']))