#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 09:31:47 2020

@author: tina
"""
import os
import json
# griftpark=['2-100807', '2-100808', '2-100809', '2-100810', '2-100811', '2-100812', 
#           '2-100813', '2-100814', '2-100815', '2-100816', '2-100817', '2-100818',
#           'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11',
#           'S11', 'S12', 'S13', 'S14']
# gosse = ['2-643084','2-643096', '2-643100','2-643162','2-643173','2-643174',
#         '2-643175','2-643182','2-643184','2-643186','2-643187','2-643251',
#         '2-643252','2-643253','2-643254','2-643255','2-643257','2-643258',
#         '2-643262','2-643263','2-643264','2-643265','2-643266','2-643267']
# water = [i[0:-4] for i in os.listdir('/hosts/linuxhome/phage/mgx/people/luc/RAT_output/') if i.endswith('.log')]
# water.sort()

env=['W19-1','W19-2','W19-3','W19-4','W19-5','W19-6',
     'W22-1','W22-2','W22-3','W22-4','W22-5','W22-6',
     'W23-1','W23-2','W23-3','W23-4','W23-5','W23-6']

# def merge_tax_tables(path_to_folders, output_folder):
    # path_to_folders is a list of folders that contain the tables you would
    # like to merge. This script goes through all the folders and merges the
    # tables per rank for all samples.
    
    # list of ranks
levels=['superkingdom','phylum','class','order','family','genus','species']        
#    levels=['species']    

def merge_tax_tables(path_to_folders, output_folder):
    smp_list=env 
    
    for level in levels:
    
        merge_dict={}
        for smp in smp_list:
            infile=[f for f in os.listdir(path_to_folders[0]+smp) if f.endswith('RAT.complete.abundance.txt')]
            for f in infile:
                table=open(path_to_folders[0]+smp+'/'+infile[0], 'r').read().strip().split('\n')[2:]
                ranks=[]
                for line in table:
                    taxon, fraction=line.split('\t')[0].strip().replace('*','').split(';'),float(line.split('\t')[2])
                    if len(line.split('\t'))>4:
                        ranks=line.split('\t')[4].strip().split(';')
                    
                    
                    if level in ranks:
                        r=ranks.index(level)
                        taxon=';'.join(taxon[0:r+1])
                    else:
                        taxon=';'.join(taxon)
                    
                    if taxon=='unclassified' or taxon=='1;131567':
                            taxon='1'
                    if not taxon in merge_dict:
                        merge_dict[taxon]={}
                    if not smp.split('.')[0] in merge_dict[taxon]:
                        merge_dict[taxon][smp.split('.')[0]]=0
                    merge_dict[taxon][smp.split('.')[0]]+=fraction
    
                with open(output_folder+'wageningen_{}.merged_read_frac.csv'.format(level), 'w+') as op:
                    for s in smp_list:
                        op.write(',{}'.format(s))
                    
                    op.write('\n')
                    for taxon in merge_dict:
                        if len(merge_dict[taxon])>=1:
                            op.write(taxon.strip('*'))
                            for s in smp_list:
                                if s in merge_dict[taxon]:
                                    op.write(',{}'.format(merge_dict[taxon][s]))
                                else:
                                    op.write(',0')
                            op.write('\n')





if __name__=='__main__':
    merge_tax_tables(['/hosts/linuxhome/phage/mgx/people/tina/wageningen/2018/'], 
    '/home/tina/Documents/RAT/profiles_wageningen/taxa_tables/')


            
                
