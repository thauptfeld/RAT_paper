#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 09:31:47 2020

@author: tina
"""
import os

def merge_tax_tables(path_to_folders, output_folder):
    # path_to_folders is a list of folders that contain the tables you would
    # like to merge. This script goes through all the folders and merges the
    # tables per rank for all samples.
    
    # list of ranks
    levels=['superkingdom','phylum','class','order','family','genus','species', 'complete']        
#    levels=['species']    
    smp_list=['2-643084','2-643096', '2-643100','2-643162','2-643173','2-643174',
              '2-643175','2-643182','2-643184','2-643186','2-643187','2-643251',
              '2-643252','2-643253','2-643254','2-643255','2-643257','2-643258',
              '2-643262','2-643263','2-643264','2-643265','2-643266','2-643267']
    
    
    for level in levels:

        merge_dict={}
        for folder in path_to_folders:
            samples=[f for f in os.listdir(folder) if f.endswith(level + '.abundance.txt')]
            for smp in samples:
                table=open(folder+smp, 'r').read().strip().split('\n')[2:]
                
                for line in table:
                    taxon, fraction=line.split('\t')[0],float(line.split('\t')[1])
                    if taxon=='unclassified' or taxon=='1;131567':
                            taxon='1'
                    if not taxon in merge_dict:
                        merge_dict[taxon]={}
                    if not smp.split('.')[0] in merge_dict[taxon]:
                        merge_dict[taxon][smp.split('.')[0]]=0
                    merge_dict[taxon][smp.split('.')[0]]+=fraction

    
            with open(output_folder+'gosse_{}.merged_reads_no_singletons.csv'.format(level), 'w+') as op:
                op.write('# lineage')
                for s in smp_list:
                    op.write(',{}'.format(s))
                op.write('\n')
                for taxon in merge_dict:
                    # taxa that aren't singletons
                    if len(merge_dict[taxon])>1:
                        op.write(taxon.strip('*'))
                        for s in smp_list:
                            if s in merge_dict[taxon]:
                                op.write(',{}'.format(merge_dict[taxon][s]))
                            else:
                                op.write(',0')
                        op.write('\n')





if __name__=='__main__':
    uber=merge_tax_tables(['/home/tina/tmp/tmptmp/gosse_mapping_quality_2/RAT_out/'], 
    '/home/tina/tmp/tmptmp/gosse_mapping_quality_2/')


            
                
