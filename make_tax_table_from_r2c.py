#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 14:24:44 2021

@author: tina
"""

import sys
import datetime


def r2c_to_complete(r2c_table, output, mode, complete2rank, taxid2rank):
    taxa={}
    
    total_reads=0
    with open(r2c_table, 'r') as r:
        for line in r:
            if not line.startswith('#') and not line=='':
                total_reads+=1
                line=line.split('\t')
            
                #If the read is classified, proceed
                if len(line)>2:
                
                    # if there is a bin classification, count towards bin taxon
                    if mode=='robust':
                        if line[2]!='':
                            taxon=line[2].replace('*','')
                        elif line[2]=='' and line[3]=='':
                            taxon='unmapped'
                    # if there is no bin classification, but a contig classification
                    # count towards contig taxon
                    if mode=='contig' or line[2]=='':
                        if len(line)>3 and line[3]!='':
                            taxon=line[3].replace('*','')
                        if mode=='contig':
                            if len(line)>4 and line[3]=='':
                                if line[4]!='':
                                    taxon=line[4].replace('*','')
                                elif len(line)>5 and line[4]=='' and line[5]!='':
                                    taxon=line[5].replace('*','').strip().replace('fw: ', '').replace('rev: ', '')
                    if taxon=='':
                        taxon='unmapped'
                    
                # otherwise, the read is unclassified
                else:
                    taxon='unclassified'
                
                if taxon not in taxa:
                    taxa[taxon]=0
                taxa[taxon]+=1
            
    with open(output, 'w') as outf:
        outf.write('#taxon\tnumber of reads\tfraction of reads\ttaxon length\tlineage ranks\n')
        for taxon in taxa:
            if not taxon.startswith('un'):
                try:
                    ranks=complete2rank[taxon]
                except KeyError:
                    ranks=';'.join([taxid2rank[t] for t in taxon.split(';')])
                outf.write('{}\t{}\t{}\t\t{}\n'.format(taxon, taxa[taxon], taxa[taxon]/total_reads,ranks))
            else:
                
                outf.write('{}\t{}\t{}\n'.format(taxon, taxa[taxon], taxa[taxon]/total_reads))
                
                    
    return taxa
            

def complete_abundance_to_dict(complete_abundance_file):
    taxa={}
    with open(complete_abundance_file) as inf:
        for line in inf:
            if not line.startswith('#') and not line.startswith('un'):
                lineage=line.split('\t')[0].replace('*','')
                ranks=line.split('\t')[4].rstrip()
                taxa[lineage]=ranks
    return taxa

def import_nodes(nodes_dmp):
    
    taxid2parent = {}
    taxid2rank = {}

    with open(nodes_dmp, 'r') as f1:
        for line in f1:
            line = line.split('\t')

            taxid = line[0]
            parent = line[2]
            rank = line[4]

            taxid2parent[taxid] = parent
            taxid2rank[taxid] = rank

    return taxid2parent, taxid2rank


def timestamp():
    now = datetime.datetime.now()
    str_ = '[{0}]'.format(now.strftime('%Y-%m-%d %H:%M:%S'))

    return str_

                
                

if __name__=='__main__':
    CAT_folder='/net/mgx/linuxhome/mgx/people/bastiaan/phage-files/CAT_prepare/CAT_prepare_20190108/'
    revision='/net/phage/linuxhome/mgx/people/tina/RAT/revision/'
    print('{} Importing nodes...'.format(timestamp()))
    taxid2parent, taxid2rank = import_nodes(CAT_folder + '2019-01-08_taxonomy/'
                                            'nodes.dmp')
    taxid2parent['469587']='816'
    taxid2parent['469586']='816'
    taxid2parent['1834200']='375288'
    taxid2parent['665937']='207244'
    taxid2rank['469587']='species'
    taxid2rank['469586']='species'
    taxid2rank['1834200']='species'
    taxid2rank['665937']='species'
    
    env=sys.argv[1]
    smp=sys.argv[2]
    mode=sys.argv[3]
    
    if env=='marine':
        prefix='2018.08.15_09.49.32'
    elif env=='plant':
        prefix='2019.09.27_13.59.10'
    

    

    print('{} load complete abundance'.format(timestamp()))
    complete2rank=complete_abundance_to_dict('/net/phage/linuxhome/mgx/people/tina/CAMI_II/'
                    '{}/simulation_short_read/{}_sample_{}/RAT_NR/'
                    '{}{}.sensitive.complete.abundance.txt'.format(env,prefix,smp,env,smp))
    
    print('{} Making {} abundance table for sample {}{}'.format(timestamp(), mode, env, smp))
    r2c_mode=r2c_to_complete('/net/phage/linuxhome/mgx/people/tina/CAMI_II/'
                    '{}/simulation_short_read/{}_sample_{}/RAT_NR/'
                    '{}{}.read2classification.txt'.format(env,prefix,smp,env,smp),
                    '/net/phage/linuxhome/mgx/people/tina/CAMI_II/'
                    '{}/simulation_short_read/{}_sample_{}/RAT_NR/'
                    '{}{}.{}.complete.abundance.txt'.format(env,prefix,smp,env,smp,mode),
                    mode, complete2rank, taxid2rank)
    print('{} Done!'.format(timestamp()))
    