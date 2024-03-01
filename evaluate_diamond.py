#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 14:28:56 2020

@author: tina
"""

import os
import gzip
import sys
import json
import datetime

'''
This script is largely made up of reused code from various CAT scripts.
This particular script was made to compare RAT performance in read 
classification to diamond performance (diamond in blastx mode, top hit)

The classification idea: each read has an assigned taxid. The script finds the
official lineage of that taxid and returns it as a list. For each official 
rank, for each read, for each method (RAT, diamond), the script scores whether
on this rank the classification is correct (1), or incorrect (0). The resulting
fraction (sum of all reads on that rank/total number of reads) will be put into
a heatmap.

Details:
    import_CAMI_results: reads in the CAMI goldstart read classification and 
    returns a dicionary with taxid as key and a list of the official ranks of 
    the given taxid as values.
    
    parse_tabular_alignment: reads in the diamond file. Output is a dictionary
    called ORF2hits that contains the taxon belonging to the hit.
    
    find_LCA_for_ORF: translates the taxon id into a taxid and gets the lineage
    of said taxid.
    
    get_official_lineages: output is a dictionary that has the read ids as keys
    and a list of the official lineage ranks as value.

TO DO:
    - make a function that reads in RAT read2classification output in the same 
    way as the CAMI data
    - make a function that compares the RAT output and the diamond output to
    the official CAMI data and scores it.
    
    Also to do
    - shrink the diamond databases to only contain the top hits (unpigz should
    finish some time tomorrow)

'''




official_ranks=['superkingdom','phylum','class','order','family','genus','species']


def get_fastais2LCAtaxid_file(database_folder):
    fastaid2LCAtaxid_file = None

    if os.path.isdir(database_folder):
        for file_ in os.listdir(database_folder):
            if file_.endswith('fastaid2LCAtaxid'):
                fastaid2LCAtaxid_file = '{0}{1}'.format(
                        database_folder, file_)
    return fastaid2LCAtaxid_file


def parse_tabular_alignment(alignment_file):

    compressed = False
    if alignment_file.endswith('.gz'):
        compressed = True

        f1 = gzip.open(alignment_file, 'rb')
    else:
        f1 = open(alignment_file, 'r')

    ORF2hits = {}
    all_hits = set()

    ORF = 'first ORF'
    for line in f1:
        if compressed:
            line = line.decode('utf-8')

        line = line.rstrip().split('\t')

        if not line[0] == ORF:
            # A new ORF is reached. This is the top hit.
            ORF = line[0]
            ORF2hits[ORF] = []
            hit = line[1]
            ORF2hits[ORF].append(hit)
            all_hits.add(hit)

    f1.close()
                
    return ORF2hits, all_hits

def find_LCA_for_ORF(hits, fastaid2LCAtaxid, taxid2parent):
    list_of_lineages = []
    for hit in hits:
            
        try:
            taxid = fastaid2LCAtaxid[hit]
            lineage = find_lineage(taxid, taxid2parent)

            list_of_lineages.append(lineage)
        except:
            # The fastaid does not have an associated taxid for some reason.
            pass
        
    if len(list_of_lineages) == 0:
        return ('no taxid found ({0})'.format(';'.join([i[0] for i in hits])), [])

    overlap = set.intersection(*map(set, list_of_lineages))

    for taxid in list_of_lineages[0]:
        if taxid in overlap:
            return (taxid, lineage)



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


def import_fastaid2LCAtaxid(fastaid2LCAtaxid_file, all_hits):
    fastaid2LCAtaxid = {}

    with open(fastaid2LCAtaxid_file, 'r') as f1:
        for line in f1:
            line = line.rstrip().split('\t')

            if line[0] in all_hits:
                # Only include fastaids that are found in hits.
                fastaid2LCAtaxid[line[0]] = line[1]

    return fastaid2LCAtaxid

        
        
        

def find_lineage(taxid, taxid2parent, lineage=None):
    if lineage == None:
        lineage = list()
#    print(lineage)
#    print(type(lineage))
    lineage.append(taxid)

    if taxid2parent[taxid] == taxid:
        return lineage
    else:
        return find_lineage(taxid2parent[taxid], taxid2parent, lineage)


def get_official_lineages(ORF2hits, fastaid2LCAtaxid, taxid2parent, taxid2rank):
    official_ids={}

    for ORF in ORF2hits:
        hitlist=[ORF2hits[ORF][i].rsplit('_', 1)[0] for i in range(len(ORF2hits[ORF]))]
        (taxid, lineage) = find_LCA_for_ORF(hitlist, fastaid2LCAtaxid, 
                                            taxid2parent)

        lineage.reverse()
        lineage=[rank.rstrip('*') for rank in lineage]
        ranks=[taxid2rank[rank] for rank in lineage]
        official_lineage=len(official_ranks)*['']
        for r in range(len(official_ranks)):
            if official_ranks[r] in ranks:
                official_lineage[r]=lineage[ranks.index(official_ranks[r])]
    
        official_ids[ORF]=(official_lineage)
    return official_ids



def get_official_lineage2(read2classification, taxid2rank):
    read2official={}
    counts={'bin':0, 'contig': 0, 'contig_dm':0, 'read_dm':0}
    for read in read2classification:
        if read2classification[read]!=4*['']:
            read2official[read]=[]
            if read2classification[read][0]!='':
                lineage=read2classification[read][0].rstrip('*').split(';')
                counts['bin']+=1
            elif read2classification[read][1]!='':
                lineage=read2classification[read][1].rstrip('*').split(';')
                counts['contig']+=1
            elif read2classification[read][2]!='':
                lineage=read2classification[read][2].rstrip('*').split(';')
                counts['contig_dm']+=1
            elif read2classification[read][3]!='':
                lineage=read2classification[read][3].rstrip('*').split(';')
                counts['read_dm']+=1
            try:
                lineage=[rank.rstrip('*') for rank in lineage]
            except UnboundLocalError:
                print(read2classification[read])
                print
                sys.exit(1)
            ranks=[taxid2rank[rank] for rank in lineage]
            official_lineage=len(official_ranks)*['']
            for r in range(len(official_ranks)):
                if official_ranks[r] in ranks:
                    official_lineage[r]=lineage[ranks.index(official_ranks[r])]
            read2official[read]=official_lineage
    print('\n\n'+str(counts)+'\n\n')
    return read2official


def import_CAMI_results(CAMI_read_mapping, taxid2parent, taxid2rank):
    read2lineage={}
    with open(CAMI_read_mapping, 'r') as f1:
        for line in f1:
            # startswith # for CAMI II, startwith @ for CAMI I
            if not line.startswith('#') and not line.rstrip()=='':
                line=line.split('\t')
                readid, taxid = line[0], line[2]
                read2lineage[readid]=len(official_ranks)*['']
                lineage=find_lineage(taxid, taxid2parent)
                ranks=[taxid2rank[rank] for rank in lineage]

                for r in range(len(official_ranks)):
                    if official_ranks[r] in ranks:
                        read2lineage[readid][r]=lineage[ranks.index(official_ranks[r])]
    
    return read2lineage
    
    

def import_RAT_results(RAT_read_mapping):
    read2classification={}
    with open(RAT_read_mapping) as f1:
        read_id='first/1'
        for line in f1:
            if not line.startswith('#'):
                line=line.rstrip().split('\t')
#                print(line)
                if line[0]==read_id.split('/')[0]:
                    read_id=line[0]+'/1'
                else:
                    read_id=line[0]+'/2'
                bin_classification=''
                contig_classification=''
                contig_dm=''
                read_dm=''
                if not line[1]=='no taxid assigned':
                    bin_classification=line[2]
                    if bin_classification=='':
                        contig_classification=line[3]
                    if bin_classification=='' and contig_classification=='':
                        contig_dm=line[4]
                        if contig_dm=='':
                            read_dm=line[5].strip().replace('fw: ', '').replace('rev: ', '')
                read2classification[read_id]=[bin_classification, contig_classification,
                                              contig_dm, read_dm]
    return read2classification
                    

def compare_results(CAMI_official, compare_official):
    correct={'correct': len(official_ranks)*[0],
             'unclassified': len(official_ranks)*[0],
             'incorrect': len(official_ranks)*[0]}
    for read in compare_official:
        for i in range(len(compare_official[read])):
            if compare_official[read][i] == CAMI_official[read][i]:
                correct['correct'][i]+=1
            elif compare_official[read][i]=='':
                correct['unclassified'][i]+=1
            else:
                correct['incorrect'][i]+=1
    return correct

def compare_results_big(CAMI_official, compare_official):
    correct={}
    for read in compare_official:
        correct[read]=[]
        for i in range(len(compare_official[read])):
            if compare_official[read][i].strip() == CAMI_official[read][i]:
                correct[read].append('c')
            elif compare_official[read][i].strip()=='':
                correct[read].append('u')
            else:
                correct[read].append('i')
        add=7-len(correct[read])
        for i in range(add):
            correct[read].append('u')
    return correct


def timestamp():
    now = datetime.datetime.now()
    str_ = '[{0}]'.format(now.strftime('%Y-%m-%d %H:%M:%S'))

    return str_
    
if __name__=='__main__':

# # # Compare diamond results to official results
    
#     CAMI_folder='/net/phage/linuxhome/mgx/people/bastiaan/CAT_prepare/CAT_prepare_20190108/'
#     # benchmark='/net/phage/linuxhome/tina/RAT/benchmark/CAMI_high/'
#     # alignment_file=sys.argv[1]
#     # sample=alignment_file.split('/')[-1]
#     # smp=sample.split('_')[1].split('smp')[-1]
    
#     alignment_file='/net/phage/linuxhome/mgx/people/tina/CAT_NT/mousegut_smp13_reads_to_translated_nt.dmnd'
#     sample='smp13'
#     smp=13
    

#     taxid2parent, taxid2rank = import_nodes(CAMI_folder+'2019-01-08_taxonomy/nodes.dmp')
#     fastaid2LCAtaxid_file=get_fastais2LCAtaxid_file(CAMI_folder+'2019-01-08_CAT_database/')
    
    
#     print('load CAMI classification')
#     CAMI_r2c=import_CAMI_results('/net/phage/linuxhome/mgx/people/tina/CAMI/CAMI_II/'
#                     'mousegut/19122017_mousegut_scaffolds/read_mappings/smp{}_reads_mapping.tsv'.format(smp),
#                     taxid2parent, taxid2rank)
#     print('load diamond results')
#     read2hits, all_hits = parse_tabular_alignment(alignment_file)
#     fastaid2LCAtaxid = import_fastaid2LCAtaxid(fastaid2LCAtaxid_file, all_hits)
#     d2c_off=get_official_lineages(read2hits, fastaid2LCAtaxid, taxid2parent, taxid2rank)
    
#     print('{}: compare the diamond results to the supposed results'.format(sample))
#     correct=compare_results(CAMI_r2c, d2c_off)
#     correct['total_reads']=len(CAMI_r2c)
#     correct['unclassified_total']=len(CAMI_r2c)-len(d2c_off)
#     correct['corr_frac']=[i/(len(CAMI_r2c)) for i in correct['correct']]
    
    
#     with open('/home/tina/Documents/RAT/d2c_stats_mousegut_{}_against_nt.txt'.format(smp), 'w+') as outf:
#         outf.write(json.dumps(correct, indent=4))



#read2classification benchmark:
#    load taxonomy data
    # CAT_folder='/net/phage/linuxhome/mgx/people/bastiaan/CAT_prepare/CAT_prepare_20190108/'
    # alignment_file=sys.argv[1]
    # sample=alignment_file.split('/')[-1].split('_')[0]
    # smp=sample.split('_')[0].split('smp')[-1]
    # print(alignment_file)
    # print(sample)
    # taxid2parent, taxid2rank = import_nodes(CAT_folder+'2019-01-08_taxonomy/nodes.dmp')
    # print('{} load CAMI classification'.format(timestamp()))
    # CAMI_r2c=import_CAMI_results('/net/phage/linuxhome/mgx/people/tina/CAMI/'
                    # 'CAMI_II/mousegut/19122017_mousegut_scaffolds/read_mappings/'
                    # '{}_reads_mapping.tsv'.format(sample),
                    # taxid2parent, taxid2rank)

    
    
    # print('{} {}: Load RAT results'.format(timestamp(), sample))
    # r2c=import_RAT_results('/net/phage/linuxhome/mgx/people/tina/RAT/benchmark/'
                           # 'CAMI_II_mousegut_w_dm/smp_vs_pooled/'
                           # '{}_nobin_clean.read2classification.txt'.format(sample))
    # print('{} {}: get official lineages for RAT results'.format(timestamp(), sample))
    # r2off=get_official_lineage2(r2c, taxid2rank)

# #    free some memory
    # del r2c
    # print('{} {}: compare the RAT results to the supposed results'.format(timestamp(), sample))
    # correct=compare_results_big(CAMI_r2c, r2off)
    # correct['total_reads']=len(CAMI_r2c)
    # # correct['unclassified_total']=len(CAMI_r2c)-len(r2off)
    # # correct['corr_frac']=[i/(len(CAMI_r2c)) for i in correct['correct']]

    # with open('/net/phage/linuxhome/mgx/people/tina/RAT/benchmark/'
              # 'CAMI_II_mousegut_w_dm/evaluation/r2c_stats_mousegut_sankey_t{}.json'
              # ''.format(smp), 'w+') as outf:
        # outf.write(json.dumps(correct))
    # print('{} Done!'.format(timestamp()))

    # del r2off, CAMI_r2c







### Get lineages for a few ORFs
    CAT_folder='/net/mgx/linuxhome/mgx/people/bastiaan/phage-files/CAT_prepare/CAT_prepare_20210430/'
    alignment_file=sys.argv[1]
    taxid2parent, taxid2rank = import_nodes(CAT_folder+'CAT_taxonomy.2021-04-30/nodes.dmp')
    fastaid2LCAtaxid_file=get_fastais2LCAtaxid_file(CAMI_folder+'CAT_database.2021-04-30/')
    read2hits, all_hits = parse_tabular_alignment(alignment_file)
    all_hits2=set([i[:-2] for i in all_hits])
    fastaid2LCAtaxid = import_fastaid2LCAtaxid(fastaid2LCAtaxid_file, all_hits2)
    d2c_off=get_official_lineages(read2hits, fastaid2LCAtaxid, taxid2parent, taxid2rank)



# # # Compare diamond results to official results -> NT database
    
    
#     # benchmark='/net/phage/linuxhome/tina/RAT/benchmark/CAMI_high/'
#     alignment_file=sys.argv[1]
#     sample=alignment_file.split('/')[-1].split('_')[1]
#     smp=sample.split('smp')[-1]
    
#     # alignment_file='/net/phage/linuxhome/mgx/people/tina/CAT_NT/mousegut_smp13_reads_to_translated_nt.dmnd'
#     # sample='smp13'
#     # smp=13
    
#     print('load CAMI classification')
#     CAMI_folder='/net/phage/linuxhome/mgx/people/bastiaan/CAT_prepare/CAT_prepare_20190108/'
#     taxid2parent, taxid2rank = import_nodes(CAMI_folder+'2019-01-08_taxonomy/nodes.dmp')
#     CAMI_r2c=import_CAMI_results('/net/phage/linuxhome/mgx/people/tina/CAMI/CAMI_II/'
#                     'mousegut/19122017_mousegut_scaffolds/read_mappings/smp{}_reads_mapping.tsv'.format(smp),
#                     taxid2parent, taxid2rank)
    
    
#     print('load diamond results')
#     CAMI_folder='/net/phage/linuxhome/mgx/people/tina/CAT_NT/'
#     taxid2parent, taxid2rank = import_nodes(CAMI_folder+'20210928_CAT_nt_taxonomy/nodes.dmp')
#     fastaid2LCAtaxid_file=get_fastais2LCAtaxid_file(CAMI_folder+'20210928_CAT_nt/')
#     read2hits, all_hits = parse_tabular_alignment(alignment_file)
#     all_hits2=set([i[:-2] for i in all_hits])
#     fastaid2LCAtaxid = import_fastaid2LCAtaxid(fastaid2LCAtaxid_file, all_hits2)
#     d2c_off=get_official_lineages(read2hits, fastaid2LCAtaxid, taxid2parent, taxid2rank)
    
#     print('{}: compare the diamond results to the supposed results'.format(sample))
#     correct=compare_results(CAMI_r2c, d2c_off)
#     correct['total_reads']=len(CAMI_r2c)
#     correct['unclassified_total']=len(CAMI_r2c)-len(d2c_off)
#     correct['corr_frac']=[i/(len(CAMI_r2c)) for i in correct['correct']]
    
    
#     with open('/home/tina/Documents/RAT/d2c_stats_mousegut_{}_against_nt.txt'.format(smp), 'w+') as outf:
#         outf.write(json.dumps(correct, indent=4))
