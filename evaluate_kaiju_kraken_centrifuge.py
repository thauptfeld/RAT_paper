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
        
def find_LCA_for_taxids(taxids, taxid2parent):
    list_of_lineages = []
    for taxid in taxids:
            
        try:
            lineage = find_lineage(taxid, taxid2parent)

            list_of_lineages.append(lineage)
        except:
            # The fastaid does not have an associated taxid for some reason.
            pass
        
    if len(list_of_lineages) == 0:
        return ('no taxid found', [])

    overlap = set.intersection(*map(set, list_of_lineages))
    for taxid in list_of_lineages[0]:
        if taxid in overlap:
            return (taxid, find_lineage(taxid, taxid2parent))



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
    


def import_centrifuge_results(results_file,taxid2parent,taxid2rank):
    read2classification={}
    with open(results_file) as inf:
        for line in inf:
            if not line.startswith('readID'):
                readID=line.split()[0]
                taxID=line.split()[2]
                if readID not in read2classification:
                    read2classification[readID]=set()
                read2classification[readID].add(taxID)
    read2lca={}
    for read in read2classification:
        (taxid, lineage) = find_LCA_for_taxids(read2classification[read], 
                                            taxid2parent)
        lineage=[rank.rstrip('*') for rank in lineage]
        ranks=[taxid2rank[rank] for rank in lineage]
        official_lineage=len(official_ranks)*['']
        for r in range(len(official_ranks)):
            if official_ranks[r] in ranks:
                official_lineage[r]=lineage[ranks.index(official_ranks[r])]
        read2lca[read+'/1']=official_lineage
        read2lca[read+'/2']=official_lineage
    return read2lca
         

        
def import_kaiju_assigned(kaiju_file, taxid2parent, taxid2rank):
    # better use lca
    read2classification={}
    with open(kaiju_file) as inf:
        for line in inf:
            readid=line.split()[1]
            tax_string=line.split()[2]
            
            (taxid, lineage) = find_LCA_for_taxids([tax_string], 
                                            taxid2parent)
            lineage=[rank.rstrip('*') for rank in lineage]
            ranks=[taxid2rank[rank] for rank in lineage]
            official_lineage=len(official_ranks)*['']
            for r in range(len(official_ranks)):
                if official_ranks[r] in ranks:
                    official_lineage[r]=lineage[ranks.index(official_ranks[r])]
            read2classification[readid+'/1']=official_lineage
            read2classification[readid+'/2']=official_lineage
            
    return read2classification


def import_kaiju_lca(kaiju_file, taxid2parent, taxid2rank):
    read2classification={}
    with open(kaiju_file) as inf:
        for line in inf:
            readid=line.split()[1]
            if not line.startswith('U'):
                taxids=set(line.split()[4].rstrip(',').split(','))
                taxids.add(line.split()[2])
            else:
                taxids=['0']
            
            (taxid, lineage) = find_LCA_for_taxids(taxids, 
                                            taxid2parent)
            lineage=[rank.rstrip('*') for rank in lineage]
            ranks=[taxid2rank[rank] for rank in lineage]
            official_lineage=len(official_ranks)*['']
            for r in range(len(official_ranks)):
                if official_ranks[r] in ranks:
                    official_lineage[r]=lineage[ranks.index(official_ranks[r])]
            read2classification[readid+'/1']=official_lineage
            read2classification[readid+'/2']=official_lineage
            
    return read2classification

def import_kraken_results(kraken_file, taxid2parent, taxid2rank):
    read2classification={}
    with open(kraken_file) as inf:
        for line in inf:
            readid=line.split()[1]
            taxids=[line.split('taxid ')[1].split(')')[0]]
            (taxid, lineage) = find_LCA_for_taxids(taxids, 
                                            taxid2parent)
            lineage=[rank.rstrip('*') for rank in lineage]
            ranks=[taxid2rank[rank] for rank in lineage]
            official_lineage=len(official_ranks)*['']
            for r in range(len(official_ranks)):
                if official_ranks[r] in ranks:
                    official_lineage[r]=lineage[ranks.index(official_ranks[r])]
            read2classification[readid+'/1']=official_lineage
            read2classification[readid+'/2']=official_lineage
    return read2classification



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
                    

def compare_results(CAMI_official, compare_official,total_reads):
    correct={'correct': len(official_ranks)*[0],
             'unclassified': len(official_ranks)*[0],
             'incorrect': len(official_ranks)*[0],
             'CAMI_unclassified_here_classified': len(official_ranks)*[0]}
    for read in compare_official:
        for i in range(len(CAMI_official[read])):
            if compare_official[read][i] == CAMI_official[read][i]:
                correct['correct'][i]+=1
            elif compare_official[read][i]!='' and compare_official[read][i]!=CAMI_official[read][i]:
                correct['incorrect'][i]+=1
            elif CAMI_official[read][i]=='' and compare_official[read][i]!='':
                correct['CAMI_unclassified_here_classified'][i]+=1
    for i in range(len(compare_official[read])):
        unclassified=total_reads-correct['correct'][i]-correct['incorrect'][i]
        correct['unclassified'][i]=unclassified
    return correct


def timestamp():
    now = datetime.datetime.now()
    str_ = '[{0}]'.format(now.strftime('%Y-%m-%d %H:%M:%S'))

    return str_
    
if __name__=='__main__':

#Kraken/Kaiju/Centrifuge benchmark:
#    load taxonomy data
    CAT_folder='/net/phage/linuxhome/mgx/people/bastiaan/CAT_prepare/CAT_prepare_20190108/'
    sample_folder=sys.argv[1]
    sample=sys.argv[2]
    smp=int(sample.split('smp')[-1])

    taxid2parent, taxid2rank = import_nodes(CAT_folder+'2019-01-08_taxonomy/nodes.dmp')
    print('{} load CAMI classification'.format(timestamp()))
    CAMI_r2c=import_CAMI_results('/net/phage/linuxhome/mgx/people/tina/CAMI/'
                    'CAMI_II/mousegut/19122017_mousegut_scaffolds/read_mappings/'
                    '{}_reads_mapping.tsv'.format(sample),
                    taxid2parent, taxid2rank)
    total_reads=len(CAMI_r2c)
    
    
    print('{} {}: Load centrifuge results'.format(timestamp(), sample))
    r2c_centrifuge=import_centrifuge_results(sample_folder +
                            'centrifuge/{}.centrifuge.out'.format(sample),
                            taxid2parent, taxid2rank)
    print('{} {}: Comparing centrifuge results to CAMI'.format(timestamp(), sample))
    correct=compare_results(CAMI_r2c, r2c_centrifuge, total_reads)
    correct['total_reads']=total_reads
    with open('/net/phage/linuxhome/mgx/people/tina/RAT/benchmark/'
              'kraken_kaiju_centrifuge/r2c_stats_mousegut_centrifuge{}.json'
              ''.format(smp), 'w+') as outf:
        outf.write(json.dumps(correct, indent=4))
    del r2c_centrifuge
    
    print('{} {}: Load kaiju results'.format(timestamp(), sample))
    r2c_kaiju=import_kaiju_lca(sample_folder +
                          'kaiju/{}.kaiju.out'.format(sample),
                          taxid2parent, taxid2rank)
    print('{} {}: Comparing kaiju results to CAMI'.format(timestamp(), sample))
    correct=compare_results(CAMI_r2c, r2c_kaiju, total_reads)
    correct['total_reads']=total_reads
    with open('/net/phage/linuxhome/mgx/people/tina/RAT/benchmark/'
              'kraken_kaiju_centrifuge/r2c_stats_mousegut_kaiju{}.json'
              ''.format(smp), 'w+') as outf:
        outf.write(json.dumps(correct, indent=4))
    del r2c_kaiju
    
    print('{} {}: Load kraken2 results'.format(timestamp(), sample))
    r2c_kraken=import_kraken_results(sample_folder +
                          'kraken2/{}.kraken2.out'.format(sample),
                          taxid2parent, taxid2rank)
    print('{} {}: Comparing kraken2 results to CAMI'.format(timestamp(), sample))
    correct=compare_results(CAMI_r2c, r2c_kraken, total_reads)
    correct['total_reads']=total_reads
    with open('/net/phage/linuxhome/mgx/people/tina/RAT/benchmark/'
              'kraken_kaiju_centrifuge/r2c_stats_mousegut_kraken{}.json'
              ''.format(smp), 'w+') as outf:
        outf.write(json.dumps(correct, indent=4))
    print('{} Done!'.format(timestamp()))




