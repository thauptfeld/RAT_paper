# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 10:19:42 2023

@author: Ernestina Hauptfeld
"""

import sys
import datetime
import gzip
import os
import json

official_ranks=['superkingdom','phylum','class','order','family','genus','species']

def timestamp():
    now = datetime.datetime.now()
    str_ = '[{0}]'.format(now.strftime('%Y-%m-%d %H:%M:%S'))

    return str_

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



def parse_tabular_alignment(alignment_files):
    ORF2hits = {}
    all_hits = set()
    
    
    for alignment_file in alignment_files:
        
        compressed = False
        if alignment_file.endswith('.gz'):
            compressed = True
    
            f1 = gzip.open(alignment_file, 'rb')
        else:
            f1 = open(alignment_file, 'r')
    
        print("Parsing {}...".format(alignment_file))
    
        ORF = 'first ORF'
        for line in f1:
            if compressed:
                line = line.decode('utf-8')
    
            line = line.rstrip().split('\t')
    
            
            # A new ORF is reached. This is the top hit.
            ORF = line[0]
            if ORF not in ORF2hits:
                ORF2hits[ORF] = []
            hit = line[1]
            ORF2hits[ORF].append(hit)
            all_hits.add(hit)
    
        f1.close()
                
    return ORF2hits, all_hits



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
    
    
    
def import_CAMI_results_marine(CAMI_read_mapping, taxid2parent, taxid2rank):
    read2lineage={}
    missingkeys=0

    with open(CAMI_read_mapping, 'r') as f1:
        for line in f1:
            # startswith # for CAMI II, startwith @ for CAMI I
            if not line.startswith('#') and not line.rstrip()=='':
                line=line.split('\t')
                readid, taxid = line[0], line[2]
                read2lineage[readid]=len(official_ranks)*['']
                try:
                    lineage=find_lineage(taxid, taxid2parent)
                except KeyError:
                    lineage=['1']
                    missingkeys+=1
                ranks=[taxid2rank[rank] for rank in lineage]

                for r in range(len(official_ranks)):
                    if official_ranks[r] in ranks:
                        read2lineage[readid][r]=lineage[ranks.index(official_ranks[r])]

    return read2lineage

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
        
        
def import_kaiju_assigned(kaiju_file, taxid2parent, taxid2rank):
    # better use lca
    read2classification={}
    with open(kaiju_file) as inf:
        for line in inf:
            readid=line.split()[1]
            tax_string=line.split()[2]
            
            (taxid, lineage) = find_LCA_for_taxids([tax_string], 
                                            taxid2parent)
            lineage=[rank.replace('*','') for rank in lineage]
            ranks=[taxid2rank[rank] for rank in lineage]
            official_lineage=len(official_ranks)*['']
            for r in range(len(official_ranks)):
                if official_ranks[r] in ranks:
                    official_lineage[r]=lineage[ranks.index(official_ranks[r])]
            read2classification[readid+'/1']=official_lineage
            read2classification[readid+'/2']=official_lineage
            
    return read2classification


def compare_results(CAMI_official, compare_official,total_reads, tool='diamond'):
    correct={}
    
    for read in CAMI_official:
        correct[read]={}
        if read in compare_official:
          
            
            for i in range(7):
                if compare_official[read][i] == CAMI_official[read][i]:
                    correct[read][i]='correct'
                elif compare_official[read][i]!='' and compare_official[read][i]!=CAMI_official[read][i]:
                    correct[read][i]='incorrect'
                else:
                    correct[read][i]='unclassified'
        else:
            correct[read]["0"]='unclassified'
    return correct


def get_fastais2LCAtaxid_file(database_folder):
    fastaid2LCAtaxid_file = None

    if os.path.isdir(database_folder):
        for file_ in os.listdir(database_folder):
            if file_.endswith('fastaid2LCAtaxid'):
                fastaid2LCAtaxid_file = '{0}{1}'.format(
                        database_folder, file_)
    return fastaid2LCAtaxid_file


def import_fastaid2LCAtaxid(fastaid2LCAtaxid_file, all_hits):
    fastaid2LCAtaxid = {}

    with open(fastaid2LCAtaxid_file, 'r') as f1:
        for line in f1:
            line = line.rstrip().split('\t')

            if line[0] in all_hits:
                # Only include fastaids that are found in hits.
                fastaid2LCAtaxid[line[0]] = line[1]

    return fastaid2LCAtaxid


def get_official_lineages(ORF2hits, fastaid2LCAtaxid, taxid2parent, taxid2rank):
    official_ids={}

    for ORF in ORF2hits:
        hitlist=[ORF2hits[ORF][i] for i in range(len(ORF2hits[ORF]))]
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


if __name__=='__main__':
    
    CAT_folder='/net/mgx/linuxhome/mgx/people/bastiaan/phage-files/CAT_prepare/CAT_prepare_20190108/'
    
    
    
    sample=sys.argv[1]
    path_to_cami='/net/phage/linuxhome/dutilh-group/tina/CAMI_II/'
    taxid2parent, taxid2rank = import_nodes(CAT_folder+'2019-01-08_taxonomy/nodes.dmp')
    taxid2parent['469587']='816'
    taxid2parent['469586']='816'
    taxid2parent['1834200']='375288'
    taxid2parent['665937']='207244'
    taxid2rank['469587']='species'
    taxid2rank['469586']='species'
    taxid2rank['1834200']='species'
    taxid2rank['665937']='species'
    
    
    
    alignment_file1='/net/phage/linuxhome/mgx/people/tina/RAT/revision/marine/diamond/marine{}_R1.10000.diamond'.format(sample)
    alignment_file2='/net/phage/linuxhome/mgx/people/tina/RAT/revision/marine/diamond/marine{}_R2.10000.diamond'.format(sample)
    

    fastaid2LCAtaxid_file=get_fastais2LCAtaxid_file(CAT_folder+'2019-01-08_CAT_database/')
    
    sample_folder='/net/phage/linuxhome/dutilh-group/tina/RAT/revision/read_classifiers_new_db_marine/results/marine{}/'.format(sample)
    
    
    # cami_r2c_file=(path_to_cami+'marine/simulation_short_read/'
    #                               '2018.08.15_09.49.32_sample_{}/reads/'
    #                               'reads_mapping.10000.tsv'.format(sample))
    
    # yuge={}
    # print('{} load CAMI classification'.format(timestamp()))
    
    # CAMI_r2c=import_CAMI_results_marine(cami_r2c_file,
    #                 taxid2parent, taxid2rank)
    # total_reads=len(CAMI_r2c)

    print('{} {}: Load kaiju results'.format(timestamp(), sample))
    r2c_kaiju=import_kaiju_assigned(sample_folder +
                          'kaiju/marine{}.10000.kaiju.out'.format(sample),
                          taxid2parent, taxid2rank)
    with open('/home/tina/Documents/RAT/marine0_kaiju_lineages.json', 'w') as f1:
        f1.write(json.dumps(r2c_kaiju, indent=4))    
    # print('{} {}: Comparing kaiju results to CAMI'.format(timestamp(), sample))
    # yuge['kaiju']=compare_results(CAMI_r2c, r2c_kaiju, total_reads)
    # yuge['total_reads']=total_reads

    print('{} {}{}: Load diamond results'.format(timestamp(), 'marine', sample))
    read2hits, all_hits = parse_tabular_alignment([alignment_file1, alignment_file2])
    print('{} Importing fastaid to LCA taxid'.format(timestamp()))
    fastaid2LCAtaxid = import_fastaid2LCAtaxid(fastaid2LCAtaxid_file, all_hits)
    print('{} Getting official lineages for diamond results'.format(timestamp()))
    d2c_off=get_official_lineages(read2hits, fastaid2LCAtaxid, taxid2parent, taxid2rank)
    with open('/home/tina/Documents/RAT/marine0_diamond_lineages.json', 'w') as f1:
        f1.write(json.dumps(d2c_off, indent=4))
    # print('{} Comparing diamond results to CAMI'.format(timestamp()))
    # yuge['diamond']=compare_results(CAMI_r2c, d2c_off, total_reads)
    
    # with open('/home/tina/Documents/RAT/kaiju_diamond_correct.json', 'w') as f1:
    #     f1.write(json.dumps(yuge, indent=4))
    
    reads_to_look_at=[]
    kaiju_diamond=json.load(open('/home/tina/Documents/RAT/kaiju_diamond_correct.json'))
    for read_id in kaiju_diamond['diamond']:
        # if kaiju_diamond['diamond'][read_id]['0']=='unclassified' and kaiju_diamond['kaiju'][read_id]['0']!='unclassified':
        #     reads_to_look_at.append(read_id)
        try:
            if kaiju_diamond['kaiju'][read_id]['5']=='correct' and kaiju_diamond['diamond'][read_id]['5']=='incorrect':
                reads_to_look_at.append(read_id)
        except KeyError:
            continue
    with open('/home/tina/Documents/RAT/kaiju_diamond_correct_reads.txt', 'w') as outf:
        for r in reads_to_look_at:
          outf.write(f'{r}\n')  
        