#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 14:02:35 2021

@author: tina
"""
import json
off_ranks=['superkingdom','phylum','class','order','family','genus','species']

# Minimum is fraction, not percent!!
def turn_kaiju_into_cami(kaiju_prefix, taxid2parent, taxid2rank, 
                         classified_reads, domain_dict, output_file, minimum=0):
    tax_dict={'superkingdom': domain_dict}
    for level in off_ranks[1:]:
        tax_dict[level]={}
        with open('{}.{}.kaiju_table.tsv'.format(kaiju_prefix, level)) as inf:
            for line in inf:
                if not line.startswith('file'):
                    taxid=line.split()[3]
                    if taxid.isdigit() and taxid in taxid2parent:
                        reads=float(line.split()[2])
                        lineage=find_lineage(taxid, taxid2parent)
                        ranks=[taxid2rank[r] for r in lineage]
                        off_lineage=len(off_ranks)*['']
                        for r in range(len(off_ranks)):
                            if off_ranks[r] in ranks:
                                off_lineage[r]=lineage[ranks.index(off_ranks[r])]
                        new_lineage='|'.join(off_lineage).rstrip('|')
                        if taxid not in tax_dict[level]:
                            tax_dict[level][taxid]={'percent': 0}
                        tax_dict[level][taxid]['percent']+=reads/classified_reads
                        tax_dict[level][taxid]['taxid']=taxid
                        tax_dict[level][taxid]['lineage']=new_lineage
                        tax_dict[level][taxid]['taxon']=''
    with open(output_file, 'w+') as outf1, open(output_file.rsplit('.', 1)[0]+ '.json', 'w+') as outf2:
            outf1.write('@SampleID:\n'
                        '@Version:0.9.1\n'
                        '@Ranks:superkingdom|phylum|class|order|family|genus|species|strain\n\n'
                        '@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n')
            for rank in tax_dict:
                for taxon in tax_dict[rank]:
                    if (not taxon=='unmapped' and not taxon=='unclassified') and tax_dict[rank][taxon]['percent']>minimum:

                        outf1.write('{}\t{}\t{}\t{}\t{}\n'.format(taxon.replace('*',''),
                                                          rank, tax_dict[rank][taxon]['lineage'],
                                                          tax_dict[rank][taxon]['taxon'],
                                                          tax_dict[rank][taxon]['percent']))
            tax_dict_json={}
            for rank in tax_dict:
                tax_dict_json[rank]={}
                for taxid in tax_dict[rank]:
                    if tax_dict[rank][taxid]['percent']>minimum:
                        tax_dict_json[rank][taxid]=tax_dict[rank][taxid]
                    
            outf2.write(json.dumps(tax_dict_json, indent=4))


                        

def calculate_domain_from_phylum(kaiju_file, taxid2parent, taxid2rank, classified_reads):
    domains={}
    with open(kaiju_file) as inf:
        for line in inf:
            if not line.startswith('file'):
                taxid=line.split()[3]
                if taxid.isdigit() and taxid in taxid2parent:
                    reads=float(line.split()[2])
                    lineage=find_lineage(taxid, taxid2parent)
                    ranks=[taxid2rank[r] for r in lineage]
                    
                    sk_index=ranks.index('superkingdom')
                    parent=lineage[sk_index]
                    if not parent in domains:
                        domains[parent]={'percent': 0,
                                         'taxid': taxid,
                                         'lineage': parent,
                                         'taxon': ''}
                    domains[parent]['percent']+=reads/classified_reads
    return domains
    
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

def gen_classified_reads(kaiju_out):
    classified=0
    with open(kaiju_out) as inf:
        for line in inf:
            if line.startswith('C'):
                classified+=1
    return classified

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


if __name__=='__main__':
    # taxid2parent, taxid2rank = import_nodes('/net/phage/linuxhome/mgx/people/'
    #                                         'bastiaan/CAT_prepare/CAT_prepare'
    #                                         '_20200304/2020-03-04_taxonomy/'
    #                                         'nodes.dmp')
    taxid2parent, taxid2rank = import_nodes('/net/mgx/linuxhome/mgx/people/'
                                            'bastiaan/phage-files/CAT_prepare/CAT_prepare'
                                            '_20190108/2019-01-08_taxonomy/'
                                            'nodes.dmp')
    taxid2parent['469587']='816'
    taxid2parent['469586']='816'
    taxid2parent['1834200']='375288'
    taxid2parent['665937']='207244'
    taxid2rank['469587']='species'
    taxid2rank['469586']='species'
    taxid2rank['1834200']='species'
    taxid2rank['665937']='species'
    taxid2parent['12884']='1'
    taxid2rank['12884']='superkingdom'

    # Minimum is fraction, not percent!!
    
    
    
    # ### Groundwater
    # for i in [19,22,23]:
    #     for j in range(1,7):
    #         print('W{}-{}...'.format(i,j))
    #         print('\tCounting classified reads...')
    #         classified_reads=gen_classified_reads('/net/phage/linuxhome/mgx/people/'
    #                                               'tina/RAT/benchmark/CAMI_II_mousegut_w_dm/wageningen_kaiju_kraken_centr/'
    #                                               'W{}-{}/kaiju/W{}-{}.kaiju.out'.format(i,j,i,j))
    #         print('\tCalculating domain abundances...')
    #         domain=calculate_domain_from_phylum('/net/phage/linuxhome/mgx/people/tina/'
    #                                             'RAT/benchmark/CAMI_II_mousegut_w_dm/wageningen_kaiju_kraken_centr/'
    #                                             'W{}-{}/kaiju/W{}-{}.phylum.kaiju_table.tsv'.format(i,j,i,j), 
    #                                             taxid2parent, taxid2rank, classified_reads)
    #         print('\tWriting Kaiju output as CAMI table...')
    #         tax_dict=turn_kaiju_into_cami('/net/phage/linuxhome/mgx/people/tina/RAT/'
    #                                       'benchmark/CAMI_II_mousegut_w_dm/wageningen_kaiju_kraken_centr/W{}-{}/'
    #                                       'kaiju/W{}-{}'.format(i,j,i,j), taxid2parent, 
    #                                       taxid2rank, classified_reads, domain,
    #                                       '/net/phage/linuxhome/mgx/people/tina/RAT/'
    #                                       'benchmark/CAMI_II_mousegut_w_dm/wageningen_kaiju_kraken_centr/'
    #                                       'W{}-{}.kaiju.0.1.profile'.format(i,j), 0.001)
    
    
    # ### Mousegut
    # for i in [6,13,23,25,26,30,33,34,38,53]:
    #     print('smp{}...'.format(i))
    #     print('\tCounting classified reads...')
    #     classified_reads=gen_classified_reads('/net/phage/linuxhome/mgx/people/'
    #                                           'tina/RAT/new_read_classifier_results/'
    #                                           'smp{}/kaiju/smp{}.kaiju.out'.format(i,i))
    #     print('\tCalculating domain abundances...')
    #     domain=calculate_domain_from_phylum('/net/phage/linuxhome/mgx/people/tina/'
    #                                         'RAT/new_read_classifier_results/'
    #                                         'smp{}/kaiju/smp{}.phylum.kaiju_table.tsv'.format(i,i), 
    #                                         taxid2parent, taxid2rank, classified_reads)
    #     print('\tWriting Kaiju output as CAMI table...')
    #     tax_dict=turn_kaiju_into_cami('/net/phage/linuxhome/mgx/people/tina/RAT/'
    #                                   'new_read_classifier_results//smp{}/'
    #                                   'kaiju/smp{}'.format(i,i), taxid2parent, 
    #                                   taxid2rank, classified_reads, domain,
    #                                   '/net/phage/linuxhome/mgx/people/tina/RAT/'
    #                                   'new_read_classifier_results/smp{}/'
    #                                   '20230321.smp{}.kaiju.0.001.profile'.format(i,i), 0.00001)
        # print('\tWriting Kaiju output as CAMI table...')
        # tax_dict=turn_kaiju_into_cami('/net/phage/linuxhome/mgx/people/tina/RAT/'
        #                               'new_read_classifier_results//smp{}/'
        #                               'kaiju/smp{}'.format(i,i), taxid2parent, 
        #                               taxid2rank, classified_reads, domain,
        #                               '/net/phage/linuxhome/mgx/people/tina/RAT/'
        #                               'new_read_classifier_results/smp{}/'
        #                               '20230130.smp{}.kaiju.nomin.profile'.format(i,i), 0)
    
    
    ### Marine
    
    path_to_results=('/net/phage/linuxhome/dutilh-group/tina/RAT/revision/'
                      'read_classifiers_new_db_marine/results/')
    
    for i in range(0,10):
        print('marine{}...'.format(i))
        print('\tCounting classified reads...')
        classified_reads=gen_classified_reads(path_to_results+'marine{}/kaiju/'
                                              'marine{}.kaiju.out'.format(i,i))
        print('\tCalculating domain abundances...')
        domain=calculate_domain_from_phylum(path_to_results+'marine{}/kaiju/'
                                            'marine{}.phylum.kaiju_table.tsv'.format(i,i), 
                                            taxid2parent, taxid2rank, classified_reads)
        print('\tWriting Kaiju output as CAMI table...')
        tax_dict=turn_kaiju_into_cami(path_to_results+'marine{}/kaiju/marine{}'.format(i,i), 
                                      taxid2parent, 
                                      taxid2rank, classified_reads, domain,
                                      path_to_results+
                                      '20230927.marine{}.kaiju.0.001.profile'.format(i,i), 0.00001)
        
    
    
    ### Plant
    path_to_results=('/net/phage/linuxhome/dutilh-group/tina/RAT/revision/'
                     'read_classifiers_new_db_plant/results/')
    
    for i in [2,3,5,8,10,12,14,15,18,19]:
        print('plant{}...'.format(i))
        print('\tCounting classified reads...')
        classified_reads=gen_classified_reads(path_to_results+'plant{}/kaiju/'
                                              'plant{}.kaiju.out'.format(i,i))
        print('\tCalculating domain abundances...')
        domain=calculate_domain_from_phylum(path_to_results+'plant{}/kaiju/'
                                            'plant{}.phylum.kaiju_table.tsv'.format(i,i), 
                                            taxid2parent, taxid2rank, classified_reads)
        print('\tWriting Kaiju output as CAMI table...')
        tax_dict=turn_kaiju_into_cami(path_to_results+'plant{}/kaiju/plant{}'.format(i,i), 
                                      taxid2parent, 
                                      taxid2rank, classified_reads, domain,
                                      path_to_results+
                                      '20230927.plant{}.kaiju.0.001.profile'.format(i,i), 0.00001)