#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 09:55:05 2021

@author: tina
"""

"""
Turn a kraken-style report into a CAMI-style taxonomic profile
"""

import json
import sys

off_ranks=['superkingdom','phylum','class','order','family','genus','species']
off_CAMI=['D','P','C','O','F','G','S']

# Minimum is fraction, not percent!!
def turn_kraken_into_cami(report_file, taxid2parent, taxid2rank, output_file, minimum=0):
    tax_dict={}
    most_recent={'D': '', 'P': '', 'C': '', 'O': '', 'F': '', 'G': ''}
    with open(report_file) as report:
        for line in report:
            if line.split()[4]=='1':
                classified_reads=float(line.split()[1])
                # print(classified_reads)
            
            not_off_rank=line.split()[3]
            if not_off_rank in most_recent:
                taxid=line.split()[4]
                most_recent[not_off_rank]=taxid
            
            if not_off_rank in off_CAMI:
                line=line.split(maxsplit=5)
                rank=off_ranks[off_CAMI.index(not_off_rank)]
                reads=float(line[1])
                taxid=line[4]
                while taxid not in taxid2parent:
                    try:
                        not_off_rank=off_CAMI[0:off_CAMI.index(not_off_rank)][-1]
                        taxid=most_recent[not_off_rank]
                    except IndexError:
                        print(line)
                        print(not_off_rank)
                        sys.exit(1)
                    
                lineage=find_lineage(taxid, taxid2parent)
                
                off_lineage=len(off_ranks)*['']
                ranks=[taxid2rank[r] for r in lineage]
                for r in range(len(off_ranks)):
                    if off_ranks[r] in ranks:
                        off_lineage[r]=lineage[ranks.index(off_ranks[r])]
                new_lineage='|'.join(off_lineage).rstrip('|')
                if rank not in tax_dict:
                    tax_dict[rank]={}
                if taxid not in tax_dict[rank]:
                    tax_dict[rank][taxid]={'percent': 0}
                tax_dict[rank][taxid]['percent']+=reads/classified_reads
                tax_dict[rank][taxid]['taxid']=taxid
                tax_dict[rank][taxid]['lineage']=new_lineage
                tax_dict[rank][taxid]['taxon']=''
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



def check_problem(kraken_file, taxid2parent):
    outd=set()
    with open(kraken_file) as inf:
        for line in inf:
            taxid=line.split()[4]
            rank=line.split()[3]
            if (not taxid in taxid2parent) and rank=='F':
                outd.add(taxid)
    return outd



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



        
        
if __name__=='__main__':
    # taxid2parent, taxid2rank = import_nodes('/net/phage/linuxhome/mgx/people/'
    #                                         'bastiaan/CAT_prepare/CAT_prepare'
    #                                         '_20200304/2020-03-04_taxonomy/'
    #                                         'nodes.dmp')
    taxid2parent, taxid2rank = import_nodes('/net/mgx/linuxhome/mgx/people/'
                                            'bastiaan/phage-files/CAT_prepare/CAT_prepare'
                                            '_20190108/2019-01-08_taxonomy/'
                                            'nodes.dmp')
    # taxid2parent['12884']='1'
    # taxid2rank['12884']='superkingdom'
    # taxid2parent['469587']='816'
    # taxid2parent['469586']='816'
    # taxid2parent['1834200']='375288'
    # taxid2parent['665937']='207244'
    # taxid2rank['469587']='species'
    # taxid2rank['469586']='species'
    # taxid2rank['1834200']='species'
    # taxid2rank['665937']='species'
    # for i in [6,13,23,25,26,30,33,34,38,53]:
    #     turn_kraken_into_cami('/net/phage/linuxhome/mgx/people/tina/RAT/'
    #                           'new_read_classifier_results/smp{}/'
    #                           'centrifuge/smp{}.centrifuge.kreport.txt'.format(i,i), 
    #                           taxid2parent, taxid2rank, 
    #                           '/net/phage/linuxhome/mgx/people/tina/RAT/'
    #                           'new_read_classifier_results/'
    #                           '20230321.smp{}.centrifuge.0.001.profile'.format(i), 0.00001)
    #     # turn_kraken_into_cami('/net/phage/linuxhome/mgx/people/tina/RAT/'
    #     #                       'new_read_classifier_results/smp{}/'
    #     #                       'centrifuge/smp{}.centrifuge.kreport.txt'.format(i,i), 
    #     #                       taxid2parent, taxid2rank, 
    #     #                       '/net/phage/linuxhome/mgx/people/tina/RAT/'
    #     #                       'new_read_classifier_results/'
    #     #                       '20230130.smp{}.centrifuge.nomin.profile'.format(i), 0)


    #     turn_kraken_into_cami('/net/phage/linuxhome/mgx/people/tina/RAT/'
    #                           'new_read_classifier_results/smp{}/'
    #                           'kraken2/smp{}.kraken2_report_bracken_species.out'.format(i,i), 
    #                           taxid2parent, taxid2rank, 
    #                           '/net/phage/linuxhome/mgx/people/tina/RAT/'
    #                           'new_read_classifier_results/'
    #                           '20230321.smp{}.bracken.0.001.profile'.format(i), 0.00001)
    #     # turn_kraken_into_cami('/net/phage/linuxhome/mgx/people/tina/RAT/'
    #     #                       'new_read_classifier_results/smp{}/'
    #     #                       'kraken2/smp{}.kraken2_report_bracken_species.out'.format(i,i), 
    #     #                       taxid2parent, taxid2rank, 
    #     #                       '/net/phage/linuxhome/mgx/people/tina/RAT/'
    #     #                       'new_read_classifier_results/'
    #     #                       '20230130.smp{}.bracken.nomin.profile'.format(i), 0)
        
        
    #     turn_kraken_into_cami('/net/phage/linuxhome/mgx/people/tina/RAT/'
    #                           'new_read_classifier_results/smp{}/'
    #                           'kraken2/smp{}.kraken2_report.out'.format(i,i), 
    #                           taxid2parent, taxid2rank, 
    #                           '/net/phage/linuxhome/mgx/people/tina/RAT/'
    #                           'new_read_classifier_results/'
    #                           '20230321.smp{}.kraken2.0.001.profile'.format(i), 0.00001)
    #     # turn_kraken_into_cami('/net/phage/linuxhome/mgx/people/tina/RAT/'
    #     #                       'new_read_classifier_results/smp{}/'
    #     #                       'kraken2/smp{}.kraken2_report.out'.format(i,i), 
    #     #                       taxid2parent, taxid2rank, 
    #     #                       '/net/phage/linuxhome/mgx/people/tina/RAT/'
    #     #                       'new_read_classifier_results/'
    #     #                       '20230130.smp{}.kraken2.nomin.profile'.format(i), 0)




        
# # Minimum is fraction, not percent!!
#     for i in [19,22,23]:
#         for j in range(1,7):
#             print('W{}-{}'.format(i,j))
#             turn_kraken_into_cami('/net/phage/linuxhome/mgx/people/tina/RAT/'
#                                   'benchmark/CAMI_II_mousegut_w_dm/wageningen_kaiju_kraken_centr/'
#                                   'W{}-{}/centrifuge/'
#                                   'W{}-{}.centrifuge.kreport.txt'.format(i,j,i,j), 
#                                   taxid2parent, taxid2rank, 
#                                   '/net/phage/linuxhome/mgx/people/tina/RAT/'
#                                   'benchmark/CAMI_II_mousegut_w_dm/wageningen_kaiju_kraken_centr/'
#                                   'W{}-{}.centrifuge.0.1.profile'.format(i,j), 0.001)
#             turn_kraken_into_cami('/net/phage/linuxhome/mgx/people/tina/RAT/'
#                                   'benchmark/CAMI_II_mousegut_w_dm/wageningen_kaiju_kraken_centr/'
#                                   'W{}-{}/kraken2/'
#                                   'W{}-{}.kraken2_report.out'.format(i,j,i,j), 
#                                   taxid2parent, taxid2rank, 
#                                   '/net/phage/linuxhome/mgx/people/tina/RAT/'
#                                   'benchmark/CAMI_II_mousegut_w_dm/wageningen_kaiju_kraken_centr/'
#                                   'W{}-{}.kraken.0.1.profile'.format(i,j), 0.001)
#             turn_kraken_into_cami('/net/phage/linuxhome/mgx/people/tina/RAT/'
#                                   'benchmark/CAMI_II_mousegut_w_dm/wageningen_kaiju_kraken_centr/'
#                                   'W{}-{}/kraken2/'
#                                   'W{}-{}.kraken2_report_bracken_species.out'.format(i,j,i,j), 
#                                   taxid2parent, taxid2rank, 
#                                   '/net/phage/linuxhome/mgx/people/tina/RAT/'
#                                   'benchmark/CAMI_II_mousegut_w_dm/wageningen_kaiju_kraken_centr/'
#                                   'W{}-{}.bracken.0.1.profile'.format(i,j), 0.001)
        
        
    path_to_results=('/net/phage/linuxhome/dutilh-group/tina/RAT/revision/'
                     'read_classifiers_new_db_plant/results/')
    for i in [2,3,5,8,10,12,14,15,18,19]:
        print('Sample plant{}...'.format(i))
        print('\tCentrifuge')
        turn_kraken_into_cami(path_to_results+'plant{}/'
                              'centrifuge/plant{}.centrifuge.kreport.txt'.format(i,i), 
                              taxid2parent, taxid2rank, 
                              path_to_results+
                              '20230927.plant{}.centrifuge.0.001.profile'.format(i), 0.00001)
        # turn_kraken_into_cami('/net/phage/linuxhome/mgx/people/tina/RAT/'
        #                       'new_read_classifier_results/smp{}/'
        #                       'centrifuge/smp{}.centrifuge.kreport.txt'.format(i,i), 
        #                       taxid2parent, taxid2rank, 
        #                       '/net/phage/linuxhome/mgx/people/tina/RAT/'
        #                       'new_read_classifier_results/'
        #                       '20230130.smp{}.centrifuge.nomin.profile'.format(i), 0)

        print('\tBracken')
        turn_kraken_into_cami(path_to_results+'plant{}/'
                              'kraken2/plant{}.kraken2_report_bracken_species.out'.format(i,i), 
                              taxid2parent, taxid2rank, 
                              path_to_results+
                              '20230927.plant{}.bracken.0.001.profile'.format(i), 0.00001)
        # turn_kraken_into_cami('/net/phage/linuxhome/mgx/people/tina/RAT/'
        #                       'new_read_classifier_results/smp{}/'
        #                       'kraken2/smp{}.kraken2_report_bracken_species.out'.format(i,i), 
        #                       taxid2parent, taxid2rank, 
        #                       '/net/phage/linuxhome/mgx/people/tina/RAT/'
        #                       'new_read_classifier_results/'
        #                       '20230130.smp{}.bracken.nomin.profile'.format(i), 0)
        
        print('\tKraken')
        turn_kraken_into_cami(path_to_results+'plant{}/'
                              'kraken2/plant{}.kraken2_report.out'.format(i,i), 
                              taxid2parent, taxid2rank, 
                              path_to_results+
                              '20230927.plant{}.kraken2.0.001.profile'.format(i), 0.00001)
        # turn_kraken_into_cami('/net/phage/linuxhome/mgx/people/tina/RAT/'
        #                       'new_read_classifier_results/smp{}/'
        #                       'kraken2/smp{}.kraken2_report.out'.format(i,i), 
        #                       taxid2parent, taxid2rank, 
        #                       '/net/phage/linuxhome/mgx/people/tina/RAT/'
        #                       'new_read_classifier_results/'
        #                       '20230130.smp{}.kraken2.nomin.profile'.format(i), 0)
        
        
        
        
    path_to_results=('/net/phage/linuxhome/dutilh-group/tina/RAT/revision/'
                     'read_classifiers_new_db_marine/results/')
    for i in range(0,10):
        print('Sample marine{}...'.format(i))
        print('\tCentrifuge')
        turn_kraken_into_cami(path_to_results+'marine{}/'
                              'centrifuge/marine{}.centrifuge.kreport.txt'.format(i,i), 
                              taxid2parent, taxid2rank, 
                              path_to_results+
                              '20230927.marine{}.centrifuge.0.001.profile'.format(i), 0.00001)
        # turn_kraken_into_cami('/net/phage/linuxhome/mgx/people/tina/RAT/'
        #                       'new_read_classifier_results/smp{}/'
        #                       'centrifuge/smp{}.centrifuge.kreport.txt'.format(i,i), 
        #                       taxid2parent, taxid2rank, 
        #                       '/net/phage/linuxhome/mgx/people/tina/RAT/'
        #                       'new_read_classifier_results/'
        #                       '20230130.smp{}.centrifuge.nomin.profile'.format(i), 0)

        print('\tBracken')
        turn_kraken_into_cami(path_to_results+'marine{}/'
                              'kraken2/marine{}.kraken2_report_bracken_species.out'.format(i,i), 
                              taxid2parent, taxid2rank, 
                              path_to_results+
                              '20230927.marine{}.bracken.0.001.profile'.format(i), 0.00001)
        # turn_kraken_into_cami('/net/phage/linuxhome/mgx/people/tina/RAT/'
        #                       'new_read_classifier_results/smp{}/'
        #                       'kraken2/smp{}.kraken2_report_bracken_species.out'.format(i,i), 
        #                       taxid2parent, taxid2rank, 
        #                       '/net/phage/linuxhome/mgx/people/tina/RAT/'
        #                       'new_read_classifier_results/'
        #                       '20230130.smp{}.bracken.nomin.profile'.format(i), 0)
        
        print('\tKraken')
        turn_kraken_into_cami(path_to_results+'marine{}/'
                              'kraken2/marine{}.kraken2_report.out'.format(i,i), 
                              taxid2parent, taxid2rank, 
                              path_to_results+
                              '20230927.marine{}.kraken2.0.001.profile'.format(i), 0.00001)
        # turn_kraken_into_cami('/net/phage/linuxhome/mgx/people/tina/RAT/'
        #                       'new_read_classifier_results/smp{}/'
        #                       'kraken2/smp{}.kraken2_report.out'.format(i,i), 
        #                       taxid2parent, taxid2rank, 
        #                       '/net/phage/linuxhome/mgx/people/tina/RAT/'
        #                       'new_read_classifier_results/'
        #                       '20230130.smp{}.kraken2.nomin.profile'.format(i), 0)