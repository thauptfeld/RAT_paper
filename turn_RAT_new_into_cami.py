#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 15:37:29 2020

@author: tina
"""
import json
import sys

# lineage       number of reads fraction of reads       taxon length       lineage ranks

off_ranks=['superkingdom','phylum','class','order','family','genus','species']

def summarize_RAT_for_CAMI(rat_table, output_file, minimum=0):
    complete_dict={}
    classified_read_frac=1
    with open(rat_table) as rat:
        for line in rat:
            pieces=line.rstrip().replace('*','').split('\t')
            if not line.startswith('#') and not line.startswith('un'):
                lineage=pieces[0].split(';')
                frac_reads=float(pieces[2])
                # ranks is piece 4 for new RAT tables, piece 5 for old RAT tables
                ranks=pieces[5].split(';')

                off_lineage=len(off_ranks)*['']
                for r in ranks:
                    if r in off_ranks:
                        off_lineage[off_ranks.index(r)]=lineage[ranks.index(r)]
                dict_key='|'.join(off_lineage)
                if dict_key not in complete_dict:
                    complete_dict[dict_key]={'percent': 0}
                complete_dict[dict_key]['percent']+=frac_reads/classified_read_frac
            elif line.startswith('un'):
                read_frac=float(line.split()[2])
                classified_read_frac-=read_frac
    
    tax_dict={}       
    for rank in off_ranks:
        tax_dict[rank]={}
        end=off_ranks.index(rank)+1
        for taxon in complete_dict:
            lineage=taxon.rstrip('|').split('|')
            if len(lineage)>=len(off_ranks[0:end]) and lineage[end-1]!='':
                rank_lineage='|'.join(lineage[0:end])
                taxid=lineage[end-1]
                if not taxid in tax_dict[rank]:
                    tax_dict[rank][taxid]={'percent': 0,
                                           'lineage': rank_lineage,
                                           'taxon': ''}
                tax_dict[rank][taxid]['percent']+=complete_dict[taxon]['percent']
                    
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
            outf2.write(json.dumps(tax_dict, indent=4))
    



if __name__=='__main__':
    # for i in [6,13,23,25,26,30,33,34,38,53]:
    #     summarize_RAT_for_CAMI('/net/phage/linuxhome/mgx/people/tina/RAT/'
    #                             'benchmark/CAMI_II_mousegut_w_dm/wbin/'
    #                             'smp{}_clean.complete.abundance.txt'
    #                             ''.format(i), 
    #                             '/net/phage/linuxhome/mgx/people/tina/RAT/'
    #                             'benchmark/CAMI_II_mousegut_w_dm/profiles_all/'
    #                             'smp{}.wbin.0.001.profile'.format(i), 0.001)
    # for i in [6,13,23,25,26,30,33,34,38,53]:
    #     summarize_RAT_for_CAMI('/net/phage/linuxhome/mgx/people/tina/RAT/'
    #                             'benchmark/CAMI_II_mousegut/abundance_files/'
    #                             'smp{}.RAT.complete.abundance.txt'
    #                             ''.format(i), 
    #                             '/net/phage/linuxhome/mgx/people/tina/RAT/'
    #                             'benchmark/CAMI_II_mousegut_w_dm/profiles_all/'
    #                             'smp{}.robust.0.001.profile'.format(i), 0.001)
    
    
    # for i in ['W19-1','W19-2','W19-3','W19-4','W19-5','W19-6',
    #           'W22-1','W22-2','W22-3','W22-4','W22-5','W22-6',
    #           'W23-1','W23-2','W23-3','W23-4','W23-5','W23-6']:
        # summarize_RAT_for_CAMI('/net/phage/linuxhome/mgx/people/tina/wageningen/'
        #                        'RAT_robust/{}.complete.abundance.txt'
        #                        ''.format(i[1:]),
        #                        '/net/phage/linuxhome/mgx/people/tina/RAT/'
        #                         'benchmark/CAMI_II_mousegut_w_dm/profiles_wageningen/'
        #                         '{}.robust.profile'.format(i))

    
    CAT_folder='/net/mgx/linuxhome/mgx/people/bastiaan/phage-files/CAT_prepare/CAT_prepare_20190108/'
    env=sys.argv[1]
    smp=sys.argv[2]
    mode=sys.argv[3]
    cut_off=float(sys.argv[4])
    
    
    if env=='marine':
        prefix='2018.08.15_09.49.32'
    elif env=='plant':
        prefix='2019.09.27_13.59.10'
        
    path_to_tina='/net/phage/linuxhome/mgx/people/tina/'
    abundance=(path_to_tina + 'CAMI_II/{}/simulation_short_read/{}_sample_{}/'
               'RAT_NR/{}{}.{}.complete.abundance.txt'
               ''.format(env,prefix,smp,env,smp,mode))
    outfile=(path_to_tina + 'RAT/revision/{}/profiles/'
             '{}{}.{}.{}.profile'.format(env,env,smp,mode,cut_off))
    
    summarize_RAT_for_CAMI(abundance,outfile,cut_off)
    
    
    
    