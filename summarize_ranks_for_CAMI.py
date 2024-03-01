#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 15:37:29 2020

@author: tina
"""
import json

# lineage       number of reads fraction of reads       taxon length    corrected fraction      lineage ranks

off_ranks=['superkingdom','phylum','class','order','family','genus','species']

def summarize_RAT_for_CAMI(rat_table, level):
    tax_dict={}
    with open(rat_table) as rat:
        for line in rat:
            pieces=line.split('\t')
            if not line.startswith('#') and not line.startswith('un'):
                lineage, num_reads, frac_reads, tax_length, ranks = pieces[0].split(';'), int(pieces[1]), float(pieces[2]), int(pieces[3]), pieces[5].strip().split(';')
#                print(ranks)
                new_lineage=level*['']
                for i in range(0,level):
                    if off_ranks[i] in ranks:
                        j=ranks.index(off_ranks[i])
                        new_lineage[i]=lineage[j].rstrip('*')
                new_lineage='|'.join(new_lineage).rstrip('|')
                
                
                if len(new_lineage.split('|'))==level:
                    taxid=new_lineage.split('|')[-1]
                    if not taxid in tax_dict:
                        tax_dict[taxid]=[num_reads, 
                                frac_reads, 
                                tax_length, 
                                new_lineage]  
                    else:
                        tax_dict[taxid][0]+=num_reads
                        tax_dict[taxid][1]+=frac_reads
                        tax_dict[taxid][2]+=tax_length
                        
                else:
                    tax_dict['unclassified'][0]+=num_reads
                    tax_dict['unclassified'][1]+=frac_reads
                    tax_dict['unclassified'][2]+=tax_length
                
            elif line.startswith('un'):
                tax_dict[pieces[0]]=[int(pieces[1]), float(pieces[2])]
                if line.startswith('unclassified'):
                    tax_dict['unclassified']+=[int(pieces[3])]
                    
    return tax_dict
    
    
def write_CAMI_profile(profile_dict, output_file):
    with open(output_file, 'w+') as outf1, open(output_file.rsplit('.', 1)[0]+ '.json', 'w+') as outf2:
        outf1.write('@SampleID:\n'
                   '@Version:0.9.1\n'
                   '@Ranks:superkingdom|phylum|class|order|family|genus|species|strain\n\n'
                   '@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n')
        for rank in profile_dict:
            for taxon in profile_dict[rank]:
                if not taxon=='unmapped' and not taxon=='unclassified':
                    outf1.write('{}\t{}\t{}\t\t{}\n'.format(taxon.rstrip('*'),
                                                     rank, profile_dict[rank][taxon][3],
                                                     profile_dict[rank][taxon][1]))
        outf2.write(json.dumps(profile_dict, indent=4))
            



if __name__=='__main__':
    # for n in [1,2,3,4,5]:
    #     output_dict={}
    #     prefix='/net/phage/linuxhome/tina/RAT/benchmark/CAMI_high/RH_S00{}'.format(n)
    #     for i in range(1,8):
    #         output_dict[off_ranks[i-1]]=summarize_RAT_for_CAMI('{}.{}.abundance.txt'.format(prefix, off_ranks[i-1]), i)
    #     write_CAMI_profile(output_dict, prefix + '.profile')
    #     with open(prefix+'.json') as outf:
    #         json.dumps(output_dict, indent=4)
    
    for n in range(0,64):
        print(n)
        output_dict={}
        prefix='/net/phage/linuxhome/tina/RAT/benchmark/CAMI_II_mousegut/smp{}.RAT'.format(n)
        for i in range(1,8):
            output_dict[off_ranks[i-1]]=summarize_RAT_for_CAMI('{}.{}.abundance.txt'.format(prefix, off_ranks[i-1]), i)
        write_CAMI_profile(output_dict, prefix + '.profile')
        with open(prefix+'.json') as outf:
            json.dumps(output_dict, indent=4)