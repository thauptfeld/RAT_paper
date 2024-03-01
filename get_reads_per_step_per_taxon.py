# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 21:19:23 2021

@author: acer
"""

import json
import copy

samples=['W19-1','W19-2','W19-3','W19-4','W19-5','W19-6',
          'W22-1','W22-2','W22-3','W22-4','W22-5','W22-6',
          'W23-1','W23-2','W23-3','W23-4','W23-5','W23-6']

samples=['W23-2']
standard={'total_reads':0,
        'bin_reads':0,
        'contig_reads':0,
        'contig_dm':0,
        'read_dm':0,
        }
output={}

for sample in samples:
    print(sample)
    output[sample]={'total_reads': 40191520}
    path='/net/phage/linuxhome/mgx/people/tina/wageningen/2018/{}/'.format(sample)
    with open(path+'{}.RAT.read2classification.txt'.format(sample)) as inf:
        
        for line in inf:
            if not line.startswith('#'):
                line=line.rstrip().split('\t')
                if len(line)>2:
                    if line[2]!='':
                        
                        # Assigned by bin
                        taxon=line[2].replace('*','').rstrip()
                        if taxon not in output[sample]:
                            output[sample][taxon]=copy.deepcopy(standard)
                        output[sample][taxon]['bin_reads']+=1
                        output[sample][taxon]['total_reads']+=1
                        
                        #?
                        if len(line)>3 and line[3]=='':
                            output[sample]['weird']+=1
                            
                    # Assigned by contig
                    elif line[2]=='' and line[3]!='':
                        taxon=line[3].replace('*','').rstrip()
                        if taxon not in output[sample]:
                            output[sample][taxon]=copy.deepcopy(standard)
                        output[sample][taxon]['contig_reads']+=1
                        output[sample][taxon]['total_reads']+=1
                        
                    # Assigned by contig mapping
                    elif line[2]=='' and line[3]=='' and line[4]!='':
                        taxon=line[4].replace('*','').rstrip()
                        if taxon not in output[sample]:
                            output[sample][taxon]=copy.deepcopy(standard)
                        output[sample][taxon]['contig_dm']+=1
                        output[sample][taxon]['total_reads']+=1
                        
                    # Assigned by read mapping
                    elif line[2]=='' and line[3]=='' and line[4]=='' and line[5]!='':
                        taxon=line[5].replace('*','').rstrip().replace('fw: ', '').replace('rev: ', '')
                        if taxon not in output[sample]:
                            output[sample][taxon]=copy.deepcopy(standard)
                        output[sample][taxon]['read_dm']+=1
                        output[sample][taxon]['total_reads']+=1
                
                    
with open('/home/tina/Documents/RAT/W23-2_reads_per_step_per_taxon.json', 'w') as outf:
    outf.write(json.dumps(output, indent=4))
    
    



samples=json.load(open('/home/tina/Documents/RAT/W23-2_reads_per_step_per_taxon.json'))
with open('../benchmark/W23-2_reads_per_step_per_taxon_min0.001.tsv', 'w') as outf:
    outf.write('# lineage\tsample\tbin\tcontig\tcontig_dm\tread_dm\ttotal_reads\n')
    for s in samples:
        for t in samples[s]:
            if t not in ['total_reads', 'unclassified'] and samples[s][t]['total_reads']>=402:
                outf.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(t,s,
                                                  samples[s][t]['bin_reads'],
                                                  samples[s][t]['contig_reads'],
                                                  samples[s][t]['contig_dm'],
                                                  samples[s][t]['read_dm'],
                                                  samples[s][t]['total_reads']))