# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 14:45:04 2023

@author: acer
"""

import gzip
from Bio import SeqIO


path_to_prophages='/net/phage/linuxhome/mgx/people/tina/prophages_Rob/'
predictions=open(path_to_prophages+'vibrio_cholerae_predictions.txt').read().strip().split()

l=len(predictions)

with open(path_to_prophages+'vibrio_prophages.fasta', 'w') as outf:
    i=0
    for p in predictions:
        i+=1
        gbk_id=p[4:-2]
        gca_path='phispy/GCA/{0}/{1}/{2}/'.format(gbk_id[:3],
                                                 gbk_id[3:6],
                                                 gbk_id[6:]
                                                 )
        gbk=gzip.open(path_to_prophages+gca_path+p+'_VOGS_phage.gbk.gz', 'rt')
        gbk_records=SeqIO.parse(gbk, 'genbank')
        
        for record in gbk_records:
            rec_id=record.id
            rec_seq=str(record.seq)
            
            outf.write('>{0}.{1}\n{2}\n'.format(p,
                                                rec_id,
                                                rec_seq))
            
        if i%1000==0:
            print('{0}/{1} records written.'.format(i,l))

