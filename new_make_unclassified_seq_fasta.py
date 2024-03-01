# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 12:40:11 2024

@author: acer
"""

import sys
import shared

def make_unclassified_seq_fasta(seq_fasta, unclassified_seq_ids, 
                                unclassified_seq_fasta, f_format, f_mode, 
                                suffix=''):
    fasta_dict={}
    if f_format == 'fasta':
        with shared.optionally_compressed_handle(seq_fasta, 'r') as f1: 
            for line in f1:
                if line.startswith('>'):
                    header = line.rstrip().split(' ')[0].lstrip('>')
                    if header in fasta_dict:
                        sys.exit('Duplicate fasta headers in the file!') # This should be prettified.
                    fasta_dict.setdefault(header, '')
                else:
                    fasta_dict[header] += line.rstrip()
                    
                    
    elif f_format=='fastq':
        with shared.optionally_compressed_handle(seq_fasta, 'r') as f1: 
            for n, line in enumerate(f1):
                if n % 4 == 0:
                    if not line.startswith('@'):
                        sys.exit('Unknown seuqence identifier symbol in line {n}!') # This should be prettified.
                    header = line.rstrip().split(' ')[0].lstrip('@') 
                    
                    if header in fasta_dict:
                        sys.exit('Duplicate fasta headers in the file!')
    else:
        sys.exit('Unknown read file format!')
    
    
    with open(unclassified_seq_fasta, f_mode) as outf: 
        # print(unclassified_seq_ids) 
        suffices=['/1', '/2', '_1', '_2']
        for seq in unclassified_seq_ids: 
            if seq in fasta_dict:
                seq=seq
            else:
                for s in suffices:
                    if seq+s in fasta_dict:
                        seq+=s
    # print(seq)
    outf.write('>{0}{1}\n{2}\n'.format(fasta_dict[seq].id, suffix,
                                       fasta_dict[seq].seq))
        
