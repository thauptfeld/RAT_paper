# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 18:42:18 2023

@author: acer
"""

import os

def parse_and_write_fasta(input_file, output_file):
    # Open the input FASTA file
    with open(input_file, 'r') as f:
        # Initialize an empty dictionary to store the sequences
        sequences = {}
        # Initialize a variable to store the current sequence name
        current_name = None
        # Iterate over each line in the file
        for line in f:
            # If the line starts with a ">" character, it's a header line
            if line.startswith(">"):
                # Extract the name of the sequence from the header
                current_name = line[1:].strip()
                # Initialize an empty string to store the sequence
                sequences[current_name] = ""
            else:
                # If the line doesn't start with a ">", it's part of the sequence
                sequences[current_name] += line.strip()

    # Open the output FASTA file
    with open(output_file, 'w') as f:
        # Iterate over each sequence in the dictionary
        for name in sequences:
            # Write the header line to the output file
            f.write(">" + name + "\n")
            # Write the sequence to the output file, 60 characters per line
            for i in len(sequences[name]):
                if i%60==59:
                    f.write()
                    
                    
def parse_and_write_fasta1(input_file, output_file):
    i=1
    # Open the input FASTA file
    with open(input_file, 'r') as f, open(output_file, 'w') as g:
        # Initialize a variable to store the current sequence name
        current_name = None
        # Iterate over each line in the file
        for line in f:
            # If the line starts with a ">" character, it's a header line
            if line.startswith(">"):
                # Extract the name of the sequence from the header
                current_name = line[1:].strip()
                g.write('>{}_{}\n'.format(i,current_name))
                i+=1
                
                
            else:
                sequence=line.strip()
                times=int(len(sequence)/60)+1
                for t in range(times):
                    g.write(sequence[0:60]+'\n')
                    sequence=sequence[60:]
                
                
if __name__=='__main__':
    parse_and_write_fasta1('/net/phage/linuxhome/mgx/people/tina/RAT/kaiju_nr/2019-01-08.nr.kaiju', 
                          '/net/phage/linuxhome/mgx/people/tina/RAT/kaiju_nr/2019-01-08.nr.kaiju_numbered')