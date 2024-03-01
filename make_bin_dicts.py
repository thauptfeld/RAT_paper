# # -*- coding: utf-8 -*-
# """
# Created on Wed Aug 10 10:09:44 2022

# @author: acer
# """

import os
import sys
import json

def process_bin_folder(path_to_folder, bin_suffix):
    """
    takes the path to a folder containing metagenomic bins and picks out all
    nucleotide fasta files. Returns a dictionary that contains the bin ids as
    keys and a set of all contig ids that belong in that bin as values. 
    Important: Only contains the binned sequences. Every contig that is not
    in the output dictionary is unbinned.
    @Tina: Don't we have a function from the BAT workflow that does this?
    """
    bin_dict={}
    bins=os.listdir(path_to_folder)
    
    for b in bins:
        # grabs all fasta files and gets the bin id
        if b.split('.')[-1]==bin_suffix.split('.')[-1] and not b.rsplit('.', 1)[0].endswith('unbinned'):
            bin_id=b.strip()
            bin_dict[bin_id]=list()
            
            # open fasta file and split it into contigs
            contigs=open(path_to_folder+b).read().split('>')[1:]
            
            # add all contig ids to set
            for contig in contigs:
                bin_dict[bin_id].append(contig.split()[0])
                
    return bin_dict

if __name__=='__main__':
    bin_folder = sys.argv[1]
    bin_suffix = sys.argv[2]
    outfile = sys.argv[3]
    
    condict=process_bin_folder(bin_folder, bin_suffix)
    json.dump(condict, open('/net/phage/linuxhome/mgx/people/tina/wageningen/2018/{}'.format(outfile), 'w'))