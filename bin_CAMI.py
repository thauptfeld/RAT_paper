# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 16:25:21 2020

@author: acer
"""

def bin_CAMI_contigs(assembly_fasta, mapping_file, output_folder):
    bin2contig={}
    
# def make_bin2contig(mapping_file):
    bin2contig={}
    # Make bin2contig dictionary
    print('processing mapping file...')
    with open(mapping_file, 'r') as m:
        for line in m:
            if not line.startswith('#') and not line.strip()=='':
                l=line.split()
                if l[1] not in bin2contig:
                    bin2contig[l[1]]=set()
                bin2contig[l[1]].add(l[0])
    # return bin2contig

# def make_contig2seq(assembly_fasta):
    # make contig2sequence dictionary
    contig2seq={}
    first=True
    print('processing assembly...')
    with open(assembly_fasta, 'r') as a:
        contig=''
        for line in a:
            if not line.startswith('>'):
                contig+=line
            else:
                if first:
                    contig+=line
                    contig_id=line.rstrip()[1:]
                    first=False
                else:
                    contig2seq[contig_id]=contig
                    contig=line
                    contig_id=line.rstrip()[1:]
        contig2seq[contig_id]=contig

# def write_bins(bin2contig, contig2seq):
    # write bin files
    print('writing bin files...')
    for b in bin2contig:
        with open('{}/{}.fna'.format(output_folder, b), 'w') as outf:
            for c in bin2contig[b]:
                outf.write(contig2seq[c])
    
    return

if __name__=='__main__':
    bin_CAMI_contigs('/net/phage/linuxhome/tina/CAMI/CAMI_II/mousegut/19122017_mousegut_scaffolds/anonymous_gsa_pooled.fasta', 
                      '/net/phage/linuxhome/tina/CAMI/CAMI_II/mousegut/19122017_mousegut_scaffolds/gsa_pooled_mapping.tsv',
                      '/net/phage/linuxhome/tina/CAMI/CAMI_II/mousegut/19122017_mousegut_scaffolds/bins')