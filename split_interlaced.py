# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 14:23:22 2020

@author: acer
"""

def split_interlaced_reads(read_file):
    with open(read_file, 'r') as reads, open(read_file.split('.fq')[0]+'_R1.fq', 'w') as outf1, open(read_file.split('.fq')[0]+'_R2.fq', 'w') as outf2:
        read=''
        i=1
        r=1
        for line in reads:
            i+=1
            read+=line
            if i%4==1:
                if r%2==1:
                    outf1.write(read)
                else:
                    outf2.write(read)
                read=''
                r+=1

    return

if __name__=='__main__':
    split_interlaced_reads('/net/phage/linuxhome/tina/CAMI/CAMI_I/toy_low/S_S001__genomes_30__insert_180_reads_anonymous.fq')