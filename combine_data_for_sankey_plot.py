# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 16:28:44 2022

@author: acer
"""

import json
import sys



if __name__=='__main__':
    
    # sample=sys.argv[1]
    sample=13
    print("loading read annotations for sample {}".format(sample))
    read_annotations=json.load(open('/net/phage/linuxhome/mgx/people/tina/tmp/diamond_dicts/'
                                    'd2c_lineages_{}_ncbi.txt'.format(sample)))
    
    print("loading contig and bin annotations")
    r2c_file=open('/net/phage/linuxhome/mgx/people/tina/RAT/benchmark/CAMI_II_mousegut_w_dm/'
             'wmetabat/smp{}_metabat_clean.read2classification.txt'.format(sample)).read().strip().split('\n')[2:]
    
    total_reads=len(r2c_file)
    
    r2c={'contig': {},
         'bin': {}}
    reads=set()
    
    print('making dict')
    for r in r2c_file:
        r=r.split('\t')
        reads.add(r[0])
        
        if len(r)>2:
            try:
                r_id=r[0]
                r_contig=r[3]
                r_bin=r[2]
                
                if r_contig!='':
                    if r_contig.split(';')[2]=='2':
                        r2c['contig'][r_id]=r_contig.split(';')[3]
                        if r_contig.split(';')[3]=='1783272':
                            r2c['contig'][r_id]=r_contig.split(';')[4]
                if r_bin!='':
                    try:
                        if r_bin.split(';')[2]=='2':
                            r2c['bin'][r_id]=r_bin.split(';')[3]
                            if r_bin.split(';')[3]=='1783272':
                                r2c['bin'][r_id]=r_bin.split(';')[4]
                    except IndexError:
                        if len(r_bin.split(';'))<=2:
                            r2c['bin'][r_id]='binned_unclassified'
                            
                
            except IndexError:
                continue

    del r2c_file
    
    sankey={}
    
    print('Making Sankey data...')
    
    for read in reads:
        r=''
        c=''
        b=''
        if read+'/1' in read_annotations:
            r=read_annotations[read+'/1'][1]
            if r=='1783272':
                r=read_annotations[read+'/1'][2]
        if read in r2c['contig']:
            c=r2c['contig'][read]
        if read in r2c['bin']:
            b=r2c['bin'][read]
            
        sank_id=';'.join([r,c,b])
        if sank_id not in sankey:
            sankey[sank_id]=0
        sankey[sank_id]+=1
        
    
    with open('/home/tina/Documents/RAT/sankey_new_data_w_binned.json', 'w') as outf:
        outf.write(json.dumps(sankey))
        