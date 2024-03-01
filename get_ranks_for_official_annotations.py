# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 16:26:28 2023

@author: Ernestina Hauptfeld
"""
import datetime
import sys

official_ranks=['superkingdom','phylum','class','order','family','genus','species']

def import_CAMI_results(CAMI_read_mapping, taxid2rank, outf):
    # ranks={}
    reads_order=[]
    reads_fam=[]
    reads_gen=[]
    with open(CAMI_read_mapping, 'r') as f1:
        for line in f1:
            # startswith # for CAMI II, startwith @ for CAMI I
            if not line.startswith('#') and not line.rstrip()=='':
                line=line.split('\t')
                taxid = line[2]
                
                rank=taxid2rank[taxid]
                # if rank not in ranks:
                #     ranks[rank]=0
                # ranks[rank]+=1
                
                if rank=='order':
                    reads_order.append(line[0])
                elif rank=='family':
                    reads_fam.append(line[0])
                elif rank=='genus':
                    reads_gen.append(line[0])
        
        print('{} Writing reads...'.format(timestamp()))
        with open(outf.format('order'), 'w') as f2:
            for r in reads_order:
                f2.write('{}\n'.format(r))
        with open(outf.format('family'), 'w') as f2:
            for r in reads_fam:
                f2.write('{}\n'.format(r))                    
        with open(outf.format('genus'), 'w') as f2:
            for r in reads_gen:
                f2.write('{}\n'.format(r))
    # print(ranks)


def get_reads_from_table(infile, read_list, taxid_column, readid_column, taxid2rank, outfile):
    
    with open(infile) as inf, open(outfile, 'w') as outf:
        for line in inf:
            if not line.startswith('readID'):
                line=line.split('\t')
                
                readid=line[readid_column]
                if readid in read_list:
                    taxid=line[taxid_column].strip()
                    if 'taxid' in taxid:
                        taxid=taxid.split('taxid ')[1][:-1]
                    
                    if taxid=='0':
                        rank='unclassified'
                    else:
                        rank=taxid2rank[taxid]
                    
                    outf.write('{}\t{}\t{}\n'.format(readid, taxid, rank))
                
    
def get_official_annotations_for_subset(infile, read_list, outfile):
    
    with open(infile) as inf, open(outfile, 'w') as outf:
        for line in inf:
            readid=line.split()[0]
            if readid in read_list:
                    
                outf.write(line)    
    


def import_nodes(nodes_dmp):
    
    taxid2parent = {}
    taxid2rank = {}

    with open(nodes_dmp, 'r') as f1:
        for line in f1:
            line = line.split('\t')

            taxid = line[0]
            parent = line[2]
            rank = line[4]

            taxid2parent[taxid] = parent
            taxid2rank[taxid] = rank

    return taxid2parent, taxid2rank




def timestamp():
    now = datetime.datetime.now()
    str_ = '[{0}]'.format(now.strftime('%Y-%m-%d %H:%M:%S'))

    return str_


if __name__=='__main__':
    env=sys.argv[1]
    smp=sys.argv[2]

    
    CAT_folder='/net/mgx/linuxhome/mgx/people/bastiaan/phage-files/CAT_prepare/CAT_prepare_20190108/'
    # print('{} Importing nodes...'.format(timestamp()))
    # taxid2parent, taxid2rank = import_nodes(CAT_folder+'2019-01-08_taxonomy/nodes.dmp')
    # print('{} Done! Making dictionary...'.format(timestamp()))
    
    
    results='/net/phage/linuxhome/dutilh-group/tina/RAT/revision/read_classifiers_new_db_marine/results/marine1/'
    
    
    
    
    if env=='marine':
        prefix='2018.08.15_09.49.32'
    elif env=='plant':
        prefix='2019.09.27_13.59.10'
    
    path_to_cami='/net/phage/linuxhome/dutilh-group/tina/CAMI_II/'
    path_to_rev='/net/phage/linuxhome/dutilh-group/tina/RAT/revision/'
    marine_folder=path_to_cami+'{}/simulation_short_read/{}_sample_{}/reads/'.format(env,prefix,smp)
    

    outf=path_to_rev+'{}/reads_no_species/{}{}'.format(env,env,smp)+'_{}.txt'
        
    # print('{} Getting reads for {}{} that are supposed to be unknown at species rank'.format(timestamp(),env,smp))
    # import_CAMI_results(marine_folder + 'reads_mapping.tsv',taxid2rank,outf)
    # print('{} Done!'.format(timestamp()))
    
    
    
    print('{} Getting reads for {}{} that are supposed to be unknown at species rank'.format(timestamp(),env,smp))
    official=marine_folder+'reads_mapping.tsv'
    family=path_to_rev+'{}/reads_no_species/{}{}_family.txt'.format(env,env,smp)
    genus=path_to_rev+'{}/reads_no_species/{}{}_genus.txt'.format(env,env,smp)
    unknowns=set(open(family).read().strip().split()+open(genus).read().strip().split())
    outfile=marine_folder+'unknowns_reads_mapping.tsv'
    get_official_annotations_for_subset(official, unknowns, outfile)
    print('{} Done!'.format(timestamp()))
    
    
    # reads_=open('/home/tina/tmp/reads_supposed_genus_marine1.txt').read().strip().split()
    # reads=set([r.split('/')[0] for r in reads_])
    
    
    
    
    # kaiju=results+'kaiju/marine1.kaiju.out'
    # kaiju_out='/home/tina/tmp/TMP_kaiju_reads_supposed_genus_marine1.txt'
    
    # get_reads_from_table(kaiju, reads, 2, 1, taxid2rank, kaiju_out)
    # print('kaiju done...')
    
    # kraken=results+'kraken2/marine1.kraken2.out'
    # kraken_out='/home/tina/tmp/TMP_kraken_reads_supposed_genus_marine1.txt'
    # get_reads_from_table(kraken, reads, 2, 1, taxid2rank, kraken_out)  
    # print('kraken done...')
    
    # centrifuge=results + 'centrifuge/marine1.centrifuge.out'
    # centrifuge_out='/home/tina/tmp/TMP_centrifuge_reads_supposed_genus_marine1.txt'
    # get_reads_from_table(centrifuge, reads, 2, 0, taxid2rank, centrifuge_out)
    # print('centrifuge done...')
    