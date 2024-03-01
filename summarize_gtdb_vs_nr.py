# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 10:09:32 2023

@author: Ernestina Hauptfeld
"""



### NR
# contig  classification  reason  lineage  lineage scores  superkingdom   phylum  class   order   family  genus   species
# bin     classification  reason  lineage  lineage scores  superkingdom   phylum  class   order   family  genus   species


def summarize_c2c_nr(nr_table,sample, outf):
    na=['no support', 'NA', '']
    bins=0
    assigned=0
    assigned_sk=0
    assigned_p=0
    assigned_c=0
    assigned_o=0
    assigned_f=0
    assigned_g=0
    assigned_sp=0
    with open(nr_table, 'r') as c2c:
        for c in c2c:
            if not c.startswith('#'):
                bins+=1
                
                if c.split('\t')[1] in ['taxid assigned', 'classified']:
                    line=c.strip().split('\t')
                    assigned+=1
                    if line[5] not in na:
                        assigned_sk+=1
                    if line[6] not in na:
                        assigned_p+=1
                    if line[7] not in na:
                        assigned_c+=1
                    if line[8] not in na:
                        assigned_o+=1
                    if line[9] not in na:
                        assigned_f+=1
                    if line[10] not in na:
                        assigned_g+=1
                    if line[11] not in na:
                        assigned_sp+=1
    
    outf.write(f'{sample},nr,{bins},{assigned},{assigned_sk},{assigned_p},'
               f'{assigned_c},{assigned_o},{assigned_f},{assigned_g},'
               f'{assigned_sp}\n')
    
    
    return



### GTDB
# contig   classification    reason    lineage    lineage scores
# bin      classification    reason    lineage    lineage scores
def summarize_c2c_gtdb(gtdb_table, sample, outf):
    bins=0
    assigned=0
    assigned_sk=0
    assigned_p=0
    assigned_c=0
    assigned_o=0
    assigned_f=0
    assigned_g=0
    assigned_sp=0
    
    with open(gtdb_table, 'r') as c2c:
        for c in c2c:
            if not c.startswith('#'):
                bins+=1
                
                if c.split('\t')[1] =='taxid assigned':
                    assigned+=1
                    lineage=c.split('\t')[3]
                    if 'd__' in lineage:
                        assigned_sk+=1
                    if 'p__' in lineage:
                        assigned_p+=1
                    if 'c__' in lineage:
                        assigned_c+=1
                    if 'o__' in lineage:
                        assigned_o+=1
                    if 'f__' in lineage:
                        assigned_f+=1
                    if 'g__' in lineage:
                        assigned_g+=1
                    if 's__' in lineage:
                        assigned_sp+=1
    
    outf.write(f'{sample},gtdb,{bins},{assigned},{assigned_sk},{assigned_p},'
               f'{assigned_c},{assigned_o},{assigned_f},{assigned_g},'
               f'{assigned_sp}\n')
    
    return


if __name__=='__main__':
    path_to_rat='/net/phage/linuxhome/dutilh-group/tina/RAT/'
    # path_to_cami='/net/phage/linuxhome/dutilh-group/tina/CAMI_II/'
    # plants=[2,3,5,8,10,12,14,15,18,19]
    
    # m_nr_bin='marine/simulation_short_read/2018.08.15_09.49.32_sample_{}/RAT_NR/marine{}.BAT.bin2classification.names.txt'
    # p_nr_bin='plant/simulation_short_read/2019.09.27_13.59.10_sample_{}/RAT_NR/plant{}.BAT.bin2classification.names.txt'
    # m_nr_contig='marine/simulation_short_read/2018.08.15_09.49.32_sample_{}/RAT_NR/marine{}.CAT.contig2classification.names.txt'
    # p_nr_contig='plant/simulation_short_read/2019.09.27_13.59.10_sample_{}/RAT_NR/plant{}.CAT.contig2classification.names.txt'
    
    # m_gtdb_bin='marine/simulation_short_read/2018.08.15_09.49.32_sample_{}/RAT_GTDB/marine{}.BAT.bin2classification.txt'
    # p_gtdb_bin='plant/simulation_short_read/2019.09.27_13.59.10_sample_{}/RAT_GTDB/plant{}.BAT.bin2classification.txt'
    # m_gtdb_contig='marine/simulation_short_read/2018.08.15_09.49.32_sample_{}/RAT_GTDB/marine{}.CAT.contig2classification.txt'
    # p_gtdb_contig='plant/simulation_short_read/2019.09.27_13.59.10_sample_{}/RAT_GTDB/plant{}.CAT.contig2classification.txt'
    
    # smp_nr_bin='benchmark/mousegut_robust/gsa_pooled.BAT.bin2classification.names.txt'
    # smp_gtdb_bin='revision/mousegut_GTDB/smp13.GTDB.BAT.bin2classifcation.txt'
    # smp_nr_contig='benchmark/mousegut_robust/gsa_pooled.CAT.contig2classification.names.txt'
    # smp_gtdb_contig='revision/mousegut_GTDB/smp13.GTDB.CAT.contig2classifcation.txt'
    
    # out_bins=path_to_rat+'revision/summary_gtdb_vs_nr_bins.csv'
    # out_contigs=path_to_rat+'revision/summary_gtdb_vs_nr_contigs.csv'

    
    
    # marine_c_nr=[path_to_cami+m_nr_contig.format(i,i) for i in range(0,10)]
    # plant_c_nr=[path_to_cami+p_nr_contig.format(i,i) for i in plants]
    
    # marine_b_nr=[path_to_cami+m_nr_bin.format(i,i) for i in range(0,10)]
    # plant_b_nr=[path_to_cami+p_nr_bin.format(i,i) for i in plants]
    
    # marine_c_gtdb=[path_to_cami+m_gtdb_contig.format(i,i) for i in range(0,10)]
    # plant_c_gtdb=[path_to_cami+p_gtdb_contig.format(i,i) for i in plants]
    
    # marine_b_gtdb=[path_to_cami+m_gtdb_bin.format(i,i) for i in range(0,10)]
    # plant_b_gtdb=[path_to_cami+p_gtdb_bin.format(i,i) for i in plants]
    
    # mouse_c_nr=[path_to_rat+smp_nr_contig]
    # mouse_c_gtdb=[path_to_rat+smp_gtdb_contig]
    
    # mouse_b_nr=[path_to_rat+smp_nr_bin]
    # mouse_b_gtdb=[path_to_rat+smp_gtdb_bin]
    
    
    # all_contig_gtdb=marine_c_gtdb+plant_c_gtdb+mouse_c_gtdb
    # all_contig_nr=marine_c_nr+plant_c_nr+mouse_c_nr
    
    # all_bin_gtdb=marine_b_gtdb+plant_b_gtdb+mouse_b_gtdb
    # all_bin_nr=marine_b_nr+plant_b_nr+mouse_b_nr
    
    # with open(out_bins, 'w') as outf:
    #     for f in all_bin_nr:
    #         sample=f.split('/')[-1].split('.')[0]
    #         summarize_c2c_nr(f, sample, outf)
    #     for f in all_bin_gtdb:
    #         sample=f.split('/')[-1].split('.')[0]
    #         summarize_c2c_gtdb(f, sample, outf)
            
    # with open(out_contigs, 'w') as outf:
    #     for f in all_contig_nr:
    #         sample=f.split('/')[-1].split('.')[0]
    #         summarize_c2c_nr(f, sample, outf)
    #     for f in all_contig_gtdb:
    #         sample=f.split('/')[-1].split('.')[0]
    #         summarize_c2c_gtdb(f, sample, outf)
            
    
    ### GTDB vs NR for biological example
    path_to_w_gtdb='/net/phage/linuxhome/mgx/people/tina/RAT/groundwater_RAT_sensitive/2018/GTDB/'
    path_to_w_nr='/net/phage/linuxhome/mgx/people/tina/RAT/groundwater_RAT_robust/wells/well{}/RAT/'
    
    w_gtdb_bin=path_to_w_gtdb+'W{}-{}/W{}-{}.GTDB.BAT.bin2classification.txt'
    w_gtdb_contig=path_to_w_gtdb+'W{}-{}/W{}-{}.GTDB.CAT.contig2classification.txt'
    w_nr_bin=path_to_w_nr+'{}-{}.BAT.bin2classification.names.txt'
    w_nr_contig=path_to_w_nr+'{}-{}.CAT.contig2classification.names.txt'
    
    out_bins=path_to_rat+'revision/groundwater_summary_gtdb_vs_nr_bins.csv'
    out_contigs=path_to_rat+'revision/groundwater_summary_gtdb_vs_nr_contigs.csv'
    
    gtdb_bin=[w_gtdb_bin.format(well,depth,well,depth) for well in [19,22,23] for depth in range(1,7)]
    nr_bin=[w_nr_bin.format(well,well,depth) for well in [19,22,23] for depth in range(1,7)]
    
    gtdb_contig=[w_gtdb_contig.format(well,depth,well,depth) for well in [19,22,23] for depth in range(1,7)]
    nr_contig=[w_nr_contig.format(well,well,depth) for well in [19,22,23] for depth in range(1,7)]
    
    with open(out_bins, 'w') as outf:
        for f in nr_bin:
            # print(f)
            sample='W'+f.split('/')[-1].split('.')[0]
            summarize_c2c_nr(f, sample, outf)
        for f in gtdb_bin:
            # print(f)
            sample=f.split('/')[-1].split('.')[0]
            summarize_c2c_gtdb(f, sample, outf)
            
    with open(out_contigs, 'w') as outf:
        for f in nr_contig:
            sample='W'+f.split('/')[-1].split('.')[0]
            summarize_c2c_nr(f, sample, outf)
        for f in gtdb_contig:
            sample=f.split('/')[-1].split('.')[0]
            summarize_c2c_gtdb(f, sample, outf)