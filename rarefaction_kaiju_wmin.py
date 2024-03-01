# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 10:46:54 2022

@author: tina
"""
import random
import sys


def do_rarefaction(input_file, list_of_percentages, iterations, path_to_out,profile):
    """
    Parameters
    ----------
    input_file : string
        Path to the kaiju.names.out file.
    list_of_percentages : list
        list of percent values that you want to do rarefactions for.
    iterations: integer
        number of iterations for rarefaction for each percent value
    
    Returns
    -------
    None.
    Writes output for each iteration in the format: sample.kaiju.rf{percent}.min_abundance.iteration.csv

    """
    read_list=open(input_file).read().strip().split('\n')
    abundant_taxa=get_abundant_taxa(profile)
    smp=input_file.split('/')[-1].split('.')[0]
    for p in list_of_percentages:
        print('{}%...'.format(p))
        sampled_dict=subsample(read_list, p, iterations)
        for i in sampled_dict:
            taxa=count_up_taxa(sampled_dict[i],abundant_taxa)
            # out1='{0}.kaiju.rf{1}.it{2}.csv'.format(smp, p, i)
            out2='{0}.kaiju.rf{1}.min{2}_new.it{3}.csv'.format(smp, p, '0.001', i)
            
            # open('{0}{1}/{2}'.format(path_to_out,p,out1), 'w') as outf1,
            with open('{0}{1}/{2}'.format('./',p,out2), 'w') as outf2:
                for taxon in taxa:
                    # outf1.write('{},{},{}\n'.format(taxon,
                    #                         taxa[taxon]['num_reads'],
                    #                         taxa[taxon]['frac_reads']))
                    # if taxa[taxon]['frac_reads']>0.00001:
                    outf2.write('{},{},{}\n'.format(taxon,
                                        taxa[taxon]['num_reads'],
                                        taxa[taxon]['frac_reads']))
            
    
    
    return


def subsample(read_list, percent, i):
    """
    Takes a read dictionary as input and subsamples p percent of reads it i times.
    Returns a dictionary with iterations as keys that are subsets of read_dict.
    """
    sample_dict={}
    for iteration in range(i):
        it=str(iteration+1).rjust(3,'0')
        number=int((percent/100)*len(read_list))
        
        sample_dict[it]=random.sample(read_list, number)
        
    return sample_dict




def count_up_taxa(read_list, abundant_taxa):
    weirdos=0
    total=len(read_list)
    tax_dict={'u': {'num_reads':0,
                    'frac_reads':0}}
    for r in read_list:
        if r.startswith('U'):
            tax_dict['u']['num_reads']+=1
            tax_dict['u']['frac_reads']+=1/total
        else:
            try:
                taxon=r.split('\t')[7]
                if taxon in abundant_taxa:
                    if taxon not in tax_dict:
                        tax_dict[taxon]={'num_reads':0,
                                         'frac_reads':0}
                    tax_dict[taxon]['num_reads']+=1
                    tax_dict[taxon]['frac_reads']+=1/total
                else:
                    tax_dict['u']['num_reads']+=1
                    tax_dict['u']['frac_reads']+=1/total
            except IndexError:
                weirdos+=1
    
    
    return tax_dict


def get_abundant_taxa(profile):
    abundant_taxa=set()
    with open(profile) as p:
        for line in p:
            taxon=line.split(',')[0]
            abundant_taxa.add(taxon)
    return abundant_taxa
            


if __name__=='__main__':
    path_to_out=('/hosts/linuxhome/phage/mgx/people/tina/RAT/benchmark/CAMI_II_mousegut_w_dm/'
                 'rarefaction/')
    
    # wells=['19', '22', '23']
    # mousegut=['smp6','smp13','smp23','smp25','smp26','smp30','smp33','smp38','smp53']
    
    i=sys.argv[1]
    
    # for i in mousegut:
        # for j in range(1,7):
            # sample='W{0}-{1}'.format(i,j)
    sample=i
    print('Working on sample {}...'.format(sample))

    input_file=('/hosts/linuxhome/phage/mgx/people/tina/RAT/benchmark/'
                'CAMI_II_mousegut_w_dm/kraken_kaiju_centrifuge/{0}/'
                'kaiju/{1}.kaiju.names.out'.format(sample,sample))
    profile=('/hosts/linuxhome/phage/mgx/people/tina/RAT/benchmark/'
             'CAMI_II_mousegut_w_dm/kraken_kaiju_centrifuge/'
             '{0}.kaiju.manually_counted.0.001.csv'.format(sample))
    print(input_file)
    list_of_percentages=[1,5,10,20,30,40,50,60,70,80,90]
    iterations=10
    do_rarefaction(input_file, list_of_percentages, iterations, path_to_out, profile)