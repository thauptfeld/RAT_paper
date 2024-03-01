# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 10:46:54 2022

@author: tina
"""
import random
import sys

def do_rarefaction(input_file, list_of_percentages, iterations, path_to_out):
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
    read_list=open(input_file).read().strip().split('\n')[1:]
    smp=input_file.split('/')[-1].split('.')[0]
    for p in list_of_percentages:
        print('{}%...'.format(p))
        sampled_dict=subsample(read_list, p, iterations)
        for i in sampled_dict:
            taxa=count_up_taxa(sampled_dict[i])
            out1='{0}.cami.rf{1}.it{2}.csv'.format(smp, p, i)
            
            with open('{0}{1}/{2}'.format('./',p,out1), 'w') as outf1:
                outf1.write('{0}'.format(len(taxa)))

    
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




def count_up_taxa(read_list):
    tax_dict=set()
    for r in read_list:
        tax_dict.add(r.split()[1])
        
    return tax_dict



if __name__=='__main__':
    path_to_out=('/hosts/linuxhome/phage/mgx/people/tina/RAT/benchmark/CAMI_II_mousegut_w_dm/'
                 'rarefaction/')
    # wells=['19', '22', '23']
    mousegut=['smp6','smp13','smp23','smp25','smp26','smp30','smp33','smp38','smp53']
    
    # for i in mousegut:
        # for j in range(1,7):
            # sample='W{0}-{1}'.format(i,j)
    sample=sys.argv[1]
    print('Working on sample {}...'.format(sample))

    input_file=('/hosts/linuxhome/phage/mgx/people/tina/CAMI/CAMI_II/mousegut/'
                '19122017_mousegut_scaffolds/read_mappings/'
                '{0}_reads_mapping.tsv'.format(sample))
    print(input_file)
    list_of_percentages=[1,5,10,20,30,40,50,60,70,80,90]
    iterations=10
    do_rarefaction(input_file, list_of_percentages, iterations, path_to_out)