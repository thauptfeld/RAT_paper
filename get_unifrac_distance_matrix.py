# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 14:43:05 2021

@author: acer
"""


import subprocess
import json
import os

def make_distance_matrix(list_of_files,list_of_CAMI):
    distance_matrix={}
    for file_path in list_of_CAMI:
        f=file_path.split('/')[-1]
        distance_matrix[f]={}
        print(f)
        for tool_file in list_of_files:
            t=tool_file.split('/')[-1]
            if t==f:
                distance_matrix[f][t]=0

            elif t in distance_matrix and f in distance_matrix[t]:
                distance_matrix[f][t]=distance_matrix[t][f]
            
            else:
                print('\t'+t)
                distance_matrix[f][t]=50.0
                # print('\n\n'+file_path)
                # print('\n\n'+tool_file)
                cmd=['/home/tina/Documents/RAT/scripts/EMDUnifrac_from_CAMI.py', 
                     '-g', file_path,
                     '-r', tool_file]
                proc=subprocess.Popen(cmd, stdout=subprocess.PIPE)
                for p in proc.stdout:
                    if p.rstrip().isdigit:
                        distance_matrix[f][t]=float(p)
        with open('/home/tina/Documents/RAT/revision_unifrac/20231020_unifrac_distance.json', 'w') as outf:
            outf.write(json.dumps(distance_matrix, indent=4))
        
def make_small_distance_matrix(env, list_of_samples, mode):
    distances={}
    
    for s in list_of_samples:
        print(f'{env}{s} {mode}')
        cami_pro=f'/home/tina/Documents/RAT/revision_reference/gs_{env}.noplasmids_{s}.profile'
        # rat_pro=f'/home/tina/Documents/RAT/revision_profiles/{env}{s}.{mode}.1e-05.profile'
        tool_pro=f'/home/tina/Documents/RAT/revision_profiles/20230927.{env}{s}.{mode}.0.001.profile'
        
        cmd=['/home/tina/Documents/RAT/scripts/EMDUnifrac_from_CAMI.py', 
             '-g', cami_pro,
             '-r', tool_pro
             # '-r', rat_pro
             ]
        proc=subprocess.Popen(cmd, stdout=subprocess.PIPE)
        for p in proc.stdout:
            if p.rstrip().isdigit:
                distances[s]=float(p)
        with open(f'/home/tina/Documents/RAT/revision_unifrac/20231020_unifrac_{env}_{mode}.json', 'w') as outf:
            outf.write(json.dumps(distances, indent=4))
            
if __name__=='__main__':
    marine=[0,1,2,3,4,5,6,7,8,9]
    plant=[2,3,5,8,10,12,14,15,18,19]
    
    # make_small_distance_matrix('marine', marine, 'sensitive')
    # make_small_distance_matrix('marine', marine, 'contig')
    # make_small_distance_matrix('marine', marine, 'robust')
    
    # make_small_distance_matrix('plant', plant, 'sensitive')
    # make_small_distance_matrix('plant', plant, 'contig')
    # make_small_distance_matrix('plant', plant, 'robust')
    
    make_small_distance_matrix('marine', marine, 'kaiju')
    make_small_distance_matrix('marine', marine, 'kraken2')
    make_small_distance_matrix('marine', marine, 'bracken')
    make_small_distance_matrix('marine', marine, 'centrifuge')
    make_small_distance_matrix('plant', plant, 'kaiju')
    make_small_distance_matrix('plant', plant, 'kraken2')
    make_small_distance_matrix('plant', plant, 'bracken')
    make_small_distance_matrix('plant', plant, 'centrifuge')
    
