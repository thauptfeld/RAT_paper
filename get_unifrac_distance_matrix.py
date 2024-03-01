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
        with open('/home/tina/Documents/RAT/202301_profiles_CAMI/20230130_unifrac_distance.json', 'w') as outf:
            outf.write(json.dumps(distance_matrix, indent=4))
        

if __name__=='__main__':
    path_to_CAMI='/home/tina/Documents/RAT/profiles_CAMI/'
    list_of_CAMI=[path_to_CAMI+f for f in os.listdir(path_to_CAMI)]
    path_to_files='/home/tina/Documents/RAT/202301_profiles_CAMI/'
    list_of_files=[path_to_files+f for f in os.listdir(path_to_files)]
    make_distance_matrix(list_of_files,list_of_CAMI)
    