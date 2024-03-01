# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 10:26:32 2023

@author: haup0007
"""


import sys
import numpy
import json

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, numpy.integer):
            return int(obj)
        if isinstance(obj, numpy.floating):
            return float(obj)
        if isinstance(obj, numpy.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

paths={'marine': '/net/phage/linuxhome/mgx/people/tina/RAT/revision/read_classifiers_new_db_marine/results/',
       'plant': '/net/phage/linuxhome/mgx/people/tina/RAT/revision/read_classifiers_new_db_plant/results/',
       'smp': '/net/phage/linuxhome/mgx/people/tina/RAT/new_read_classifier_results/',
       'gw': '/net/phage/linuxhome/mgx/people/tina/RAT/revision/centrifuge_on_groundwater/results/'
       }

env=sys.argv[1]
sample=sys.argv[2]

if env=='gw':
    env2=''
else:
    env2=env
centrifuge=paths[env]+'{}{}/centrifuge/{}{}.centrifuge.out'.format(env2,sample,env2,sample)


counts={}
all_counts=[]

print('Collecting hit lengths...')

with open(centrifuge) as inf:
    current=''
    for line in inf:
        if not line.startswith('read'):
            line=line.split()
            
            read=line[0]
            if read!=current:
                current=read
                hit=line[5]
                
                all_counts.append(int(hit))
                if hit not in counts:
                    counts[hit]=0
                counts[hit]+=1
                
print('Calculating distributions...')
stats={}
hit_lengths=numpy.array(all_counts)    
stats['avg']=numpy.mean(hit_lengths)
stats['std']=numpy.std(hit_lengths)
stats['minimum']=numpy.min(hit_lengths)
stats['maximum']=numpy.max(hit_lengths)
stats['quart1']=numpy.percentile(hit_lengths, 25)
stats['quart2']=numpy.percentile(hit_lengths, 50)
stats['quart3']=numpy.percentile(hit_lengths, 75)


print('Writing output...')
with open('/net/phage/linuxhome/mgx/people/tina/RAT/revision/check_datasets/{}{}.centrifuge.hitlengths.txt'.format(env2,sample),'w') as outf:
    outf.write('Stats:\n')
    outf.write(json.dumps(stats, indent=4, cls=NpEncoder))
    outf.write('\n\nCounts:\n')
    outf.write(json.dumps(counts, indent=4))