# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 09:47:46 2023

@author: Ernestina Hauptfeld
"""


import os
import json

def append_tool_json(json_file, outf, tool, env):
    tool_table=json.load(open(json_file))
    
    for sample in tool_table:
        value=float(tool_table[sample])
        if not (env=='plant' and sample in ['8','15']):
            outf.write(f'{tool},{env}{sample},{value}\n')
    
    return
    

if __name__=='__main__':
    
    plant=['../unifrac/'+f for f in os.listdir('../unifrac') if 'plant' in f]
    marine=['../unifrac/'+f for f in os.listdir('../unifrac') if 'marine' in f]
    
    with open('../evaluations/20231028_plant_unifrac.csv', 'w') as outf:
        outf.write('X,variable,value\n')
        for p in plant:
            tool=p.split('_')[-1][:-5]
            
            append_tool_json(p, outf, tool, 'plant')
            
    with open('../evaluations/20231028_marine_unifrac.csv', 'w') as outf:
        outf.write('X,variable,value\n')
        for m in marine:
            tool=m.split('_')[-1][:-5]
            
            append_tool_json(m, outf, tool, 'marine')