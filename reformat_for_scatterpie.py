# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 09:35:16 2023

@author: acer
"""

data_dict={}

c='Correctly classified'
i='Incorrectly classified'
u='Unclassified'

with open('../evaluations/20231020.CAMI2_marine.csv') as inf:
    for line in inf:
        correct, rank, value, tool, sample=line.strip().split(',')
        if sample not in data_dict:
            data_dict[sample]={}
        if rank not in data_dict[sample]:
            data_dict[sample][rank]={}
        if tool not in data_dict[sample][rank]:
            data_dict[sample][rank][tool]={}
        data_dict[sample][rank][tool][correct]=float(value)
        
with open('../evaluations/20231024.CAMI2_marine_scatterpie.csv', 'w') as outf:
    outf.write('tool,rank,sample,TPR,A,B,C,D\n')
    for sample in data_dict:
        for rank in data_dict[sample]:
            for tool in data_dict[sample][rank]:
                tool_dict=data_dict[sample][rank][tool]
                TPR=tool_dict[c]/(tool_dict[c]+tool_dict[i])
                
                B=tool_dict[i]+tool_dict[u]
                A=1-B
                
                D=tool_dict[u]
                C=1-D
                outf.write(f'{tool},{rank},{sample},{TPR},{A},{B},{C},{D}\n')
                
                
data_dict={}
with open('../evaluations/20231020_CAMI2_plant.csv') as inf:
    for line in inf:
        correct, rank, value, tool, sample=line.strip().split(',')
        if sample not in data_dict:
            data_dict[sample]={}
        if rank not in data_dict[sample]:
            data_dict[sample][rank]={}
        if tool not in data_dict[sample][rank]:
            data_dict[sample][rank][tool]={}
        data_dict[sample][rank][tool][correct]=float(value)
        
with open('../evaluations/20231024.CAMI2_plant_scatterpie.csv', 'w') as outf:
    outf.write('tool,rank,sample,TPR,A,B,C,D\n')
    for sample in data_dict:
        for rank in data_dict[sample]:
            for tool in data_dict[sample][rank]:
                tool_dict=data_dict[sample][rank][tool]
                TPR=tool_dict[c]/(tool_dict[c]+tool_dict[i])
                B=tool_dict[i]+tool_dict[u]
                A=1-B
                
                D=tool_dict[u]
                C=1-D
                outf.write(f'{tool},{rank},{sample},{TPR},{A},{B},{C},{D}\n')
                
                
# with open('../evaluations/20230127_CAMI2_mousegut.csv') as inf:
#     for line in inf:
#         correct, rank, value, tool, sample=line.strip().split(',')
#         if sample not in data_dict:
#             data_dict[sample]={}
#         if rank not in data_dict[sample]:
#             data_dict[sample][rank]={}
#         if tool not in data_dict[sample][rank]:
#             data_dict[sample][rank][tool]={}
#         data_dict[sample][rank][tool][correct]=float(value)
        
# with open('../evaluations/20231027.CAMI2_mousegut_scatterpie.csv', 'w') as outf:
#     outf.write('tool,rank,sample,TPR,A,B,C,D\n')
#     for sample in data_dict:
#         for rank in data_dict[sample]:
#             for tool in data_dict[sample][rank]:
#                 tool_dict=data_dict[sample][rank][tool]
#                 TPR=tool_dict[c]/(tool_dict[c]+tool_dict[i])
#                 B=tool_dict[i]+tool_dict[u]
#                 A=1-B
                
#                 D=tool_dict[u]
#                 C=1-D
#                 outf.write(f'{tool},{rank},{sample},{TPR},{A},{B},{C},{D}\n')
                
        