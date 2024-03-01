# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 18:53:18 2023

@author: acer
"""
import json
import csv

def write_nested_dict_to_csv(nested_dict, file_path):
    # Get the keys of the first dictionary (which will be used as column names)
    column_names = list(nested_dict.keys())

    # Get the keys of the nested dictionary (which will be used as row names)
    row_names = list(nested_dict[column_names[0]].keys())

    # Open the CSV file for writing
    with open(file_path, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)

        # Write the header row
        writer.writerow([''] + column_names)

        # Write each row of data
        for row_name in row_names:
            row_data = [row_name]
            for col_name in column_names:
                row_data.append(nested_dict[col_name][row_name])
            writer.writerow(row_data)
            
            
if __name__=='__main__':
    in_dict=json.load(open('../RAT_wbin_reads_per_step.json'))
    write_nested_dict_to_csv(in_dict, '../RAT_wbin_reads_per_step.csv')