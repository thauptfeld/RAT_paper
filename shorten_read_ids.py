# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 11:44:49 2023

@author: haup0007
"""


import sys

def shorten_line(line):
    return '{}\n'.format(line[:-3])

def shorten_every_fourth_line(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    modified_lines = [line if i % 4 != 0 else shorten_line(line) for i, line in enumerate(lines)]

    with open(filename, 'w') as file:
        file.writelines(modified_lines)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    shorten_every_fourth_line(filename)