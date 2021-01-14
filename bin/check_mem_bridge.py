#!/usr/bin/env python3
import argparse
import pickle
import pandas as pd
import sys
from collections import Counter
from Bio.PDB.Polypeptide import three_to_one

parser = argparse.ArgumentParser()

parser.add_argument("mems", type=str, help="Membrane pickle")
parser.add_argument("threeline", type=str, help="3line in file")
parser.add_argument("bridges", type=str, help="Bridge pickle")

args = parser.parse_args()
membranes = pickle.load(open(args.mems, 'rb'))
bridges = pickle.load(open(args.bridges, 'rb'))

threeline_data = {}
with open(args.threeline, 'r') as TMHandle:
    pdb_id = ''
    rowNum = 0
    fa = ''
    for line in TMHandle:
        if line[0] == '>':
            # pid = line[1:].strip().split('|')[1]
            if '|' in line:
                pdb_id = line[1:].strip().split('|')[1]
            else:
                pdb_id = line[1:].strip()
            # pid = line[1:].strip()# .split('|')[2]
        elif rowNum % 3 == 1:
            fa = line.strip()
        elif rowNum % 3 == 2:
            if not "MISSING" in pdb_id:
                threeline_data[pdb_id] = [fa, line.strip()]
        else:
            print("You should not be here...")
        rowNum += 1
print("Number of proteins: {}".format(len(threeline_data)))

print(bridges.keys())
print("Number of total bridges: {}".format(len(bridges["all"])))
print("Number of membrane bridges: {}".format(len(bridges["mems"])))
print("Number of local bridges: {}".format(len(bridges["local"])))
r = 0
w = 0
index = 0
for key, bridge in bridges["local"].items():
    if key in threeline_data:
        seq = threeline_data[key][0]
        for b in bridge:
            first_aa = three_to_one(b[1])
            first_ix = b[0]
            second_aa = three_to_one(b[3])
            second_ix = b[2]
            if first_ix > len(seq):
                index += 1
                continue
            if first_aa == seq[first_ix-1]:
                r += 1
            else:
                print(key)
                print(seq)
                print(b)
                w += 1

            # print(seq[first_ix-1], first_aa)
            # print(seq[second_ix-1], second_aa)
print(r, w, index)

    # print(key, bridge)
