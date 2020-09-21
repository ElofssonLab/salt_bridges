#!/usr/bin/env python3
""" Extract membrane regions from 3line topology file"""

import argparse
import pickle

parser = argparse.ArgumentParser()

parser.add_argument("in_file", type=argparse.FileType('r'), help="3line infile")
parser.add_argument("out_file", type=str, help="Output file")

args = parser.parse_args()

data = {}
helicies = {}
membranes = {}
pid = ''
seq = ''
currentLine = 1

for line in args.in_file:
    modu = currentLine % 3
    if line[0] == '>' and modu == 1:
        if '|' in line:
            pid = line[1:].strip().split('|')[1]
        else:
            pid = line[1:].strip()
    elif modu == 2:
        seq = line.strip()
    elif modu == 0:
        topo = line.strip()
        curr_letter = topo[0]
        start_pos = 0
        helices = []
        for i, c in enumerate(topo):
            if c != curr_letter:
                if curr_letter in 'MH':  # Also include H to have compatability wityh dssp globular topology
                    memSeq = seq[start_pos:i]
                    if 'X' not in memSeq and\
                       'U' not in memSeq and\
                       'Z' not in memSeq and\
                       'O' not in memSeq:
                        helices.append([start_pos, memSeq])
                start_pos = i
                curr_letter = c

        if len(helices) > 0 or not "MISSING" in pid:
            membranes[pid] = helices
    currentLine += 1

with open(args.out_file, 'wb') as saveFile:
    pickle.dump(membranes, saveFile)
