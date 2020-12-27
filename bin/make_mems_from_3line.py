#!/usr/bin/env python3
import argparse
import pickle

parser = argparse.ArgumentParser()

parser.add_argument("in_file", type=argparse.FileType('r'), help="3line infile")
parser.add_argument("out_file", type=str, help="Output file")
parser.add_argument("-t", "--tolerance", type=bool, default=True, help="Be tolerant for dssp fault?")

args = parser.parse_args()

if args.tolerance:
    accept_letters = 'MmH'
else:
    accept_letters = 'MH'

data = {}
helicies = {}
membranes = {}
pid = ''
seq = ''
currentLine = 1
mem_len = 17  # At least 17 in length

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
        # print(pid)
        topo = line.strip()
        # print(topo)
        curr_letter = topo[0]
        start_pos = 0
        helices = []
        for i, c in enumerate(topo):
            if c != curr_letter:
                if curr_letter in accept_letters:  # Also include H to have compatability wityh dssp globular topology
                    memSeq = seq[start_pos:i]
                    if len(memSeq) >= mem_len and 'X' not in memSeq and\
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
