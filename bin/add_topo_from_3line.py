#!/usr/bin/env python3
import sys

fasta_file = sys.argv[1]
threeline_file = sys.argv[2]

ids = []
with open(fasta_file, 'r') as fasta_handle:
    for line in fasta_handle:
        if line[0] == '>':
            ids.append(line.strip())

with open(threeline_file, 'r') as threeline_handle:
    pid = ''
    for i, line in enumerate(threeline_handle):
        mod = i % 3
        if mod == 0:
            pid = line.strip()
        elif mod == 1:
            seq = line.strip()
        else:
            if pid in ids:
                print(pid)
                print(seq)
                print(line.strip())
