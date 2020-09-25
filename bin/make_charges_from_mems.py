#!/usr/bin/env python3
"""Membrane charge caluclations."""
# import matplotlib.pyplot as plt
import numpy as np
import pickle
import argparse
import sys

data = {}
positive = 'KR'
negative = 'DE'
charged = 'KRDE'
chargedplus = 'KRDEH'
positiveplus = 'RHK'

parser = argparse.ArgumentParser()

parser.add_argument("helices", type=str, help="In membrane/helix pickle")
parser.add_argument("out_pickle", type=str, help="Out pickle")

args = parser.parse_args()

helicies = pickle.load(open(args.helices, 'rb'))
totalmemstring = ''
totalmembranes = 0
samehitcounter = [0, 0, 0, 0, 0, 0, 0]
hitcounter = [0, 0, 0, 0, 0, 0, 0]
misscounter = [0, 0, 0, 0, 0, 0, 0]
totalcharges = []
pairs = [[], [], [], [], [], [], []]
mem_length = 17
aas = 'ACDEFGHIKLMNPQRSTVWY'
aaHits = np.zeros([7, 20, 20])
aaPairs = np.zeros(7)
# Dim 1, first and second in pair
aaCount = np.zeros([2, 7, 20])
memLens = []
prolinMems = []
pPositions = []
for key, membranes in helicies.items():
    for start, mem in membranes:
        # Uncomment next line for only mid mem
        if len(mem) < mem_length:
            continue
        midMem = mem[5:-5]
        # midMem = mem
        memLen = len(midMem)

        for place, aa in enumerate(midMem):
            # Now count the for each gap
            for i in range(1, 8):
                if memLen <= place + i:
                    break
                if aa in aas:
                    first = aas.index(aa)
                    # print(midMem[place+i])
                    second_aa = midMem[place + i]
                    if second_aa in aas:
                        aaPairs[i - 1] += 1
                        second = aas.index(midMem[place + i])
                        aaHits[i - 1][first][second] += 1
                        # Count both first and second AA for each gap sized i
                        aaCount[0][i - 1][first] += 1
                        aaCount[1][i - 1][second] += 1

with open(args.out_pickle, 'wb') as dataPickle:
    pickle.dump([aaCount, aaPairs, aaHits], dataPickle)
