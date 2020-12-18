#!/usr/bin/env python3
"""Membrane charge caluclations."""
# import matplotlib.pyplot as plt
# import numpy as np
# import math
import pickle
import dgCalc
import saltBridges
import parser
import os.path
import argparse
import urllib.request
from urllib.error import HTTPError
import sys

parser = argparse.ArgumentParser()

parser.add_argument("mems_pickle", type=str, help="Pickle file of membranes")
parser.add_argument("threeline", type=str, help="3line in file")
parser.add_argument("-b", "--bridges", type=int, default=2, help="Required connections per bridge")
parser.add_argument("-t", "--tolerant", type=bool, default=False, help="Use tolerant membranes (also include m, not just M)")

args = parser.parse_args()
mem_length = 17

stats = {'aas':0, 'chr':0, 'pos':0, 'neg':0, 'pospair':0, 'negpair':0, 'chargedpair':0, 'opppair':0}
con_req = args.bridges  # How many saltbridge connections to make a bridge? Default 2
aaMap = {'ARG': 'R',
         'HIS': 'H',
         'LYS': 'K',
         'ASP': 'D',
         'GLU': 'E',
         'SER': 'S',
         'THR': 'T',
         'ASN': 'N',
         'GLN': 'Q',
         'CYS': 'C',
         'SEC': 'U',
         'GLY': 'G',
         'PRO': 'P',
         'ALA': 'A',
         'VAL': 'V',
         'ILE': 'I',
         'LEU': 'L',
         'MET': 'M',
         'PHE': 'F',
         'TYR': 'Y',
         'TRP': 'W'}

pdbURL = "https://files.rcsb.org/download/"
positive = 'KR'
negative = 'DE'
charged = 'KRDE'
chargedplus = 'KRDEH'
positiveplus = 'RHK'
helicies = pickle.load(open(args.mems_pickle,
                            'rb'))


def countMems(topoStr):
    num = 0
    currTopo = topoStr[0]
    # print(currTopo)
    for _, topo in enumerate(topoStr):
        # print(currTopo, topo)
        if topo != currTopo:
            if currTopo in 'Mm':
                # print("In loop")
                num += 1
            currTopo = topo
    return num


def whatMem(topoStr, aaIndex, tol=False):
    if tol:
        mem_letters = 'Mm'
    else:
        mem_letters = 'M'
    memNum = 1
    currTopo = topoStr[0]
    for i, topo in enumerate(topoStr):
        if topo != currTopo:
            if currTopo in mem_letters:
                if i >= aaIndex:
                    return memNum
                memNum += 1
            currTopo = topo
    return "Can't determine membrane number"


aas = 'ACDEFGHIKLMNPQRSTVWY'
proteinsExamples = {}
proteinSet = set()
fullproteinSet = set()
longproteinSet = set()
totalMems = 0
longMems = 0
correctMems = 0

TMdata = {}
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
                TMdata[pdb_id] = [fa, line.strip()]
        else:
            print("You should not be here...")
        rowNum += 1


for key, membranes in helicies.items():
    totalMems += len(membranes)
    for mem_place, mem in membranes:

        fullproteinSet.add(key)
        if len(mem) < mem_length:
            # print("Short mem")
            continue
        longproteinSet.add(key)
        longMems += 1
        # Uncomment next line for only mid mem
        edge = 5
        midMem = mem[edge:-edge]
        fullMem = mem
        if '?' in fullMem:
            continue
        correctMems += 1
        # midMem = mem
        memLen = len(midMem)
#         if mem in TMdata[key][0]:
        for place, aa in enumerate(midMem):
            stats['aas'] += 1
            if aa in charged:
                stats['chr'] += 1
            if aa in positive:
                stats['pos'] += 1
            if aa in negative:
                stats['neg'] += 1
            for i in range(1, 8):
                if memLen <= place + i:
                    break
                # aaPairs[i - 1] += 1
                # first = aas.index(aa)
                # print(midMem[place+i])
                # second = aas.index(midMem[place + i])
                secAA = midMem[place + i]
                if aa in charged and secAA in charged:
                    deltaG = dgCalc.calc_segment_DG(fullMem)
                    # Add five is midmem
                    globalPlace = mem_place + edge
                    # globalPlace = TMdata[key][0].index(midMem)
                    # Plus one for the globalPlace and plus 1 for the place
                    # in the membrane, 0 -> 1 offset
                    proteinsExamples[key +
                                     str(globalPlace + place) +
                                     str(globalPlace + place + i)] = \
                                    [midMem,
                                     place,
                                     i,
                                     globalPlace + place,
                                     aa,
                                     secAA,
                                     deltaG]
                    proteinSet.add(key)
                    if aa in positive and secAA in positive:
                        stats['pospair'] += 1
                    if aa in negative and secAA in negative:
                        stats['negpair'] += 1
                    if (aa in positive and secAA in negative) or (aa in negative and secAA in positive):
                        stats['opppair'] += 1
                    if aa in charged and secAA in charged:
                        stats['chargedpair'] += 1
                # aaHits[i - 1][first][second] += 1
                # # Count both first and second AA for each gap sized i
                # aaCount[0][i - 1][first] += 1
                # aaCount[1][i - 1][second] += 1
print("Number of proteins: ", len(helicies.keys()))
print("Total membranes: ", totalMems)
print("Long membranes: ", longMems)
print("Correct membranes: ", correctMems)
print("Proteins that contain charged pairs: ", len(proteinSet))
print("Membrane regions with charged pairs: ", len(proteinsExamples))

sortedProt = sorted(proteinsExamples.keys())
debug = False
has_bridge = set()
has_bridge_within_7 = set()
lacks_bridge = set()
for full_pdb_id in TMdata.keys():
    pdb_id = full_pdb_id[:4]
    chain = full_pdb_id[-1]
    filename = pdb_id + '.pdb'

    filepath = 'data/pdbFiles/' + pdb_id + '.pdb'
    saltpath = 'data/bridgeFiles/' + pdb_id + chain + str(con_req) + 'Bridges.pickle'
    if not os.path.exists(filepath):
        try:
            urllib.request.urlretrieve(pdbURL + filename, filepath)
        except HTTPError as err:
            print("{} not exists, skipping...".format(filename), file=sys.stderr)
            continue
    if not os.path.exists(saltpath):
        bridges = saltBridges.calcSaltBridges(filepath, chain, 4, con_req)
        pickle.dump(bridges, open(saltpath, 'wb'))
    else:
        bridges = pickle.load(open(saltpath, 'rb'))
    if len(bridges)>0:
        has_bridge.add(full_pdb_id)
    else:
        lacks_bridge.add(full_pdb_id)
    for bridge in bridges:
        if bridge[1][3] - bridge[0][3] < 8:
            has_bridge_within_7.add(full_pdb_id)
print("Has bridges {}".format(len(has_bridge)))
print("Has bridges within 7: {}".format(len(has_bridge_within_7)))
print("Lacks bridges {}".format(len(lacks_bridge)))
print("Has mems and bridges {}".format(len(has_bridge & fullproteinSet)))
print("Has long mems and bridges {}".format(len(has_bridge & longproteinSet)))
print("Has long mems,charge and bridges {}".format(len(has_bridge & proteinSet)))
print("Has long mems,charge and bridges within 7: {}".format(len(has_bridge_within_7 & proteinSet)))
# for key in sortedProt:
#     midMem, place, i, globalPlace, aa, secAA, dG = proteinsExamples[key]
#     pureKey = key[:5]
# 
#     # if pureKey[:4] not in topcons_added_keys:
#     #     continue
#     chain = key[4]
#     filename = pureKey[:4] + '.pdb'
#     filepath = 'data/pdbFiles/' + pureKey[:4] + '.pdb'
#     saltpath = 'data/bridgeFiles/' + pureKey[:4] + chain + str(con_req) + 'Bridges.pickle'
# 
#     if not os.path.exists(filepath):
#         try:
#             urllib.request.urlretrieve(pdbURL + filename, filepath)
#         except HTTPError as err:
#             print("{} not exists, skipping...".format(filename), file=sys.stderr)
#             continue
#     if not os.path.exists(saltpath):
#         bridges = saltBridges.calcSaltBridges(filepath, pureKey[4], 4, con_req)
#         pickle.dump(bridges, open(saltpath, 'wb'))
#     else:
#         bridges = pickle.load(open(saltpath, 'rb'))
#     if len(bridges)>0:
#         has_bridge += 1
#     else:
#         lacks_bridge += 1
#     does_have_saltbridge = False
# print(has_bridge)
# print(lacks_bridge)
    # if (i == 1 or i == 3 or i == 4) and ((aa in positiveplus and secAA in negative) or (aa in negative and secAA in positiveplus)):
    #     if not os.path.exists(filepath):
    #         try:
    #             urllib.request.urlretrieve(pdbURL + filename, filepath)
    #         except HTTPError as err:
    #             print("{} not exists, skipping...".format(filename), file=sys.stderr)
    #             continue
    #     if not os.path.exists(saltpath):
    #         bridges = saltBridges.calcSaltBridges(filepath, pureKey[4], 4, con_req)
    #         pickle.dump(bridges, open(saltpath, 'wb'))
    #     else:
    #         bridges = pickle.load(open(saltpath, 'rb'))
    #     
    #     numMem = countMems(TMdata[pureKey][1])
    #     memNumber = whatMem(TMdata[pureKey][1], globalPlace, args.tolerant)
    #     span = 'multi span' if numMem > 1 else 'single span'
    #     for bridge in bridges:
    #         # print("In bridges")
    #         # print(bridge)
    #         if (bridge[0][3] == globalPlace
    #            and bridge[0][2] == chain)\
    #            or (bridge[0][3] == (globalPlace + i)
    #            and bridge[0][2] == chain):
    #             # and aaMap[bridge[0][1]] == aa:
    #             # print(bridge)
    #             # print(str(i))
    #             # if (bridge[1][3] - bridge[0][3]) == i:
    #             does_have_saltbridge = True
    #             print(pureKey[:4] + "(" + pureKey[4] + ")",
    #                   span,
    #                   aa + secAA + '-pair with gap ' + str(i),
    #                   'startindex ' + str(globalPlace) + '(' +
    #                   str(memNumber) + ' membrane region)',
    #                   'dG ' + '{0:.2f}'.format(dG),
    #                   sep=', ')
    #             if len(bridge[0][4]) > 0:
    #                 first_alt = "alt. conf {}".format(bridge[0][4])
    #             else:
    #                 first_alt = ''
    #             if len(bridge[1][4]) > 0:
    #                 second_alt = "alt. conf {}".format(bridge[1][4])
    #             else:
    #                 second_alt = ''

    #             print(bridge[0][0],
    #                   first_alt,
    #                   bridge[0][1],
    #                   bridge[0][2],
    #                   bridge[0][3],
    #                   bridge[1][0],
    #                   second_alt,
    #                   bridge[1][1],
    #                   bridge[1][2],
    #                   bridge[1][3],
    #                   "{0:.2f}".format(bridge[2]) + "Å")
    # # print(bridges)
    #     if does_have_saltbridge:
    #         print(TMdata[pureKey][0])
    #         print(TMdata[pureKey][1])
#   #   if p[9] > 5:
#   #       # 4c9g(A), multi span (6 membranes), DK-pair with gap 1,
#   #       # startindex 108(second membrane region), dG 7.86
#   #       # print(p)
