#!/usr/bin/env python3
"""Membrane charge caluclations."""
# import matplotlib.pyplot as plt
# import numpy as np
# import math
import pickle
import dgCalc
import saltBridges
# import parser
import os.path
import argparse
import urllib.request
from urllib.error import HTTPError
import sys

parser = argparse.ArgumentParser()

parser.add_argument("bridge_file", type=str, help="Bridge file")
parser.add_argument("mems_pickle", type=str, help="Pickle file of membranes")
parser.add_argument("threeline", type=str, help="3line in file")
parser.add_argument("-b", "--bridges", type=int, default=1, help="Required connections per bridge")
parser.add_argument("-g", "--gap", type=int, default=7, help="Must be within number of residues, (0 = unlimited)")
parser.add_argument("-t", "--tolerant", type=str, default="True", help="Use tolerant membranes (also include m, not just M)")
parser.add_argument("-s", "--stats", type=bool, default=False, help="Only calculate stats, do not generate list")
parser.add_argument("-a3m", "--a3m", type=bool, default=False, help="Special steps for a3m, skip bridges")

args = parser.parse_args()
mem_length = 17
debug = False
stats = {'aas':0, 'chr':0, 'pos':0, 'neg':0, 'pospair':0, 'negpair':0, 'chargedpair':0, 'opppair':0}
con_req = args.bridges  # How many saltbridge connections to make a bridge? Default 1
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

local_bridges = pickle.load(open(args.bridge_file, 'rb'))['local']

def countMems(topoStr, tol=True):
    if tol:
        mem_letters = 'Mm'
    else:
        mem_letters = 'M'
    num = 0
    currTopo = topoStr[0]
    # print(currTopo)
    for i, topo in enumerate(topoStr):
        # print(currTopo, topo)
        if topo != currTopo:
            if currTopo in mem_letters:
                # print("In loop")
                num += 1
            currTopo = topo
    return num


def whatMem(topoStr, aaIndex, tol=True):
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
membraneSet = set()
totalMems = 0
totalLongMems = 0
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
    # if key == "1SU4A":
    #     debug = True
    totalMems += len(membranes)
    for mem_place, mem in membranes:
        mem = mem.strip()
        if len(mem) < mem_length:
            # print("Short mem")
            continue
        totalLongMems += 1  # If membrane is longer than or equal 17
        # Uncomment next line for only mid mem
        edge = 5
        midMem = mem[edge:-edge]
        fullMem = mem
        if '?' in fullMem:
            continue
        # midMem = mem
        memLen = len(midMem)
        #### If not tolerant, exclude membranes with 'm', indicating DSSP is not fully correct
        if args.tolerant == "False":
            topo = TMdata[key][1]
            t_start = mem_place
            t_end = t_start + len(mem)
            topo = topo[t_start:t_end][edge:-edge]
            if 'm' in topo:
                continue
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
                    membraneSet.add(key + "-" + str(mem_place))
                    # print(membraneSet)
                    # sys.exit()
                    deltaG = dgCalc.calc_segment_DG(fullMem)
                    # Add five is midmem
                    globalPlace = mem_place + edge + 1  # Residue numbering
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
print("Total membranes > 17: ", totalLongMems)
print("Proteins that contain charged pairs: ", len(proteinSet))
print("Membrane regions with charged pairs: ", len(membraneSet))  # Number of helices with charged pairs
print("Total number of charged pairs: ", len(proteinsExamples))  # Number of total pairs
if args.stats or args.a3m:
    sys.exit()
sortedProt = sorted(proteinsExamples.keys())
debug = False
for key in sortedProt:
    midMem, place, i, globalPlace, aa, secAA, dG = proteinsExamples[key]
    pureKey = key[:5]
    # if pureKey[:4] not in topcons_added_keys:
    #     continue
    chain = key[4]
    filename = pureKey[:4] + '.pdb'
    filepath = 'data/pdbFiles/pdb' + pureKey[:4] + '.ent'
    saltpath = 'data/bridgeFiles/' + pureKey[:4] + chain + str(con_req) + 'Bridges.pickle'
    does_have_saltbridge = False
    if (i == 1 or i == 3 or i == 4) and ((aa in positiveplus and secAA in negative) or (aa in negative and secAA in positiveplus)):
        if not os.path.exists(filepath):
            try:
                urllib.request.urlretrieve(pdbURL + filename, filepath)
            except HTTPError as err:
                print("{} not exists, skipping...".format(filename), file=sys.stderr)
                continue
        if pureKey in local_bridges.keys():
            bridges = local_bridges[pureKey]
        elif not os.path.exists(saltpath):
            bridges = saltBridges.calcSaltBridges(filepath, pureKey[4], 4, con_req)
            pickle.dump(bridges, open(saltpath, 'wb'))
        else:
            bridges = pickle.load(open(saltpath, 'rb'))
        numMem = countMems(TMdata[pureKey][1])
        memNumber = whatMem(TMdata[pureKey][1], globalPlace)
        span = 'multi span' if numMem > 1 else 'single span'
        # if pureKey[:4] == "6RQP":
        #     print(globalPlace ,chain ,i, aa, secAA ,place, midMem)
        for bridge in bridges:
            # print("In bridges")
            # if pureKey[:4] == "6RQP":
            #     print(bridge)
            # sys.exit()
            # if bridge[1][3]-bridge[0][3] == i:
            # [164, 'ARG', 166, 'GLU', 'A', 3.0681937683268963]
            #print(bridge)
            if args.gap != 0 and bridge[2] - bridge[0] > args.gap:
                continue
            if (bridge[0] == globalPlace
                    and bridge[2] == (globalPlace + i) and bridge[4] == chain):
            # if ((bridge[0][3] == globalPlace
            #    and bridge[0][2] == chain)\
            #    and (bridge[1][3] == (globalPlace + i)
            #        and bridge[1][2] == chain)) or ((bridge[0][3] == globalPlace + i
            #    and bridge[0][2] == chain)\
            #    and (bridge[1][3] == (globalPlace + 2*i)
            #        and bridge[1][2] == chain)):
                # and aaMap[bridge[0][1]] == aa:
                # print(bridge)
                # print(str(i))
                # if (bridge[1][3] - bridge[0][3]) == i:
                does_have_saltbridge = True
                print(pureKey[:4] + "(" + pureKey[4] + ")",
                      span,
                      aa + secAA + '-pair with gap ' + str(i),
                      'startindex ' + str(globalPlace) + '(' +
                      str(memNumber) + ' membrane region)',
                      'dG ' + '{0:.2f}'.format(dG),
                      sep=', ')
                # if len(bridge[0][4]) > 0:
                #     first_alt = "alt. conf {}".format(bridge[0][4])
                # else:
                #     first_alt = ''
                # if len(bridge[1][4]) > 0:
                #     second_alt = "alt. conf {}".format(bridge[1][4])
                # else:
                #     second_alt = ''

                print(bridge[0],
                      bridge[1],
                      bridge[2],
                      bridge[3],
                      "{0:.2f}".format(bridge[5]) + "Å")
        # sys.exit()
    # print(bridges)
        if does_have_saltbridge:
            print(TMdata[pureKey][0])
            print(TMdata[pureKey][1])
#     if p[9] > 5:
#         # 4c9g(A), multi span (6 membranes), DK-pair with gap 1,
#         # startindex 108(second membrane region), dG 7.86
#         # print(p)
