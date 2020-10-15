#!/usr/bin/env python3
"""Membrane charge caluclations."""

# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# import numpy as np
import os
import sys
import pickle
import argparse
import random

parser = argparse.ArgumentParser()

parser.add_argument("top", type=str, help="Topological file")
parser.add_argument("a3m_folder", type=str, help="A3m file folder")
parser.add_argument("out_file", type=str, help="Out file")

args = parser.parse_args()
data = {}
positive = 'KR'
negative = 'DE'
charged = 'KRDE'
chargedplus = 'KRDEH'
positiveplus = 'RHK'
# with open('../data/red_TM_prop.fa') as faFile:
# with open('../data/red_TM_prop.top') as topFile:
random.seed(3141)
k=200
mem_part_length = 15
helicies = {}
with open(args.top) as topFile:
    name = ''
    for line in topFile:
        if line[0] == '>':
            name = line[1:].strip()
            helicies[name] = []
        else:
            topo = line.strip()

            curr_letter = topo[0]
            start_pos = 0
            for i, c in enumerate(topo):
                if c != curr_letter:
                    if curr_letter == 'M':
                        helicies[name].append((start_pos, i))
                        # for line in faFile:
                        #     if line[0] == '>':
                        #         continue
                        #     prefix = line[:start_pos]
                        #     inserts = sum(1 for c in prefix if c.islower())
                        #     print(line)
                    start_pos = i
                    curr_letter = c
# print(data)


membranes = {}
name = ''
raw_lines = []
for a3mfile in os.listdir(args.a3m_folder):
    with open(args.a3m_folder + '/' + a3mfile) as faFile:
        name = faFile.readline()[1:].strip()
        # print(name)
        # originalSequence = faFile.readline().strip()
        helixNumbers = helicies[name]
        # print(helixNumbers)
        raw_lines = []
        currentLine = 1
        for line in faFile:
            # print(line)
            if line[0] != '>':
                mems = []
                currentLine += 1
                # print(line.strip())
                for start, end in helixNumbers:
                    # print(start, end)
                    prefix = line[:start]
                    insertions = sum([1 for c in prefix if c.islower()])
                    # print(insertions)
                    while True:
                        # start += insertions
                        prefix = line[:start + insertions]
                        insertionsExt = sum(1 for c in prefix if c.islower())
                        if insertions == insertionsExt:
                            break
                        insertions = insertionsExt
                    memPart = line[start + insertions:end + insertions]
                    if memPart.count('-') == 0 and\
                       sum(1 for c in memPart if c.islower()) == 0 and\
                       'X' not in memPart and len(memPart) >= mem_part_length:
                        # membranes[name].append(memPart)
                        mems.append([start, memPart])
                    # print(start+insertions, end+insertions)
                    # print(memPart)
                # print(mems)
                # print(len(mems))
                # sys.exit()
                if len(mems) > 0:
                    raw_lines.extend(mems)
            # if currentLine > 100:
            #     break
    if len(raw_lines) > k:
        membranes[name] = random.choices(raw_lines, k=k)
    elif len(raw_lines) > 0:
        membranes[name] = raw_lines

# print(membranes)

# print(len(membranes))
with open(args.out_file, 'wb') as saveFile:
    pickle.dump(membranes, saveFile)
# with open('membraneUniref.pickle','rb') as loadFile:
    # membranes = pickle.load(loadFile)
# print(len(membranes))
# print(len(membranes['1c17M']))
#     helicies[name] = []
#     curr_letter = topo[0]
#     start_pos = 0
#     for i, c in enumerate(topo):
#         if c != curr_letter:
#             if curr_letter == 'M':
#                 helicies[name].append((start_pos,i))
#                 # for line in faFile:
#                 #     if line[0] == '>':
#                 #         continue
#                 #     prefix = line[:start_pos]
#                 #     inserts = sum(1 for c in prefix if c.islower())
#                 #     print(line)
#             start_pos = i
#             curr_letter = c
# print(helicies)
#     for line in faFile:
#         if line[0] == '>':
#             name = line[1:].strip()
#         else:
#             data[name] = [line.strip()]
# for key, topo in data.items():
#     # for test in re.finditer('[O]{1,20}M+', value[1]):
#     # print(key, value)
#     helicies[key] = []
#     curr_letter = topo[0]
#     start_pos = 0
#     for i, c in enumerate(topo):
#         if c != curr_letter:
#             if curr_letter == 'M':
#                 helicies[key].append((start_pos,i))
#                 # helicies[key].append(value[1][start_pos:i])
#             start_pos = i
#             curr_letter = c
# print(helicies)
# # print(len(helicies))
# simpledistances = []
# permdistances = []
# posnegdist = []
# negposdist = []
# diffdist = []
# rangedist = []
# rangedistplus = []
# for key, value in helicies.items():
#     for mem in value:
#         memlen = len(mem)
#         distList = []
#         distListPlus = []
#         curr = [None, 0]
#         currplus = [None, 0]
#         currlog = [None, 0]
#         # print(key, value)
#         for i, c in enumerate(mem):
#             ######### Regular charge variant #######
#             if c in charged:
#                 # If c has any charge (in regular charged group)
#                 distList.append(i)
#                 if curr[0] is None:
#                     # If first charged, set up initial
#                     curr[0] = 'pos' if c in positive else 'neg'
#                     curr[1] = i
#                 for x in range(1, 8):
#                     # Check the coming 7 residues for opposite charge
#                     if i + x >= memlen:
#                         break
#                     if c in positive:
#                         if mem[i+x] in negative:
#                             rangedist.append(x)
#                     if c in negative:
#                         if mem[i+x] in positive:
#                             rangedist.append(x)
#             ####### Charged plus variant #########
#             if c in chargedplus:
#                 # If c has any charge (in regular charged group)
#                 distListPlus.append(i)
#                 if currplus[0] is None:
#                     # If first charged, set up initial
#                     currplus[0] = 'pos' if c in positiveplus else 'neg'
#                     currplus[1] = i
#                 for x in range(1, 8):
#                     # Check the coming 7 residues for opposite charge
#                     if i + x >= memlen:
#                         break
#                     if c in positiveplus:
#                         if mem[i+x] in negative:
#                             rangedistplus.append(x)
#                     if c in negative:
#                         if mem[i+x] in positiveplus:
#                             rangedistplus.append(x)
#             if c in positive and curr[0] == 'neg':
#                 # print("in positive")
#                 negposdist.append(i-curr[1])
#                 diffdist.append(i-curr[1])
#                 curr[0], curr[1] = 'pos', i
#             elif c in negative and curr[0] == 'pos':
#                 posnegdist.append(i-curr[1])
#                 diffdist.append(i-curr[1])
#                 curr[0], curr[1] = 'neg', i
#             if c in chargedplus:
#                 if curr[0] is None:
#                     # If first charged, set up initial
#                     curr[0] = 'pos' if c in positiveplus else 'neg'
#                     curr[1] = i
#                 for x in range(1, 8):
#                     # Check the coming 7 residues for opposite charge
#                     if i + x >= memlen:
#                         break
#                     if c in positiveplus:
#                         if mem[i+x] in negative:
#                             rangedistplus.append(x)
#                     if c in negative:
#                         if mem[i+x] in positiveplus:
#                             rangedistplus.append(x)
#         # print(distList)
#         simpledistances.extend(np.diff(distList))
#         for i, place in enumerate(distList):
#             for j in distList[i+1:]:
#                 permdistances.append(j-i)
#                 # print(j-i)
# # print(posnegdist)
# plt.hist(simpledistances, bins=max(simpledistances), align='left')
# plt.xticks(np.arange(0, max(simpledistances)+1, 1.0))
# plt.title("Distance between neighbooring charged AA's")
# plt.savefig('chargessimple.png')
# plt.clf()
# plt.hist(permdistances, bins=max(permdistances), align='left')
# plt.xticks(np.arange(0, max(simpledistances)+1, 1.0))
# plt.title("Distance between all permuations of charged AA's")
# plt.savefig('chargesperm.png')
# plt.clf()
# plt.hist(posnegdist, bins=max(posnegdist), align='left')
# plt.xticks(np.arange(0, max(posnegdist)+1, 1.0))
# plt.title("Distances between positive to negative charge")
# plt.savefig('postoneg.png')
# plt.clf()
# plt.hist(negposdist, bins=max(negposdist), align='left')
# plt.xticks(np.arange(0, max(negposdist)+1, 1.0))
# plt.title("Distances between negative to positive charge")
# plt.savefig('negtopos.png')
# plt.clf()
# plt.hist(diffdist, bins=max(diffdist), align='left')
# plt.xticks(np.arange(0, max(diffdist)+1, 1.0))
# plt.title("Distance between alternating charges")
# plt.savefig('diff.png')
# plt.clf()
# plt.hist(rangedist, bins=max(rangedist), align='left')
# plt.xticks(np.arange(0, max(rangedist)+1, 1.0))
# plt.title("Distance in range 1-7 from charged, alternating charge")
# plt.savefig('range.png')
# plt.clf()
# plt.hist(rangedistplus, bins=max(rangedist), log=True)
# plt.xticks(np.arange(0, max(rangedist)+1, 1.0))
# plt.title("Distance in range 1-7 from charged, alternating charge, inkl His")
# plt.savefig('rangeHis.png')
# # print(len(permdistances))
# # print(helicies)
