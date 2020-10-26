#!/usr/bin/env python3
"""Membrane charge caluclations."""
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import math
import pickle
import random
import argparse
import collections
data = {}
positive = 'KR'
negative = 'DE'
charged = 'KRDE'
chargedplus = 'KRDEH'
positiveplus = 'RHK'
parser = argparse.ArgumentParser()

parser.add_argument("input_file", type=str, help="Input membrane file")

args = parser.parse_args()
helicies = pickle.load(open(args.input_file, 'rb'))
# helicies = pickle.load(open('helicies.pickle','rb'))
totalmemstring = ''
totalmembranes = 0
samehitcounter = [0, 0, 0, 0, 0, 0, 0,0]
hitcounter = [0, 0, 0, 0, 0, 0, 0,0]
misscounter = [0, 0, 0, 0, 0, 0, 0,0]
totalcharges = []
pairs = [[], [], [], [], [], [], [], []]
chargeCount =[[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
poschargeCount =[[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
negchargeCount =[[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
aaCount =[0,0,0,0,0,0,0,0]
for key, membranes in helicies.items():
    # print(key)
    for mem_start, mem in membranes:
        # print(mem)
        if len(mem) < 17:
            continue
        seg_trimmed = mem[5:-5]
        totalmembranes += 1
        totalmemstring += seg_trimmed
        memLen = len(seg_trimmed)
        chargesinmem = 0
        for mem_start, aa in enumerate(seg_trimmed):
            if aa in charged:
                chargesinmem += 1
            for i in range(1, 9):
                if memLen <= mem_start + i:
                    break
                pairs[i-1].append(aa+seg_trimmed[mem_start+i])
                nextAA = seg_trimmed[mem_start + i]
                aaCount[i-1] += 1
                if aa in charged:
                    chargeCount[i-1][0] += 1
                    if aa in positive:
                        poschargeCount[i-1][0] += 1
                        if nextAA in negative:
                            negchargeCount[i-1][1] += 1
                            hitcounter[i-1] = hitcounter[i-1] + 1
                            chargeCount[i-1][1] += 1
                        elif nextAA in positive:
                            poschargeCount[i-1][1] += 1
                            samehitcounter[i-1] = samehitcounter[i-1] + 1
                            chargeCount[i-1][1] += 1
                        else:
                            misscounter[i-1] = misscounter[i-1] + 1
                    if aa in negative:
                        negchargeCount[i-1][0] += 1
                        if nextAA in positive:
                            poschargeCount[i-1][1] += 1
                            hitcounter[i-1] = hitcounter[i-1] + 1
                            chargeCount[i-1][1] += 1
                        elif nextAA in negative:
                            negchargeCount[i-1][1] += 1
                            samehitcounter[i-1] = samehitcounter[i-1] + 1
                            chargeCount[i-1][1] += 1
                        else:
                            misscounter[i-1] = misscounter[i-1] + 1
                else:
                    if nextAA in charged:
                        chargeCount[i-1][1] += 1
                        if nextAA in positive:
                            poschargeCount[i-1][1] += 1
                        elif nextAA in negative:
                            negchargeCount[i-1][1] += 1

        totalcharges.append(chargesinmem)

totalhits = np.array(hitcounter) + np.array(samehitcounter)
freq = collections.Counter(totalcharges)
logOdds = []
logOddsOpp = []
logOddsSame = []
ci = []
ciOpp = []
ciSame = []
for i in range(8):
    ###### For All charges ######
    a = totalhits[i]                # Number of pairs, there is a charge
    b = len(pairs[i])       # Number of total observed pairs
    c = chargeCount[i][0]*chargeCount[i][1]             # How many baseline pairs for this distance?
    d = aaCount[i]**2       # Number of non-hits for baseline
    odds = (a/b)/(c/d)
    logOdds.append(math.log(odds))
    ci.append((math.log(odds)+1.96*math.sqrt(1/a+1/b+1/c+1/d), math.log(odds)-1.96*math.sqrt(1/a+1/b+1/c+1/d)))

    ###### For Opp charges ######
    a = hitcounter[i]                 # Number of opposite pairs, only opposite
    b = len(pairs[i])       # Number of total observed pairs
    c = poschargeCount[i][0]*negchargeCount[i][1] + negchargeCount[i][0]*poschargeCount[i][1]  # Only calculate opposite pairs
    d = aaCount[i]**2       # Number of non-hits for baseline
    odds = (a/b)/(c/d)
    logOddsOpp.append(math.log(odds))
    ciOpp.append((math.log(odds)+1.96*math.sqrt(1/a+1/b+1/c+1/d), math.log(odds)-1.96*math.sqrt(1/a+1/b+1/c+1/d)))

    ###### For same charges ######
    a = samehitcounter[i]                 # Number of same pairs, only same
    b = len(pairs[i])       # Number of total observed pairs
    c = poschargeCount[i][0]*poschargeCount[i][1] + negchargeCount[i][0]*negchargeCount[i][1]  # Only calculate opposite pairs
    d = aaCount[i]**2       # Number of non-hits for baseline
    # Quickfix for topcons step 8, pseudo count, does not impact anything else
    if a == 0:
        a = 1
    if c == 0:
        c = 1
    odds = (a/b)/(c/d)
    logOddsSame.append(math.log(odds))
    ciSame.append((math.log(odds)+1.96*math.sqrt(1/a+1/b+1/c+1/d), math.log(odds)-1.96*math.sqrt(1/a+1/b+1/c+1/d)))

f, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 10)) # , sharex=True)
###### Totalt ######
y_r = [logOdds[i] - ci[i][0] for i in range(len(ci))]
threshold = 0
values = np.array(logOdds)
x = range(1, len(values)+1)
above_threshold = np.maximum(values - threshold, 0)
below_threshold = np.minimum(values, threshold)
ax1.bar(x, logOdds, yerr=y_r, color=['g', 'r', 'g', 'g', 'r', 'r', 'g', 'r'], alpha=0.8, align='center')
ax1.set_ylabel('All charges')
ax1.axhline(xmax=8, color='black')

###### Opp ######
y_r = [logOddsOpp[i] - ciOpp[i][0] for i in range(len(ciOpp))]
threshold = 0
values = np.array(logOddsOpp)
x = range(1, len(values)+1)
above_threshold = np.maximum(values - threshold, 0)
below_threshold = np.minimum(values, threshold)
ax2.bar(x, logOddsOpp, yerr=y_r, color=['g', 'r', 'g', 'g', 'r', 'r', 'g', 'r'], alpha=0.8, align='center')
ax2.set_ylabel('Opposite charges')
ax2.axhline(xmax=8, color='black')

###### Same ######
y_r = [logOddsSame[i] - ciSame[i][0] for i in range(len(ciSame))]
threshold = 0
values = np.array(logOddsSame)
x = range(1, len(values)+1)
above_threshold = np.maximum(values - threshold, 0)
below_threshold = np.minimum(values, threshold)
ax3.bar(x, logOddsSame, yerr=y_r, color=['g', 'r', 'g', 'g', 'r', 'r', 'g', 'r'], alpha=0.8, align='center')
ax3.set_ylabel('Same charges')
ax3.set_xlabel('Step')
ax3.axhline(xmax=8, color='black')

ax1.set_xticks([])
ax1.xaxis.set_tick_params(length=0)
ax2.xaxis.set_tick_params(length=0)
ax2.set_xticks([])
ax3.set_xticks([1, 2, 3, 4, 5, 6, 7, 8])

degree_cmap = mpl.colors.ListedColormap(mpl.cm.get_cmap('viridis_r').colors + mpl.cm.get_cmap('viridis').colors)
for ax in [ax1, ax2, ax3]:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    for i in range(8):
        color_index = (i*100 + 100) % 360
        ax.patches[i].set_facecolor(degree_cmap(color_index/360))

name = args.input_file.split('/')[-1].split('.')[0]
for ext in [".png", ".svg"]:
    plt.savefig('images/' + name + '_logodds' + ext)

