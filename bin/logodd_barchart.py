#!/usr/bin/env python3
"""Membrane charge caluclations."""
import matplotlib
matplotlib.use('Agg')
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
# with open('../data/red_TM_prop.fa') as faFile:
#     name = ''
#     for line in faFile:
#         if line[0] == '>':
#             name = line[1:].strip()
#         else:
#             data[name] = [line.strip()]
# with open('../data/red_TM_prop.top') as topFile:
#     name = ''
#     for line in topFile:
#         if line[0] == '>':
#             name = line[1:].strip()
#         else:
#             data[name].append(line.strip())
# helicies = {}
# for key, value in data.items():
#     # for test in re.finditer('[O]{1,20}M+', value[1]):
#     # print(key, value)
#     helicies[key] = []
#     curr_letter = value[1][0]
#     start_pos = 0
#     for i, c in enumerate(value[1]):
#         if c != curr_letter:
#             if curr_letter == 'M':
#                 helicies[key].append(value[0][start_pos:i])
#                 # helicies[key].append(value[1][start_pos:i])
#             start_pos = i
#             curr_letter = c
# pickle.dump(helicies,open("helicies.pickle","wb"))
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
    for place, mem in membranes:
        # print(mem)
        if len(mem) < 17:
            continue
        mem = mem[5:-5]
        totalmembranes += 1
        totalmemstring += mem
        memLen = len(mem)
        chargesinmem = 0
        for place, aa in enumerate(mem):
            if aa in charged:
                chargesinmem += 1
            for i in range(1, 9):
                if memLen <= place + i:
                    break
                pairs[i-1].append(aa+mem[place+i])
                nextAA = mem[place + i]
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

# seqs = random.sample(totalmemstring,int(len(totalmemstring)/20))
# loops = 1000
# examples = []
# loopnr = 1
# for randomshuffles in range(loops):
#     print('Running loop {}...'.format(loopnr))
#     loopnr+=1
#     basememstring = totalmemstring
#     random.seed(2017+randomshuffles)
#     basememstring = ''.join(random.sample(basememstring,len(basememstring)))
#     # print(basememstring, len(basememstring))
#     basepairs = [[],[],[],[],[],[],[]]
#     for step in range(0, len(basememstring), 20):
#         # print(step, basememstring[step:step+20])
#         mem = basememstring[step:step+20]
#         for place, aa in enumerate(mem):
#             for i in range(1,8):
#                 if len(mem) <= place + i:
#                     break
#                 # print(place, i)
#                 basepairs[i-1].append(aa+mem[place+i])
#     basehits = [[],[],[],[],[],[],[]]
#     for dist, pairlist in enumerate(basepairs):
#         temphits = 0
#         # print(pair)
#         for basepair in pairlist:
#             if basepair[0] in chargedplus and basepair[1] in chargedplus:
#                 # print("hit")
#                 temphits += 1
#         # print(temphits)
#         basehits[dist].append(temphits)
#     # print(basehits)
#     examples.append(basehits)
# # print(np.array(examples))
# print(np.mean(np.array(examples),axis=0))
# print(basehits)
# print(basepairs)
# Baseline using 1000 loops of shuffling and measuring how many pairs we find
# baselinehits = [48,46,43,41,38,36,33] # For original
# baselinehits = [3123, 2961, 2796, 2629, 2468, 2301, 2137]  # For uniref aligned
totalhits = np.array(hitcounter) + np.array(samehitcounter)
print("Samehitcounter: {}".format(samehitcounter))
print("Hit counter:    {}".format(hitcounter))
print("Total hits:     {}".format(totalhits))
print("Misscounter:    {}".format(misscounter))
# print("Baseline:       {}".format(baselinehits))
print("Total charges:  {}".format(sum(totalcharges)))
freq = collections.Counter(totalcharges)
print(freq)
print("Total membranes: {}".format(totalmembranes))
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
    print("Odds ratio for step {}: {:.2f}".format(i+1, (odds)))
    print("Log odds ratio for step {}: {:.2f}".format(i+1, math.log(odds)))
    logOdds.append(math.log(odds))
    print("Upper 95% CI: {:.2f}".format((math.log(odds)+1.96*math.sqrt(1/a+1/b+1/c+1/d))))
    print("Lower 95% CI: {:.2f}".format((math.log(odds)-1.96*math.sqrt(1/a+1/b+1/c+1/d))))
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
# print(logOdds)
# print(ci)
# print(logOddsSame)
# print(ciSame)
# print(logOddsOpp)
# print(ciOpp)

f, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
###### Totalt ######
y_r = [logOdds[i] - ci[i][0] for i in range(len(ci))]
threshold = 0
values = np.array(logOdds)
x = range(1, len(values)+1)
# plt.bar(range(len(logOdds)), logOdds, yerr=y_r, alpha=0.2, align='center')
# x = range(1, len(values)+1)
above_threshold = np.maximum(values - threshold, 0)
below_threshold = np.minimum(values, threshold)
# plt.bar(x, below_threshold, 0.35, color="g", yerr=y_r)
# plt.bar(x, above_threshold, 0.35, color="r", yerr=y_r, bottom=below_threshold)
ax1.bar(x, logOdds, yerr=y_r, color=['g', 'r', 'g', 'g', 'r', 'r', 'g', 'r'], alpha=0.8, align='center')
# plt.xticks(range(len(logOdds)), [str(i) for i in range(1,len(logOdds)+1)])
ax1.set_ylabel('All charges')
# ax1.set_xlabel('Step distance from first charged residue')
ax1.set_title('Log odds ratios')
ax1.axhline(xmax=8, color='black')
ax1.grid(True)
###### Opp ######
y_r = [logOddsOpp[i] - ciOpp[i][0] for i in range(len(ciOpp))]
threshold = 0
values = np.array(logOddsOpp)
x = range(1, len(values)+1)
# plt.bar(range(len(logOdds)), logOdds, yerr=y_r, alpha=0.2, align='center')
# x = range(1, len(values)+1)
above_threshold = np.maximum(values - threshold, 0)
below_threshold = np.minimum(values, threshold)
# plt.bar(x, below_threshold, 0.35, color="g", yerr=y_r)
# plt.bar(x, above_threshold, 0.35, color="r", yerr=y_r, bottom=below_threshold)
ax2.bar(x, logOddsOpp, yerr=y_r, color=['g', 'r', 'g', 'g', 'r', 'r', 'g', 'r'], alpha=0.8, align='center')
# plt.xticks(range(len(logOdds)), [str(i) for i in range(1,len(logOdds)+1)])
ax2.set_ylabel('Opposite charges')
# ax2.set_xlabel('Step distance from first charged residue')
# ax2.set_title('Log odds ratio for opposite pairs')
ax2.axhline(xmax=8, color='black')
ax2.grid(True)
###### Totalt ######
y_r = [logOddsSame[i] - ciSame[i][0] for i in range(len(ciSame))]
threshold = 0
values = np.array(logOddsSame)
x = range(1, len(values)+1)
# plt.bar(range(len(logOdds)), logOdds, yerr=y_r, alpha=0.2, align='center')
# x = range(1, len(values)+1)
above_threshold = np.maximum(values - threshold, 0)
below_threshold = np.minimum(values, threshold)
# plt.bar(x, below_threshold, 0.35, color="g", yerr=y_r)
# plt.bar(x, above_threshold, 0.35, color="r", yerr=y_r, bottom=below_threshold)
ax3.bar(x, logOddsSame, yerr=y_r, color=['g', 'r', 'g', 'g', 'r', 'r', 'g', 'r'], alpha=0.8, align='center')
# plt.xticks(range(len(logOdds)), [str(i) for i in range(1,len(logOdds)+1)])
ax3.set_ylabel('Same chages')
ax3.set_xlabel('Step distance from first charged residue')
# ax3.set_title('Log odds ratio for same pairs')
ax3.axhline(xmax=8, color='black')
ax3.grid(True)
# plt.axis('off')
ax1.set_xticks([])
ax2.set_xticks([])

for ax in [ax1, ax2, ax3]:
    ax.grid(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.patches[0].set_facecolor("#b4a06d")
    ax.patches[1].set_facecolor("#d3605b")
    ax.patches[2].set_facecolor("#9bbe77")
    ax.patches[3].set_facecolor("#8dcc7c")
    ax.patches[4].set_facecolor("#c97961")
    ax.patches[5].set_facecolor("#bf8e67")
    ax.patches[6].set_facecolor("#7dd980")
    ax.patches[7].set_facecolor("#a8b072")

name = args.input_file.split('/')[-1].split('.')[0]
plt.savefig('images/' + name + '_logodds.png')

# ###### Build above and below graph
# allthreshold = chargedbaseline**2
# values = np.array(chargedpairs)
# x = range(1, len(values)+1)
# above_threshold = np.maximum(values - allthreshold, 0)
# below_threshold = np.minimum(values, allthreshold)
# # fig, ax = plt.subplots()
# axarr[0, 0].bar(x, below_threshold, 0.35, color="g")
# axarr[0, 0].bar(x, above_threshold, 0.35, color="r", bottom=below_threshold)
# axarr[0, 0].set_ylim([0, 0.0065])
# axarr[0, 0].plot([0.,7],[allthreshold, allthreshold], 'k--')
# axarr[0, 0].set_title("Above or below expected pair, all charged")
# # fig.savefig("pairsAllCharged.png")


# # print("Same charge hit: ", samehitcounter)
# # print("Opposite charge hit: ", hitcounter)
# # print("Non charged: :", misscounter)
# totalresidues = len(totalmemstring)
# aaRatios = {}
# normsamehitcounter = [i/sum(samehitcounter) for i in samehitcounter]
# normhitcounter = [i/sum(hitcounter) for i in hitcounter]
# normmisscounter = [i/sum(misscounter) for i in misscounter]
# for c in 'ACDEFGHIKLMNPQRSTVWY':
#     aaRatios[c] = totalmemstring.count(c)/totalresidues
# chargedbaseline = aaRatios['K']+aaRatios['R']+aaRatios['D']+aaRatios['E']+aaRatios['H']
# negbaseline = aaRatios['D']+aaRatios['E']
# posbaseline = aaRatios['K']+aaRatios['R']+aaRatios['H']
# # print(pairs)
# print("Charged baseline: ", aaRatios['K']+aaRatios['R']+aaRatios['D']+aaRatios['E']+aaRatios['H'])
# print(hitcounter[2]/misscounter[2])
# print(hitcounter[4]/misscounter[4])
# print("Odds ratio 5: ", ((hitcounter[2]/misscounter[2])/(hitcounter[4]/misscounter[4])))
# print("Log odds 5: ", math.log((hitcounter[2]/misscounter[2])/(hitcounter[4]/misscounter[4])))
# print("Odds ratio 6: ", ((hitcounter[2]/misscounter[2])/(hitcounter[5]/misscounter[5])))
# print("Log odds 6: ", math.log((hitcounter[2]/misscounter[2])/(hitcounter[5]/misscounter[5])))
# chargedpairs= []
# diffpairs = []
# samepairs = []
# posnegpairs = []
# negpospairs = []
# pospospairs = []
# negnegpairs = []
# for dist in range(7):
#     tempcharge = 0
#     tempdiff = 0
#     tempsame = 0
#     tempposneg = 0
#     tempnegpos = 0
#     temppospos = 0
#     tempnegneg = 0
#     for pair in pairs[dist]:
#         first, second = pair
#         if first in chargedplus and second in chargedplus:
#             tempcharge += 1
#         if (first in positiveplus and second in negative) or\
#             (first in negative and second in positiveplus):
#                 tempdiff += 1
#         if (first in positiveplus and second in positiveplus) or\
#             (first in negative and second in negative):
#                 tempsame += 1
#         if (first in positiveplus and second in negative):
#                 tempposneg += 1
#         if (first in negative and second in positiveplus):
#                 tempnegpos += 1
#         if (first in positiveplus and second in positiveplus):
#                 temppospos += 1
#         if (first in negative and second in negative):
#                 tempnegneg += 1
#     chargedpairs.append(tempcharge/len(pairs[dist]))
#     diffpairs.append(tempdiff/len(pairs[dist]))
#     samepairs.append(tempsame/len(pairs[dist]))
#     posnegpairs.append(tempposneg/len(pairs[dist]))
#     negpospairs.append(tempnegpos/len(pairs[dist]))
#     pospospairs.append(temppospos/len(pairs[dist]))
#     negnegpairs.append(tempnegneg/len(pairs[dist]))
# # for i, pair in enumerate(chargedpairs):
# #     print("Hitpair position {}: {}".format(i+1, pair*100))
# # print("Expected hitpar: ", chargedbaseline**2*100)
# # for i in range(len(pairs)):
# #     print("Length of pair {}: {}".format(i+1, len(pairs[i])))
# fig, axarr = plt.subplots(3, 3, figsize=(25,15))
# ###### Build above and below graph
# allthreshold = chargedbaseline**2
# values = np.array(chargedpairs)
# x = range(1, len(values)+1)
# above_threshold = np.maximum(values - allthreshold, 0)
# below_threshold = np.minimum(values, allthreshold)
# # fig, ax = plt.subplots()
# axarr[0, 0].bar(x, below_threshold, 0.35, color="g")
# axarr[0, 0].bar(x, above_threshold, 0.35, color="r", bottom=below_threshold)
# axarr[0, 0].set_ylim([0, 0.0065])
# axarr[0, 0].plot([0.,7],[allthreshold, allthreshold], 'k--')
# axarr[0, 0].set_title("Above or below expected pair, all charged")
# # fig.savefig("pairsAllCharged.png")
# # plt.clf()
# #######
# samethreshold = posbaseline**2 + negbaseline**2
# values = np.array(samepairs)
# x = range(1, len(values)+1)
# above_threshold = np.maximum(values - samethreshold, 0)
# below_threshold = np.minimum(values, samethreshold)
# # fig, ax = plt.subplots()
# axarr[0, 1].bar(x, below_threshold, 0.35, color="g")
# axarr[0, 1].bar(x, above_threshold, 0.35, color="r", bottom=below_threshold)
# axarr[0, 1].set_ylim([0, 0.0045])
# axarr[0, 1].plot([0.,7],[samethreshold, samethreshold], 'k--')
# axarr[0, 1].set_title("Above or below expected pair, same charged")
# # fig.savefig("pairsSameCharge.png")
# # plt.clf()
# ##############
# diffthreshold = posbaseline*negbaseline + negbaseline*posbaseline
# values = np.array(diffpairs)
# x = range(1, len(values)+1)
# above_threshold = np.maximum(values - diffthreshold, 0)
# below_threshold = np.minimum(values, diffthreshold)
# # fig, ax = plt.subplots()
# axarr[0, 2].bar(x, below_threshold, 0.35, color="g")
# axarr[0, 2].bar(x, above_threshold, 0.35, color="r", bottom=below_threshold)
# axarr[0, 2].set_ylim([0, 0.0045])
# axarr[0, 2].plot([0.,7],[diffthreshold, diffthreshold], 'k--')
# axarr[0, 2].set_title("Above or below expected pair, different charged")
# # fig.savefig("pairsDiffCharge.png")
# # plt.clf()
# ##################
# posnegthreshold = posbaseline*negbaseline
# values = np.array(posnegpairs)
# x = range(1, len(values)+1)
# above_threshold = np.maximum(values - posnegthreshold, 0)
# below_threshold = np.minimum(values, posnegthreshold)
# # fig, ax = plt.subplots()
# axarr[1, 0].bar(x, below_threshold, 0.35, color="g")
# axarr[1, 0].bar(x, above_threshold, 0.35, color="r", bottom=below_threshold)
# axarr[1, 0].set_ylim([0, 0.0025])
# axarr[1, 0].plot([0.,7],[posnegthreshold, posnegthreshold], 'k--')
# axarr[1, 0].set_title("Above or below expected pair, pos to neg ")
# # fig.savefig("pairsPosNeg.png")
# # plt.clf()
# ##################
# negposthreshold = posbaseline*negbaseline
# values = np.array(negpospairs)
# x = range(1, len(values)+1)
# above_threshold = np.maximum(values - negposthreshold, 0)
# below_threshold = np.minimum(values, negposthreshold)
# # fig, ax = plt.subplots()
# axarr[1, 1].bar(x, below_threshold, 0.35, color="g")
# axarr[1, 1].bar(x, above_threshold, 0.35, color="r", bottom=below_threshold)
# axarr[1, 1].set_ylim([0, 0.0025])
# axarr[1, 1].plot([0.,7],[negposthreshold, negposthreshold], 'k--')
# axarr[1, 1].set_title("Above or below expected pair, neg to pos")
# # fig.savefig("pairsNegPos.png")
# # plt.clf()
# ##################
# posposthreshold = posbaseline**2
# values = np.array(pospospairs)
# x = range(1, len(values)+1)
# above_threshold = np.maximum(values - posposthreshold, 0)
# below_threshold = np.minimum(values, posposthreshold)
# # fig, ax = plt.subplots()
# axarr[1, 2].bar(x, below_threshold, 0.35, color="g")
# axarr[1, 2].bar(x, above_threshold, 0.35, color="r", bottom=below_threshold)
# axarr[1, 2].set_ylim([0, 0.0025])
# axarr[1, 2].plot([0.,7],[posposthreshold, posposthreshold], 'k--')
# axarr[1, 2].set_title("Above or below expected pair, pos to pos")
# # fig.savefig("pairsPosPos.png")
# # plt.clf()
# ##################
# negnegthreshold = negbaseline**2
# values = np.array(negnegpairs)
# x = range(1, len(values)+1)
# above_threshold = np.maximum(values - negnegthreshold, 0)
# below_threshold = np.minimum(values, negnegthreshold)
# # fig, ax = plt.subplots()
# axarr[2, 0].bar(x, below_threshold, 0.35, color="g")
# axarr[2, 0].bar(x, above_threshold, 0.35, color="r", bottom=below_threshold)
# axarr[2, 0].set_ylim([0, 0.0025])
# axarr[2, 0].plot([0.,7],[negnegthreshold, negnegthreshold], 'k--')
# axarr[2, 0].set_title("Above or below expected pair, neg to neg")
# # fig.savefig("pairsNegNeg.png")
# # plt.clf()
# ##################
# fig.subplots_adjust(hspace=0.3)
# fig.savefig("unirefpairs.png")
# plt.clf()
# # plt.bar(range(1,len(normsamehitcounter)+1),normsamehitcounter)
# # plt.title("Same hit counter")
# # plt.savefig('normsamehitcounter.png')
# # plt.clf()
# # plt.bar(range(1,len(normhitcounter)+1),normhitcounter)
# # plt.title("Hit counter")
# # plt.savefig('normhitcounter.png')
# # plt.clf()
# # plt.bar(range(1,len(normmisscounter)+1),normmisscounter)
# # plt.title("Miss counter")
# # # plt.savefig('normmisscounter.png')
