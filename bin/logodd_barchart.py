#!/usr/bin/env python3
"""Membrane charge caluclations."""
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import math
import pickle
import random
import collections
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("input_file", type=str, help="Input membrane pickle")

args = parser.parse_args()

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
helicies = pickle.load(open(args.input_file, 'rb'))
# helicies = pickle.load(open('helicies.pickle','rb'))
mem_length = 17
totalmemstring = ''
totalmembranes = 0
samehitcounter = [0, 0, 0, 0, 0, 0, 0]
hitcounter = [0, 0, 0, 0, 0, 0, 0]
misscounter = [0, 0, 0, 0, 0, 0, 0]
totalcharges = []
pairs = [[], [], [], [], [], [], []]
aas = 'ACDEFGHIKLMNPQRSTVWY'
aaHits = np.zeros([7, 20, 20])
aaPairs = np.zeros(7)
# Dim 1, first and second in pair
aaCount = np.zeros([2, 7, 20])
memLens = []
prolinMems = []                                                                                                                                                                                                                                                                         
pPositions = []
for key, membranes in helicies.items():
    for startid, mem in membranes:

        # print(mem)
        # mem = mems[1]
        # Uncomment next line for only mid mem
        if len(mem) < mem_length:
            continue
        midMem = mem[5:-5]
        # midMem = mem
        memLen = len(midMem)
        # The follow two lines for Globular with mem length > 15
        # if memLen < 15:
        #     continue

#         # middle = memLen//2
#         # for pos, aa in enumerate(midMem):
#         #     if aa == 'P':
#         #         pPositions.append(pos-middle)
#         memLens.append(memLen)
#         if 'P' in midMem:
#             prolinMems.append(memLen)
#         # print(memLen)
# f, axarr = plt.subplots(1, 1, figsize=(12, 8))
# plt.hist(memLens, alpha=0.5, bins=40, label='Membrane lengths')
# plt.hist(prolinMems, alpha=0.5, bins=40, label='Contains Prolin')
# plt.legend(loc='upper right')
# # ax.title.set_text("Amino Acid Count MidMem")
# # ax = sns.heatmap(normdf, ax=axarr[1], annot=True, fmt='.1%')
# # ax.title.set_text("Amino Acid Count MidMem Norm")
# plt.tight_layout()
# f.savefig('images/histogramOfProlinMemLensTM.png')
        for place, aa in enumerate(midMem):
            # print(aa)
            # aaCount[aas.index(aa)] += 1 # Old way of counting AA
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

logOdds = np.zeros([7])
ci = []
logOdds_opp = np.zeros([7])
ci_opp = []
logOdds_same = np.zeros([7])
ci_same = []
# print([aas.index(c) for c in chargedplus])
aaIndex = [aas.index(c) for c in positiveplus]
# print(sum(aaHits[0][chargesPlusIndex][chargesPlusIndex]))
# for ind, aa in enumerate(aasOri):
#    print(aa, aaCount[ind])
for i in range(7):
    # for first in range(20):
    #     for second in range(20):
    a = 0
    b_temp = aaPairs[i]
    c_first = 0
    c_second = 0
    c = 0
    d_temp = aaPairs[i]**2

    a_opp = 0
    b_temp_opp = 0
    c_first_opp = 0
    c_second_opp = 0
    d_temp_opp = 0

    a_same = 0
    b_temp_same = 0
    c_first_same = 0
    c_second_same = 0
    d_temp_same = 0
    
    for first in [aas.index(c) for c in charged]:
        for second in [aas.index(c) for c in charged]:
            a += aaHits[i][first][second]                                                                                                                                                                                                                                               
            # b_temp += aaPairs[i]
            # b = aaPairs[i] - a
            # print(a, b, a/b)
            c_first += aaCount[0][i][first]
            c_second += aaCount[1][i][second]
            c += aaCount[0][i][first] * aaCount[1][i][second]
            # Need square to handle pairs, first and second
            # d = sum(aaCount)**2 # Old version
            # d is aaPairs for this gap times 2 for each of the two
            # amino acids as each pair contains two AAs.
            # d_temp += aaPairs[i]
            # d = aaPairs[i]**2 - c
            # print(a, b, c, d)
            # if (first in positive and second in negative) or (first in negative and second in positive):
            #     a_opp += aaHits[i][first][second]                                                                                                                                                                                                                                  
            #     b_temp_opp += aaPairs[i]
            #     c_first_opp += aaCount[0][i][first]
            #     c_second_opp += aaCount[1][i][second]
            #     d_temp_opp += aaPairs[i]
            # if (first in positive and second in positive) or (first in negative and second in negative):
            #     a_same += aaHits[i][first][second]                                                                                                                                                                                                                                  
            #     b_temp_same += aaPairs[i]
            #     c_first_same += aaCount[0][i][first]
            #     c_second_same += aaCount[1][i][second]
            #     d_temp_same += aaPairs[i]

    b = (b_temp - a)
    # c = c_first * c_second
    d = (d_temp-c)
    odds = (a / b) / (c / d)

    # b_opp = b_temp_opp - a_opp
    # c_opp = c_first_opp*c_second_opp
    # d_opp = d_temp_opp**2-c_opp
    # odds_opp = (a_opp / b_opp) / (c_opp / d_opp)

    # b_same = b_temp_same - a_same
    # c_same = c_first_same*c_second_same
    # d_same = d_temp_same**2-c_same
    # odds_same = (a_same / b_same) / (c_same / d_same)
#     # print("Odds ratio for gap {}: {:.2f}".format(i+1, (odds)))
# print("Log odds ratio for gap {}: {:.2f}".format(i+1, math.log(odds)))
# logOdds.append(math.log(odds))
    if odds == 0.0:
        logOddsValue = -math.inf
        print("Odds 0 at gap: ", i + 1)
        print("A-D: ", a, b, c, d)
        print("First and second AA: ", aasOri[first], aasOri[second])
        # print(aasOri[first], aasOri[second], a, b, c, d, odds)
    else:
        logOddsValue = math.log(odds)
    logOdds[i] = logOddsValue
    ci.append((math.log(odds)+1.96*math.sqrt(1/a+1/b+1/c+1/d), math.log(odds)-1.96*math.sqrt(1/a+1/b+1/c+1/d)))
print(logOdds)
print(ci)
# print(aaCount)
# print(aaPairs)
# print(aaHits)
# sys.exit()
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
# totalhits = np.array(hitcounter) + np.array(samehitcounter)
# print("Samehitcounter: {}".format(samehitcounter))
# print("Hit counter:    {}".format(hitcounter))
# print("Total hits:     {}".format(totalhits))
# print("Misscounter:    {}".format(misscounter))
# print("Baseline:       {}".format(baselinehits))
# print("Total charges:  {}".format(sum(totalcharges)))
# freq = collections.Counter(totalcharges)
# print(freq)
# print("Total membranes: {}".format(totalmembranes))
# logOdds = []
# ci = []
# for i in range(7):
#     a = totalhits[i]                # Number of pairs, there is a charge
#     b = sum(totalcharges) - a       # Number of charges that are not pair for this distance
#     c = baselinehits[i]             # How many baseline pairs for this distance?
#     d = sum(totalcharges) - c       # Number of non-hits for baseline
#     odds = (a/b)/(c/d)
#     # print("Odds ratio for gap {}: {:.2f}".format(i+1, (odds)))
#     print("Log odds ratio for gap {}: {:.2f}".format(i+1, math.log(odds)))
#     logOdds.append(math.log(odds))
#     print("Upper 95% CI: {:.2f}".format((math.log(odds)+1.96*math.sqrt(1/a+1/b+1/c+1/d))))
#     print("Lower 95% CI: {:.2f}".format((math.log(odds)-1.96*math.sqrt(1/a+1/b+1/c+1/d))))
#     ci.append((math.log(odds)+1.96*math.sqrt(1/a+1/b+1/c+1/d), math.log(odds)-1.96*math.sqrt(1/a+1/b+1/c+1/d)))
# print(logOdds)
# print(ci)

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
plt.bar(x, logOdds, yerr=y_r, color=['g', 'r', 'g', 'g', 'r', 'r', 'r'], alpha=0.8, align='center')
# plt.xticks(range(len(logOdds)), [str(i) for i in range(1,len(logOdds)+1)])
plt.ylabel('Log odds ratio')
plt.xlabel('Gap distance from first charged residue')
plt.title('Log odds ratio from a charged residue')
plt.axhline(xmax=8, color='black')
plt.grid(True)
plt.show()
# plt.savefig("logoddsthreshold2.png")

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
