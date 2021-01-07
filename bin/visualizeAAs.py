#!/usr/bin/env python3
"""Membrane charge caluclations."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import seaborn as sns
import pickle
import argparse
import scipy
import sys
# import random
# import collections
plt.rcParams.update({'font.size': 14})
parser = argparse.ArgumentParser()

parser.add_argument("charges", type=str, help="In charge pickle")
parser.add_argument("out", type=str, help="Out image name")
parser.add_argument("title", type=str, help="Image title")

args = parser.parse_args()

data = {}
positive = 'KR'
negative = 'DE'
charged = 'KRDE'
chargedplus = 'KRDEH'
positiveplus = 'RHK'
# aas = 'ACDEFGHIKLMNPQRSTVWY'
aasOri = 'ACDEFGHIKLMNPQRSTVWY'
# Regular AA group, [nonpolar, polar, positive, negative]
aas = 'GAVCPLIMWFSTYNQHKRDE'
# Engleman order
aas = 'FMILVCWATGSPYHQNEKDR'
# in_pickle = 'data/pdbtm/pdbtm_redundant_alpha_struct_scop_reduced_culled_charge.pickle'
in_pickle = args.charges
# out_image = in_pickle.split('/')[-1].split('.')[0] + "_AA_vis.svg"
out_image = args.out
aaCount, aaPairs, aaHits = pickle.load(open(in_pickle, 'rb'))
    # open('data/chargeDataGlobularFromPDBRed50TrimmedLen15v2.pickle',
    #      'rb'))
log_odds_file = args.out + ".txt"
log_odds_text_list = ["Separation\tFirst\tSecond\tlogOdds\terror\tp"]
logOdds = np.zeros([10, 20, 20])
ci = []
# print([aas.index(c) for c in chargedplus])
aaIndex = [aas.index(c) for c in positiveplus]
# print(sum(aaHits[0][chargesPlusIndex][chargesPlusIndex]))
# for ind, aa in enumerate(aasOri):
#    print(aa, aaCount[ind])

mask = np.zeros([10, 20, 20])
for i in range(8):
    for first in range(20):
        for second in range(20):
            a = aaHits[i][first][second]
            b = aaPairs[i]  #  - a  # Not odds
            # print(a, b, a/b)
            c = aaCount[0][i][first] * aaCount[1][i][second]
            # Need square to handle pairs, first and second
            # d = sum(aaCount)**2 # Old version
            # d is aaPairs for this step times 2 for each of the two
            # amino acids as each pair contains two AAs.
            d = aaPairs[i]**2  # - c  # same as above, not odds
            # print(a, b, c, d)
            odds = (a / b) / (c / d)
            if odds == 0.0:
                logOddsValue = -math.inf
                print("Odds 0 at step: ", i + 1)
                print("A-D: ", a, b, c, d)
                print("First and second AA: ", aasOri[first], aasOri[second])
                # print(aasOri[first], aasOri[second], a, b, c, d, odds)
            else:
                logOddsValue = math.log(odds)
            logOdds[i][aas.index(aasOri[first])][aas.index(aasOri[second])]\
                = logOddsValue

            SE = math.sqrt(1/a+1/b+1/c+1/d)
            z = abs(logOddsValue/SE)  ## Two sided
            # pm = math.exp(-0.717*z-0.416*z**2)*20*20*8*2 ## two sided
            # p = (1-scipy.stats.norm.cdf(abs(z)))*2
            psf = (scipy.stats.norm.sf(abs(z)))*2*20*20*8  ## Make simple multiple hypothesis correction

            # if i == 0:
            mask[i][aas.index(aasOri[first])][aas.index(aasOri[second])]= 1 if psf > 0.05 else 0
            # print("{}-{} {}: logOdds: {:.3f} SE: {:.3f} p: {:.2e}".format(aasOri[first], aasOri[second],i+1, logOddsValue, SE, psf))
            log_odds_text_list.append("{}\t{}\t{}\t{:.3f}\t{:.3f}\t{:.2e}".format(i+1, aasOri[first], aasOri[second],logOddsValue, SE, psf))
            ci.append((math.log(odds)+1.96*math.sqrt(1/a+1/b+1/c+1/d), math.log(odds)-1.96*math.sqrt(1/a+1/b+1/c+1/d), psf))
            # if psf < 0.05:
            #     print("{}-{} {}: logOdds: {:.3f} SE: {:.3f} p: {:.2e}".format(aas[first], aas[second],i+1, logOddsValue, SE, psf))
#     # print("Odds ratio for step {}: {:.2f}".format(i+1, (odds)))
    # print("Log odds ratio for step {}: {:.2f}".format(i+1, math.log(odds)))
    # logOdds.append(math.log(odds))
    # print("Upper 95% CI: {:.2f}".format((math.log(odds)
    # +1.96*math.sqrt(1/a+1/b+1/c+1/d))))
    # print("Lower 95% CI: {:.2f}".format((math.log(odds)
    # -1.96*math.sqrt(1/a+1/b+1/c+1/d))))
    # ci.append((math.log(odds)+1.96*math.sqrt(1/a+1/b+1/c+1/d),
    # math.log(odds)-1.96*math.sqrt(1/a+1/b+1/c+1/d)))
# print(statsData[7][18][13])
# f, axarr = plt.subplots(3, 3, figsize=(16, 14))
prefix_dot = out_image.rfind('.')
out_log_odds = out_image[:prefix_dot] + "_logOdds.txt"

with open(out_log_odds, 'w') as log_out_handle:
    log_out_handle.write('\n'.join(log_odds_text_list))

for o in range(2):
    f, axarr = plt.subplots(2, 4, figsize=(25, 12))
# axarr[-1, -1].axis('off')
# axarr[-1, -2].axis('off')
# plt.suptitle("TMs alpha helices, trimmed, pdb50, len > 15, V2")
    plt.suptitle(args.title)
    if o:
        out_image = out_image[:prefix_dot] + "_mask" + out_image[prefix_dot]
    for i in range(8):
        # print(i)
        # if i == 3 or i==7:
        #     show_cbar=True
        # else:
        #     show_cbar=False
        show_cbar = False
        df = pd.DataFrame(logOdds[i],
                          index=[c for c in aas],
                          columns=[c for c in aas])
        if o:
            ax = sns.heatmap(df,
                             ax=axarr[i // 4, i % 4],
                             vmin=-2,
                             vmax=2,
                             cmap="coolwarm",
                             center=0, cbar=show_cbar,mask=mask[i])
        else:
            ax = sns.heatmap(df,
                             ax=axarr[i // 4, i % 4],
                             vmin=-2,
                             vmax=2,
                             cmap="coolwarm",
                             center=0, cbar=show_cbar,)
        ax.title.set_text("Step " + str(i + 1))
        for x_tick in ax.get_xticklabels():
            if x_tick.get_text() in charged:
                x_tick.set_weight("bold")
                x_tick.set_fontsize(14)
        for y_tick in ax.get_yticklabels():
            if y_tick.get_text() in charged:
                y_tick.set_weight("bold")
                y_tick.set_fontsize(14)
    plt.tight_layout(rect=[0,0,24/25,1])
# plt.subplots_adjust(top=0.93)
    cb_ax = f.add_axes([24/25, 0.044, 0.01, 0.9])
    if o:
        sns.heatmap(df,
                    ax=ax,
                    vmin=-2,
                    vmax=2,
                    cmap="coolwarm",
                    center=0, cbar=True, mask=mask[i], cbar_ax=cb_ax)
    else:
        sns.heatmap(df,
                    ax=ax,
                    vmin=-2,
                    vmax=2,
                    cmap="coolwarm",
                    center=0, cbar=True, cbar_ax=cb_ax)
    f.savefig(out_image)
    plt.clf()
# f.savefig('images/'
#           + 'logOddsGlobularTrimmed'
#           + 'HeatmapGroupOrderFromPDBRed50Len15v2.svg')
# # print(logOdds)
# # print(ci)
# y_r = [logOdds[i] - ci[i][0] for i in range(len(ci))]
# threshold = 0
# values = np.array(logOdds)
# x = range(1, len(values)+1)
# # plt.bar(range(len(logOdds)), logOdds, yerr=y_r, alpha=0.2, align='center')
# # x = range(1, len(values)+1)
# above_threshold = np.maximum(values - threshold, 0)
# below_threshold = np.minimum(values, threshold)
# # plt.bar(x, below_threshold, 0.35, color="g", yerr=y_r)
# # plt.bar(x, above_threshold, 0.35, color="r",
#   yerr=y_r, bottom=below_threshold)
# plt.bar(x, logOdds, yerr=y_r,
#   color=['g' if c > 0 else 'r' for c in logOdds], alpha=0.8, align='center')
# # plt.xticks(range(len(logOdds)), [str(i) for i in range(1,len(logOdds)+1)])
# plt.ylabel('Log odds ratio')
# plt.xlabel('Step distance from first charged residue')
# plt.title('Log odds ratio from a charged residue')
# plt.axhline(xmax=8, color='black')
# plt.grid(True)
# plt.savefig("positiveplus.svg")
# # ###### Build above and below graph
# # allthreshold = chargedbaseline**2
# # values = np.array(chargedpairs)
# # x = range(1, len(values)+1)
# # above_threshold = np.maximum(values - allthreshold, 0)
# # below_threshold = np.minimum(values, allthreshold)
# # # fig, ax = plt.subplots()
# # axarr[0, 0].bar(x, below_threshold, 0.35, color="g")
# # axarr[0, 0].bar(x, above_threshold, 0.35,
#   color="r", bottom=below_threshold)
# # axarr[0, 0].set_ylim([0, 0.0065])
# # axarr[0, 0].plot([0.,7],[allthreshold, allthreshold], 'k--')
# # axarr[0, 0].set_title("Above or below expected pair, all charged")
# # # fig.savefig("pairsAllCharged.svg")
# # # print("Same charge hit: ", samehitcounter)
# # # print("Opposite charge hit: ", hitcounter)
# # # print("Non charged: :", misscounter)
# # totalresidues = len(totalmemstring)
# # aaRatios = {}
# # normsamehitcounter = [i/sum(samehitcounter) for i in samehitcounter]
# # normhitcounter = [i/sum(hitcounter) for i in hitcounter]
# # normmisscounter = [i/sum(misscounter) for i in misscounter]
# # for c in 'ACDEFGHIKLMNPQRSTVWY':
# #     aaRatios[c] = totalmemstring.count(c)/totalresidues
# # chargedbaseline = aaRatios['K']+aaRatios['R']
#   +aaRatios['D']+aaRatios['E']+aaRatios['H']
# # negbaseline = aaRatios['D']+aaRatios['E']
# # posbaseline = aaRatios['K']+aaRatios['R']+aaRatios['H']
# # # print(pairs)
# # print("Charged baseline: ", aaRatios['K']+aaRatios['R']
#   +aaRatios['D']+aaRatios['E']+aaRatios['H'])
# # print(hitcounter[2]/misscounter[2])
# # print(hitcounter[4]/misscounter[4])
# # print("Odds ratio 5: ",
#   ((hitcounter[2]/misscounter[2])/(hitcounter[4]/misscounter[4])))
# # print("Log odds 5: ",
#    math.log((hitcounter[2]/misscounter[2])/(hitcounter[4]/misscounter[4])))
# # print("Odds ratio 6: ",
#   ((hitcounter[2]/misscounter[2])/(hitcounter[5]/misscounter[5])))
# # print("Log odds 6: ",
#   math.log((hitcounter[2]/misscounter[2])/(hitcounter[5]/misscounter[5])))
# # chargedpairs= []
# # diffpairs = []
# # samepairs = []
# # posnegpairs = []
# # negpospairs = []
# # pospospairs = []
# # negnegpairs = []
# # for dist in range(7):
# #     tempcharge = 0
# #     tempdiff = 0
# #     tempsame = 0
# #     tempposneg = 0
# #     tempnegpos = 0
# #     temppospos = 0
# #     tempnegneg = 0
# #     for pair in pairs[dist]:
# #         first, second = pair
# #         if first in chargedplus and second in chargedplus:
# #             tempcharge += 1
# #         if (first in positiveplus and second in negative) or\
# #             (first in negative and second in positiveplus):
# #                 tempdiff += 1
# #         if (first in positiveplus and second in positiveplus) or\
# #             (first in negative and second in negative):
# #                 tempsame += 1
# #         if (first in positiveplus and second in negative):
# #                 tempposneg += 1
# #         if (first in negative and second in positiveplus):
# #                 tempnegpos += 1
# #         if (first in positiveplus and second in positiveplus):
# #                 temppospos += 1
# #         if (first in negative and second in negative):
# #                 tempnegneg += 1
# #     chargedpairs.append(tempcharge/len(pairs[dist]))
# #     diffpairs.append(tempdiff/len(pairs[dist]))
# #     samepairs.append(tempsame/len(pairs[dist]))
# #     posnegpairs.append(tempposneg/len(pairs[dist]))
# #     negpospairs.append(tempnegpos/len(pairs[dist]))
# #     pospospairs.append(temppospos/len(pairs[dist]))
# #     negnegpairs.append(tempnegneg/len(pairs[dist]))
# # # for i, pair in enumerate(chargedpairs):
# # #     print("Hitpair position {}: {}".format(i+1, pair*100))
# # # print("Expected hitpar: ", chargedbaseline**2*100)
# # # for i in range(len(pairs)):
# # #     print("Length of pair {}: {}".format(i+1, len(pairs[i])))
# # fig, axarr = plt.subplots(3, 3, figsize=(25,15))
# # ###### Build above and below graph
# # allthreshold = chargedbaseline**2
# # values = np.array(chargedpairs)
# # x = range(1, len(values)+1)
# # above_threshold = np.maximum(values - allthreshold, 0)
# # below_threshold = np.minimum(values, allthreshold)
# # # fig, ax = plt.subplots()
# # axarr[0, 0].bar(x, below_threshold, 0.35, color="g")
# # axarr[0, 0].bar(x, above_threshold, 0.35, color="r",
#   bottom=below_threshold)
# # axarr[0, 0].set_ylim([0, 0.0065])
# # axarr[0, 0].plot([0.,7],[allthreshold, allthreshold], 'k--')
# # axarr[0, 0].set_title("Above or below expected pair, all charged")
# # # fig.savefig("pairsAllCharged.svg")
# # # plt.clf()
# # #######
# # samethreshold = posbaseline**2 + negbaseline**2
# # values = np.array(samepairs)
# # x = range(1, len(values)+1)
# # above_threshold = np.maximum(values - samethreshold, 0)
# # below_threshold = np.minimum(values, samethreshold)
# # # fig, ax = plt.subplots()
# # axarr[0, 1].bar(x, below_threshold, 0.35, color="g")
# # axarr[0, 1].bar(x, above_threshold, 0.35, color="r",
#   bottom=below_threshold)
# # axarr[0, 1].set_ylim([0, 0.0045])
# # axarr[0, 1].plot([0.,7],[samethreshold, samethreshold], 'k--')
# # axarr[0, 1].set_title("Above or below expected pair, same charged")
# # # fig.savefig("pairsSameCharge.svg")
# # # plt.clf()
# # ##############
# # diffthreshold = posbaseline*negbaseline + negbaseline*posbaseline
# # values = np.array(diffpairs)
# # x = range(1, len(values)+1)
# # above_threshold = np.maximum(values - diffthreshold, 0)
# # below_threshold = np.minimum(values, diffthreshold)
# # # fig, ax = plt.subplots()
# # axarr[0, 2].bar(x, below_threshold, 0.35, color="g")
# # axarr[0, 2].bar(x, above_threshold, 0.35, color="r",
#   bottom=below_threshold)
# # axarr[0, 2].set_ylim([0, 0.0045])
# # axarr[0, 2].plot([0.,7],[diffthreshold, diffthreshold], 'k--')
# # axarr[0, 2].set_title("Above or below expected pair, different charged")
# # # fig.savefig("pairsDiffCharge.svg")
# # # plt.clf()
# # ##################
# # posnegthreshold = posbaseline*negbaseline
# # values = np.array(posnegpairs)
# # x = range(1, len(values)+1)
# # above_threshold = np.maximum(values - posnegthreshold, 0)
# # below_threshold = np.minimum(values, posnegthreshold)
# # # fig, ax = plt.subplots()
# # axarr[1, 0].bar(x, below_threshold, 0.35, color="g")
# # axarr[1, 0].bar(x, above_threshold, 0.35, color="r",
#   bottom=below_threshold)
# # axarr[1, 0].set_ylim([0, 0.0025])
# # axarr[1, 0].plot([0.,7],[posnegthreshold, posnegthreshold], 'k--')
# # axarr[1, 0].set_title("Above or below expected pair, pos to neg ")
# # # fig.savefig("pairsPosNeg.svg")
# # # plt.clf()
# # ##################
# # negposthreshold = posbaseline*negbaseline
# # values = np.array(negpospairs)
# # x = range(1, len(values)+1)
# # above_threshold = np.maximum(values - negposthreshold, 0)
# # below_threshold = np.minimum(values, negposthreshold)
# # # fig, ax = plt.subplots()
# # axarr[1, 1].bar(x, below_threshold, 0.35, color="g")
# # axarr[1, 1].bar(x, above_threshold, 0.35, color="r",
#   bottom=below_threshold)
# # axarr[1, 1].set_ylim([0, 0.0025])
# # axarr[1, 1].plot([0.,7],[negposthreshold, negposthreshold], 'k--')
# # axarr[1, 1].set_title("Above or below expected pair, neg to pos")
# # # fig.savefig("pairsNegPos.svg")
# # # plt.clf()
# # ##################
# # posposthreshold = posbaseline**2
# # values = np.array(pospospairs)
# # x = range(1, len(values)+1)
# # above_threshold = np.maximum(values - posposthreshold, 0)
# # below_threshold = np.minimum(values, posposthreshold)
# # # fig, ax = plt.subplots()
# # axarr[1, 2].bar(x, below_threshold, 0.35, color="g")
# # axarr[1, 2].bar(x, above_threshold, 0.35, color="r",
#   bottom=below_threshold)
# # axarr[1, 2].set_ylim([0, 0.0025])
# # axarr[1, 2].plot([0.,7],[posposthreshold, posposthreshold], 'k--')
# # axarr[1, 2].set_title("Above or below expected pair, pos to pos")
# # # fig.savefig("pairsPosPos.svg")
# # # plt.clf()
# # ##################
# # negnegthreshold = negbaseline**2
# # values = np.array(negnegpairs)
# # x = range(1, len(values)+1)
# # above_threshold = np.maximum(values - negnegthreshold, 0)
# # below_threshold = np.minimum(values, negnegthreshold)
# # # fig, ax = plt.subplots()
# # axarr[2, 0].bar(x, below_threshold, 0.35, color="g")
# # axarr[2, 0].bar(x, above_threshold, 0.35, color="r",
#   bottom=below_threshold)
# # axarr[2, 0].set_ylim([0, 0.0025])
# # axarr[2, 0].plot([0.,7],[negnegthreshold, negnegthreshold], 'k--')
# # axarr[2, 0].set_title("Above or below expected pair, neg to neg")
# # # fig.savefig("pairsNegNeg.svg")
# # # plt.clf()
# # ##################
# # fig.subplots_adjust(hspace=0.3)
# # fig.savefig("unirefpairs.svg")
# # plt.clf()
# # # plt.bar(range(1,len(normsamehitcounter)+1),normsamehitcounter)
# # # plt.title("Same hit counter")
# # # plt.savefig('normsamehitcounter.svg')
# # # plt.clf()
# # # plt.bar(range(1,len(normhitcounter)+1),normhitcounter)
# # # plt.title("Hit counter")
# # # plt.savefig('normhitcounter.svg')
# # # plt.clf()
# # # plt.bar(range(1,len(normmisscounter)+1),normmisscounter)
# # # plt.title("Miss counter")
# # # # plt.savefig('normmisscounter.svg')
