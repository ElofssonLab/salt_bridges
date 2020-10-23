#!/usr/bin/env python3
"""Membrane charge caluclations."""
import matplotlib.pyplot as plt
import numpy as np
import cmocean
import pandas as pd
import math
import seaborn as sns
import pickle
import argparse
import sys
# import random
# import collections
parser = argparse.ArgumentParser()

parser.add_argument("glob", type=str, help="In glob charges pickle")
parser.add_argument("mems", type=str, help="In membrane charges pickle")
parser.add_argument("out", type=str, help="Out image name")
parser.add_argument("title", type=str, help="Image title")

args = parser.parse_args()

data = {}
positive = 'KR'
negative = 'DE'
charged = 'KRDE'
chargedplus = 'KRDEH'
positiveplus = 'RHK'
# keep = [c for c in 'HQNEKDR']
keep = [c for c in 'EKDR']
keep_len = len(keep)
# aas = 'ACDEFGHIKLMNPQRSTVWY'
aasOri = 'ACDEFGHIKLMNPQRSTVWY'
# Regular AA group, [nonpolar, polar, positive, negative]
# aas = 'GAVCPLIMWFSTYNQHKRDE'
# Hessa order
# aas = 'ILFVCMAWTYGSNHPQREKD'
# Engleman order
aas = 'FMILVCWATGSPYHQNEKDR'
prefix = ["Globular"] * 20 + ["Membranes"]*20
new_cols = []
for i, p in enumerate(prefix):
    new_cols.append(p + " " + aas[i%20])
# in_pickle = 'data/pdbtm/pdbtm_redundant_alpha_struct_scop_reduced_culled_charge.pickle'
in_glob_pickle = args.glob
in_mem_pickle = args.mems
# out_image = in_pickle.split('/')[-1].split('.')[0] + "_AA_vis.svg"
# out_image = args.out
glob_aaCount, glob_aaPairs, glob_aaHits = pickle.load(open(in_glob_pickle, 'rb'))
mem_aaCount, mem_aaPairs, mem_aaHits = pickle.load(open(in_mem_pickle, 'rb'))
    # open('data/chargeDataGlobularFromPDBRed50TrimmedLen15v2.pickle',
    #      'rb'))

globlogOdds = np.zeros([7, 20, 20])
memlogOdds = np.zeros([7, 20, 20])
ci = []
# print([aas.index(c) for c in chargedplus])
aaIndex = [aas.index(c) for c in positiveplus]
# print(sum(aaHits[0][chargesPlusIndex][chargesPlusIndex]))
# for ind, aa in enumerate(aasOri):
#    print(aa, aaCount[ind])
for i in range(7):
    for first in range(20):
        for second in range(20):
            a = glob_aaHits[i][first][second]
            b = glob_aaPairs[i] - a
            # print(a, b, a/b)
            c = glob_aaCount[0][i][first] * glob_aaCount[1][i][second]
            # Need square to handle pairs, first and second
            # d = sum(aaCount)**2 # Old version
            # d is aaPairs for this step times 2 for each of the two
            # amino acids as each pair contains two AAs.
            d = glob_aaPairs[i]**2 - c
            # print(a, b, c, d)
            odds = (a / b) / (c / d)
            if odds == 0.0:
                logOddsValue = -math.inf
                print("Odds 0 at step: ", i + 1)
                print("A-D: ", a, b, c, d)
                print("First and second AA: ", aasOri[first], aasOri[second])
                # print(aasOri[first], aasOri[second], a, b, c, d, odds)
            else:
                # std_err = math.sqrt(1/a + 1/b + 1/c + 1/d)
                logOddsValue = math.log(odds)
            globlogOdds[i][aas.index(aasOri[first])][aas.index(aasOri[second])]\
                = logOddsValue
for i in range(7):
    for first in range(20):
        for second in range(20):
            a = mem_aaHits[i][first][second]
            b = mem_aaPairs[i] - a
            # print(a, b, a/b)
            c = mem_aaCount[0][i][first] * mem_aaCount[1][i][second]
            # Need square to handle pairs, first and second
            # d = sum(aaCount)**2 # Old version
            # d is aaPairs for this step times 2 for each of the two
            # amino acids as each pair contains two AAs.
            d = mem_aaPairs[i]**2 - c
            odds = (a / b) / (c / d)
            if odds == 0.0:
                logOddsValue = -math.inf
                print("Odds 0 at step: ", i + 1)
                print("A-D: ", a, b, c, d)
                print("First and second AA: ", aasOri[first], aasOri[second])
                # print(aasOri[first], aasOri[second], a, b, c, d, odds)
            else:
                logOddsValue = math.log(odds)
                # std_err = math.sqrt(1/a + 1/b + 1/c + 1/d)
            memlogOdds[i][aas.index(aasOri[first])][aas.index(aasOri[second])]\
                = logOddsValue

# print(df)
# df.columns = new_cols
# used_types = ["Globular", "Membrane"]
# type_pal = sns.husl_palette(2, s=.45)
# type_lut = dict(zip(map(str, used_types), type_pal))
# types = df.columns.get_level_values("Type")
# print(types)
# used_types = ["Globular", "Membrane"]
# type_pal = sns.husl_palette(2, s=.45)
# type_lut = dict(zip(map(str, used_types), type_pal))
# type_colors = pd.Series(types, index=df.columns).map(type_lut)
# 
# df = df[:, 1:]
# methods = ['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward' ]
methods = ['complete', 'ward' ]
for method in methods:
# method='complete'

# plt.savefig(method + "_" + out_image)
# f, axarr = plt.subplots(3, 3, figsize=(16, 14))
# axarr[-1, -1].axis('off')
# axarr[-1, -2].axis('off')
# plt.suptitle("TMs alpha helices, trimmed, pdb50, len > 15, V2")
# plt.suptitle(args.title)
# for method in methods:
    df = pd.DataFrame()
    for i in range(7):
        df1 = pd.DataFrame(globlogOdds[i],
                          index=["Globular " + str(c) for c in aas],
                          columns=[c for c in aas])
        G_keep = ["Globular " + str(c) for c in keep]
        df1 = df1.loc[G_keep][keep]
        # df1.insert(loc=0,column= "Type",  value="Globular")
        df2 = pd.DataFrame(memlogOdds[i],
                          index=["Membranes " + str(c) for c in aas],
                          columns=[c for c in aas])
        # df2.insert( loc=0, column= "Type", value="Membranes")
        M_keep = ["Membranes " + str(c) for c in keep]
        df2 = df2.loc[M_keep][keep]
        df_temp = pd.concat([df1, df2])
        df_temp.columns = ["Step " + str(i+1) + "-" + c for c in df_temp.columns]
        df = pd.concat([df, df_temp], axis=1)
    df.insert(loc=0, column="Type", value=["Globular"]*keep_len + ["Membranes"]*keep_len)
# for i in range(7):
# i = 3
    plt.suptitle(args.title +  " sep " + str(i+1))
        # df1 = pd.DataFrame(globlogOdds[i],
        #                   index=["Globular " + str(c) for c in aas],
        #                   columns=[c for c in aas])
        # df1.insert(loc=0,column= "Type",  value="Globular")
        # df2 = pd.DataFrame(memlogOdds[i],
        #                   index=["Membranes " + str(c) for c in aas],
        #                   columns=[c for c in aas])
        # df2.insert( loc=0, column= "Type", value="Membranes")
        # df = pd.concat([df1, df2])
    types =df.pop("Type") 
# print(df.columns)
# sys.exit()
    step_list = [item for i in range(1, 8) for item in [i]*keep_len]
    steps = pd.Series(step_list, name="Steps", index=df.columns)
    step_c_list = [item for i in range(1, 8) for item in [(((100*i)%360)/360)]*keep_len]
    steps_c = pd.Series(step_c_list, name="C_cols", index=df.columns)
# type_list = pd.Series(["Globular"] * 20 + ["Membranes"]*20)
    lut = dict(zip(types.unique(), sns.husl_palette(2, s=.45)))
# lut_cols = dict(zip(gaps.unique(), sns.color_palette()))
# lut_cols = dict(zip(gaps.unique(), ["Crimson", "DarkSlateBlue", "Crimson", "Crimson", "DarkSlateBlue",  "DarkSlateBlue", "Crimson"] ))
    lut_cols = dict(zip(steps.unique(), ["#b4a06d", "#d3605b", "#9bbe77", "#8dcc7c", "#c97961",  "#bf8e67", "#7dd980"] ))
# cycl_cols = {v:cmocean.cm.phase(v) for v in gaps_c.unique()}
# cycl_cols = {v:plt.get_cmap("twilight")(v) for v in gaps_c.unique()}
    row_colors = types.map(lut)
# col_colors = gaps_c.map(cycl_cols)
    col_colors = steps.map(lut_cols)
    chart = sns.clustermap(df, row_cluster=True, col_cluster=True, figsize=(30, 10), cmap="coolwarm", center=0,
                   method=method, row_colors=row_colors, col_colors=col_colors, cbar_pos=None, vmin=-2, vmax=2,
                   linewidth=.75,
                   dendrogram_ratio=(0.05, 0.15),
                   colors_ratio=(0.005, 0.02),
                   )
    plt.setp(chart.ax_heatmap.xaxis.get_majorticklabels(), rotation=45, ha="right")
    plt.setp(chart.ax_row_colors.xaxis.get_majorticklabels(), rotation=45, ha="right")
    chart.savefig(args.out.replace(".svg", method + ".svg"))
##################################
# plt.figure(i+2)
# ax = sns.clustermap(df, ax=axarr[i // 3, i % 3], row_cluster=True, col_cluster=False, cmap="coolwarm", center=0, method=method, row_colors=row_colors)
    # df = pd.DataFrame(logOdds[i],
    #                   index=[c for c in aas],
    #                   columns=[c for c in aas])
    # ax = sns.heatmap(df,
    #                  ax=axarr[i // 3, i % 3],
    #                  vmin=-2,
    #                  vmax=2,
    #                  cmap="coolwarm",
    #                  center=0)
    # ax.title.set_text("Step " + str(i + 1))
# plt.tight_layout()
# plt.subplots_adjust(top=0.93)
# plt.show()
#########################################################################
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
