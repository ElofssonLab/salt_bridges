#!/usr/bin/env python3
"""Membrane charge caluclations."""
import matplotlib as mpl
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
plt.rcParams.update({'font.size': 14})
parser = argparse.ArgumentParser()

# parser.add_argument("glob", type=str, help="In glob charges pickle")
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
# prefix = ["Globular"] * 20 + ["Membranes"]*20
prefix = ["Membranes"]*20
new_cols = []
for i, p in enumerate(prefix):
    new_cols.append(p + " " + aas[i%20])
# in_pickle = 'data/pdbtm/pdbtm_redundant_alpha_struct_scop_reduced_culled_charge.pickle'
# in_glob_pickle = args.glob
in_mem_pickle = args.mems
# out_image = in_pickle.split('/')[-1].split('.')[0] + "_AA_vis.svg"
# out_image = args.out
# glob_aaCount, glob_aaPairs, glob_aaHits = pickle.load(open(in_glob_pickle, 'rb'))
mem_aaCount, mem_aaPairs, mem_aaHits = pickle.load(open(in_mem_pickle, 'rb'))
    # open('data/chargeDataGlobularFromPDBRed50TrimmedLen15v2.pickle',
    #      'rb'))

# globlogOdds = np.zeros([7, 20, 20])
memlogOdds = np.zeros([8, 20, 20])
ci = []
# print([aas.index(c) for c in chargedplus])
aaIndex = [aas.index(c) for c in positiveplus]
# print(sum(aaHits[0][chargesPlusIndex][chargesPlusIndex]))
# for ind, aa in enumerate(aasOri):
#    print(aa, aaCount[ind])
# for i in range(7):
#     for first in range(20):
#         for second in range(20):
#             a = glob_aaHits[i][first][second]
#             b = glob_aaPairs[i] - a
#             # print(a, b, a/b)
#             c = glob_aaCount[0][i][first] * glob_aaCount[1][i][second]
#             # Need square to handle pairs, first and second
#             # d = sum(aaCount)**2 # Old version
#             # d is aaPairs for this step times 2 for each of the two
#             # amino acids as each pair contains two AAs.
#             d = glob_aaPairs[i]**2 - c
#             # print(a, b, c, d)
#             odds = (a / b) / (c / d)
#             if odds == 0.0:
#                 logOddsValue = -math.inf
#                 print("Odds 0 at step: ", i + 1)
#                 print("A-D: ", a, b, c, d)
#                 print("First and second AA: ", aasOri[first], aasOri[second])
#                 # print(aasOri[first], aasOri[second], a, b, c, d, odds)
#             else:
#                 # std_err = math.sqrt(1/a + 1/b + 1/c + 1/d)
#                 logOddsValue = math.log(odds)
#             globlogOdds[i][aas.index(aasOri[first])][aas.index(aasOri[second])]\
#                 = logOddsValue
for i in range(8):
    for first in range(20):
        for second in range(20):
            a = mem_aaHits[i][first][second]
            b = mem_aaPairs[i]  # - a, not odds
            # print(a, b, a/b)
            c = mem_aaCount[0][i][first] * mem_aaCount[1][i][second]
            # Need square to handle pairs, first and second
            # d = sum(aaCount)**2 # Old version
            # d is aaPairs for this step times 2 for each of the two
            # amino acids as each pair contains two AAs.
            d = mem_aaPairs[i]**2  #  - c, not odds
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
# for method in methods:
method='complete'

# plt.savefig(method + "_" + out_image)
# f, axarr = plt.subplots(3, 3, figsize=(16, 14))
# axarr[-1, -1].axis('off')
# axarr[-1, -2].axis('off')
# plt.suptitle("TMs alpha helices, trimmed, pdb50, len > 15, V2")
# plt.suptitle(args.title)
# for method in methods:
# df = pd.DataFrame()
df_full = pd.DataFrame()
df = pd.DataFrame()
for i in range(8):
    # df1 = pd.DataFrame(globlogOdds[i],
    #                   index=["Globular " + str(c) for c in aas],
    #                   columns=[c for c in aas])
    # G_keep = ["Globular " + str(c) for c in keep]
    # df1 = df1.loc[G_keep][keep]
    # df1.insert(loc=0,column= "Type",  value="Globular")
    df_full_temp = pd.DataFrame(memlogOdds[i],
                      index=["Membranes " + str(c) for c in aas],
                      columns=[c for c in aas])
    # df2.insert( loc=0, column= "Type", value="Membranes")
    M_keep = ["Membranes " + str(c) for c in keep]
    df_temp = df_full_temp.loc[M_keep][keep]
    df_full_temp.columns = ["Step " + str(i+1) + "-" + c for c in df_full_temp.columns]
    df_temp.columns = ["Step " + str(i+1) + "-" + c for c in df_temp.columns]
    # df_temp = pd.concat([df1, df2])
    # df2.columns = ["Step " + str(i+1) + "-" + c for c in df2.columns]
    df = pd.concat([df, df_temp], axis=1)
    df_full = pd.concat([df_full, df_full_temp], axis=1)
df.insert(loc=0, column="Type", value=["Membranes"]*keep_len)
df_full.insert(loc=0, column="Type", value=["Membranes"]*20)
# print(df)
# Actual graphing
# plt.suptitle(args.title +  " sep " + str(i+1))
types =df.pop("Type") 
step_list = [item for i in range(1, 9) for item in [i]*keep_len]
steps = pd.Series(step_list, name="Steps", index=df.columns)
step_c_list = [item for i in range(1, 9) for item in [(((100*i)%360)/360)]*keep_len]
steps_c = pd.Series(step_c_list, name="C_cols", index=df.columns)
# type_list = pd.Series(["Globular"] * 20 + ["Membranes"]*20)
lut = dict(zip(types.unique(), sns.husl_palette(2, s=.45)))
###################################################
## Gradient from #6BE585 to #DD3E54 using 19 bins (for 0-180 degrees) at https://colordesigner.io/gradient-generator  #####
###################################################
## 0    #6be585
## 10   #74df83
## 20   #7dd980
## 30   #85d37e
## 40   #8dcc7c
## 50   #94c579
## 60   #9bbe77
## 70   #a1b774
## 80   #a8b072
## 90   #aea86f
## 100  #b4a06d
## 110  #b9976a
## 120  #bf8e67
## 130  #c48464
## 140  #c97961
## 150  #ce6d5e
## 160  #d3605b
## 170  #d85157
## 180  #dd3e54

degree_cmap = mpl.colors.ListedColormap(mpl.cm.get_cmap('viridis_r').colors + mpl.cm.get_cmap('viridis').colors)
#                                    100        200         300         40          140     240             340     80
step_index = [x%360 for x in range(100, 900, 100)]
lut_cols = dict(zip(steps.unique(), [degree_cmap(s/360) for s in step_index] ))
row_colors = types.map(lut)
col_colors = steps.map(lut_cols)
chart = sns.clustermap(df, row_cluster=True, col_cluster=True, figsize=(30, 6), cmap="coolwarm", center=0,
               method=method, col_colors=col_colors, cbar_pos=None, vmin=-2, vmax=2,
               linewidth=.75,
               dendrogram_ratio=(0.05, 0.15),
               colors_ratio=(0.005, 0.02),
               )
plt.setp(chart.ax_heatmap.xaxis.get_majorticklabels(), rotation=45, ha="right")
name = args.out
for ext in [".png", ".svg"]:
    chart.savefig(name + ext)

##### Make the same with full data
plt.clf()
# plt.suptitle(args.title +  " sep " + str(i+1))
types =df_full.pop("Type") 
step_list = [item for i in range(1, 9) for item in [i]*20]
steps = pd.Series(step_list, name="Steps", index=df_full.columns)
step_c_list = [item for i in range(1, 9) for item in [(((100*i)%360)/360)]*20]
steps_c = pd.Series(step_c_list, name="C_cols", index=df_full.columns)
# type_list = pd.Series(["Globular"] * 20 + ["Membranes"]*20)
lut = dict(zip(types.unique(), sns.husl_palette(2, s=.45)))
degree_cmap = mpl.colors.ListedColormap(mpl.cm.get_cmap('viridis_r').colors + mpl.cm.get_cmap('viridis').colors)
#                                    100        200         300         40          140     240             340     80
step_index = [x%360 for x in range(100, 900, 100)]
lut_cols = dict(zip(steps.unique(), [degree_cmap(s/360) for s in step_index] ))
row_colors = types.map(lut)
col_colors = steps.map(lut_cols)
chart = sns.clustermap(df_full, row_cluster=True, col_cluster=False, figsize=(40, 10), cmap="coolwarm", center=0,
               method=method, col_colors=col_colors, cbar_pos=None, vmin=-2, vmax=2,
               linewidth=.75,
               dendrogram_ratio=(0.05, 0.15),
               colors_ratio=(0.005, 0.02),
               )
plt.setp(chart.ax_heatmap.xaxis.get_majorticklabels(), rotation=45, ha="right")
name = args.out + "_full"
for ext in [".png", ".svg"]:
    chart.savefig(name + ext)
##################################
