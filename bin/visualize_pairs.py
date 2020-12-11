#!/usr/bin/env python3
"""Membrane charge caluclations."""
import numpy as np
import sys
import pickle
import argparse
# import dgCalc
# import saltBridges
import math
import os.path
# import urllib.request
import collections
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import scipy
# from matplotlib.sankey import Sankey

plt.rcParams.update({'font.size': 14})
parser = argparse.ArgumentParser()

parser.add_argument("pickle_file", type=str, help="Pickle with membranes")

args = parser.parse_args()

positive = 'KR'
negative = 'DE'
charged = 'KRDE'
chargedplus = 'KRDEH'
positiveplus = 'RHK'
helices = pickle.load(open(args.pickle_file, 'rb'))
name = args.pickle_file.split('/')[-1].split('.')[0] + "_pairs"
number_of_segs = 0
edge = 5  # Number of residues to trim from start and end
seg_aas = ""
seg_aas_trimmed = ""

proteins_with_charges_full = set()
proteins_with_charges_full_sanH = set()
proteins_with_charges_trimmed = set()
proteins_with_charges_trimmed_sanH = set()

proteins_with_multi_charges_full = set()
proteins_with_multi_charges_full_sanH = set()
proteins_with_multi_charges_trimmed = set()
proteins_with_multi_charges_trimmed_sanH = set()

segs_with_charges_full = 0
segs_with_charges_full_sanH = 0
segs_with_charges_trimmed = 0
segs_with_charges_trimmed_sanH = 0

segs_with_multi_charges_full = list()
segs_with_multi_charges_full_sanH = list()
segs_with_multi_charges_trimmed = list()
segs_with_multi_charges_trimmed_sanH = list()

AAS = 'ACDEFGHIKLMNPQRSTVWY'
aa_pairs_full = np.zeros((len(AAS), len(AAS)), dtype=int)
aa_pairs_trimmed = np.zeros((len(AAS), len(AAS)), dtype=int)
charged_steps_trimmed = [[] for i in range(10)]
charged_steps_trimmed_sanH = [[] for i in range(10)]
segs_charged_steps_trimmed = [set() for i in range(10)]
segs_charged_steps_trimmed_sanH = [set() for i in range(10)]
prot_charged_steps_trimmed = [set() for i in range(10)]
prot_charged_steps_trimmed_sanH = [set() for i in range(10)]
opp_steps_trimmed = [[] for i in range(10)]
opp_steps_trimmed_sanH = [[] for i in range(10)]
segs_opp_steps_trimmed = [set() for i in range(10)]
segs_opp_steps_trimmed_sanH = [set() for i in range(10)]
prot_opp_steps_trimmed = [set() for i in range(10)]
prot_opp_steps_trimmed_sanH = [set() for i in range(10)]
same_steps_trimmed = [[] for i in range(10)]
same_steps_trimmed_sanH = [[] for i in range(10)]
segs_same_steps_trimmed = [set() for i in range(10)]
segs_same_steps_trimmed_sanH = [set() for i in range(10)]
prot_same_steps_trimmed = [set() for i in range(10)]
prot_same_steps_trimmed_sanH = [set() for i in range(10)]

protein_set = set()
protein_set_with_charges = set()
total_num_proteins = len(helices.items())
seg_charges = []
seg_charges_sanH = []
seg_charges_trimmed = []
seg_charges_sanH_trimmed = []
total_segments = 0
total_segments_with_charges = 0
same_pair_ids = [[],[],[],[],[],[],[],[],[],[]]
opp_pair_ids = [[],[],[],[],[],[],[],[],[],[]]
#######################################################
#### This is from logodd_barchart.py, could be better refactored
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
#######################################################
for key, membranes in helices.items():
    # Total number of membrane segments
    # print(key, membranes)
    for mem_num, mem_data in enumerate(membranes):
        mem_start, mem = mem_data
        # Only use transmembrane regions that are longer than 17 residues as per Baezo-Delgado 2012
        total_segments += 1
        if len(mem) < 17:
            continue
        protein_set.add(key)
        number_of_segs += 1
        # Create seg_full for the full segment and seg_trimmed
        # for the trimmed segment (5 caps)
        # Also fill the seg_aa* variables with all the amino acids
        seg_trimmed = mem[edge:-edge]
        seg_full = mem
        seg_aas_trimmed += seg_trimmed
        seg_aas += seg_full

        # Calculate number of charges for protein and segments
        # This is only used for number statistics
        # Proteins that contain charges in membrane segments
        if len(re.findall('[' + chargedplus + ']', seg_full)) > 0:
            proteins_with_charges_full.add(key)
            segs_with_charges_full += 1
            seg_charges.append(len(re.findall('[' + chargedplus + ']', seg_full)))
        # Proteins that contain charges(minus H) in membrane segments
        if len(re.findall('[' + charged + ']', seg_full)) > 0:
            proteins_with_charges_full_sanH.add(key)
            segs_with_charges_full_sanH += 1
            seg_charges_sanH.append(len(re.findall('[' + charged + ']', seg_full)))

        # Proteins that contain charges in membrane segments, trimmed
        if len(re.findall('[' + chargedplus + ']', seg_trimmed)) > 0:
            proteins_with_charges_trimmed.add(key)
            segs_with_charges_trimmed += 1
            seg_charges_trimmed.append(len(re.findall('[' + chargedplus + ']', seg_trimmed)))
        # Proteins that contain charges(minus H) in membrane segments, trimmed
        if len(re.findall('[' + charged + ']', seg_trimmed)) > 0:
            proteins_with_charges_trimmed_sanH.add(key)
            segs_with_charges_trimmed_sanH += 1
            seg_charges_sanH_trimmed.append(len(re.findall('[' + charged + ']', seg_trimmed)))

        # Multi charges in membrane segments
        if len(re.findall('[' + chargedplus + ']', seg_full)) > 1:
            proteins_with_multi_charges_full.add(key)
            segs_with_multi_charges_full.append(len(re.findall(
                                                '[' + chargedplus + ']',
                                                seg_full)))
        # Multi charges(minus H) in membrane segments
        if len(re.findall('[' + charged + ']', seg_full)) > 1:
            proteins_with_multi_charges_full_sanH.add(key)
            segs_with_multi_charges_full_sanH.append(len(re.findall(
                                                     '[' + charged + ']',
                                                     seg_full)))

        # Multi charges in membrane segments, trimmed
        if len(re.findall('[' + chargedplus + ']', seg_trimmed)) > 1:
            proteins_with_multi_charges_trimmed.add(key)
            segs_with_multi_charges_trimmed.append(len(re.findall(
                                                   '[' + chargedplus + ']',
                                                   seg_trimmed)))
        # Multi charges(minus H) in membrane segments, trimmed
        if len(re.findall('[' + charged + ']', seg_trimmed)) > 1:
            proteins_with_multi_charges_trimmed_sanH.add(key)
            segs_with_multi_charges_trimmed_sanH.append(len(re.findall(
                                                       '[' + charged + ']',
                                                       seg_trimmed)))

        # Use the full segments for full pairs
        seg_full_len = len(seg_full)
        for place, first_aa in enumerate(seg_full):
            for i in range(1, 11):
                if seg_full_len <= place + i:
                    break

                if first_aa not in AAS:
                    continue
                first_aa_index = AAS.index(first_aa)
                second_aa = seg_full[place + i]

                if second_aa not in AAS:
                    continue
                second_aa_index = AAS.index(second_aa)
                aa_pairs_full[first_aa_index][second_aa_index] += 1

        # Use the trimmed segments for the rest
        seg_trimmed_len = len(seg_trimmed)
        for place, first_aa in enumerate(seg_trimmed):
            for i in range(1, 11):
                if seg_trimmed_len <= place + i:
                    break

                if first_aa not in AAS:
                    continue
                first_aa_index = AAS.index(first_aa)
                second_aa = seg_trimmed[place + i]
                if second_aa not in AAS:
                    continue
                second_aa_index = AAS.index(second_aa)

                aa_pairs_trimmed[first_aa_index][second_aa_index] += 1

                if first_aa in chargedplus and second_aa in chargedplus:
                    # Any charges pairs
                    charged_steps_trimmed[i-1].append(first_aa + second_aa)
                    segs_charged_steps_trimmed[i-1].add(key + str(mem_num))
                    prot_charged_steps_trimmed[i-1].add(key)
                    if first_aa != "H" and second_aa != "H":
                        charged_steps_trimmed_sanH[i-1].append(first_aa +
                                                              second_aa)
                        segs_charged_steps_trimmed_sanH[i-1].add(key +
                                                                str(mem_num))
                        prot_charged_steps_trimmed_sanH[i-1].add(key)
                    # Opposite charges pairs
                    if ((first_aa in positiveplus and
                         second_aa in negative) or
                        (first_aa in negative and
                         second_aa in positiveplus)):
                            opp_steps_trimmed[i-1].append(first_aa + second_aa)
                            segs_opp_steps_trimmed[i-1].add(key + str(mem_num))
                            prot_opp_steps_trimmed[i-1].add(key)
                            if first_aa != 'H' and second_aa != 'H':
                                if first_aa in negative and second_aa in positive:
                                    opp_pair_ids[i-1].append(key + '\tnegpos')
                                else:
                                    opp_pair_ids[i-1].append(key + '\tposneg')
                                opp_steps_trimmed_sanH[i-1].append(
                                                       first_aa + second_aa)
                                segs_opp_steps_trimmed_sanH[i-1].\
                                    add(key + str(mem_num))
                                prot_opp_steps_trimmed_sanH[i-1].add(key)

                    # Same charge pairs
                    if ((first_aa in positiveplus and
                         second_aa in positiveplus) or
                        (first_aa in negative and
                         second_aa in negative)):
                            same_steps_trimmed[i-1].append(first_aa +
                                                          second_aa)
                            segs_same_steps_trimmed[i-1].add(key +
                                                            str(mem_num))
                            prot_same_steps_trimmed[i-1].add(key)
                            if first_aa != 'H' and second_aa != 'H':
                                if first_aa in negative:
                                    same_pair_ids[i-1].append(key + '\tneg')
                                else:
                                    same_pair_ids[i-1].append(key + '\tpos')
                                same_steps_trimmed_sanH[i-1].append(
                                                       first_aa + second_aa)
                                segs_same_steps_trimmed_sanH[i-1].\
                                    add(key + str(mem_num))
                                prot_same_steps_trimmed_sanH[i-1].add(key)

            ######################################
            # Below is cut in from logodd_barchart.py, could be refactored with the above 
            ######################################
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

num_prots = len(protein_set)

################################################
### From logodd_barchart.py, could be refactored better
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
    logOdd = math.log(odds)
    logOdds.append(logOdd)
    SE = math.sqrt(1/a+1/b+1/c+1/d)
    z = abs(logOdd/SE)  ## Two sided
    pm = math.exp(-0.717*z-0.416*z**2)*20*20*8*2 ## two sided
    p = (1-scipy.stats.norm.cdf(abs(z)))*2
    psf = (scipy.stats.norm.sf(abs(z)))*2
    ci.append((math.log(odds)+1.96*math.sqrt(1/a+1/b+1/c+1/d), math.log(odds)-1.96*math.sqrt(1/a+1/b+1/c+1/d), p))
    print("Total {}: logOdds: {} SE: {} p: {:.5f}".format(i+1, logOdd, SE, psf))
    ###### For Opp charges ######
    a = hitcounter[i]                 # Number of opposite pairs, only opposite
    b = len(pairs[i])       # Number of total observed pairs
    c = poschargeCount[i][0]*negchargeCount[i][1] + negchargeCount[i][0]*poschargeCount[i][1]  # Only calculate opposite pairs
    d = aaCount[i]**2       # Number of non-hits for baseline
    odds = (a/b)/(c/d)
    logOdd = math.log(odds)
    logOddsOpp.append(logOdd)
    SE = math.sqrt(1/a+1/b+1/c+1/d)
    z = abs(logOdd/SE)  ## Two sided
    pm = math.exp(-0.717*z-0.416*z**2)*20*20*8*2 ## two sided
    p = (1-scipy.stats.norm.cdf(abs(z)))*2
    psf = (scipy.stats.norm.sf(abs(z)))*2
    ciOpp.append((math.log(odds)+1.96*math.sqrt(1/a+1/b+1/c+1/d), math.log(odds)-1.96*math.sqrt(1/a+1/b+1/c+1/d), p))
    print("Opp {}: logOdds: {} SE: {} p: {:.5f}".format(i+1, logOdd, SE, psf))

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
    logOdd = math.log(odds)
    logOddsSame.append(logOdd)
    SE = math.sqrt(1/a+1/b+1/c+1/d)
    z = abs(logOdd/SE)  ## Two sided
    pm = math.exp(-0.717*z-0.416*z**2)*20*20*8*2 ## two sided
    p = (1-scipy.stats.norm.cdf(abs(z)))*2
    psf = (scipy.stats.norm.sf(abs(z)))*2
    ciSame.append((math.log(odds)+1.96*math.sqrt(1/a+1/b+1/c+1/d), math.log(odds)-1.96*math.sqrt(1/a+1/b+1/c+1/d), p))
    print("Same {}: logOdds: {} SE: {} p: {:.5f}".format(i+1, logOdd, SE, psf))
################################################
##########
# Graphs #
##########
#### First make the top graph
charge_hist_data = list(np.concatenate(
                        [[i+1]*len(pos)
                            for i, pos in enumerate(charged_steps_trimmed)]))
charge_hist_data_sanH = list(np.concatenate(
                        [[i+1]*len(pos)
                            for i, pos in
                            enumerate(charged_steps_trimmed_sanH)]))
opp_hist_data = list(np.concatenate(
                        [[i+1]*len(pos)
                            for i, pos in enumerate(opp_steps_trimmed)]))
opp_hist_data_sanH = list(np.concatenate(
                        [[i+1]*len(pos)
                            for i, pos in enumerate(opp_steps_trimmed_sanH)]))
same_hist_data = list(np.concatenate(
                        [[i+1]*len(pos)
                            for i, pos in enumerate(same_steps_trimmed)]))
same_hist_data_sanH = list(np.concatenate(
                        [[i+1]*len(pos)
                            for i, pos in enumerate(same_steps_trimmed_sanH)]))
df_charge = pd.DataFrame(charge_hist_data_sanH, columns=["Step"])
df_charge["Type"] = "Charged"
df_same = pd.DataFrame(same_hist_data_sanH, columns=["Step"])
df_same["Type"] = "Same"
df_opp = pd.DataFrame(opp_hist_data_sanH, columns=["Step"])
df_opp["Type"] = "Opp"
df = pd.concat([df_same, df_opp], ignore_index=True)
df = df[df["Step"] < 9]
df["Step"] = "i+" + df["Step"].astype(int).astype(str)

nodes = [0, 0.25, 0.5, 0.75, 1]
colors = ["#FDE725FF", "#440154FF", "#FDE725FF"]
# Regular red and green
# colors = ["#6BE585", "#DD3E54", "#6BE585"]
# degree_cmap = LinearSegmentedColormap.from_list("", list(zip(nodes, colors)))
# degree_cmap = LinearSegmentedColormap.from_list("", list(zip(nodes, colors)))
degree_cmap = mpl.colors.ListedColormap(mpl.cm.get_cmap('viridis_r').colors + mpl.cm.get_cmap('viridis').colors)
# print(cmap)
# grid = plt.GridSpec(3,6, wspace=0.4, hspace=0.1)
sns.set_theme(style="white", context="paper")
# sns.set(fontsize=14)
f, axes = plt.subplots(5,1,figsize=(7, 12), gridspec_kw={"height_ratios":[64,1,16,16,16]})
h = sns.histplot(df, x="Step", color="grey", hue="Type", discrete=True, multiple="stack", shrink=.8, ax=axes[0])
axes[0].get_legend().remove()
hatches = {0:"///", 1:"\\\\\\", 2:"|||"}
for i in range(2):
    for j in range(8):
        color_index = (j*100 + 100) % 360
        axes[0].patches[j+i*8].set_facecolor(degree_cmap(color_index/360))
        axes[0].patches[j+i*8].set_hatch(hatches[i])
sns.despine(bottom=True, left=True)
##### Here we start the second graph
################################################ 
# f, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 10)) # , sharex=True)
###### Totalt ######
axes[1].set_visible(False)
y_r = [logOdds[i] - ci[i][0] for i in range(len(ci))]
threshold = 0
values = np.array(logOdds)
x = ['i+' + str(x) for x in range(1, len(values)+1)]
above_threshold = np.maximum(values - threshold, 0)
below_threshold = np.minimum(values, threshold)
axes[2].bar(x, logOdds, yerr=y_r, color=['g', 'r', 'g', 'g', 'r', 'r', 'g', 'r'], alpha=0.8, align='center')
axes[2].set_ylabel('All charges')
axes[2].set_yticks([-1,-0.5, 0, 0.5, 1, 1.5, 2])
axes[2].axhline(xmax=8, color='black')

###### Opp ######
y_r = [logOddsOpp[i] - ciOpp[i][0] for i in range(len(ciOpp))]
threshold = 0
values = np.array(logOddsOpp)
x = ['i+' + str(x) for x in range(1, len(values)+1)]
above_threshold = np.maximum(values - threshold, 0)
below_threshold = np.minimum(values, threshold)
axes[3].bar(x, logOddsOpp, yerr=y_r, color=['g', 'r', 'g', 'g', 'r', 'r', 'g', 'r'], alpha=0.8, align='center')
axes[3].set_ylabel('Opposite charges')
axes[3].set_yticks([-1,-0.5, 0, 0.5, 1, 1.5, 2])
axes[3].axhline(xmax=8, color='black')

###### Same ######
y_r = [logOddsSame[i] - ciSame[i][0] for i in range(len(ciSame))]
threshold = 0
values = np.array(logOddsSame)
x = ['i+' + str(x) for x in range(1, len(values)+1)]
above_threshold = np.maximum(values - threshold, 0)
below_threshold = np.minimum(values, threshold)
axes[4].bar(x, logOddsSame, yerr=y_r, color=['g', 'r', 'g', 'g', 'r', 'r', 'g', 'r'], alpha=0.8, align='center')
axes[4].set_ylabel('Same charges')
axes[4].set_yticks([-1,-0.5, 0, 0.5, 1, 1.5, 2])
axes[4].set_xlabel('Step')
axes[4].axhline(xmax=8, color='black')

axes[2].set_xticks([])
axes[2].xaxis.set_tick_params(length=0)
axes[3].xaxis.set_tick_params(length=0)
axes[3].set_xticks([])
# axes[4].set_xticks([1, 2, 3, 4, 5, 6, 7, 8])

degree_cmap = mpl.colors.ListedColormap(mpl.cm.get_cmap('viridis_r').colors + mpl.cm.get_cmap('viridis').colors)
for ax in [axes[2], axes[3], axes[4]]:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    for i in range(8):
        color_index = (i*100 + 100) % 360
        ax.patches[i].set_facecolor(degree_cmap(color_index/360))
###########################################
plt.tight_layout()
####################################
# #### Add the color wheel ###########
degrees = np.linspace(0,2*np.pi, 18, endpoint=False).tolist()
theta = [2*np.pi*((x%360)/360) for x in range(0, 1701, 100)]
order = [0, 11, 4, 15, 8, 1, 12, 5, 16, 9, 2, 13, 6, 17, 10, 3, 14, 7]
c_ax = f.add_axes([0.60, 0.68, 0.40,0.40],projection='polar')
c_ax.set_theta_offset(np.pi/2)
c_ax.set_theta_direction(-1)
# norm = mpl.colors.Normalize(0.0, 2*np.pi)
# cb = mpl.colorbar.ColorbarBase(c_ax, cmap=degree_cmap,
#                                    norm=norm,
#                                    orientation='horizontal')
# cb.set_ticks(np.linspace(0,2*np.pi, 14, endpoint=False).tolist())
# cb.set_ticklabels(["", " ", "Step 4", " ", "Step 8", "Step 1", " ", "Step 5", " ", " ", "Step 2", " ", "Step 6", " ", " ", "Step 3", " ", "Step 7"])
# cb.ax.plot([0,0],[0,1])
# cb.ax.set_title("First charged residue\n\n")
# cb.ax.xaxis.set_tick_params(pad=16, length=0)
# cb.outline.set_visible(False)                                 
# c_ax.set_rlim([-2,1])
####################################
c_ax.plot([x for _,x in sorted(zip(order, degrees))][:9],len(degrees[:9])*[0.9], zorder=1, linestyle='--', color='lightgray')
first = c_ax.scatter(theta[:9], len(theta[:9])*[0.9], c=theta[:9], s=300, cmap=degree_cmap, zorder=2, alpha=0.5,edgecolor='black')
second = c_ax.scatter(theta[9:], len(theta[9:])*[0.9], s=300, color="lightgray", zorder=3,edgecolor='black')
c_ax.xaxis.set_tick_params(pad=16, length=0)
c_ax.set_xticks(degrees)
c_ax.get_yaxis().set_visible(False)
# ax.set_xticklabels(["", "i+11", "i+4", "i+15", "i+8", "i+1", "i+12", "i+5", "i+16",
#                    "i+9", "i+2", "i+13", "i+6", "i+17", "i+10", "i+3", "i+14", "i+7"])
c_ax.get_xaxis().set_visible(False)
c_ax.grid(False)
c_ax.spines['polar'].set_visible(False)
c_ax.set_rlim([0,1.1])
for l,x,y in sorted(zip(order, degrees, len(degrees)*[0.9]))[:9]:
    if l == 0:
        label = "i"
    else:
        label = "+{}".format(l)
    c_ax.annotate(label, (x,y), ha="center", va="center",weight='bold')
#### adding a) and b) text
plt.text(-0.05, 0.95, "a)", fontsize=16, fontdict={"weight":'bold'}, transform=axes[0].transAxes)
plt.text(-0.05, 0.10, "b)", fontsize=16, fontdict={"weight":'bold'}, transform=axes[1].transAxes)
# plt.savefig("test.png")
for ext in [".png", ".svg"]:
    plt.savefig('images/' + name + ext)
