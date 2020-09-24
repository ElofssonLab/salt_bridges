#!/usr/bin/env python3
"""Membrane charge caluclations."""
import numpy as np
import sys
import pickle
import argparse
# import dgCalc
# import saltBridges
import os.path
# import urllib.request
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
# from matplotlib.sankey import Sankey

parser = argparse.ArgumentParser()

parser.add_argument("pickle_file", type=str, help="Pickle with membranes")
# parser.add_argument("threeline_file", type=str, help="Corresponding 3line files")

args = parser.parse_args()

# pickle_file_helices = 'data/topcons/TMs.pickle'
# full_3line_file = 'data/topcons/TMs.3line'
# pdbURL = "https://files.rcsb.org/download/"
positive = 'KR'
negative = 'DE'
charged = 'KRDE'
chargedplus = 'KRDEH'
positiveplus = 'RHK'
helices = pickle.load(open(args.pickle_file, 'rb'))
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
charged_gaps_trimmed = [[] for i in range(10)]
charged_gaps_trimmed_sanH = [[] for i in range(10)]
segs_charged_gaps_trimmed = [set() for i in range(10)]
segs_charged_gaps_trimmed_sanH = [set() for i in range(10)]
prot_charged_gaps_trimmed = [set() for i in range(10)]
prot_charged_gaps_trimmed_sanH = [set() for i in range(10)]
opp_gaps_trimmed = [[] for i in range(10)]
opp_gaps_trimmed_sanH = [[] for i in range(10)]
segs_opp_gaps_trimmed = [set() for i in range(10)]
segs_opp_gaps_trimmed_sanH = [set() for i in range(10)]
prot_opp_gaps_trimmed = [set() for i in range(10)]
prot_opp_gaps_trimmed_sanH = [set() for i in range(10)]
same_gaps_trimmed = [[] for i in range(10)]
same_gaps_trimmed_sanH = [[] for i in range(10)]
segs_same_gaps_trimmed = [set() for i in range(10)]
segs_same_gaps_trimmed_sanH = [set() for i in range(10)]
prot_same_gaps_trimmed = [set() for i in range(10)]
prot_same_gaps_trimmed_sanH = [set() for i in range(10)]
# proteinsExamples = {}
# proteinsExamplesSansH = 0
# proteinSet = set()
# proteinChargedSet = set()
# proteinSetAll = set()
# proteinChargedSetAll = set()
# proteinSetTrimmed = set()
# proteinChargedSetTrimmed = set()
# proteinSetAllTrimmed = set()
# proteinChargedSetAllTrimmed = set()
# mems_with_two_or_more_charges = []
# mems_with_two_or_more_charges_trimmed = []
# mems_with_two_or_more_charges_plus = []
# mems_with_two_or_more_charges_plus_trimmed = []
# totalMems = 0
# TMdata = {}
# mems_aa = ''
# mems_aa_trimmed = ''
# opposite_gaps_trimmed = [[],
#                          [],
#                          [],
#                          [],
#                          [],
#                          [],
#                          []]
# charge_gaps_trimmed = [[],
#                        [],
#                        [],
#                        [],
#                        [],
#                        [],
#                        []]
# opposite_gaps_trimmed_noH = [[],
#                              [],
#                              [],
#                              [],
#                              [],
#                              [],
#                              []]
# 
# with open(full_3line_file, 'r') as TMHandle:
#     pid = ''
#     rowNum = 0
#     fa = ''
#     for line in TMHandle:
#         if line[0] == '>':
#             pid = line[1:].strip()
#         elif rowNum % 3 == 1:
#             fa = line.strip()
#         elif rowNum % 3 == 2:
#             TMdata[pid] = [fa, line.strip()]
#         else:
#             print("You should not be here...")
#         rowNum += 1
# testlist = []
protein_set = set()
total_num_proteins = len(helices.items())
total_segments = 0
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
        # Proteins that contain charges(minus H) in membrane segments
        if len(re.findall('[' + charged + ']', seg_full)) > 0:
            proteins_with_charges_full_sanH.add(key)
            segs_with_charges_full_sanH += 1

        # Proteins that contain charges in membrane segments, trimmed
        if len(re.findall('[' + chargedplus + ']', seg_trimmed)) > 0:
            proteins_with_charges_trimmed.add(key)
            segs_with_charges_trimmed += 1
        # Proteins that contain charges(minus H) in membrane segments, trimmed
        if len(re.findall('[' + charged + ']', seg_trimmed)) > 0:
            proteins_with_charges_trimmed_sanH.add(key)
            segs_with_charges_trimmed_sanH += 1

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
                    charged_gaps_trimmed[i-1].append(first_aa + second_aa)
                    segs_charged_gaps_trimmed[i-1].add(key + str(mem_num))
                    prot_charged_gaps_trimmed[i-1].add(key)
                    if first_aa != "H" and second_aa != "H":
                        charged_gaps_trimmed_sanH[i-1].append(first_aa +
                                                              second_aa)
                        segs_charged_gaps_trimmed_sanH[i-1].add(key +
                                                                str(mem_num))
                        prot_charged_gaps_trimmed_sanH[i-1].add(key)
                    # Opposite charges pairs
                    if ((first_aa in positiveplus and
                         second_aa in negative) or
                        (first_aa in negative and
                         second_aa in positiveplus)):
                            opp_gaps_trimmed[i-1].append(first_aa + second_aa)
                            segs_opp_gaps_trimmed[i-1].add(key + str(mem_num))
                            prot_opp_gaps_trimmed[i-1].add(key)
                            if first_aa != 'H' and second_aa != 'H':
                                opp_gaps_trimmed_sanH[i-1].append(
                                                       first_aa + second_aa)
                                segs_opp_gaps_trimmed_sanH[i-1].\
                                    add(key + str(mem_num))
                                prot_opp_gaps_trimmed_sanH[i-1].add(key)

                    # Same charge pairs
                    if ((first_aa in positiveplus and
                         second_aa in positiveplus) or
                        (first_aa in negative and
                         second_aa in negative)):
                            same_gaps_trimmed[i-1].append(first_aa +
                                                          second_aa)
                            segs_same_gaps_trimmed[i-1].add(key +
                                                            str(mem_num))
                            prot_same_gaps_trimmed[i-1].add(key)
                            if first_aa != 'H' and second_aa != 'H':
                                same_gaps_trimmed_sanH[i-1].append(
                                                       first_aa + second_aa)
                                segs_same_gaps_trimmed_sanH[i-1].\
                                    add(key + str(mem_num))
                                prot_same_gaps_trimmed_sanH[i-1].add(key)

num_prots = len(protein_set)
##########
# Graphs #
##########

charge_hist_data = list(np.concatenate(
                        [[i+1]*len(pos)
                            for i, pos in enumerate(charged_gaps_trimmed)]))
charge_hist_data_sanH = list(np.concatenate(
                        [[i+1]*len(pos)
                            for i, pos in
                            enumerate(charged_gaps_trimmed_sanH)]))
opp_hist_data = list(np.concatenate(
                        [[i+1]*len(pos)
                            for i, pos in enumerate(opp_gaps_trimmed)]))
opp_hist_data_sanH = list(np.concatenate(
                        [[i+1]*len(pos)
                            for i, pos in enumerate(opp_gaps_trimmed_sanH)]))
same_hist_data = list(np.concatenate(
                        [[i+1]*len(pos)
                            for i, pos in enumerate(same_gaps_trimmed)]))
same_hist_data_sanH = list(np.concatenate(
                        [[i+1]*len(pos)
                            for i, pos in enumerate(same_gaps_trimmed_sanH)]))
df_charge = pd.DataFrame(charge_hist_data_sanH, columns=["Gap"])
df_charge["Type"] = "Charged"
df_same = pd.DataFrame(same_hist_data_sanH, columns=["Gap"])
df_same["Type"] = "Same"
df_opp = pd.DataFrame(opp_hist_data_sanH, columns=["Gap"])
df_opp["Type"] = "Opp"
df = pd.concat([df_charge, df_same, df_opp], ignore_index=True)
df = df[df["Gap"] < 8]

name = args.pickle_file.split('/')[-1].split('.')[0]
sns.set_theme(style="white", context="talk")
f, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
f.suptitle(name + ", distance between charges")
sns.histplot(df[df["Type"]=="Charged"], x="Gap", color="grey", discrete=True, ax=ax1)
ax1.axhline(0, color="k", clip_on=False)
ax1.set_ylabel("Charged")
sns.histplot(df[df["Type"]=="Opp"], x="Gap", color="grey", discrete=True, ax=ax2)
ax2.axhline(0, color="k", clip_on=False)
ax2.set_ylabel("Opp")
sns.histplot(df[df["Type"]=="Same"], x="Gap", color="grey", discrete=True, ax=ax3)
ax3.axhline(0, color="k", clip_on=False)
ax3.set_ylabel("Same")
for ax in [ax1, ax2, ax3]:
    ax.patches[0].set_facecolor("crimson")
    ax.patches[2].set_facecolor("orange")
    ax.patches[3].set_facecolor("yellow")
sns.despine(bottom=True)
plt.tight_layout(h_pad=2)
plt.savefig('images/' + name + '.png')
# sns_plot = sns.displot(df, x="Gap", binwidth=1, discrete=True, hue="Type")
# sns_plot.fig.set_figwidth(15)
# sns_plot.fig.set_figheight(10)
# sns_plot.savefig("images/test.png")

# f, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7,5), sharex=True)
# # # Charged pairs
# # plt.hist(charge_hist_data, bins=np.arange(12)-0.5)
# # plt.xticks(range(1, 11))
# # plt.yticks(range(0, 200, 5))
# # plt.xlim([0, 11])
# # plt.ylim([0, 200])
# # plt.title("Distance between charged AA's")
# # plt.savefig('images/charged_pairs.svg')
# # plt.clf()
# # Charged pairs sanH
# plt.hist(charge_hist_data_sanH, bins=np.arange(12)-0.5)
# plt.xticks(range(1, 11))
# plt.yticks(range(0, 400, 5))
# plt.xlim([0, 11])
# plt.ylim([0, 400])
# plt.title("Distance between charged AA's san H")
# plt.savefig('images/charged_pairs_sanH.svg')
# plt.clf()
# # # Opposite pairs
# # plt.hist(opp_hist_data, bins=np.arange(12)-0.5)
# # plt.xticks(range(1, 11))
# # plt.yticks(range(0, 200, 5))
# # plt.xlim([0, 11])
# # plt.ylim([0, 200])
# # plt.title("Distance between opposite charged AA's")
# # plt.savefig('images/opp_pairs.svg')
# # plt.clf()
# # Opposite pairs sanH
# plt.hist(opp_hist_data_sanH, bins=np.arange(12)-0.5)
# plt.xticks(range(1, 11))
# plt.yticks(range(0, 400, 5))
# plt.xlim([0, 11])
# plt.ylim([0, 400])
# plt.title("Distance between opposite charged AA's san H")
# plt.savefig('images/opp_pairs_sanH.svg')
# plt.clf()
# # # Same charged pairs
# # plt.hist(same_hist_data, bins=np.arange(12)-0.5)
# # plt.xticks(range(1, 11))
# # plt.yticks(range(0, 200, 5))
# # plt.xlim([0, 11])
# # plt.ylim([0, 200])
# # plt.title("Distance between same charged AA's")
# # plt.savefig('images/same_pairs.svg')
# # plt.clf()
# # Same charged pairs sanH
# plt.hist(same_hist_data_sanH, bins=np.arange(12)-0.5)
# plt.xticks(range(1, 11))
# plt.yticks(range(0, 400, 5))
# plt.xlim([0, 11])
# plt.ylim([0, 400])
# plt.title("Distance between same charged AA's san H")
# plt.savefig('images/same_pairs_sanH.svg')
# plt.clf()
# 
# sys.exit()
#
#                     deltaG = dgCalc.calc_segment_DG(fullMem)
#                     # Add five is midmem
#                     globalPlace = mem[0] + edge
#                     # globalPlace = TMdata[key][0].index(midMem)
#                     # Plus one for the globalPlace and plus 1 for the place
#                     # in the membrane, 0 -> 1 offset
#                     proteinsExamples[key +
#                                      str(globalPlace + place) +
#                                      str(globalPlace + place + i)] = \
#                                     [midMem,
#                                      place,
#                                      i,
#                                      globalPlace + place,
#                                      aa,
#                                      secAA,
#                                      deltaG]
#                     proteinSetAll.add(key)
#                     if aa != 'H' and secAA != 'H':
#                         proteinSet.add(key)
#                         proteinsExamplesSansH += 1
#                 # aaHits[i - 1][first][second] += 1
#                 # # Count both first and second AA for each gap sized i
#                 # aaCount[0][i - 1][first] += 1
#                 # aaCount[1][i - 1][second] += 1


print("###################################################")
print("## Numbers within parenthesis are inclusive of H ##")
print("## Trimmed means the transmembrane segment has   ##")
print("## been trimmed to not include the first and     ##")
print("## last 5 residues to prevent boundary charges.  ##")
print("###################################################")
print()
print()
print("Total number of proteins:".ljust(70), total_num_proteins)
print("Total number of transmembrane segments: ".ljust(70), total_segments)
print()
print("Number of proteins(with transmembrane regions >17):".ljust(70), num_prots)
print("Number of transmembrane segments(min 17 residues): ".ljust(70), number_of_segs)
print()
print("Amino acids in trans segments: ".ljust(70), len(seg_aas))
print("Amino acids in trans segments, trimmed: ".ljust(70),
      len(seg_aas_trimmed))
print()
print("Proteins with charged amino acids(in trans seg):".ljust(70),
      len(proteins_with_charges_full_sanH),
      '(', len(proteins_with_charges_full), ')')
print("Proteins with charged amino acids(in trans seg), trimmed:"
      .ljust(70),
      len(proteins_with_charges_trimmed_sanH),
      '(', len(proteins_with_charges_trimmed), ')')
print("Segments with charged amino acids:".ljust(70),
      segs_with_charges_full_sanH,
      '(', segs_with_charges_full, ')')
print("Segments with charged amino acids, trimmed:".ljust(70),
      segs_with_charges_trimmed_sanH,
      '(', segs_with_charges_trimmed, ')')
print()
print("Proteins with multi charged amino acids(in trans seg):".ljust(70),
      len(proteins_with_multi_charges_full_sanH),
      '(', len(proteins_with_multi_charges_full), ')')
print("Proteins with multi charged amino acids(in trans seg), trimmed:"
      .ljust(70),
      len(proteins_with_multi_charges_trimmed_sanH),
      '(', len(proteins_with_multi_charges_trimmed), ')')
print()
print("Segments with multi charged amino acids:".ljust(70),
      len(segs_with_multi_charges_full_sanH),
      '(', len(segs_with_multi_charges_full), ')')
print("Number of charges and count:")
print(pd.Series(segs_with_multi_charges_full_sanH).value_counts().to_string())
print("Number of charges and count with H:")
print(pd.Series(segs_with_multi_charges_full).value_counts().to_string())
print()
print("Segments with multi charged amino acids, trimmed:".ljust(70),
      len(segs_with_multi_charges_trimmed_sanH),
      '(', len(segs_with_multi_charges_trimmed), ')')
print("Number of charges and count:")
print(pd.Series(segs_with_multi_charges_trimmed_sanH).
      value_counts().to_string())
print("Number of charges and count with H:")
print(pd.Series(segs_with_multi_charges_trimmed).value_counts().to_string())
print()
print("Charged amino acids in trans segments:".ljust(70),
      len(re.findall('[' + charged + ']', seg_aas)),
      '(', len(re.findall('[' + chargedplus + ']', seg_aas)), ')')
print(pd.Series(re.findall('[' + chargedplus + ']',
                seg_aas)).value_counts().to_string())
print("Charged amino acids in trans segments, trimmed:".ljust(70),
      len(re.findall('[' + charged + ']', seg_aas_trimmed)),
      '(', len(re.findall('[' + chargedplus + ']', seg_aas_trimmed)), ')')
print(pd.Series(re.findall('[' + chargedplus + ']',
      seg_aas_trimmed)).value_counts().to_string())
print()
print("##################################################")
print("# Everything below is using the trimmed segments #")
print("##################################################")
print()
#######################
# Proteins with pairs #
#######################
print("Proteins with charged amino acids(pairs) in pos:")
for i in range(1, 11):
    print("i+{}: ".format(i),
          len(prot_charged_gaps_trimmed_sanH[i-1]),
          '(', len(prot_charged_gaps_trimmed[i-1]), ')')
print()
print("Proteins with opposite charged pairs:")
for i in range(1, 11):
    print("i+{}: ".format(i),
          len(prot_opp_gaps_trimmed_sanH[i-1]),
          '(', len(prot_opp_gaps_trimmed[i-1]), ')')
print()
print("Proteins with same charged pairs:")
for i in range(1, 11):
    print("i+{}: ".format(i),
          len(prot_same_gaps_trimmed_sanH[i-1]),
          '(', len(prot_same_gaps_trimmed[i-1]), ')')
print()
#######################
# Segments with pairs #
#######################
print("Segments with charged amino acids(pairs) in pos:")
for i in range(1, 11):
    print("i+{}: ".format(i),
          len(segs_charged_gaps_trimmed_sanH[i-1]),
          '(', len(segs_charged_gaps_trimmed[i-1]), ')')
print()
print("Segments with opposite charged pairs:")
for i in range(1, 11):
    print("i+{}: ".format(i),
          len(segs_opp_gaps_trimmed_sanH[i-1]),
          '(', len(segs_opp_gaps_trimmed[i-1]), ')')
print()
print("Segments with same charged pairs:")
for i in range(1, 11):
    print("i+{}: ".format(i),
          len(segs_same_gaps_trimmed_sanH[i-1]),
          '(', len(segs_same_gaps_trimmed[i-1]), ')')
print()
####################################
# Charged amino acids in positions #
####################################
print("Charged amino acids(pairs) in pos:")
for i in range(1, 11):
    print("i+{}: ".format(i),
          len(charged_gaps_trimmed_sanH[i-1]),
          '(', len(charged_gaps_trimmed[i-1]), ')')
    print(pd.Series(charged_gaps_trimmed[i-1]).value_counts().to_string())
print()
print("Opposite charged pairs:")
for i in range(1, 11):
    print("i+{}: ".format(i),
          len(opp_gaps_trimmed_sanH[i-1]),
          '(', len(opp_gaps_trimmed[i-1]), ')')
    print(pd.Series(opp_gaps_trimmed[i-1]).value_counts().to_string())
print()
print("Same charged pairs:")
for i in range(1, 11):
    print("i+{}: ".format(i),
          len(same_gaps_trimmed_sanH[i-1]),
          '(', len(same_gaps_trimmed[i-1]), ')')
    print(pd.Series(same_gaps_trimmed[i-1]).value_counts().to_string())
# print("i+2: ", len(''.join(charged_gaps_trimmed[1]).replace('H', '')), '(', len(''.join(charged_gaps_trimmed[1])), ')')
# print(pd.Series(charged_gaps_trimmed[1]).value_counts().to_string())
# print("i+3: ", len(''.join(charged_gaps_trimmed[2]).replace('H', '')), '(', len(''.join(charged_gaps_trimmed[2])), ')')
# print(pd.Series(charged_gaps_trimmed[2]).value_counts().to_string())
# print("i+4: ", len(''.join(charged_gaps_trimmed[3]).replace('H', '')), '(', len(''.join(charged_gaps_trimmed[3])), ')')
# print(pd.Series(charged_gaps_trimmed[3]).value_counts().to_string())
# print("i+5: ", len(''.join(charged_gaps_trimmed[4]).replace('H', '')), '(', len(''.join(charged_gaps_trimmed[4])), ')')
# print(pd.Series(charged_gaps_trimmed[4]).value_counts().to_string())
# print("i+6: ", len(''.join(charged_gaps_trimmed[5]).replace('H', '')), '(', len(''.join(charged_gaps_trimmed[5])), ')')
# print(pd.Series(charged_gaps_trimmed[5]).value_counts().to_string())
# print("i+7: ", len(''.join(charged_gaps_trimmed[6]).replace('H', '')), '(', len(''.join(charged_gaps_trimmed[6])), ')')
# print(pd.Series(charged_gaps_trimmed[6]).value_counts().to_string())
print()

# ### Opposite charged pairs ###
# print("Opposite charged pairs")
# print("i+1: ", len(opposite_gaps_trimmed_noH[0]), '(', len(opposite_gaps_trimmed[0]), ')')
# print(pd.Series(opposite_gaps_trimmed[0]).value_counts().to_string())
# print("i+2: ", len(opposite_gaps_trimmed_noH[1]), '(', len(opposite_gaps_trimmed[1]), ')')
# print(pd.Series(opposite_gaps_trimmed[1]).value_counts().to_string())
# print("i+3: ", len(opposite_gaps_trimmed_noH[2]), '(', len(opposite_gaps_trimmed[2]), ')')
# print(pd.Series(opposite_gaps_trimmed[2]).value_counts().to_string())
# print("i+4: ", len(opposite_gaps_trimmed_noH[3]), '(', len(opposite_gaps_trimmed[3]), ')')
# print(pd.Series(opposite_gaps_trimmed[3]).value_counts().to_string())
# print("i+5: ", len(opposite_gaps_trimmed_noH[4]), '(', len(opposite_gaps_trimmed[4]), ')')
# print(pd.Series(opposite_gaps_trimmed[4]).value_counts().to_string())
# print("i+6: ", len(opposite_gaps_trimmed_noH[5]), '(', len(opposite_gaps_trimmed[5]), ')')
# print(pd.Series(opposite_gaps_trimmed[5]).value_counts().to_string())
# print("i+7: ", len(opposite_gaps_trimmed_noH[6]), '(', len(opposite_gaps_trimmed[6]), ')')
# print(pd.Series(opposite_gaps_trimmed[6]).value_counts().to_string())
# 

# flows = [# len(mems_aa_trimmed),
#          len(re.findall('[' + charged + ']', mems_aa_trimmed)),
#          len(''.join(charge_gaps_trimmed[0]).replace('H', '')) +
#          len(''.join(charge_gaps_trimmed[2]).replace('H', '')) +
#          len(''.join(charge_gaps_trimmed[3]).replace('H', '')),
#          len(opposite_gaps_trimmed_noH[0]) +
#          len(opposite_gaps_trimmed_noH[2]) +
#          len(opposite_gaps_trimmed_noH[3])]
# labels = [# "Amino acids in transmemgrane segments",
#           "Charged amino acids in segments",
#           "Charged pairs in position 1,3 and 4",
#           "Opposite charged pairs in position 1, 3 and 4"]
# colors = [# "#FF4000",
#           "#FF8000", "#FFBF00", "#FFFF00"]
# # print(flows)
# fig = plt.figure()#figsize=(8, 14))
# ax = fig.add_subplot(1, 1, 1, xticks=[], yticks=[])
# sankey = Sankey(ax=ax, scale=0.03, offset=-.1)
# for input_numbers, output_numbers, label, prior, color in zip(flows[:-1],
#                                                               flows[1:],
#                                                               labels,
#                                                               [None, 0],
#                                                               colors):
#     # print(prior)
#     if prior != 0:
#         sankey.add(flows=[input_numbers, -output_numbers,
#                           output_numbers - input_numbers],
#                    orientations=[0, 0, 1],
#                    patchlabel=label,
#                    labels=['', None, 'Excluded'],
#                    prior=prior,
#                    connect=(1, 0),
#                    pathlengths=[0, 0, 2],
#                    trunklength=5.,
#                    rotation=-90,
#                    facecolor=color)
#     else:
#         sankey.add(flows=[input_numbers, -output_numbers,
#                           output_numbers - input_numbers],
#                    orientations=[0, 0, 1],
#                    patchlabel=label,
#                    labels=['', labels[-1], 'Excluded'],
#                    prior=prior,
#                    connect=(1, 0),
#                    pathlengths=[0, 0, 10],
#                    trunklength=5.,
#                    rotation=-90,
#                    facecolor=color)
#          
# diagrams = sankey.finish()
# for diagram in diagrams:
#     diagram.text.set_fontweight('bold')
#     diagram.text.set_fontsize('10')
#     for text in diagram.texts:
#         text.set_fontsize('10')
# ylim = plt.ylim()
# plt.ylim(ylim[0]*1.05, ylim[1])
# plt.title("Charge # flows")
# plt.tight_layout()
# plt.savefig("sankey_charges.png")
# # plt.show()
