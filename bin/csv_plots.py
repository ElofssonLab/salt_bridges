#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys
import argparse
import matplotlib.pyplot as plt
from collections import Counter

parser = argparse.ArgumentParser()

parser.add_argument("aa_csv", type=str, help="AA csv file")
parser.add_argument("pair_csv", type=str, help="Pairs csv file")

args = parser.parse_args()
CHARGEDwithH = "DEKRH"
CHARGED = "DEKR"
POLAR = "DEKRHNPQ"
POS = "KR"
NEG = "DE"

aa_data = pd.read_csv(args.aa_csv, delimiter=',')
pair_data = pd.read_csv(args.pair_csv, delimiter=',')
prefix = args.aa_csv.split('/')[-1].split('.')[0]
stats_data = []
stats_file = "stats/" + prefix + "_csv_stats.txt"
# print(prefix)
# sys.exit()
num_pairs = pair_data["PID"].count()
num_pairs_noH_local = pair_data[(pair_data["Res1"] != 'H') & (pair_data["Res2"] != 'H') & (pair_data["Step"] < 7) & (pair_data["Pair type"]=="Opp") & (pair_data["Local saltbridge"].notna())]
# print(num_pairs_noH_local)

pair_gap = []
num_pair_gap = 0
for k, row in num_pairs_noH_local.iterrows():
    for s in row["Local saltbridge"].split(';'):
        if row["Step"] != 5:
            pair_gap.append(row["Step"])
            num_pair_gap += 1
pair_counter = Counter(pair_gap)
fig, ax = plt.subplots()
# print(np.sum(np.array(list(pair_counter.values()))/num_pair_gap))
ax.bar(pair_counter.keys(), np.array(list(pair_counter.values()))/num_pair_gap)
# ax.hist(gaps, bins=[-4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5], align="mid", rwidth=0.9) 
ax.set_xticks([0,1,2,3,4])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel("Salt bridge separation")
ax.set_ylabel("Fraction of all local salt bridges")
# plt.show()
plt.savefig("images/{}_pair_bridge_gap.png".format(prefix))
fig.clear()
# sys.exit()
# plt.savefig("saltbridge_gap.png")
# num_pairs_noH_count = pair_data[(pair_data["Res1"] != 'H') & (pair_data["Res2"] != 'H') & (pair_data["Step"] < 7)]["Local saltbridge"].count()
# num_pairs_noH = pair_data[(pair_data["Res1"] != 'H') & (pair_data["Res2"] != 'H') & (pair_data["Local saltbridge"].notna())].groupby("Step")["PID"].count()
num_pairs_noH_count = pair_data[(pair_data["Res1"] != 'H') & (pair_data["Res2"] != 'H') & (pair_data["Local saltbridge"].notna()) & (pair_data["Step"] < 7)]["PID"].count()

stats_data.append("Pairs: {}".format(num_pairs))
stats_data.append("Local saltbridges: {}".format(num_pairs_noH_count))
num_residues = aa_data['PID'].count()
stats_data.append("Number of residues: {}".format(num_residues))
polar_aa_data = aa_data[aa_data["Res"].isin(list(POLAR))]
charged_aa_data = aa_data[aa_data["Res"].isin(list(CHARGED))]

localbridge_data = charged_aa_data[charged_aa_data["localbridge"].notna()]
### Look at salt bridge distribution
saltbridge_place = []
saltbridge_gap = []
for k, row in localbridge_data.iterrows():
    local_place = row["ResN"]-row["ResCoreStart"]
    for s in row["localbridge"].split(';'):
        saltbridge_place.append(local_place)
        saltbridge_gap.append(int(s[1:])-row["ResN"])

fig, ax = plt.subplots()
ax.scatter(saltbridge_place,saltbridge_gap)
ax.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel("Place in core mem")
ax.set_ylabel("Relative bridge gap")
plt.savefig("images/{}_bridge_scatter.png".format(prefix))
# plt.show()
fig.clear()
cols = ["-7","-6","-5","-4","-3","-2","-1","1","2","3","4","5","6","7"]
# print(charged_aa_data)
def opp(row):
    res = row["Res"]
    df = row[cols].apply(lambda x: x if (res in POS and str(x) in NEG) or (res in NEG and str(x) in POS) else np.nan)
    return df

def same(row):
    res = row["Res"]
    df = row[cols].apply(lambda x: x if (res in POS and str(x) in POS) or (res in NEG and str(x) in NEG) else np.nan)
    return df

charged_aa_data_sanH = charged_aa_data.copy(deep=True)
charged_aa_data_sanH[cols] = charged_aa_data_sanH[cols].applymap(lambda x: np.nan if x == 'H' else x)
# print(charged_aa_data_sanH)
charged_aa_data_sanH_opp = charged_aa_data_sanH.copy(deep=True)
charged_aa_data_sanH_same = charged_aa_data_sanH.copy(deep=True)
charged_aa_data_sanH_opp[cols] = charged_aa_data_sanH.apply(opp,axis=1)
charged_aa_data_sanH_same[cols] = charged_aa_data_sanH.apply(same,axis=1)
opp_series = charged_aa_data_sanH_opp[cols].count()
same_series = charged_aa_data_sanH_same[cols].count()
# charged_aa_data_sanH_opp = charged_aa_data_sanH.apply(lambda x: x if x.name in cols and ((charged_aa_data_sanH.loc[x.index]["Res"].values[0] in POS and x in NEG) (charged_aa_data_sanH.loc[x.index]["Res"].values[0] in NEG and x in POS)) else np.nan)
# print(charged_aa_data_sanH_opp)
index = [-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7]
df = pd.DataFrame({"Same":same_series.tolist(), "Opp":opp_series.tolist()}, index=index)

fig, ax = plt.subplots()
df.plot.bar(ax=ax)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel("Pairing distances")
ax.set_ylabel("Count")
plt.savefig("images/{}_charged_aa_pairings.png".format(prefix))
fig.clear()
# plt.show()
stats_data.append("Number of polarresidues: {} ({:.1%})".format(polar_aa_data['PID'].count(), polar_aa_data['PID'].count()/num_residues))
stats_data.append("Number of charged residues: {} ({:.1%})".format(charged_aa_data['PID'].count(),charged_aa_data['PID'].count()/num_residues))

same_charge34 = charged_aa_data[(charged_aa_data["Res"].isin(list(POS)) &
                            (charged_aa_data["-3"].isin(list(POS)) | (charged_aa_data["-4"].isin(list(POS))) |
                            charged_aa_data["3"].isin(list(POS)) | (charged_aa_data["4"].isin(list(POS))))) |
                            (charged_aa_data["Res"].isin(list(NEG)) &
                            (charged_aa_data["-3"].isin(list(NEG)) | (charged_aa_data["-4"].isin(list(NEG))) |
                            charged_aa_data["3"].isin(list(NEG)) | (charged_aa_data["4"].isin(list(NEG)))))]
opp_charge34 = charged_aa_data[(charged_aa_data["Res"].isin(list(POS)) &
                            (charged_aa_data["-3"].isin(list(NEG)) | (charged_aa_data["-4"].isin(list(NEG))) |
                            charged_aa_data["3"].isin(list(NEG)) | (charged_aa_data["4"].isin(list(NEG))))) |
                            (charged_aa_data["Res"].isin(list(NEG)) &
                            (charged_aa_data["-3"].isin(list(POS)) | (charged_aa_data["-4"].isin(list(POS))) |
                            charged_aa_data["3"].isin(list(POS)) | (charged_aa_data["4"].isin(list(POS)))))]
same_charge134 = charged_aa_data[(charged_aa_data["Res"].isin(list(POS)) &
                            (charged_aa_data["-3"].isin(list(POS)) | (charged_aa_data["-4"].isin(list(POS))) |
                            charged_aa_data["3"].isin(list(POS)) | (charged_aa_data["4"].isin(list(POS))) |
                            charged_aa_data["-1"].isin(list(POS)) | (charged_aa_data["1"].isin(list(POS))))) |
                            (charged_aa_data["Res"].isin(list(NEG)) &
                            (charged_aa_data["-3"].isin(list(NEG)) | (charged_aa_data["-4"].isin(list(NEG))) |
                            (charged_aa_data["-1"].isin(list(NEG)) | (charged_aa_data["1"].isin(list(NEG))) |
                            charged_aa_data["3"].isin(list(NEG)) | (charged_aa_data["4"].isin(list(NEG))))))]
opp_charge134 = charged_aa_data[(charged_aa_data["Res"].isin(list(POS)) &
                            (charged_aa_data["-3"].isin(list(NEG)) | (charged_aa_data["-4"].isin(list(NEG))) |
                            charged_aa_data["3"].isin(list(NEG)) | (charged_aa_data["4"].isin(list(NEG))) |
                            charged_aa_data["-1"].isin(list(NEG)) | (charged_aa_data["1"].isin(list(NEG))))) |
                            (charged_aa_data["Res"].isin(list(NEG)) &
                            (charged_aa_data["-3"].isin(list(POS)) | (charged_aa_data["-4"].isin(list(POS))) |
                            (charged_aa_data["-1"].isin(list(POS)) | (charged_aa_data["1"].isin(list(POS))) |
                            charged_aa_data["3"].isin(list(POS)) | (charged_aa_data["4"].isin(list(POS))))))]    
opp_charge_sanH = opp_charge34[(opp_charge34["saltbridge"].notna())]
opp_charge_sanH_local = opp_charge34[(opp_charge34["localbridge"].notna())]
opp_charge_sanH = opp_charge_sanH[opp_charge_sanH["saltbridge"].str.count('H') <= opp_charge_sanH["saltbridge"].str.count(';')]
opp_charge_sanH_local = opp_charge_sanH_local[opp_charge_sanH_local["localbridge"].str.count('H') <= opp_charge_sanH_local["localbridge"].str.count(';')]
opp_charge_sanH_134 = opp_charge134[(opp_charge134["saltbridge"].notna())]
opp_charge_sanH_local_134 = opp_charge134[(opp_charge134["localbridge"].notna())]
opp_charge_sanH_134 = opp_charge_sanH_134[opp_charge_sanH_134["saltbridge"].str.count('H') <= opp_charge_sanH_134["saltbridge"].str.count(';')]
opp_charge_sanH_local_134 = opp_charge_sanH_local_134[opp_charge_sanH_local_134["localbridge"].str.count('H') <= opp_charge_sanH_local_134["localbridge"].str.count(';')]
### Local Saltbridge separation 

has_local = []
for k, r in opp_charge_sanH_local_134["localbridge"].str.split(";").iteritems():
    has_local.append(k)

gaps = []
num_local_bridges = 0
for k in has_local:
    # print(k)
    row = opp_charge_sanH_local_134.loc[k] 
    for l in row["localbridge"].split(';'):
        if l[0] != 'H' and int(l[1:])-row["ResN"] < 7:
            gaps.append(int(l[1:])-row["ResN"]) 
            num_local_bridges += 1
gap_counts = Counter(gaps)
# print(gaps)
# print(Counter(gaps))
#### Figure for saltbridge gap ####
fig, ax = plt.subplots()
# print(np.sum(np.array(list(gap_counts.values()))/num_local_bridges))
ax.bar(gap_counts.keys(), np.array(list(gap_counts.values()))/num_local_bridges)
# ax.hist(gaps, bins=[-4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5], align="mid", rwidth=0.9) 
ax.set_xticks([-4,-3,-2,-1, 0,1,2,3,4])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel("Salt bridge separation")
ax.set_ylabel("Fraction of all local salt bridges")
# plt.show()
# sys.exit()
plt.savefig("images/{}_saltbridge_gap.png".format(prefix))
################################

stats_data.append("Number of charged residues with opp 3 and 4: {}".format(opp_charge34["PID"].count()))
stats_data.append("Number of charged residues with same 3 and 4: {}".format(same_charge34["PID"].count()))
stats_data.append("Number of charged residues in pair(1,3,4): {} ({:.1%})".format(same_charge134["PID"].count()+opp_charge134["PID"].count(),(same_charge134["PID"].count()+opp_charge134["PID"].count())/num_residues))
stats_data.append("Number of charged residues with opp 1, 3 and 4: {} ({:.1%})".format(opp_charge134["PID"].count(),opp_charge134["PID"].count()/num_residues))
stats_data.append("Number of charged residues with same 1, 3 and 4: {} ({:.1%})".format(same_charge134["PID"].count(),same_charge134["PID"].count()/num_residues))
stats_data.append("34 with saltbridges: {}".format(opp_charge_sanH["PID"].count()))
stats_data.append("34 with local saltbridges: {}".format(opp_charge_sanH_local["PID"].count()))
stats_data.append("134 with saltbridges: {}".format(opp_charge_sanH_134["PID"].count()))
stats_data.append("134 with local saltbridges: {}".format(opp_charge_sanH_local_134["PID"].count()))
with open(stats_file, 'w') as stats_handle:
    stats_handle.write('\n'.join(stats_data))
# print(pair_data[(pair_data["Res1"] != 'H') & (pair_data["Res2"] != 'H') & (pair_data["Local saltbridge"].notna())])
# print(opp_charge_sanH_local_134.groupby("PID")["PID"].count())
# print(charged_aa_data[charged_aa_data["PID"] == "6Y7FA"])
