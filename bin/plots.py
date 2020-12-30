#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
from collections import Counter

CHARGEDwithH = "DEKRH"
CHARGED = "DEKR"
POLAR = "DEKRHNPQ"
POS = "KR"
NEG = "DE"

aa_data = pd.read_csv("pdbtm.clust.aas.csv",delimiter=',')
pair_data = pd.read_csv("pdbtm.clust.pairs.csv",delimiter=',')

num_pairs = pair_data["PID"].count()
num_pairs_noH = pair_data[(pair_data["Res1"] != 'H') & (pair_data["Res2"] != 'H') & (pair_data["Step"] < 7)]["Local saltbridge"].count()
# num_pairs_noH = pair_data[(pair_data["Res1"] != 'H') & (pair_data["Res2"] != 'H') & (pair_data["Local saltbridge"].notna())].groupby("Step")["PID"].count()
num_pairs_noH = pair_data[(pair_data["Res1"] != 'H') & (pair_data["Res2"] != 'H') & (pair_data["Local saltbridge"].notna()) & (pair_data["Step"] < 7)]["PID"].count()

print("Pairs: {}".format(num_pairs))
print("Local saltbridges: {}".format(num_pairs_noH))
num_residues = aa_data['PID'].count()
print("Number of residues: {}".format(num_residues))
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
plt.savefig("bridge_scatter.png")
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
plt.savefig("charged_aa_pairings.png")
fig.clear()
# plt.show()
print("Number of polarresidues: {} ({:.1%})".format(polar_aa_data['PID'].count(), polar_aa_data['PID'].count()/num_residues))
print("Number of charged residues: {} ({:.1%})".format(charged_aa_data['PID'].count(),charged_aa_data['PID'].count()/num_residues))

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
for k in has_local:
    # print(k)
    row = opp_charge_sanH_local_134.loc[k] 
    for l in row["localbridge"].split(';'):
        if l[0] != 'H' and int(l[1:])-row["ResN"] < 7:
            gaps.append(int(l[1:])-row["ResN"])  
    
# print(gaps)
# print(Counter(gaps))
#### Figure for saltbridge gap ####
fig, ax = plt.subplots()
ax.hist(gaps, bins=[-4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5], align="mid", rwidth=0.9) 
ax.set_xticks([-4,-3,-2,-1,0,1,2,3,4])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel("Salt bridge separation")
ax.set_ylabel("Count")
# plt.show()
plt.savefig("saltbridge_gap.png")
################################

print("Number of charged residues with opp 3 and 4: {}".format(opp_charge34["PID"].count()))
print("Number of charged residues with same 3 and 4: {}".format(same_charge34["PID"].count()))
print("Number of charged residues in pair(1,3,4): {} ({:.1%})".format(same_charge134["PID"].count()+opp_charge134["PID"].count(),(same_charge134["PID"].count()+opp_charge134["PID"].count())/num_residues))
print("Number of charged residues with opp 1, 3 and 4: {} ({:.1%})".format(opp_charge134["PID"].count(),opp_charge134["PID"].count()/num_residues))
print("Number of charged residues with same 1, 3 and 4: {} ({:.1%})".format(same_charge134["PID"].count(),same_charge134["PID"].count()/num_residues))
print("34 with saltbridges: {}".format(opp_charge_sanH["PID"].count()))
print("34 with local saltbridges: {}".format(opp_charge_sanH_local["PID"].count()))
print("134 with saltbridges: {}".format(opp_charge_sanH_134["PID"].count()))
print("134 with local saltbridges: {}".format(opp_charge_sanH_local_134["PID"].count()))
# print(pair_data[(pair_data["Res1"] != 'H') & (pair_data["Res2"] != 'H') & (pair_data["Local saltbridge"].notna())])
# print(opp_charge_sanH_local_134.groupby("PID")["PID"].count())
# print(charged_aa_data[charged_aa_data["PID"] == "6Y7FA"])
