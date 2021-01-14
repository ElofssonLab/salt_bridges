#!/usr/bin/env python3
import argparse
import pickle
import pandas as pd
import sys
from collections import Counter

parser = argparse.ArgumentParser()

parser.add_argument("pair_csv", type=str, help="Pairs csv file")
parser.add_argument("mems", type=str, help="Membrane pickle")
parser.add_argument("bridges", type=str, help="Bridge pickle")

args = parser.parse_args()

pair_data = pd.read_csv(args.pair_csv, delimiter=',')
# No Histidine
pair_data["R_type"] = pair_data["Res1"] + pair_data["Res2"]
pair_data_sanH = pair_data[((pair_data["Res1"] != 'H') & (pair_data["Res2"] != 'H'))]    
pair_data_sanH_salt = pair_data_sanH[pair_data_sanH["Local saltbridge"].notna()]
pair_data_salt = pair_data[pair_data["Local saltbridge"].notna()]
pair_data_sanH_chr_types = pair_data_sanH[["Step","R_type"]]
pair_data_chr_types = pair_data[["Step","R_type"]]
pair_data_sanH_bridge_types = pair_data_sanH_salt[["Step","R_type"]]
pair_data_bridge_types = pair_data_salt[["Step","R_type"]]

pair_data_sanH_chr_count = pair_data_sanH_chr_types.value_counts(subset=["Step", "R_type"])
pair_data_chr_count = pair_data_chr_types.value_counts(subset=["Step", "R_type"])
pair_data_sanH_bridge_count = pair_data_sanH_bridge_types.value_counts(subset=["Step", "R_type"])
pair_data_bridge_count = pair_data_bridge_types.value_counts(subset=["Step", "R_type"])
# print(pair_data_bridge_count)
# sys.exit()

charges = {}
bridges = {}
pairs = ["DK", "KD", "ER", "RE", "RK", "KR", "EK", "KE", "DR", "RD"]

print("Charges")
for row in pair_data_sanH_chr_count.iteritems():
    charges[str(row[0][0]) + "-" + row[0][1]] = row[1]
    if row[0][0] in [1,3,4,6] and row[0][1] in pairs:  # + ["RR"]:  #["DK", "EK", "RR"]:
        print(row)
print("Charges RR")
for row in pair_data_sanH_chr_count.iteritems():
    charges[str(row[0][0]) + "-" + row[0][1]] = row[1]
    if row[0][1] == "RR":  # + ["RR"]:  #["DK", "EK", "RR"]:
        print(row)
print("Charges with H")
for row in pair_data_chr_count.iteritems():
    # charges[str(row[0][0]) + "-" + row[0][1]] = row[1]
    if row[0][0] in [1,3,4,6] and "H" in row[0][1]:
        print(row)
print("Bridges")
for row in pair_data_sanH_bridge_count.iteritems():
    bridges[str(row[0][0]) + "-" + row[0][1]] = row[1]
    if row[0][0] in [1,3,4]:  # and row[0][1] in ["DK", "EK", "RR"]:
        print(row)
# # print(pair_data_chr_count["Step"])
# print(charges)
# print(bridges)

print("Fractions")
for i in ["1", "3", "4"]:
    for p in pairs:
    # for p in ["DK", "EK"]:
        if i + "-" + p in bridges:
            frac = bridges[i + "-" + p] / charges[i + "-" + p]
            print(i, p, "{:.2%}".format(frac))

RR_mids = []
RR_bridges = []
pickle_bridges = pickle.load(open(args.bridges,
                            'rb'))

print("RR at 6")
mems = pickle.load(open(args.mems, 'rb'))
for key, membranes in mems.items():
    local_bridges = []
    mems_bridges = []
    if key in pickle_bridges['local']:
        local_bridges = pickle_bridges['local'][key]
    if key in pickle_bridges['mems']:
        mems_bridges = pickle_bridges['mems'][key]
    # if len(local_bridges) > 0:
    #     print(key)
    #     print(local_bridges)
    #     print(membranes)
    # continue
            
    for mem_place, mem in membranes:

        if len(mem) < 17:
            continue
        edge = 5
        midMem = mem[edge:-edge]
        fullMem = mem
        # Only use mems with all known amino acids
        if '?' in fullMem:
            continue
        memLen = len(midMem)
        fullmemLen = len(fullMem)

        mem_bridge = []
        local_bridge = []
        for mb in mems_bridges:
            ss_first = mb[0]
            ss_second = mb[2]
            if ss_first >= (mem_place+1) and ss_first <= (mem_place + 1 + memLen) or\
                    ss_second >= (mem_place + 1) and ss_second <= (mem_place + 1 + memLen):
                        mem_bridge.append(mb)
        # for lb in local_bridges:
        #     ss_first = lb[0]
        #     ss_second = lb[2]
        #     if ss_first > (mem_place+5) and ss_first <= (mem_place + fullmemLen - 5) and\
        #             ss_second > (mem_place+5) and ss_second <= (mem_place + fullmemLen - 5):
        #                 local_bridge.append(lb)
        for place, first_aa in enumerate(midMem):
            if first_aa != 'R':
                continue
            if memLen <= place + 6:
                break
            second_aa = midMem[place + 6]
            if second_aa != 'R':
                continue
            # print(midMem, mem_place, place)        
            pot_places = [mem_place + place + 5 + 1, mem_place + place + 8 + 1, mem_place + place + 11 + 1]
            # print(pot_places)
            # print(mem_bridge)
            # print(key, mem_place, place, midMem[place:place + 6 + 1])
            # print(midMem)
            # print(mem)
            m_bridges = 0
            for mb in mem_bridge:
                ss_first = mb[0]
                ss_second = mb[2]
                if ss_first in pot_places or ss_second in pot_places:
                    # print("Has bridge", mb)
                    m_bridges += 1
            RR_mids.append(midMem[place+3])
            RR_bridges.append([midMem[place+3],m_bridges])

            # print(midMem[place+3])
            # print(mem_bridge)
            # sys.exit()
# print(RR_mids)
print(Counter(RR_mids))
print(RR_bridges)
            
        # if len(mem_bridge) > 0:     
        #     print(mem_bridge)
        #     print(local_bridge)
        #     print(mem)
        #     sys.exit()
