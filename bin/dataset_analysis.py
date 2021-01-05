#!/usr/bin/env python3
"""Membrane charge caluclations."""
# import matplotlib.pyplot as plt
import numpy as np
# import math
import pickle
import dgCalc
import saltBridges
# import parser
import os.path
import argparse
import urllib.request
from urllib.error import HTTPError
import sys

parser = argparse.ArgumentParser()

parser.add_argument("mems_pickle", type=str, help="Pickle file of membranes")
parser.add_argument("bridge_pickle", type=str, help="Pickle file of bridges")
parser.add_argument("threeline", type=str, help="3line in file")

# parser.add_argument("-b", "--bridges", type=int, default=1, help="Required connections per bridge")
# parser.add_argument("-t", "--tolerant", type=bool, default=False, help="Use tolerant membranes (also include m, not just M)")

args = parser.parse_args()
mem_length = 17
aas = 'ACDEFGHIKLMNPQRSTVWY'
CHARGED = "DEKRH"
POS = "KRH"
NEG = "DE"
aaHits = np.zeros([10, 20, 20])
aaPairs = np.zeros(10)
# Dim 1, first and second in pair
aaCount = np.zeros([2, 10, 20])

stats = {'aas':0, 'chr':0, 'pos':0, 'neg':0, 'pospair':0, 'negpair':0, 'chargedpair':0, 'opppair':0, 'polar':0}
# con_req = args.bridges  # How many saltbridge connections to make a bridge? Default 1
aaMap = {'ARG': 'R',
         'HIS': 'H',
         'LYS': 'K',
         'ASP': 'D',
         'GLU': 'E',
         'SER': 'S',
         'THR': 'T',
         'ASN': 'N',
         'GLN': 'Q',
         'CYS': 'C',
         'SEC': 'U',
         'GLY': 'G',
         'PRO': 'P',
         'ALA': 'A',
         'VAL': 'V',
         'ILE': 'I',
         'LEU': 'L',
         'MET': 'M',
         'PHE': 'F',
         'TYR': 'Y',
         'TRP': 'W'}

pdbURL = "https://files.rcsb.org/download/"
positive = 'KR'
negative = 'DE'
charged = 'KRDE'
chargedplus = 'KRDEH'
positiveplus = 'RHK'
helicies = pickle.load(open(args.mems_pickle,
                            'rb'))
bridges = pickle.load(open(args.bridge_pickle,
                            'rb'))
charges_file = args.mems_pickle.replace("mems", "charges")
mems_index = args.mems_pickle.index("mems")
csv_file_pairs = args.mems_pickle[:mems_index] + "pairs.csv"
csv_file_aas = args.mems_pickle[:mems_index] + "aas.csv"

# charged_pair_res = set()
# 
# def countMems(topoStr):
#     num = 0
#     currTopo = topoStr[0]
#     # print(currTopo)
#     for _, topo in enumerate(topoStr):
#         # print(currTopo, topo)
#         if topo != currTopo:
#             if currTopo in 'Mm':
#                 # print("In loop")
#                 num += 1
#             currTopo = topo
#     return num
# 
# 
# def whatMem(topoStr, aaIndex, tol=False):
#     if tol:
#         mem_letters = 'Mm'
#     else:
#         mem_letters = 'M'
#     memNum = 1
#     currTopo = topoStr[0]
#     for i, topo in enumerate(topoStr):
#         if topo != currTopo:
#             if currTopo in mem_letters:
#                 if i >= aaIndex:
#                     return memNum
#                 memNum += 1
#             currTopo = topo
#     return "Can't determine membrane number"
# 
# 
# aas = 'ACDEFGHIKLMNPQRSTVWY'
# polaraas="DEKRHNPQ"
# proteinsExamples = {}
# proteinSet = set()
# fullproteinSet = set()
# longproteinSet = set()
# totalMems = 0
# longMems = 0
# correctMems = 0
# 
TMdata = {}
with open(args.threeline, 'r') as TMHandle:
    pdb_id = ''
    rowNum = 0
    fa = ''
    for line in TMHandle:
        if line[0] == '>':
            # pid = line[1:].strip().split('|')[1]
            if '|' in line:
                pdb_id = line[1:].strip().split('|')[1]
            else:
                pdb_id = line[1:].strip()
            # pid = line[1:].strip()# .split('|')[2]
        elif rowNum % 3 == 1:
            fa = line.strip()
        elif rowNum % 3 == 2:
            if not "MISSING" in pdb_id:
                TMdata[pdb_id] = [fa, line.strip()]
        else:
            print("You should not be here...")
        rowNum += 1


proteins_with_mems= set()
proteins_with_correct_mems = set()
local_bridge_list = set()
any_saltbridge_prot = set()
true_local_saltbridge_prot = set()
correct_mems = 0
mem_bridge_count = 0
local_bridge_count = 0
tot_num_mem_bridges = 0
tot_num_mem_proteins = 0
tot_num_local_bridges = 0
tot_num_local_proteins = 0
csv_text_pairs = ["PID,Pair type,Core start,Core end,ResN1,ResN2," +
            "Step,Res1,Res2,Any saltbridge,Local saltbridge"]
csv_text_aas = ["PID,ResN,Res,ResCoreStart,ResCoreEnd,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,saltbridge,localbridge"]
for k, m in bridges["mems"].items():
    tot_num_mem_proteins += 1
    tot_num_mem_bridges += len(m) 
for k, m in bridges["local"].items():
    tot_num_local_proteins += 1
    tot_num_local_bridges += len(m) 
print("Membrane bridges: ", tot_num_mem_bridges)
print("Membrane proteins: ", tot_num_mem_proteins)
print("Local bridges: ", tot_num_local_bridges)
print("Local proteins: ", tot_num_local_proteins)
for key, membranes in helicies.items():
    local_bridges = []
    mems_bridges = []
    proteins_with_mems.add(key)
    if key in bridges['local']:
        local_bridges = bridges['local'][key]
    if key in bridges['mems']:
        mems_bridges = bridges['mems'][key]
    full_seq = TMdata[key][0]
    for mem_place, mem in membranes:
        # Only use mems of length 17 or more
        # print(mem_place, mem)
        if len(mem) < mem_length:
            continue
        # Cut the 5 caps off on each side
        edge = 5
        midMem = mem[edge:-edge]
        fullMem = mem
        #  ext_mem = full_seq[mem_place - 2: mem_place + len(mem) + 2]
        # print(ext_mem)
        # print('  ' + fullMem + '  ')
        # print('       ' + midMem + '       ')
        # sys.exit()
        # Only use mems with all known amino acids
        if '?' in fullMem:
            continue
        proteins_with_correct_mems.add(key)
        correct_mems += 1
        memLen = len(midMem)
        fullmemLen = len(fullMem)
        mem_bridge = []
        local_bridge = []
        for mb in mems_bridges:
            ss_first = mb[0]
            ss_second = mb[2]
            if ss_first > (mem_place+5) and ss_first <= (mem_place + fullmemLen-5) or\
                    ss_second > (mem_place+5) and ss_second <= (mem_place + fullmemLen-5):
                        any_saltbridge_prot.add(key)  # Make sure to add for any protein
                        mem_bridge.append(mb)
        for lb in local_bridges:
            ss_first = lb[0]
            ss_second = lb[2]
            if ss_first > (mem_place+5) and ss_first <= (mem_place + fullmemLen - 5) and\
                    ss_second > (mem_place+5) and ss_second <= (mem_place + fullmemLen - 5):
                        any_saltbridge_prot.add(key)  # Make sure to add for any protein
                        true_local_saltbridge_prot.add(key)  # Proteins with local salt bridge
                        local_bridge.append(lb)
        # Run through each mem and save the pairs
        for place, aa in enumerate(midMem):
            # Charge CSV-file
            if aa in CHARGED:
                aa_place = mem_place + 1 + 5 + place  # Residue numbering
                charge_csv_line = [key, str(aa_place), aa, str(mem_place+1+5),str(mem_place+len(fullMem)-5)]
                # print("       " + aa + "       ")
                #  ext_mem = full_seq[mem_place - 2: mem_place + len(mem) + 2]
                # ext_mem = full_seq[mem_place - 2 + place: mem_place + len(midMem) + 5]
                ext_mem = full_seq[aa_place-7-1:aa_place + 7]  # Remember, aa_place is residue numbering
                # print("       " + aa + "       |")
                # print(ext_mem, len(ext_mem))
                # sys.exit()
                # Check all pairings -7 to 7 incl.
                # if key == "6Y7FA":
                #     print(ext_mem)
                for i in range(-7,8):
                    if i == 0:
                        pass
                    else:
                        if not i >= len(ext_mem)-7:
                            # if key == "6Y7FA":
                            #     print(aa_place, i, ext_mem[i+7])
                            paired_aa = ext_mem[i+7]
                            if paired_aa in CHARGED:
                                charge_csv_line.append(paired_aa)
                            else:
                                charge_csv_line.append("")
                        else:
                            charge_csv_line.append("")
                # Check bridges
                charge_mem_bridge_text = ''
                charge_local_bridge_text = ''
                # stop = False
                if len(mem_bridge)>0:
                    for mb in mem_bridge:
                        ss_first = mb[0]
                        ss_second = mb[2]
                        if ss_first == aa_place or ss_second == aa_place:
                            if ss_first == aa_place:
                                paired_res = ss_second
                                matched_res = aaMap[mb[3]]
                            else:
                                paired_res = ss_first
                                matched_res = aaMap[mb[1]]

                            if len(charge_mem_bridge_text) > 0 and abs(ss_second - ss_first) > 7:
                                charge_mem_bridge_text += ";"+matched_res+str(paired_res)
                            elif abs(ss_second - ss_first) > 7:
                                # print(key, paired_res, full_seq, len(full_seq))
                                charge_mem_bridge_text = matched_res+str(paired_res)
                            elif abs(ss_second - ss_first) < 8:
                                if len(charge_local_bridge_text) > 0:
                                    charge_local_bridge_text += ';' + matched_res+str(paired_res)
                                else:
                                    charge_local_bridge_text = matched_res+str(paired_res)
                            else:
                                print("something wrong", key, mem_bridge)
                                sys.exit()
                                # stop = True
                charge_csv_line.append(charge_mem_bridge_text)
                charge_csv_line.append(charge_local_bridge_text)
            else:
                aa_place = mem_place + 1 + 5 + place  # Residue numbering
                charge_csv_line = [key, str(aa_place), aa, str(mem_place+1+5),str(mem_place+len(fullMem)-5)]
                # print("       " + aa + "       ")
                ext_mem = full_seq[mem_place - 2 + place: mem_place + 5 + 7 + 1]
                # print(ext_mem, len(ext_mem))
                # Check all pairings -7 to 7 incl.
                for i in range(-7,8):
                    if i == 0:
                        pass
                    else:
                        if not i >= len(ext_mem)-7:
                            paired_aa = ext_mem[i+7]
                            if paired_aa in CHARGED:
                                charge_csv_line.append(paired_aa)
                            else:
                                charge_csv_line.append("")
                        else:
                            charge_csv_line.append("")

                charge_csv_line.append('')  # mem bridge
                charge_csv_line.append('')  # local bridge



                # if len(local_bridge)>0:
                #     for lb in local_bridge:
                #         ss_first = lb[0]
                #         ss_second = lb[2]
                #         if ss_first == (mem_place+5+1+place) and ss_second == (mem_place + 5 + 1 + place + i):
                #             has_local_bridge = 1
                #             local_bridges_text = str(ss_first)+"-"+str(ss_second)
                # if stop:
            csv_text_aas.append(','.join(charge_csv_line))

                #     sys.exit()
            # Pair CSV-file
            # Now count the for each gap
            for i in range(1, 11):
                if memLen <= place + i:
                    break
                if aa in aas:
                    first = aas.index(aa)
                    second_aa = midMem[place + i]
                    if second_aa in aas:
                        aaPairs[i - 1] += 1
                        second = aas.index(midMem[place + i])
                        aaHits[i - 1][first][second] += 1
                        # Count both first and second AA for each gap sized i
                        aaCount[0][i - 1][first] += 1
                        aaCount[1][i - 1][second] += 1
                        
                        # Time for the charged pairs
                        if aa in CHARGED and second_aa in CHARGED:
                            has_mem_bridge = 0
                            mem_bridges_text = ""
                            local_bridges_text = ""
                            has_mem_bridge = 0
                            has_local_bridge = 0
                            if (aa in POS and second_aa in POS) or (aa in NEG and second_aa in NEG):
                                # Same charge
                                pair_type = "Same"
                                if len(mem_bridge)>0:
                                    for mb in mem_bridge:
                                        ss_first = mb[0]
                                        aa_first = aaMap[mb[1]]
                                        ss_second = mb[2]
                                        aa_second = aaMap[mb[3]]
                                        if ss_first == (mem_place+5+1+place) or ss_first == (mem_place + 5 + 1 + place + i) or\
                                                ss_second == (mem_place+5+1+place) or ss_second == (mem_place + 5 + 1 + place + i) :
                                                    has_mem_bridge += 1
                                                    if len(mem_bridges_text) > 0:
                                                        mem_bridges_text += ";"+aa_first+str(ss_first)+"-"+aa_second+str(ss_second)
                                                    else:
                                                        mem_bridges_text = aa_first+str(ss_first)+"-"+aa_second+str(ss_second)
                                mem_bridge_count += has_mem_bridge
                                # if len(local_bridge)>0:
                                #     for lb in local_bridge:
                                #         ss_first = lb[0]
                                #         ss_second = lb[2]
                                #         if ss_first == (mem_place+5+1+place) and ss_second == (mem_place + 5 + 1 + place + i):
                                #             has_local_bridge = 1
                                #             local_bridges_text = str(ss_first)+"-"+str(ss_second)
                                #         
                                # # print(key, "SAME", i, has_mem_bridge, has_local_bridge)
                                # if has_local_bridge:
                                #     # Should be impossible to have a local bridge between same charges
                                #     print("SAME local saltbridge?", local_bridges_text)
                                #     print(aa, second_aa)
                                #     print(ss_first, ss_second)
                                #     print(aa_first, aa_second)
                                # local_bridge_count += has_local_bridge
                            else:
                                # Opp charge
                                pair_type = "Opp"
                                if len(mem_bridge)>0:
                                    for mb in mem_bridge:
                                        ss_first = mb[0]
                                        aa_first = aaMap[mb[1]]
                                        ss_second = mb[2]
                                        aa_second = aaMap[mb[3]]
                                        if ss_first == (mem_place+5+1+place) or ss_first == (mem_place + 5 + 1 + place + i) or\
                                                ss_second == (mem_place+5+1+place) or ss_second == (mem_place + 5 + 1 + place + i) :
                                                    has_mem_bridge += 1
                                                    if len(mem_bridges_text) > 0:
                                                        mem_bridges_text += ";"+aa_first+str(ss_first)+"-"+aa_second+str(ss_second)
                                                    else:
                                                        mem_bridges_text = aa_first+str(ss_first)+"-"+aa_second+str(ss_second)
                                has_local_bridge = 0
                                if len(local_bridge)>0:
                                    for lb in local_bridge:
                                        ss_first = lb[0]
                                        aa_first = aaMap[lb[1]]
                                        ss_second = lb[2]
                                        aa_second = aaMap[lb[3]]
                                        if ss_first == (mem_place+5+1+place) and ss_second == (mem_place + 5 + 1 + place + i):
                                            has_local_bridge += 1
                                            local_bridge_list.add(key)
                                            local_bridges_text = aa_first + str(ss_first)+"-"+aa_second + str(ss_second)
                                        
                                # print(key, "OPP", midMem, mem_place + 5, i, has_mem_bridge, has_local_bridge, local_bridge)
                                # if has_mem_bridge > 1:
                                #     print(key, mem_bridge)
                                mem_bridge_count += has_mem_bridge
                                local_bridge_count += has_local_bridge
                            # charge_csv_line = [key, str(aa_place), aa, str(mem_place+1+5),str(mem_place+len(fullMem)-5)]
                            csv_text_pairs.append(','.join([key, pair_type, str(mem_place +1 + 5), str(mem_place+len(fullMem)-5),
                                                      str(mem_place + 1+ 5 + place), str(mem_place + 1 + 5 + place + i),
                                                      str(i), aa, second_aa, mem_bridges_text, local_bridges_text]))


with open(charges_file, 'wb') as dataPickle:
    pickle.dump([aaCount, aaPairs, aaHits], dataPickle)
with open(csv_file_pairs, 'w') as csv_handle:
    csv_handle.write('\n'.join(csv_text_pairs))
with open(csv_file_aas, 'w') as csv_handle:
    csv_handle.write('\n'.join(csv_text_aas))
print("Proteins with mems: {}".format(len(proteins_with_mems)))
print("Proteins with correct mems: {}".format(len(proteins_with_correct_mems)))
print("Proteins with any salt bridge: {}".format(len(any_saltbridge_prot)))
print("Proteins with local salt bridge: {}".format(len(true_local_saltbridge_prot)))
print("Number of correct mems: {}".format(correct_mems))
print("Mem bridges: {}".format(mem_bridge_count))
print("Local bridges: {}".format(local_bridge_count))
# #         if mem in TMdata[key][0]:
#         for place, aa in enumerate(midMem):
#             stats['aas'] += 1
#             if aa in polaraas:
#                 stats['polar'] += 1
#             if aa in charged:
#                 stats['chr'] += 1
#             if aa in positive:
#                 stats['pos'] += 1
#             if aa in negative:
#                 stats['neg'] += 1
#             for i in range(1, 8):
#                 if memLen <= place + i:
#                     break
#                 # aaPairs[i - 1] += 1
#                 # first = aas.index(aa)
#                 # print(midMem[place+i])
#                 # second = aas.index(midMem[place + i])
#                 secAA = midMem[place + i]
#                 if aa in charged and secAA in charged:
#                     deltaG = dgCalc.calc_segment_DG(fullMem)
#                     pair_found = True
#                     # Add five is midmem
#                     globalPlace = mem_place + edge
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
#                     proteinSet.add(key)
#                     if aa in positive and secAA in positive:
#                         stats['pospair'] += 1
#                     if aa in negative and secAA in negative:
#                         stats['negpair'] += 1
#                     if (aa in positive and secAA in negative) or (aa in negative and secAA in positive):
#                         stats['opppair'] += 1
#                     if aa in charged and secAA in charged:
#                         stats['chargedpair'] += 1
#                     charged_pair_res.add(key + '_' + str(globalPlace))
#                     charged_pair_res.add(key + '_' + str(globalPlace+i))
#                 # aaHits[i - 1][first][second] += 1
#                 # # Count both first and second AA for each gap sized i
#                 # aaCount[0][i - 1][first] += 1
#                 # aaCount[1][i - 1][second] += 1
#         if pair_found == True:
#             mems_with_charge += 1
# print("Number of proteins: ", len(helicies.keys()))
# print("Total membranes: ", totalMems)
# print("Long membranes: ", longMems)
# print("Correct membranes: ", correctMems)
# print("Proteins that contain charged pairs: ", len(proteinSet))
# print("Membrane regions with charged pairs: ", mems_with_charge)
# print("Number of AAs in core region: {}".format(stats['aas']))
# print("Number of polar AAs in core region: {} ({:.1%})".format(stats['polar'],stats['polar']/stats['aas']))
# print("Number of charged AAs in core region: {} ({:.1%})".format(stats['chr'],stats['chr']/stats['aas']))
# print("Number of charged pair AAs: {} ({:.1%}) ({:.1%})".format(len(charged_pair_res), len(charged_pair_res)/stats['aas'], len(charged_pair_res)/stats['chr']))
# print("Number of charged pairs: {}".format(stats["chargedpair"]))
# print("Number of opp pairs: {} ({:.1%})".format(stats["opppair"], stats["opppair"]/stats['chargedpair']))
# # print(charged_pair_res)
# 
# sortedProt = sorted(proteinsExamples.keys())
# debug = False
# has_bridge = set()
# has_bridge_within_7 = set()
# lacks_bridge = set()
# for full_pdb_id in TMdata.keys():
#     pdb_id = full_pdb_id[:4]
#     chain = full_pdb_id[-1]
#     filename = pdb_id + '.pdb'
# 
#     filepath = 'data/pdbFiles/' + pdb_id + '.pdb'
#     saltpath = 'data/bridgeFiles/' + pdb_id + chain + str(con_req) + 'Bridges.pickle'
#     if not os.path.exists(filepath):
#         try:
#             urllib.request.urlretrieve(pdbURL + filename, filepath)
#         except HTTPError as err:
#             print("{} not exists, skipping...".format(filename), file=sys.stderr)
#             continue
#     if not os.path.exists(saltpath):
#         bridges = saltBridges.calcSaltBridges(filepath, chain, 4, con_req)
#         pickle.dump(bridges, open(saltpath, 'wb'))
#     else:
#         bridges = pickle.load(open(saltpath, 'rb'))
#     if len(bridges)>0:
#         has_bridge.add(full_pdb_id)
#     else:
#         lacks_bridge.add(full_pdb_id)
#     for bridge in bridges:
#         if bridge[1][3] - bridge[0][3] < 8:
#             has_bridge_within_7.add(full_pdb_id)
# print("Has bridges {}".format(len(has_bridge)))
# print("Has bridges within 7: {}".format(len(has_bridge_within_7)))
# print("Lacks bridges {}".format(len(lacks_bridge)))
# print("Has mems and bridges {}".format(len(has_bridge & fullproteinSet)))
# print("Has long mems and bridges {}".format(len(has_bridge & longproteinSet)))
# print("Has long mems,charge and bridges {}".format(len(has_bridge & proteinSet)))
# print("Has long mems,charge and bridges within 7: {}".format(len(has_bridge_within_7 & proteinSet)))
# true_saltbridge_prot = set()
# true_saltbridge = 1
# for key in sortedProt:
#     midMem, place, i, globalPlace, aa, secAA, dG = proteinsExamples[key]
#     pureKey = key[:5]
# 
#     # if pureKey[:4] not in topcons_added_keys:
#     #     continue
#     chain = key[4]
#     filename = pureKey[:4] + '.pdb'
#     filepath = 'data/pdbFiles/' + pureKey[:4] + '.pdb'
#     saltpath = 'data/bridgeFiles/' + pureKey[:4] + chain + str(con_req) + 'Bridges.pickle'
# 
#     if not os.path.exists(filepath):
#         try:
#             urllib.request.urlretrieve(pdbURL + filename, filepath)
#         except HTTPError as err:
#             print("{} not exists, skipping...".format(filename), file=sys.stderr)
#             continue
#     if not os.path.exists(saltpath):
#         bridges = saltBridges.calcSaltBridges(filepath, pureKey[4], 4, con_req)
#         pickle.dump(bridges, open(saltpath, 'wb'))
#     else:
#         bridges = pickle.load(open(saltpath, 'rb'))
#     # if len(bridges)>0:
#     #     has_bridge += 1
#     # else:
#     #     lacks_bridge += 1
#     does_have_saltbridge = False
# # print(has_bridge)
# # print(lacks_bridge)
#     if (i < 8) and ((aa in positiveplus and secAA in negative) or (aa in negative and secAA in positiveplus)):
#         if not os.path.exists(filepath):
#             try:
#                 urllib.request.urlretrieve(pdbURL + filename, filepath)
#             except HTTPError as err:
#                 print("{} not exists, skipping...".format(filename), file=sys.stderr)
#                 continue
#         if not os.path.exists(saltpath):
#             bridges = saltBridges.calcSaltBridges(filepath, pureKey[4], 4, con_req)
#             pickle.dump(bridges, open(saltpath, 'wb'))
#         else:
#             bridges = pickle.load(open(saltpath, 'rb'))
#         
#         numMem = countMems(TMdata[pureKey][1])
#         memNumber = whatMem(TMdata[pureKey][1], globalPlace, args.tolerant)
#         span = 'multi span' if numMem > 1 else 'single span'
#         for bridge in bridges:
#             # print("In bridges")
#             # print(bridge)
#             if (bridge[0][3] == globalPlace
#                and bridge[0][2] == chain)\
#                and (bridge[1][3] == (globalPlace + i)
#                and bridge[1][2] == chain):
#                 # and aaMap[bridge[0][1]] == aa:
#                 # print(bridge)
#                 # print(str(i))
#                 # if (bridge[1][3] - bridge[0][3]) == i:
#                 true_saltbridge += 1
#                 true_saltbridge_prot.add(pureKey[:5])
#                 does_have_saltbridge = True
#                 # print(pureKey[:4] + "(" + pureKey[4] + ")",
#                 #       span,
#                 #       aa + secAA + '-pair with gap ' + str(i),
#                 #       'startindex ' + str(globalPlace) + '(' +
#                 #       str(memNumber) + ' membrane region)',
#                 #       'dG ' + '{0:.2f}'.format(dG),
#                 #       sep=', ')
#                 # if len(bridge[0][4]) > 0:
#                 #     first_alt = "alt. conf {}".format(bridge[0][4])
#                 # else:
#                 #     first_alt = ''
#                 # if len(bridge[1][4]) > 0:
#                 #     second_alt = "alt. conf {}".format(bridge[1][4])
#                 # else:
#                 #     second_alt = ''
# 
#                 # print(bridge[0][0],
#                 #       first_alt,
#                 #       bridge[0][1],
#                 #       bridge[0][2],
#                 #       bridge[0][3],
#                 #       bridge[1][0],
#                 #       second_alt,
#                 #       bridge[1][1],
#                 #       bridge[1][2],
#                 #       bridge[1][3],
#                 #       "{0:.2f}".format(bridge[2]) + "Å")
#     # print(bridges)
#         # if does_have_saltbridge:
#         #     print(TMdata[pureKey][0])
#         #     print(TMdata[pureKey][1])
#       # if p[9] > 5:
#           # 4c9g(A), multi span (6 membranes), DK-pair with gap 1,
#           # startindex 108(second membrane region), dG 7.86
#           # print(p)
# print("Proteins with true salt bridge: {}".format(len(true_saltbridge_prot)))
# print("Number of true salt bridge: {}".format(true_saltbridge))
# print("Number of true salt bridge per charged pair: {:.1%}".format(true_saltbridge/stats["chargedpair"]))
# print("Number of true salt bridge per opp pair: {:.1%}".format(true_saltbridge/stats["opppair"]))
# # print(true_saltbridge_prot)
