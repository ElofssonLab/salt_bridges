#!/usr/bin/env python3
import sys
import os
# import xml.etree.ElementTree as etree
# import glob
# import pandas as pd
import saltBridges
import re
import pickle
from Bio.PDB import *
from Bio.PDB.Polypeptide import three_to_one
import argparse

# if len(sys.argv) != 4:
#     print("Usage: " + __file__ + " <pdb list> <xml of proteins> <out file>")
#     sys.exit()
# 
# # input_folder = sys.argv[1]
# list_file = sys.argv[1]
# xml_file = sys.argv[2]
# out_file = sys.argv[3]
# bridge_file = sys.argv[3][:-6]+ "_bridges.pickle"
# 
parser = argparse.ArgumentParser()

parser.add_argument("mem_pickle", type=str, help="Membrane pickle")
# parser.add_argument("pdb_dir", type=str, help="PDB file directory")

parser.add_argument("-b", "--bridges", type=int, default=1, help="Required connections per bridge")
parser.add_argument("-t", "--tolerant", type=bool, default=False, help="Use tolerant membranes (also include m, not just M)")

args = parser.parse_args()

bridge_file = args.mem_pickle[:-7]+ "_bridges.pickle"
helices = pickle.load(open(args.mem_pickle, 'rb'))
min_mem_len = 17
gap_offset = 5
# debug = False
# struct_dict = {}
all_bridge_dict = {}
mem_bridge_dict = {}
local_bridge_dict = {}

# TMdata = {}
# with open(args.threeline, 'r') as TMHandle:
#     pdb_id = ''
#     rowNum = 0
#     fa = ''
#     for line in TMHandle:
#         if line[0] == '>':
#             # pid = line[1:].strip().split('|')[1]
#             if '|' in line:
#                 pdb_id = line[1:].strip().split('|')[1]
#             else:
#                 pdb_id = line[1:].strip()
#             # pid = line[1:].strip()# .split('|')[2]
#         elif rowNum % 3 == 1:
#             fa = line.strip()
#         elif rowNum % 3 == 2:
#             if not "MISSING" in pdb_id:
#                 TMdata[pdb_id] = [fa, line.strip()]
#         else:
#             print("You should not be here...")
#         rowNum += 1

# print(TMdata)
# with open(list_file) as list_handle:
#     list_handle.readline()  # Header line
#     for line in list_handle:
#         pdb_id = line.strip().split()[0]
#         pdb = pdb_id[:-1]
#         chain = pdb_id[-1]
#         if pdb in struct_dict:
#             struct_dict[pdb].append(chain)
#         else:
#             struct_dict[pdb] = [chain]
# # struct_list = [p.lower() for p in raw_list]
# # print(len(chain_list))
# # print(len(pdb_list))
# # print(chain_list[:5])
# # print(pdb_list[:5])
# # sys.exit()
# # print(pdb_list[:10])
# namespace = "{http://pdbtm.enzim.hu}"
# raw_data = []
# pdb_dir = 'data/pdbFiles'
# # for sub_folder in glob.glob(input_folder + "/*"):
# #     for input_file in glob.glob(sub_folder + "/*.xml"):
# # print(xml_file)
# tree = etree.parse(xml_file)
# root = tree.getroot()
# # print(root)
# # sys.exit()
# out_text = ""
pdb_parser = PDBParser(QUIET=True)
pdblist = PDBList(pdb="data/pdbFiles", verbose=False)
ppBuilder = Polypeptide.PPBuilder(4)
# count = 0
tot = len(helices)
i = 0
for key, membranes in helices.items():
    i += 1
    keep_running = True
    pdb_seq_offset = {}
    # Total number of membrane segments
    # print(key, membranes)
    pdb_id = key[:4].upper()
    chain_id = key[4]

    pdb_file = pdblist.retrieve_pdb_file(pdb_id, pdir="data/pdbFiles", file_format="pdb")
    if not os.path.isfile(pdb_file):
        # print("No file {}".format(pdb_file))
        continue
    pdb_struct = pdb_parser.get_structure(pdb_id, pdb_file)
    struct_type = pdb_struct.header["structure_method"]
    full_chain_id = pdb_id + chain_id
    pdb_seq_offset[full_chain_id] = 0
    with open(pdb_file) as pdb_handle:
        for line in pdb_handle:
            if line.startswith("DBREF"):
                parts = line.split()[:4]
                if parts[0] in ["DBREF", "DBREF1"]:
                    pdb_seq_offset[parts[1] + parts[2]] = int(re.sub('[^0-9]','', parts[3])) - 1  # The pdb seq starts numbering @
    # print(pdb_seq_offset)
    # sys.exit()
    # Only use well established methods
    if not struct_type in ["x-ray diffraction", "electron microscopy","solution nmr"]:
        continue
    pp_list = ppBuilder.build_peptides(pdb_struct[0][chain_id])

    bridges = saltBridges.calcSaltBridges(pdb_file, chain_id, 4)
    all_bridge_dict[str(pdb_id) + chain_id] = bridges

    for mem_num, mem_data in enumerate(membranes):
        if not keep_running:
            break
        # print(key, pdb_id, chain_id, mem_num, mem_data)
        # sys.exit()
        mem = mem_data[1] 
        mem_len = len(mem)
        if mem_len < min_mem_len:  # Only use membranes longer than 17
            continue    
        global_place = mem_data[0] + pdb_seq_offset[full_chain_id]
        mem_start = global_place + 1 
        mem_end = mem_start + mem_len
        ### bridge = [resi1, res1, resi2, res2, chain, dist]
        # print(mem)
        # print(mem_start, mem_end)
        for bridge in bridges:
            save_bridge = [bridge[0], bridge[1], bridge[2], bridge[3], bridge[4], bridge[5]]
            # print(mem[bridge[0]-mem_start], mem[bridge[2]-mem_start])
            # print(save_bridge)
            if (bridge[0] >= mem_start and bridge[0] <= mem_end) or (bridge[2] >= mem_start and bridge[2] <= mem_end):
                # print("Mem bridge")
                s = save_bridge[0]
                first_aa = three_to_one(save_bridge[1])
                e = save_bridge[2]
                second_aa = three_to_one(save_bridge[3])
                seq_first = int(save_bridge[0])-mem_start
                seq_second = int(save_bridge[2])-mem_start
                # print("**********************")
                # print(bridge)
                # print(mem_start, mem_end)
                if seq_first > -1 and seq_first < mem_len:
                    if mem[seq_first] != first_aa:
                        print("Membrane bridge out of sync {}".format(full_chain_id))
                        keep_running = False
                        break
                        # print("ERROR! Membrane first")
                        # print(mem[seq_first], first_aa)
                        # print(full_chain_id, bridge, mem, mem_start)
                        # sys.exit()
                if seq_second > -1 and seq_second < mem_len:
                    if mem[seq_second] != second_aa:
                        print("Membrane bridge out of sync {}".format(full_chain_id))
                        keep_running = False
                        break
                        # print("ERROR! Membrane second")
                        # print(mem[seq_second], second_aa)
                        # print(full_chain_id, bridge, mem, mem_start)
                        # sys.exit()
                # print(s,e, first_aa, second_aa, seq_first, seq_second)
                if full_chain_id in mem_bridge_dict:
                    if save_bridge not in mem_bridge_dict[str(pdb_id) + chain_id]:
                        mem_bridge_dict[str(pdb_id) + chain_id].append(save_bridge)
                else:
                    mem_bridge_dict[str(pdb_id) + chain_id] = [save_bridge]


            if mem_end - mem_start >= 17:  # Only use core parts if saving local bridges
                if (bridge[0] >= mem_start+gap_offset and bridge[0] <= mem_end-gap_offset) and (bridge[2] >= mem_start+gap_offset and bridge[2] <= mem_end-gap_offset):
                    # print("Local bridge")
                    # print("**********************")
                    # print(bridge)
                    # print(mem_start, mem_end)
                    s = save_bridge[0]
                    first_aa = three_to_one(save_bridge[1])
                    e = save_bridge[2]
                    second_aa = three_to_one(save_bridge[3])
                    seq_first = int(save_bridge[0])-mem_start
                    seq_second = int(save_bridge[2])-mem_start
                    # print(mem_data)
                    # print(full_chain_id, seq_first, seq_second)
                    if (mem[seq_first] != first_aa) or (mem[seq_second] != second_aa):
                        print("Local bridge out of sync {}".format(full_chain_id))
                        keep_running = False
                        break
                        # print("ERROR! Local bridge")
                        # print(full_chain_id, bridge, mem, mem_start)
                        # sys.exit()
                    if full_chain_id in local_bridge_dict:
                        local_bridge_dict[str(pdb_id) + chain_id].append(save_bridge)
                    else:
                        local_bridge_dict[str(pdb_id) + chain_id] = [save_bridge]

    print("Done {}/{}...".format(i, tot))
        # sys.exit()
               # stop = True

#                         T = 'M'
#                         m_len = pdb_end - pdb_start + 1
#                         if m_len > 10:
#                             # ss_string = ''
#                             for ss in range(pdb_start + 5, pdb_end + 1 - 5):
#                                 try:
#                                     # print(ss, dssp[chain_id, (' ', ss, ' ')][2])
#                                     # ss_string += dssp[chain_id, (' ', ss, ' ')][2]
#                                     # if not dssp[chain_id, (' ', ss, ' ')][2] in 'H -':
#                                     # If not the core of the longer membranes are ONLY helix, skip them but retain the protein and other membranes
#                                     if not dssp[chain_id, (' ', ss, ' ')][2] in 'H':
#                                         # print("DSSP not correct", chain_id, ss)
#                                         T = 'm'
#                                         # save_protein = False 
#                                         # break
#                                 except:
#                                     print("dssp fail", pdb_id, chain_id, ss)
#                                     save_protein = False
#                                     break
#                             # print(pdb_id, ss_string)
#                         seq += raw_seq[reg_start-1:reg_end]
#                     else:
#                         T = 'i'
#                         seq += raw_seq[reg_start-1:reg_end]
#                     topo += (reg_end - reg_start + 1)*T
#                 if not save_protein:
#                     print("Not saving protein {}".format(full_chain_id))
#                     continue
#                 if len(str(seq)) != len(topo):
#                     print("Seq and topo not matching, {}".format(pdb_id))
#                     print(str(seq))
#                     print(topo)
#                     # sys.exit()
#                     continue
# 
#                 out_text += ">" + pdb_id + chain_id + '\n'
#                 out_text += str(seq) + '\n'  # [first_seq-1:end_seq] + '\n'
#                 out_text += topo + '\n'

# with open(out_file, 'w') as out_handle:
#     out_handle.write(out_text)
pickle.dump({'all':all_bridge_dict, 'mems':mem_bridge_dict, 'local':local_bridge_dict},open(bridge_file,'wb'))
