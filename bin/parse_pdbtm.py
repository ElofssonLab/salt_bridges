#!/usr/bin/env python3
import sys
import os
import xml.etree.ElementTree as etree
import glob
import pandas as pd
from Bio.PDB import *

if len(sys.argv) != 4:
    print("Usage: " + __file__ + " <pdb list> <xml of proteins> <out 3lin>")
    sys.exit()

# input_folder = sys.argv[1]
list_file = sys.argv[1]
xml_file = sys.argv[2]
out_file = sys.argv[3]

struct_dict = {}
with open(list_file) as list_handle:
    list_handle.readline()  # Header line
    for line in list_handle:
        pdb_id = line.strip().split()[0]
        pdb = pdb_id[:-1]
        chain = pdb_id[-1]
        if pdb in struct_dict:
            struct_dict[pdb].append(chain)
        else:
            struct_dict[pdb] = [chain]
# struct_list = [p.lower() for p in raw_list]
# print(len(chain_list))
# print(len(pdb_list))
# print(chain_list[:5])
# print(pdb_list[:5])
# sys.exit()
# print(pdb_list[:10])
namespace = "{http://pdbtm.enzim.hu}"
raw_data = []
pdb_dir = 'data/pdbFiles'
# for sub_folder in glob.glob(input_folder + "/*"):
#     for input_file in glob.glob(sub_folder + "/*.xml"):
# print(xml_file)
tree = etree.parse(xml_file)
root = tree.getroot()
# print(root)
# sys.exit()
out_text = ""
parser = PDBParser(QUIET=True)
pdblist = PDBList(pdb="data/pdbFiles", verbose=False)
ppBuilder = Polypeptide.PPBuilder(4)
# count = 0
for prot in root.iter(namespace + "pdbtm"):
# for prot in root.findall(namespace + "pdbtm"):
    pdb_id = prot.attrib["ID"].upper()
    if pdb_id in struct_dict:
    # if pdb_id == "6A93":
        if prot.attrib["TMP"] != "yes":
            print("Not found TMP")
            print(prot.attrib)
            continue
        pdb_file = pdblist.retrieve_pdb_file(pdb_id, pdir="data/pdbFiles", file_format="pdb")
        # print(pdb_file)
        if not os.path.isfile(pdb_file):
            continue
        pdb_struct = parser.get_structure(pdb_id, pdb_file)
        if not "x-ray diffraction" == pdb_struct.header["structure_method"]:
            continue
        try:
            dssp = DSSP(pdb_struct[0], pdb_file)
        except:
            continue
        # print(dssp.keys()[0])
        pdb_chains = [c.get_id() for c in pdb_struct[0].get_chains()]

        for chain in prot.findall(namespace + "CHAIN"):
            num_mems = 0
            chain_id = chain.attrib["CHAINID"]
            # print("Chain {}".format(chain_id))
            if not chain_id in pdb_chains:
                continue
            pp_list = ppBuilder.build_peptides(pdb_struct[0][chain_id])
            if not len(pp_list) > 0:
                continue
            pp_iter = iter(pp_list)
            # print(pp)
            # print(pdb_start, pdb_end)
            # print(pp[0][start:end])
            # sys.exit()
            # print(chain.attrib)
            if int(chain.attrib["NUM_TM"]) > 0 and chain.attrib["TYPE"] == "alpha" and chain.attrib["CHAINID"] in struct_dict[pdb_id]:
                # count += 1
                # print(prot.attrib["ID"])
                # print(chain.attrib)
                # seq = chain.find(namespace + "SEQ").text.replace(' ', '').replace('\n', '').strip()
                topo = ''
                first_seq = None
                pp = next(pp_iter)
                # print(pp)
                pdb_start = pp[0].get_id()[1]
                pdb_end = pp[-1].get_id()[1]
                seq = pp.get_sequence()
                # print(seq)
                # Set prev_end low to mitigate when the start sequence is negative in the pdb/region
                prev_end = -999
                for region in chain.findall(namespace + "REGION"):
                    # if not first_seq:
                    #     first_seq = int(region.attrib["pdb_beg"])
                    # end_seq = int(region.attrib["pdb_end"])
                    # print(region.attrib)
                    save_region = True
                    reg_type = region.attrib["type"]
                    if reg_type == 'H':
                        T = 'M'
                    else:
                        T = 'i'
                    start = int(region.attrib["pdb_beg"])
                    end = int(region.attrib["pdb_end"])
                    # print(pdb_start, pdb_end)
                    # print(start, end )
                    # If the current peptide ends before the current region starts, check for next peptide
                    if pdb_end < start:
                        try:
                            pp = next(pp_iter)
                            # print(pp)   
                        except StopIteration:
                            break
                        pdb_start = pp[0].get_id()[1]
                        pdb_end = pp[-1].get_id()[1]
                        seq += pp.get_sequence()
                    # if region starts before the peptide starts, start from the peptide
                    if start < pdb_start:
                        start = pdb_start
                    # if the current peptide start before the previous region ended, fill up the last part of that region
                    if pdb_start < prev_end:
                        # print(pdb_id)
                        add_topo = prev_T*(prev_end - pdb_start+1)
                        topo += add_topo
                        # print(add_topo, len(add_topo))
                    # If the region ends later than the peptide, stop at the end of the peptide
                    # but save the region type and end if the next peptide cover this part
                    if end > pdb_end:
                        prev_end = end
                        prev_T = T
                        end = pdb_end
                    # Length of the current segment
                    part_num = end - start + 1
                    if T == 'M':
                        for ss in range(start, end + 1):
                            try:
                                if not 'H' == dssp[chain_id, (' ', ss, ' ')][2]:
                                    save_region = False 
                            except:
                                print("dssp fail", pdb_id, chain_id, ss)
                                save_region = False
                                break
                    topo += part_num*T
                    # print(part_num*T, len(part_num*T))
                if not save_region:
                    continue
                # print(seq, pdb_start, pdb_end)
                # count += 1
                if len(str(seq)) != len(topo):
                    print("Seq and topo no matching, {}".format(pdb_id))
                    print(str(seq))
                    print(topo)
                    continue
                out_text += ">" + pdb_id + chain_id + '\n'
                out_text += str(seq) + '\n'  # [first_seq-1:end_seq] + '\n'
                out_text += topo + '\n'
        # break
        # print(out_text)
        # sys.exit()
    # if count > 10:
    #     break
    #     print(out_text)
    #     sys.exit()

with open(out_file, 'w') as out_handle:
    out_handle.write(out_text)
                # print(">" + prot.attrib["ID"].upper() + chain.attrib["CHAINID"])
                # print(seq[first_seq-1:end_seq])
                # print(topo)
                # print(out_text)
                # sys.exit()
                    # if region.attrib["type"] == 'H':
                        # row = [input_file.split('/')[-1].split('.')[0] +
                        #             chain.attrib["CHAINID"],
                        #             seq.text.replace('\n', '').replace(' ', ''),
                        #             len(seq.text.replace('\n', '').replace(' ', '')),
                        #             int(child.attrib["seq_beg"]),
                        #             int(child.attrib["pdb_end"])]
                        # raw_data.append(row)
                        # print(raw_data)
                    # print(child.attrib["type"])

        # sys.exit()

#     if int(chain.attrib["NUM_TM"]) == 1:
#         if chain.attrib["TYPE"] != "alpha":
#             print("WRONG! " + input_file)
#             sys.exit()
#         seq = chain.find(namespace + "SEQ")
#         # print(input_file.split('/')[-1])
#         # print(chain.attrib["CHAINID"])
#         # print(seq.text.replace('\n', '').replace(' ', ''))
#         for child in chain.findall(namespace + "REGION"):
#             if child.attrib["type"] == 'H':
#                 # print(child.attrib["seq_beg"], child.attrib["pdb_end"])
#                 row = [input_file.split('/')[-1].split('.')[0] +
#                             chain.attrib["CHAINID"],
#                             seq.text.replace('\n', '').replace(' ', ''),
#                             len(seq.text.replace('\n', '').replace(' ', '')),
#                             int(child.attrib["seq_beg"]),
#                             int(child.attrib["pdb_end"])]
#                 raw_data.append(row)
#                         # print(raw_data)
#                     # print(child.attrib["type"])
#         # sys.exit()
# 
# df = pd.DataFrame(raw_data, columns=["PID",
#                                      "Sequence",
#                                      "Seq_len",
#                                      "TM_start",
#                                      "TM_end"])
# # print(df.head())
# df.to_pickle("data/pdbtm.pkl")
