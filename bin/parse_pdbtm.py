#!/usr/bin/env python3
import sys
import os
import xml.etree.ElementTree as etree
import glob
import pandas as pd
import saltBridges
from Bio.PDB import *
from Bio.PDB.Polypeptide import three_to_one

if len(sys.argv) != 4:
    print("Usage: " + __file__ + " <pdb list> <xml of proteins> <out 3lin> <salt bridges>")
    sys.exit()

# input_folder = sys.argv[1]
list_file = sys.argv[1]
xml_file = sys.argv[2]
out_file = sys.argv[3]

debug = True
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
count = 0
for prot in root.iter(namespace + "pdbtm"):
# for prot in root.findall(namespace + "pdbtm"):
    pdb_id = prot.attrib["ID"].upper()
    if pdb_id in struct_dict:
    # if pdb_id == "5L25":
        if debug:
            print("Running {}...".format(pdb_id))
        if prot.attrib["TMP"] != "yes":
            print("Not found TMP")
            print(prot.attrib)
            continue
        pdb_file = pdblist.retrieve_pdb_file(pdb_id, pdir="data/pdbFiles", file_format="pdb")
        if not os.path.isfile(pdb_file):
            print("No file {}".format(pdb_file))
            continue
        pdb_struct = parser.get_structure(pdb_id, pdb_file)
        # if not "x-ray diffraction" == pdb_struct.header["structure_method"]:
        #     continue
        try:
            dssp = DSSP(pdb_struct[0], pdb_file)
        except:
            print("Failed dssp {}".format(pdb_id))
            continue
        # print(dssp.keys()[0])
        pdb_chains = [c.get_id() for c in pdb_struct[0].get_chains()]

        for chain in prot.findall(namespace + "CHAIN"):
            # num_mems = 0
            chain_id = chain.attrib["CHAINID"]
            if not chain_id in pdb_chains:
                continue
            if debug:
                print(chain_id)
            # pp_list = ppBuilder.build_peptides(pdb_struct[0][chain_id])
            # if not len(pp_list) > 0:
            #     continue
            # struct_end_id = pp_list[-1][-1].get_id()[1]
            # pp_iter = iter(pp_list)
            # if debug:
            #     print(pp_list)
            if int(chain.attrib["NUM_TM"]) > 0 and chain.attrib["TYPE"] == "alpha" and chain.attrib["CHAINID"] in struct_dict[pdb_id]:
                seq = ""
                topo = ""
                # pp = next(pp_iter)
                j = 0
                # while True:
                #     if pp[j].get_id()[1] < struct_end_id:
                #         struct_pdb_start = pp[j].get_id()[1]
                #         break
                #     j += 1
                # struct_pdb_end = pp[-1].get_id()[1]
                bridges = saltBridges.calcSaltBridges(pdb_file, chain, 4)
                print(bridges)
                sys.exit()
                raw_seq = chain.find(namespace + "SEQ").text.replace(' ', '').replace('\n', '')
                # res_ids = set()
                # num_ids = set(range(struct_pdb_start, struct_pdb_end + 1))
                # for res in pp:
                #     r_id = res.get_id()[1]
                #     res_ids.add(r_id)
                # peptide_missing_res = num_ids - res_ids

                save_protein = True
                for region in chain.findall(namespace + "REGION"):
                    reg_type = region.attrib["type"]
                    reg_start = int(region.attrib["seq_beg"])
                    reg_end = int(region.attrib["seq_end"])
                    if reg_type == 'H':
                        pdb_start = int(region.attrib["pdb_beg"])
                        pdb_end = int(region.attrib["pdb_end"])
                        T = 'M'
                        m_len = pdb_end - pdb_start + 1
                        if m_len > 10:
                            # ss_string = ''
                            for ss in range(pdb_start + 5, pdb_end + 1 - 5):
                                try:
                                    # print(ss, dssp[chain_id, (' ', ss, ' ')][2])
                                    # ss_string += dssp[chain_id, (' ', ss, ' ')][2]
                                    # if not dssp[chain_id, (' ', ss, ' ')][2] in 'H -':
                                    # If not the core of the longer membranes are ONLY helix, skip them but retain the protein and other membranes
                                    if not dssp[chain_id, (' ', ss, ' ')][2] in 'H':
                                        # print("DSSP not correct", chain_id, ss)
                                        T = 'm'
                                        # save_protein = False 
                                        # break
                                except:
                                    print("dssp fail", pdb_id, chain_id, ss)
                                    save_protein = False
                                    break
                            # print(pdb_id, ss_string)
                        seq += raw_seq[reg_start-1:reg_end]
                    else:
                        T = 'i'
                        seq += raw_seq[reg_start-1:reg_end]
                    topo += (reg_end - reg_start + 1)*T
                if not save_protein:
                    continue
                if len(str(seq)) != len(topo):
                    print("Seq and topo not matching, {}".format(pdb_id))
                    print(str(seq))
                    print(topo)
                    # sys.exit()
                    continue

                out_text += ">" + pdb_id + chain_id + '\n'
                out_text += str(seq) + '\n'  # [first_seq-1:end_seq] + '\n'
                out_text += topo + '\n'
                if debug:
                    print(out_text)
                    sys.exit()
                    #########################################################
                    #########################################################
                    #########################################################
                    # reg_start = int(region.attrib["pdb_beg"])
                    # reg_end = int(region.attrib["pdb_end"])
                    # # print("Region {} {}".format(reg_start, reg_end))
                    # # print(reg_start, reg_end)
                    # while True:
                    #     # If the region is before the start of the peptide, go for the next region
                    #     if reg_end < pdb_start:
                    #         break
                    #     # Normal behavior
                    #     ########################
                    #     if reg_start < pdb_start:
                    #         start = pdb_start
                    #     elif reg_start > struct_end_id and reg_end > struct_end_id:
                    #         # All done, no more pdb
                    #         break
                    #     elif reg_start > struct_end_id and reg_end:
                    #         # Special if first seq is dbxref, ie 1001 or the like
                    #         start = pdb_start
                    #     else:
                    #         start = reg_start

                    #     if reg_end > pdb_end:
                    #         end = pdb_end
                    #     else:
                    #         end = reg_end
                    #     ########################
                    #     # if debug: 
                    #     #     print(start, end)

                    #     curr_missing = {i for i in peptide_missing_res if i >= start and i <= end}
                    #     if debug:
                    #         if len(curr_missing) > 0:
                    #             print("Missing ", curr_missing)
                    #     part_num = end - start + 1 - len(curr_missing)
                    #     # if T == 'M':
                    #     #     reg_dssp = ''
                    #     #     # print("Starting dssp with ",start, end)
                    #     #     # for t in range(start, end+1):
                    #     #     #     print("DSSP for {}".format(t))
                    #     #     #     print(dssp[chain_id, (' ', t, ' ')][2])
                    #     #     for ss in set(range(start, end + 1)) - curr_missing:
                    #     #         try:
                    #     #             # print(dssp[chain_id, (' ', ss, ' ')][2])
                    #     #             # print('M', start, end+1)
                    #     #             # reg_dssp += dssp[chain_id, (' ', ss, ' ')][2]
                    #     #             if not 'H' == dssp[chain_id, (' ', ss, ' ')][2]:
                    #     #                 print("DSSP not correct")
                    #     #                 save_protein = False 
                    #     #                 # break
                    #     #         except:
                    #     #             print("dssp fail", pdb_id, chain_id, ss)
                    #     #             save_protein = False
                    #     #             break
                    #     #     # if (len(reg_dssp) - reg_dssp.count('H') - reg_dssp.count('G') - reg_dssp.count('I')) > 1 and len(reg_dssp) > 10:
                    #     #     #     if debug:
                    #     #     #         print("DSSP not correct {}".format(pdb_id))
                    #     #     #         print(reg_dssp)
                    #     #     #     # Just skip this membrane?
                    #     #     #     # break

                    #     #     # #     # sys.exit()
                    #     #     #     save_protein = False
                    #     #     #     break
                    #     if debug:
                    #         print(part_num*T, start, end)
                    #     topo += part_num*T
                    #     # Region end is before the peptide end, break and grab next region
                    #     if reg_end < pdb_end:
                    #         break
                    #     # Otherwise grab the next peptide
                    #     else:
                    #         try:
                    #             pp = next(pp_iter)
                    #         except StopIteration:
                    #             break
                    #         pdb_start = pp[0].get_id()[1]
                    #         pdb_end = pp[-1].get_id()[1]
                    #         for res in pp:
                    #             if not res.get_id()[1] > struct_end_id:
                    #                 seq += three_to_one(res.get_resname())
                    #         res_ids = set()
                    #         num_ids = set(range(pdb_start, pdb_end + 1))
                    #         for res in pp:
                    #             r_id = res.get_id()[1]
                    #             res_ids.add(r_id)
                    #         peptide_missing_res = num_ids - res_ids
                    #         
                    # # if not save_protein:
                    # #     break

                    #########################################################
                    #########################################################
                    #########################################################

                   # print(part_num*T, len(part_num*T))

                #     ######################################3
                #         
                #     # if region starts before the peptide starts, start from the peptide
                #     if start < pdb_start:
                #         start = pdb_start
                #     # if the current peptide start before the previous region ended, fill up the last part of that region
                #     if pdb_start < prev_end:
                #         # print(pdb_id)
                #         if prev_T == 'M':
                #             for ss in range(pdb_start, prev_end + 1):
                #                 try:
                #                     if not 'H' == dssp[chain_id, (' ', ss, ' ')][2]:
                #                         save_region = False 
                #                         break
                #                 except:
                #                     print("dssp fail", pdb_id, chain_id, ss)
                #                     save_region = False
                #                     break
                #         curr_missing = [i for i in peptide_missing_res if i >= pdb_start and i <= prev_end]
                #         add_topo = prev_T*(prev_end - pdb_start + 1 - len(curr_missing))
                #         if debug:
                #             print(add_topo, "peptide start before", pdb_start, prev_end)
                #         topo += add_topo
                #         # We have already filled out the missing residues, "reset" prev_end
                #         prev_end = -999
                #         # print(add_topo, len(add_topo))
                #     # If the current peptide ends before the current region starts, check for next peptide
                #     while pdb_end < start:
                #         try:
                #             pp = next(pp_iter)
                #             # print(pp)   
                #         except StopIteration:
                #             break
                #         pdb_start = pp[0].get_id()[1]
                #         pdb_end = pp[-1].get_id()[1]
                #         res_ids = set()
                #         num_ids = set(range(pdb_start, pdb_end + 1))
                #         for res in pp:
                #             r_id = res.get_id()[1]
                #             res_ids.add(r_id)
                #         peptide_missing_res = num_ids - res_ids
                #         if debug:
                #             if len(peptide_missing_res) > 0:
                #                 print(peptide_missing_res)
                #         # print("PDB", pdb_start, pdb_end)
                #         # if len(pp.get_sequence()) != (pdb_end - pdb_start + 1):
                #         #     print("Lengths not matchin")
                #         seq += pp.get_sequence()
                #         if debug:
                #             print(pp.get_sequence())
                #     # If the region ends later than the peptide, stop at the end of the peptide
                #     # but save the region type and end if the next peptide cover this part
                #     if end > pdb_end:
                #         prev_end = end
                #         prev_T = T
                #         end = pdb_end
                #     # Length of the current segment
                #     curr_missing = [i for i in peptide_missing_res if i >= start and i <= end]
                #     part_num = end - start + 1 - len(curr_missing)
                #     if T == 'M':
                #         for ss in range(start, end + 1):
                #             try:
                #                 if not 'H' == dssp[chain_id, (' ', ss, ' ')][2]:
                #                     save_region = False 
                #             except:
                #                 print("dssp fail", pdb_id, chain_id, ss)
                #                 save_region = False
                #                 break
                #     if debug:
                #         print(part_num*T, "Standard", start, end)
                #     topo += part_num*T
                #     # print(part_num*T, len(part_num*T))
                # if not save_region:
                #     continue
                # # print(seq, pdb_start, pdb_end)
                # # count += 1
                # if len(str(seq)) != len(topo):
                #     print("Seq and topo no matching, {}".format(pdb_id))
                #     print(str(seq))
                #     print(topo)
                #     continue
                # out_text += ">" + pdb_id + chain_id + '\n'
                # out_text += str(seq) + '\n'  # [first_seq-1:end_seq] + '\n'
                # out_text += topo + '\n'
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
