#!/usr/bin/env python3
import sys
import xml.etree.ElementTree as etree
import glob
import pandas as pd

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
        pdb = pdb_id[:-1].lower()
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
# for sub_folder in glob.glob(input_folder + "/*"):
#     for input_file in glob.glob(sub_folder + "/*.xml"):
# print(xml_file)
tree = etree.parse(xml_file)
root = tree.getroot()
# print(root)
# sys.exit()
out_text = ""
for prot in root.findall(namespace + "pdbtm"):
    # print(prot.attrib["ID"])
    if prot.attrib["ID"] in struct_dict:
        # print(prot.attrib["TMP"])
        if prot.attrib["TMP"] != "yes":
            print("Not found TMP")
            print(prot.attrib)
            sys.exit()
        for chain in prot.findall(namespace + "CHAIN"):
            # print(chain.attrib)
            if int(chain.attrib["NUM_TM"]) > 0 and chain.attrib["TYPE"] == "alpha" and chain.attrib["CHAINID"] in struct_dict[prot.attrib["ID"]]:
                # print(prot.attrib["ID"])
                # print(chain.attrib)
                seq = chain.find(namespace + "SEQ").text.replace(' ', '').replace('\n', '').strip()
                # print(seq)
                topo = ''
                first_seq = None
                for region in chain.findall(namespace + "REGION"):
                    if not first_seq:
                        first_seq = int(region.attrib["seq_beg"])
                    end_seq = int(region.attrib["seq_end"])
                    # print(region.attrib)
                    reg_type = region.attrib["type"]
                    if reg_type == 'H':
                        T = 'M'
                    else:
                        T = 'i'
                    start = int(region.attrib["seq_beg"])
                    stop = int(region.attrib["seq_end"])
                    part_num = stop - start + 1
                    topo += part_num*T
                    # print(reg_type)
                out_text += ">" + prot.attrib["ID"].upper() + chain.attrib["CHAINID"] + '\n'
                out_text += seq[first_seq-1:end_seq] + '\n'
                out_text += topo + '\n'

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
                        #             int(child.attrib["seq_end"])]
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
#                 # print(child.attrib["seq_beg"], child.attrib["seq_end"])
#                 row = [input_file.split('/')[-1].split('.')[0] +
#                             chain.attrib["CHAINID"],
#                             seq.text.replace('\n', '').replace(' ', ''),
#                             len(seq.text.replace('\n', '').replace(' ', '')),
#                             int(child.attrib["seq_beg"]),
#                             int(child.attrib["seq_end"])]
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
