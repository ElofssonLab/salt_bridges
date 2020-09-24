#!/usr/bin/env python3
###################################################################
# This file creates 3line files from OPM proteins.
#
# The resulting 3line has uniprot id and pdb id
###################################################################

import json
import requests
import re
import sys
import os
import argparse
import gzip

page_size = 3000
base_url = "https://lomize-group-opm.herokuapp.com/"
class_url = "classtypes/{!s}/primary_structures?pageSize=" +\
             str(page_size)

parser = argparse.ArgumentParser()
parser.add_argument("input_file", type=argparse.FileType('r'), help="Input json")
parser.add_argument("map_file", type=str, help="pdb_chain_uniprot.tsv.gz")
parser.add_argument("output_file", type=argparse.FileType('w'), help="Output 3line file")
args = parser.parse_args()

uniprot_captured = {}
primary_structure_url = "primary_structures/{}"
membrane_url = "membranes/{!s}"
pdbmap = {}
out_text = ""

with gzip.open(args.map_file, 'rb') as pdbmap_file:
    pdbmap_file.readline()
    pdbmap_file.readline()
    for line in pdbmap_file:
        parts = str(line, 'utf-8').split("\t")
        pdbmap[parts[0].upper() + parts[1]] = [parts[2], int(parts[7])]

# with open("data/opm/primary_structures_{}.json".format(args.prot_class)) as class_handle:
class_list = json.load(args.input_file)

top_list = []
full_entries = {}
tot = len(class_list["objects"])
print("Running 1/{}".format(tot))
for i, entry in enumerate(class_list["objects"]):
    debug = True
    if debug and (i+1)%50 == 0:
        print("Running {}/{}".format(i+1, tot))
    chain_segments = dict()
    struct_id = entry["id"]
    structure = json.loads(
                     requests.get(
                              base_url +
                              primary_structure_url.format(struct_id)).text)
    if structure["pdbid"] == "1su4":
        print("Got 1su4...")
    subunit_segments = int(structure["subunit_segments"])
    if not subunit_segments > 0:
        continue
    pdbid_base = structure["pdbid"].upper()
    # struct_type = structure["resolution"]
    topology_start = "I" if structure["topology_show_in"] is True\
                     else "O"
    mem_name = structure["membrane_name_cache"]

    if mem_name in ["Secreted", "Viral", "Undefined"]:
        continue
    for subunit in structure["subunits"]:
        mem_segments = []
        opm_mem_segments = subunit["segment"]
        sp_stop = 0
        chain = subunit["protein_letter"]
        patt = r"(?P<start>[0-9]{1,4})\s*-\s*(?P<end>[0-9]{1,4})"
        segs = re.findall(patt, subunit["segment"])

        if pdbid_base + chain in pdbmap:
            uniprot = pdbmap[pdbid_base + chain][0]
            if uniprot in uniprot_captured:
                header = '>' + uniprot + '|' + pdbid_base + chain
                out_text += header + '\n'
                out_text += uniprot_captured[uniprot][1] + '\n'
                out_text += uniprot_captured[uniprot][2] + '\n'
                continue
            try:
                raw_fasta = requests.get("https://www.uniprot.org/uniprot/{}.fasta".format(uniprot)).text.strip()
                if not len(raw_fasta) > 0:
                    continue
                fasta = raw_fasta.split("\n")
                info_text = requests.get("https://www.uniprot.org/uniprot/{}.txt".format(uniprot)).text

                for line in info_text.split("\n"):
                    if line.startswith("FT   SIGNAL"):
                        # print(line.strip())
                        if '.' in line[20:].strip():
                            # print(line[20:28])
                            sp_stop = int(line[20:].strip().replace('..', '.').split('.')[-1])
                        else:
                            sp_stop = int(line[20:].strip())
                    elif line.startswith("FT   TRANSMEM"):
                        mem_segments.append([int(num) for num in line[13:].strip().replace('..', '.').split('.')])
                        continue

                header = '>' + uniprot + '|' + pdbid_base + chain
                fa = ''.join(fasta[1:])
                topo = ""
                curr_topo = topology_start
                curr_pos = 0
                for mem in mem_segments:
                    if curr_pos == 0 and sp_stop > 0:
                        topo += 'S' * sp_stop
                        curr_pos = sp_stop
                    start, stop = mem
                    topo += curr_topo * (start - curr_pos - 1)
                    topo += 'M' * (stop - start + 1)
                    curr_pos = stop
                    curr_topo = 'I' if curr_topo == 'O' else 'O'
                topo += curr_topo * (len(fa) - curr_pos)
                if len(fa) == len(topo):
                    out_text += header + '\n'
                    out_text += fa + '\n'
                    out_text += topo + '\n'
                else:
                    print("Unequal fasta and topo:")
                    print(header)
                    print(fa)
                    print(topo)
                uniprot_captured[uniprot] = [header, fa, topo]
            except requests.exceptions.RequestException as e:
                print("Request error for {}".format(uniprot))
                continue
            except requests.exceptions.ConnectionError as e:
                print("Connection error for {}".format(uniprot))
                continue

args.output_file.write(out_text.strip())
