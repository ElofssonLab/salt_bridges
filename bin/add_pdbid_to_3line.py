#!/usr/bin/env python3
import argparse
import gzip
import sys

parser = argparse.ArgumentParser()

parser.add_argument("input_file", type=argparse.FileType('r'), help="Input 3line file")
parser.add_argument("map_file", type=str, help="pdb_chain_uniprot.tsv.gz")
parser.add_argument("output_file", type=argparse.FileType('w'), help="Output 3line file")

args = parser.parse_args()
uniprot_to_pdb_map = {}

with gzip.open(args.map_file, 'rb') as map_handle:
    # Read the two header lines
    map_handle.readline()
    map_handle.readline()

    for line in map_handle.readlines():
        pdb_id, chain, uniprot_id = str(line, 'utf-8').strip().split()[:3]
        full_pdb_id = pdb_id.upper() + chain

        if not uniprot_id in uniprot_to_pdb_map:
            uniprot_to_pdb_map[uniprot_id] = [full_pdb_id]
        else:
            uniprot_to_pdb_map[uniprot_id].append(full_pdb_id)

out_text = ""
for i, line in enumerate(args.input_file):
    if line.startswith('>'):
        uniprot_id = line[1:].strip()
        if uniprot_id in uniprot_to_pdb_map:
            pdb_id = uniprot_to_pdb_map[uniprot_id]
        else:
            pdb_id = False
    elif i % 3 == 1:
        seq = line.strip()
    elif i % 3 == 2:
        topo = line.strip()
        if pdb_id:
            out_text += ">" + uniprot_id + "|" + pdb_id[0] + '\n'
            out_text += seq + '\n'
            out_text += topo + '\n'
args.output_file.write(out_text.strip())
