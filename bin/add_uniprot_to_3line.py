#!/usr/bin/env python3
import argparse
import gzip
import sys

parser = argparse.ArgumentParser()

parser.add_argument("input_file", type=argparse.FileType('r'), help="Input 3line file")
parser.add_argument("map_file", type=str, help="pdb_chain_uniprot.tsv.gz")
parser.add_argument("output_file", type=argparse.FileType('w'), help="Output 3line file")

args = parser.parse_args()
pdb_to_uniprot_map = {}

with gzip.open(args.map_file, 'rb') as map_handle:
    # Read the two header lines
    map_handle.readline()
    map_handle.readline()

    for line in map_handle.readlines():
        pdb_id, chain, uniprot_id = str(line, 'utf-8').strip().split()[:3]
        full_pdb_id = pdb_id.upper() + chain

        if not full_pdb_id in pdb_to_uniprot_map:
            pdb_to_uniprot_map[full_pdb_id] = [uniprot_id]
        else:
            pdb_to_uniprot_map[full_pdb_id].append(uniprot_id)

out_text = ""
for i, line in enumerate(args.input_file):
    if line.startswith('>'):
        pdb_id = line[1:].strip().upper()
        if pdb_id in pdb_to_uniprot_map:
            uniprot_id = pdb_to_uniprot_map[pdb_id]
        else:
            uniprot_id = False
    elif i % 3 == 1:
        seq = line.strip()
    elif i % 3 == 2:
        topo = line.strip()
        if uniprot_id:
            out_text += ">" + uniprot_id[0] + "|" + pdb_id + '\n'
            out_text += seq + '\n'
            out_text += topo + '\n'
args.output_file.write(out_text.strip())
