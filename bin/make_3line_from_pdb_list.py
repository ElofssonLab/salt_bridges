#!/usr/bin/env python3
import argparse
import gzip
import re
import sys

parser = argparse.ArgumentParser()

parser.add_argument("pdb_list", type=argparse.FileType('r'), help="List of pdb chain IDs")
# parser.add_argument("map_file", type=str, help="pdb_chain_uniprot.tsv.gz")
parser.add_argument("ss_file", type=str, help="ss.txt.gzip")
parser.add_argument("out_file", type=argparse.FileType('w'), help="Output file")
parser.add_argument("-t", "--helix_char", type=str, default='M', help="Topo char")

args = parser.parse_args()


def read_ss_file(ss_filename):
    """
    Reads in an ss.txt file from pdb with secondary structures.
    - H = alpha helix
    - B = residue in isolated beta-bridge
    - E = extended strand, participates in beta ladder
    - G = 3-helix (3/10 helix)
    - I = 5 helix (pi helix)
    - T = hydrogen bonded turn
    - S = bend
    """
    structures = {}
    with gzip.open(ss_filename, 'rb') as ssFile:
        pid = ''
        gotSeq = False
        seq = ''
        struct = ''
        first = True
        for line in ssFile:
            strippedLine = str(line, 'utf-8').strip('\n')
            if strippedLine.startswith('>'):
                if strippedLine.endswith('sequence'):
                    if not first:
                        if len(seq) != len(struct):
                            print("Missmatch in length for ", pid)
                            sys.exit(1)
                        # row = [pid, seq, struct]
                        if pid in structures:
                            print("Got duplicates")
                        sub_struct = re.sub("[H]", args.helix_char, struct)
                        structures[pid] = [seq, sub_struct]
                        pid = ''
                        seq = ''
                        struct = ''
                        gotSeq = False
                    else:
                        first = False
                    pid = ''.join(strippedLine[1:].split(':')[:2])
                else:
                    gotSeq = True
            elif not gotSeq:
                seq += strippedLine
            else:
                struct += re.sub("[^H]", "O", strippedLine)  # Only alpha helix, not 3/10 or pi
        # And the last structure as well
        sub_struct = re.sub("[H]", args.helix_char, struct)
        structures[pid] = [seq, sub_struct]
    return structures

seq_and_ss_structs = read_ss_file(args.ss_file)

# pdb_to_uniprot_map = {}

# with gzip.open(args.map_file, 'rb') as map_handle:
#     # Read the two header lines
#     map_handle.readline()
#     map_handle.readline()
# 
#     for line in map_handle.readlines():
#         pdb_id, chain, uniprot_id = str(line, 'utf-8').strip().split()[:3]
#         full_pdb_id = pdb_id.upper() + chain
# 
#         if full_pdb_id in pdb_to_uniprot_map:
#             if pdb_to_uniprot_map[full_pdb_id] != uniprot_id:
#                 print("Error", full_pdb_id, "already in map")
#                 sys.exit()
#         pdb_to_uniprot_map[full_pdb_id] = uniprot_id
#         # if not uniprot_id in uniprot_to_pdb_map:
#         #     uniprot_to_pdb_map[uniprot_id] = [full_pdb_id]
#         # else:
#         #     uniprot_to_pdb_map[uniprot_id].append(full_pdb_id)
# sys.exit()
mapped_ids = []

out_text = ""
args.pdb_list.readline()
for line in args.pdb_list:
    pdb_id = line.strip().split()[0]
    header = '>' + pdb_id
    if not pdb_id in seq_and_ss_structs:
        print("Lacking SS {}".format(pdb_id))
        continue
    seq, sec_struct = seq_and_ss_structs[pdb_id]
    out_text += header + '\n' + seq + '\n' + sec_struct + '\n'

args.out_file.write(out_text.strip())
