#!/usr/bin/env python3
import argparse
import gzip
import re
import sys

parser = argparse.ArgumentParser()

parser.add_argument("in_file", type=str, help="3line file")

args = parser.parse_args()


out_text = ""
with open(args.in_file, 'r') as in_handle:
    row = 0
    for line in in_handle:
        if line.startswith('>'):
            header = line
        elif row % 3 == 1:
            seq = line
        elif row % 3 == 2:
            if re.search("[^O]",line.strip()):
                out_text += header + seq + line
        row += 1

with open(args.in_file, 'w') as out_handle:
    out_handle.write(out_text)

