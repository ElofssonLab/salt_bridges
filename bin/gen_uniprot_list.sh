#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "No input given, fasta-style file needed"
	exit 0
fi

grep -e "^>" $1 | cut -c2- | cut -f1 -d'|'
