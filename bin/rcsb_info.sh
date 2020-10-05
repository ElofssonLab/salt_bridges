#!/bin/bash

echo "Pos	Charge	ID	Name";
while read row;
do
	p=$(echo $row | cut -d' ' -f1)	
	pid=$(echo $row | cut -d' ' -f2)	
	tp=$(echo $row | cut -d' ' -f3)	
	rcsb_id=$(echo $pid | cut -c-4)
	# echo $rcsb_id
	echo -n $p '	'$tp '	'$pid '	'
	curl -s https://files.rcsb.org/header/${rcsb_id}.pdb | grep TITLE | cut -c11- | tr -d '\n'
	echo
done <$1
