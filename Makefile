.PHONY: all clean
all: data/pdb_chain_uniprot.tsv.gz\
	 data/ss.txt.gz\
	 data/TOPCONS.zip\
	 data/topcons/TMs.3line\
	 data/topcons/Globular.3line

data/ss.txt.gz:
	wget -P data/ https://cdn.rcsb.org/etl/kabschSander/ss.txt.gz

data/pdb_chain_uniprot.tsv.gz: 
	wget -P data/ ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz

data/TOPCONS.zip:
	wget -O data/TOPCONS.zip http://topcons.net/static/download/TOPCONS2.0_datasets.zip

data/topcons/TMs.3line data/topcons/Globular.3line:
	mkdir -p data/topcons
	unzip -o data/TOPCONS.zip -d data/
	# Generate up lists of uniprot IDs to map against PDB and ss from RCSB
	./bin/gen_uniprot_list.sh data/TOPCONS2_datasets/Globular.3line > data/topcons/g1.txt
	./bin/gen_uniprot_list.sh data/TOPCONS2_datasets/Globular+SP.3line > data/topcons/g2.txt
	grep -h "" data/topcons/g1.txt data/topcons/g2.txt > data/topcons/globularlist.txt
	# Actually generate the 3line from the list
	./bin/make_3line_from_uniprot_list.py data/topcons/globularlist.txt data/pdb_chain_uniprot.tsv.gz data/ss.txt.gz data/topcons/Globular.3line
	# Add in PDB ids for the SP+TM.3line, keep the original topology annotation
	# Merge TM and SPTM together
	./bin/add_pdbid_to_3line.py data/TOPCONS2_datasets/SP+TM.3line data/pdb_chain_uniprot.tsv.gz data/topcons/SPTM.3line
	grep -h "" data/TOPCONS2_datasets/TM.3line data/topcons/SPTM.3line > data/topcons/TMs.3line
	# Clean up
	rm -rf data/TOPCONS2_datasets
	rm -f data/topcons/globularlist.txt
	rm -f data/topcons/memtemp
	rm -f data/topcons/g1.txt
	rm -f data/topcons/g2.txt
	rm -f data/topcons/SPTM.3line

clean:
	rm -rf data/pdb_chain_uniprot.tsv.gz
	rm -rf ss.txt.gz
	rm -rf data/TOPCONS.zip
	rm -rf data/topcons
