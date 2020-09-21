DATADIR := data
TOPCONSDIR := $(DATADIR)/topcons
3LINES := $(addprefix ${TOPCONSDIR}/,Globular.3line TMs.3line)
3LINESCLUST := $(addprefix ${TOPCONSDIR}/,Globular.clust.3line TMs.clust.3line)
FASTAS := $(addprefix ${TOPCONSDIR}/,Globular.fa TMs.fa)
CLUST := $(addprefix ${TOPCONSDIR}/,Globular.clust.fa TMs.clust.fa)

all: $(DATADIR)/pdb_chain_uniprot.tsv.gz\
	 $(DATADIR)/ss.txt.gz\
	 $(DATADIR)/TOPCONS.zip\
	 $(3LINES)\
	 $(3LINESCLUST)\
	 $(FASTAS)\
	 $(CLUST)

$(DATADIR)/ss.txt.gz:
	wget -P data/ https://cdn.rcsb.org/etl/kabschSander/ss.txt.gz

$(DATADIR)/pdb_chain_uniprot.tsv.gz: 
	wget -P data/ ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz

$(DATADIR)/TOPCONS.zip:
	wget -O data/TOPCONS.zip http://topcons.net/static/download/TOPCONS2.0_datasets.zip

$(3LINES): $(DATADIR)/pdb_chain_uniprot.tsv.gz $(DATADIR)/ss.txt.gz $(DATADIR)/TOPCONS.zip
	mkdir -p $(TOPCONSDIR)
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
	rm -f data/topcons/memtemp
	rm -f data/topcons/g1.txt
	rm -f data/topcons/g2.txt
	rm -f data/topcons/SPTM.3line

$(FASTAS) : %.fa: %.3line
	sed '3~3d' $< > $@

$(CLUST) : %.clust.fa: %.fa
	cd-hit -i $< -o $@ -c 0.5 -n 3 -T 0 -M 0 -d 0

$(3LINESCLUST) : %.clust.3line: %.clust.fa
	./bin/add_topo_from_3line.py $< $(subst .clust,,$@) > $@

.PHONY: clean deepclean
clean:
	rm -rf $(TOPCONSDIR)
deepclean:
	rm -rf $(DATADIR)/pdb_chain_uniprot.tsv.gz
	rm -rf $(DATADIR)/ss.txt.gz
	rm -rf $(DATADIR)/TOPCONS.zip
	rm -rf $(TOPCONSDIR)
