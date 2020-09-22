DATADIR := data
TOPCONSDIR := $(DATADIR)/topcons
SCAMPIDIR := $(DATADIR)/scampi
PROCDIR := $(DATADIR)/processed
3LINES := $(addprefix ${PROCDIR}/,Globular.3line TMs.3line)
3LINESSCAMPI := $(addprefix ${SCAMPIDIR}/,Globular.3line TMs.3line)
3LINESTOPCONS := $(addprefix ${TOPCONSDIR}/,Globular.3line TMs.3line)
3LINESCLUST := $(addprefix ${PROCDIR}/,Globular.clust.3line TMs.clust.3line)
FASTAS := $(addprefix ${PROCDIR}/,Globular.fa TMs.fa)
CLUST := $(addprefix ${PROCDIR}/,Globular.clust.fa TMs.clust.fa)
MEMS := $(addprefix ${PROCDIR}/,Globular.clust.mems.pickle TMs.clust.mems.pickle)
STATS := $(addprefix ${PROCDIR}/,Globular_stats.txt TMs_stats.txt)

all: $(STATS)

$(DATADIR)/ss.txt.gz:
	wget -P data/ https://cdn.rcsb.org/etl/kabschSander/ss.txt.gz

$(DATADIR)/pdb_chain_uniprot.tsv.gz: 
	wget -P data/ ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz

$(DATADIR)/TOPCONS.zip:
	wget -O data/TOPCONS.zip http://topcons.net/static/download/TOPCONS2.0_datasets.zip

$(3LINESSCAMPI): $(DATADIR)/SCAMPI.zip
	mkdir -p $(SCAMPIDIR)
	unzip -o $(DATADIR)/SCAMPI.zip -d $(DATADIR)/scampi
	# Add in PDB ids for the SP+TM.3line, keep the original topology annotation
	./bin/add_uniprot_to_3line.py $(SCAMPIDIR)/membranes.3line $(DATADIR)/pdb_chain_uniprot.tsv.gz $(SCAMPIDIR)/TMs.3line
	# Generate up lists of uniprot IDs to map against PDB and ss from RCSB
	./bin/gen_uniprot_list.sh $(SCAMPIDIR)/globular.3line > $(SCAMPIDIR)/globularlist.txt
	# Actually generate the 3line from the list, use the -t H flag to use H for alpha helix topology
	./bin/make_3line_from_uniprot_list.py $(SCAMPIDIR)/globularlist.txt $(DATADIR)/pdb_chain_uniprot.tsv.gz $(DATADIR)/ss.txt.gz $(SCAMPIDIR)/Globular.3line -t H
	# Make sure all proteins contain at least one helix
	./bin/remove_nonTMs_from_3line.py $(SCAMPIDIR)/Globular.3line
	# Clean up
	rm -f $(SCAMPIDIR)/membranes.3line
	rm -f $(SCAMPIDIR)/globularlist.txt
	rm -f $(SCAMPIDIR)/globular.3line

$(3LINESTOPCONS): $(DATADIR)/pdb_chain_uniprot.tsv.gz $(DATADIR)/ss.txt.gz $(DATADIR)/TOPCONS.zip
	mkdir -p $(TOPCONSDIR)
	unzip -o $(DATADIR)/TOPCONS.zip -d $(DATADIR)
	# Generate up lists of uniprot IDs to map against PDB and ss from RCSB
	./bin/gen_uniprot_list.sh data/TOPCONS2_datasets/Globular.3line > data/topcons/g1.txt
	./bin/gen_uniprot_list.sh data/TOPCONS2_datasets/Globular+SP.3line > data/topcons/g2.txt
	grep -h "" data/topcons/g1.txt data/topcons/g2.txt > data/topcons/globularlist.txt
	# Actually generate the 3line from the list, use the -t H flag to use H for alpha helix topology
	./bin/make_3line_from_uniprot_list.py $(TOPCONSDIR)/globularlist.txt $(DATADIR)/pdb_chain_uniprot.tsv.gz $(DATADIR)/ss.txt.gz $(TOPCONSDIR)/Globular.3line -t H
	# Make sure all proteins contain at least one helix
	./bin/remove_nonTMs_from_3line.py $(TOPCONSDIR)/Globular.3line
	# Add in PDB ids for the SP+TM.3line, keep the original topology annotation
	# Merge TM and SPTM together
	./bin/add_pdbid_to_3line.py $(DATADIR)/TOPCONS2_datasets/SP+TM.3line $(DATADIR)/pdb_chain_uniprot.tsv.gz $(TOPCONSDIR)/SPTM.3line
	grep -h "" $(DATADIR)/TOPCONS2_datasets/TM.3line $(TOPCONSDIR)/SPTM.3line > $(TOPCONSDIR)/TMs.3line
	# Clean up
	rm -rf $(DATADIR)/TOPCONS2_datasets
	rm -f $(TOPCONSDIR)/globularlist.txt
	rm -f $(TOPCONSDIR)/memtemp
	rm -f $(TOPCONSDIR)/g1.txt
	rm -f $(TOPCONSDIR)/g2.txt
	rm -f $(TOPCONSDIR)/SPTM.3line

$(3LINES): $(3LINESSCAMPI) $(3LINESTOPCONS)
	mkdir -p $(PROCDIR)
	grep -h "" $(SCAMPIDIR)/TMs.3line $(TOPCONSDIR)/TMs.3line > $(PROCDIR)/TMs.3line
	grep -h "" $(SCAMPIDIR)/Globular.3line $(TOPCONSDIR)/Globular.3line > $(PROCDIR)/Globular.3line

$(FASTAS) : %.fa: %.3line | $(3LINES)
	sed '3~3d' $< > $@

$(CLUST) : %.clust.fa: %.fa | $(FASTAS)
	cd-hit -i $< -o $@ -c 0.4 -n 2 -T 0 -M 0 -d 0

$(3LINESCLUST) : %.clust.3line: %.clust.fa | $(CLUST)
	./bin/add_topo_from_3line.py $< $(subst .clust,,$@) > $@

$(MEMS) : %.clust.mems.pickle: %.clust.3line | $(3LINESCLUST)
	./bin/make_mems_from_3line.py $< $@

$(STATS) : %_stats.txt: %.clust.mems.pickle | $(MEMS)
	./bin/statsCharges.py $< > $@

.PHONY: clean deepclean
clean:
	rm -rf $(TOPCONSDIR)
	rm -rf $(SCAMPIDIR)
	rm -rf $(PROCDIR)
deepclean:
	rm -rf $(DATADIR)/pdb_chain_uniprot.tsv.gz
	rm -rf $(DATADIR)/ss.txt.gz
	rm -rf $(DATADIR)/TOPCONS.zip
	rm -rf $(TOPCONSDIR)
	rm -rf $(SCAMPIDIR)
