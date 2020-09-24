DATADIR := data
STATDIR := stats
IMAGEDIR := images
RAWDIR := $(DATADIR)/raw
TOPCONSDIR := $(DATADIR)/topcons
SCAMPIDIR := $(DATADIR)/scampi
OPMDIR := $(DATADIR)/OPM
PDBTMDIR := $(DATADIR)/pdbtm
PROCDIR := $(DATADIR)/processed
3LINES := $(addprefix ${PROCDIR}/,Topcons_Globular.3line Topcons_TMs.3line Scampi_Globular.3line Scampi_TMs.3line pdbtm.3line opm.3line)
3LINESSCAMPI := $(addprefix ${SCAMPIDIR}/,Globular.3line TMs.3line)
3LINESTOPCONS := $(addprefix ${TOPCONSDIR}/,Globular.3line TMs.3line)
3LINESPDBTM := $(PDBTMDIR)/pdbtm.3line
3LINESOPM := $(OPMDIR)/opm.3line
3LINESCLUST := $(addprefix ${PROCDIR}/,Topcons_Globular.clust.3line Topcons_TMs.clust.3line Scampi_Globular.clust.3line Scampi_TMs.clust.3line pdbtm.clust.3line opm.clust.3line)
FASTAS := $(addprefix ${PROCDIR}/,Topcons_Globular.fa Topcons_TMs.fa Scampi_Globular.fa Scampi_TMs.fa pdbtm.fa opm.fa)
CLUST := $(addprefix ${PROCDIR}/,Topcons_Globular.clust.fa Topcons_TMs.clust.fa Scampi_Globular.clust.fa Scampi_TMs.clust.fa pdbtm.clust.fa opm.clust.fa)
MEMS := $(addprefix ${PROCDIR}/,Topcons_Globular.clust.mems.pickle Topcons_TMs.clust.mems.pickle Scampi_Globular.clust.mems.pickle Scampi_TMs.clust.mems.pickle pdbtm.clust.mems.pickle opm.clust.mems.pickle)
STATS := $(addprefix ${STATDIR}/,Topcons_Globular_stats.txt Topcons_TMs_stats.txt Scampi_Globular_stats.txt Scampi_TMs_stats.txt pdbtm_stats.txt opm_stats.txt)



all: $(STATS)

$(RAWDIR)/opm_poly.json:
	wget -O $@ https://lomize-group-opm.herokuapp.com/classtypes/1/primary_structures?pageSize=3000

$(RAWDIR)/opm_bi.json:
	wget -O $@ https://lomize-group-opm.herokuapp.com/classtypes/11/primary_structures?pageSize=3000

$(RAWDIR)/pdbtm_alpha_entries.xml:
	wget -O $@ http://pdbtm.enzim.hu/data/pdbtmalpha

$(RAWDIR)/pdbtm_non_redundant_alpha_struct.txt:
	curl http://pdbtm.enzim.hu/data/pdbtm_alpha_nr.list | tr '[:lower:]' '[:upper:]' | tr -d _ | sort | uniq  > $(RAWDIR)/pdbtm_non_redundant_alpha_struct.txt

$(3LINESPDBTM) : $(RAWDIR)/pdbtm_non_redundant_alpha_struct.txt $(RAWDIR)/pdbtm_alpha_entries.xml
	mkdir -p $(PDBTMDIR)
	./bin/parse_pdbtm.py $^ $@

$(RAWDIR)/ss.txt.gz:
	wget -P $(RAWDIR)/ https://cdn.rcsb.org/etl/kabschSander/ss.txt.gz

$(RAWDIR)/pdb_chain_uniprot.tsv.gz: 
	wget -P $(RAWDIR) ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz

$(RAWDIR)/TOPCONS.zip:
	wget -O $(RAWDIR)/TOPCONS.zip http://topcons.net/static/download/TOPCONS2.0_datasets.zip

$(3LINESOPM): $(RAWDIR)/pdb_chain_uniprot.tsv.gz $(RAWDIR)/opm_bi.json $(RAWDIR)/opm_poly.json
	mkdir -p $(OPMDIR)
	./bin/opm_uniprot_api.py $(RAWDIR)/opm_bi.json $< $(OPMDIR)/opm_bi.3line
	./bin/opm_uniprot_api.py $(RAWDIR)/opm_poly.json $< $(OPMDIR)/opm_poly.3line
	grep -h "" $(OPMDIR)/opm_bi.3line $(OPMDIR)/opm_poly.3line > $@

$(3LINESSCAMPI): $(RAWDIR)/pdb_chain_uniprot.tsv.gz $(RAWDIR)/ss.txt.gz $(RAWDIR)/SCAMPI.zip
	mkdir -p $(SCAMPIDIR)
	unzip -o $(RAWDIR)/SCAMPI.zip -d $(SCAMPIDIR)
	# Add in PDB ids for the SP+TM.3line, keep the original topology annotation
	./bin/add_uniprot_to_3line.py $(SCAMPIDIR)/membranes.3line $(RAWDIR)/pdb_chain_uniprot.tsv.gz $(SCAMPIDIR)/TMs.3line
	# Generate up lists of uniprot IDs to map against PDB and ss from RCSB
	./bin/gen_uniprot_list.sh $(SCAMPIDIR)/globular.3line > $(SCAMPIDIR)/globularlist.txt
	# Actually generate the 3line from the list, use the -t H flag to use H for alpha helix topology
	./bin/make_3line_from_uniprot_list.py $(SCAMPIDIR)/globularlist.txt $(RAWDIR)/pdb_chain_uniprot.tsv.gz $(RAWDIR)/ss.txt.gz $(SCAMPIDIR)/Globular.3line -t H
	# Make sure all proteins contain at least one helix
	./bin/remove_nonTMs_from_3line.py $(SCAMPIDIR)/Globular.3line
	# Clean up
	rm -f $(SCAMPIDIR)/membranes.3line
	rm -f $(SCAMPIDIR)/globularlist.txt
	rm -f $(SCAMPIDIR)/globular.3line

$(3LINESTOPCONS): $(RAWDIR)/pdb_chain_uniprot.tsv.gz $(RAWDIR)/ss.txt.gz $(RAWDIR)/TOPCONS.zip
	mkdir -p $(TOPCONSDIR)
	unzip -o $(RAWDIR)/TOPCONS.zip -d $(RAWDIR)
	# Generate up lists of uniprot IDs to map against PDB and ss from RCSB
	./bin/gen_uniprot_list.sh $(RAWDIR)/TOPCONS2_datasets/Globular.3line > $(TOPCONSDIR)/g1.txt
	./bin/gen_uniprot_list.sh $(RAWDIR)/TOPCONS2_datasets/Globular+SP.3line > $(TOPCONSDIR)/g2.txt
	grep -h "" $(TOPCONSDIR)/g1.txt $(TOPCONSDIR)/g2.txt > $(TOPCONSDIR)/globularlist.txt
	# Actually generate the 3line from the list, use the -t H flag to use H for alpha helix topology
	./bin/make_3line_from_uniprot_list.py $(TOPCONSDIR)/globularlist.txt $(RAWDIR)/pdb_chain_uniprot.tsv.gz $(RAWDIR)/ss.txt.gz $(TOPCONSDIR)/Globular.3line -t H
	# Make sure all proteins contain at least one helix
	./bin/remove_nonTMs_from_3line.py $(TOPCONSDIR)/Globular.3line
	# Add in PDB ids for the SP+TM.3line, keep the original topology annotation
	# Merge TM and SPTM together
	./bin/add_pdbid_to_3line.py $(RAWDIR)/TOPCONS2_datasets/SP+TM.3line $(RAWDIR)/pdb_chain_uniprot.tsv.gz $(TOPCONSDIR)/SPTM.3line
	grep -h "" $(RAWDIR)/TOPCONS2_datasets/TM.3line $(TOPCONSDIR)/SPTM.3line > $(TOPCONSDIR)/TMs.3line
	# Clean up
	rm -rf $(RAWDIR)/TOPCONS2_datasets
	rm -f $(TOPCONSDIR)/globularlist.txt
	rm -f $(TOPCONSDIR)/memtemp
	rm -f $(TOPCONSDIR)/g1.txt
	rm -f $(TOPCONSDIR)/g2.txt
	rm -f $(TOPCONSDIR)/SPTM.3line

$(3LINES): $(3LINESSCAMPI) $(3LINESTOPCONS) $(3LINESPDBTM) $(3LINESOPM)
	mkdir -p $(PROCDIR)
	cp $(SCAMPIDIR)/TMs.3line $(PROCDIR)/Scampi_TMs.3line
	cp $(SCAMPIDIR)/Globular.3line $(PROCDIR)/Scampi_Globular.3line
	cp $(TOPCONSDIR)/TMs.3line $(PROCDIR)/Topcons_TMs.3line
	cp $(TOPCONSDIR)/Globular.3line $(PROCDIR)/Topcons_Globular.3line
	cp $(3LINESPDBTM) $(PROCDIR)/pdbtm.3line
	cp $(3LINESOPM) $(PROCDIR)/opm.3line

$(FASTAS) : %.fa: %.3line | $(3LINES)
	sed '3~3d' $< > $@

$(CLUST) : %.clust.fa: %.fa | $(FASTAS)
	cd-hit -i $< -o $@ -c 0.4 -n 2 -T 0 -M 0 -d 0

$(3LINESCLUST) : %.clust.3line: %.clust.fa | $(CLUST)
	./bin/add_topo_from_3line.py $< $(subst .clust,,$@) > $@

$(MEMS) : %.clust.mems.pickle: %.clust.3line | $(3LINESCLUST)
	./bin/make_mems_from_3line.py $< $@

$(STATS) : $(STATDIR)/%_stats.txt : $(PROCDIR)/%.clust.mems.pickle | $(MEMS)
	./bin/statsCharges.py $< > $@

.PHONY: clean deepclean
clean:
	rm -rf $(TOPCONSDIR)
	rm -rf $(SCAMPIDIR)
	rm -rf $(PDBTMDIR)
	rm -rf $(PROCDIR)
	rm -r $(IMAGEDIR)/*
	rm -r $(STATS)
deepclean: clean
	rm -rf $(OPMDIR)
	rm -rf $(RAWDIR)/pdb_chain_uniprot.tsv.gz
	rm -rf $(RAWDIR)/ss.txt.gz
	rm -rf $(RAWDIR)/TOPCONS.zip
	rm -rf $(RAWDIR)/pdbtm_alpha_entries.xml
	rm -rf $(RAWDIR)/pdbtm_non_redundant_alpha_struct.txt
