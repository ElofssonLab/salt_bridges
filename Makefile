DATADIR := data
STATDIR := stats
IMAGEDIR := images
RAWDIR := $(DATADIR)/raw
# TOPCONSDIR := $(DATADIR)/topcons
# SCAMPIDIR := $(DATADIR)/scampi
# OPMDIR := $(DATADIR)/OPM
PDBTMDIR := $(DATADIR)/pdbtm
PROCDIR := $(DATADIR)/processed
3LINES := $(addprefix ${PROCDIR}/, pdbtm.3line)
3LINESRED := $(addprefix ${PROCDIR}/, pdbtm_redundant.3line)
# 3LINESSCAMPI := $(addprefix ${SCAMPIDIR}/,Globular.3line TMs.3line)
# 3LINESTOPCONS := $(addprefix ${TOPCONSDIR}/,Globular.3line TMs.3line)
3LINESPDBTM := $(addprefix ${PDBTMDIR}/, pdbtm.3line)
3LINESPDBTMRED:= $(addprefix ${PDBTMDIR}/, pdbtm_redundant.3line)
3LINESGLOB := $(PROCDIR)/scop_glob.3line
# 3LINESOPM := $(OPMDIR)/opm.3line
3LINESCLUST := $(addprefix ${PROCDIR}/, pdbtm.clust.3line scop_glob.clust.3line)
FASTAS := $(addprefix ${PROCDIR}/, pdbtm.fa scop_glob.fa)
CLUST := $(addprefix ${PROCDIR}/, pdbtm.clust.fa scop_glob.clust.fa)
MEMS := $(addprefix ${PROCDIR}/, pdbtm.clust.mems.pickle scop_glob.clust.mems.pickle)
A3MMEMS := $(addprefix ${PROCDIR}/, pdbtm_a3m.mems.pickle scop_glob_a3m.mems.pickle)
MEMSRED := $(PROCDIR)/pdbtm_redundant.mems.pickle
MEMSGLOB := $(PROCDIR)/scop_glob.mems.pickle
STATS := $(addprefix ${STATDIR}/,pdbtm_stats.txt)
# LOGODDS := $(addprefix ${IMAGEDIR}/,Topcons_TMs_logodds.svg pdbtm_logodds.svg)
LISTS := $(addprefix ${STATDIR}/,pdbtm_redundant_1_list.txt ) # pdbtm_redundant_1_tolerant_list.txt pdbtm_redundant_2_list.txt pdbtm_redundant_1_all_list.txt pdbtm_redundant_1_all_tolerant_list.txt pdbtm_redundant_2_all_list.txt)
PROTS := $(addprefix ${STATDIR}/,pdbtm_same_pairs.txt pdbtm_opp_pairs.txt)
RCSB := $(addprefix ${STATDIR}/,pdbtm_same_info.txt pdbtm_opp_info.txt)
GLOBCHARGES := $(PROCDIR)/scop_glob.charges.pickle
# MEMCHARGES := $(PROCDIR)/pdbtm.clust.charges.pickle
A3MCHARGES := $(PROCDIR)/pdbtm_a3m.charges.pickle $(PROCDIR)/scop_glob_a3m.charges.pickle
VISIMAGES := $(addprefix ${IMAGEDIR}/, pdbtm_vis.svg pdbtm_vis.png scop_glob_vis.svg scop_glob_vis.png mem_cluster.svg mem_cluster.png mem_cluster_full.svg mem_cluster_full.png pdbtm_pairs.png pdbtm_pairs.svg pdbtm_stats.png pdbtm_stats.svg)

all:  $(STATS)  $(VISIMAGES) $(LISTS) # $(RCSB)
images:  $(VISIMAGES) # $(RCSB)
# $(RAWDIR)/opm_poly.json:
# 	wget -O $@ https://lomize-group-opm.herokuapp.com/classtypes/1/primary_structures?pageSize=3000
# 
# $(RAWDIR)/opm_bi.json:
# 	wget -O $@ https://lomize-group-opm.herokuapp.com/classtypes/11/primary_structures?pageSize=3000

$(RAWDIR)/pdbtm_alpha_entries.xml:
	wget -O $@ http://pdbtm.enzim.hu/data/pdbtmalpha

$(RAWDIR)/pdbtm_non_redundant_alpha_list.txt:
	curl http://pdbtm.enzim.hu/data/pdbtm_alpha_nr.list | tr '[:lower:]' '[:upper:]' | tr -d _ | sort | uniq  > $(RAWDIR)/pdbtm_non_redundant_alpha_list.txt

$(RAWDIR)/pdbtm_redundant_alpha_list.txt:
	curl http://pdbtm.enzim.hu/data/pdbtm_alpha.list | tr '[:lower:]' '[:upper:]' | tr -d _ | sort | uniq  > $(RAWDIR)/pdbtm_redundant_alpha_list.txt

$(3LINESPDBTM) : $(RAWDIR)/pdbtm_non_redundant_alpha_list.txt $(RAWDIR)/pdbtm_alpha_entries.xml 
	mkdir -p $(PDBTMDIR)
	./bin/parse_pdbtm.py $^ $@

$(3LINESPDBTMRED) : $(RAWDIR)/pdbtm_redundant_alpha_list.txt $(RAWDIR)/pdbtm_alpha_entries.xml 
	mkdir -p $(PDBTMDIR)
	./bin/parse_pdbtm.py $^ $@
	cp $(3LINESPDBTMRED) $(PROCDIR)/
	cp $(PDBTMDIR)/pdbtm_redundant_bridges.pickle $(PROCDIR)/


$(RAWDIR)/scop_globular_alpha.txt:
	wget -O $(RAWDIR)/scop_raw.txt http://scop.mrc-lmb.cam.ac.uk/files/scop-cla-latest.txt
	grep -e "TP=1," $(RAWDIR)/scop_raw.txt | grep -e "CL=1000000" | cut -c9-14 | tr -d ' ' |sort |uniq > $(RAWDIR)/scop_globular_alpha.txt
	# rm $(RAWDIR)/scop_raw.txt

$(PROCDIR)/scop_globular_alpha_reduced.txt: $(RAWDIR)/scop_globular_alpha.txt
	comm -23 $< $(RAWDIR)/pdbtm_redundant_alpha_list.txt > $@

$(3LINESGLOB): $(PROCDIR)/scop_globular_alpha_reduced.txt $(RAWDIR)/ss.txt.gz
	mkdir -p $(PROCDIR)
	./bin/make_3line_from_pdb_list.py $< $(RAWDIR)/ss.txt.gz $@

$(RAWDIR)/ss.txt.gz:
	wget -P $(RAWDIR)/ https://cdn.rcsb.org/etl/kabschSander/ss.txt.gz

$(RAWDIR)/pdb_chain_uniprot.tsv.gz: 
	wget -P $(RAWDIR) ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz


# $(3LINESOPM): $(RAWDIR)/pdb_chain_uniprot.tsv.gz $(RAWDIR)/opm_bi.json $(RAWDIR)/opm_poly.json
# 	mkdir -p $(OPMDIR)
# 	./bin/opm_uniprot_api.py $(RAWDIR)/opm_bi.json $< $(OPMDIR)/opm_bi.3line
# 	./bin/opm_uniprot_api.py $(RAWDIR)/opm_poly.json $< $(OPMDIR)/opm_poly.3line
# 	grep -h "" $(OPMDIR)/opm_bi.3line $(OPMDIR)/opm_poly.3line > $@

# $(3LINESSCAMPI): $(RAWDIR)/pdb_chain_uniprot.tsv.gz $(RAWDIR)/ss.txt.gz $(RAWDIR)/SCAMPI.zip
# 	mkdir -p $(SCAMPIDIR)
# 	unzip -o $(RAWDIR)/SCAMPI.zip -d $(SCAMPIDIR)
# 	# Add in PDB ids for the SP+TM.3line, keep the original topology annotation
# 	./bin/add_uniprot_to_3line.py $(SCAMPIDIR)/membranes.3line $(RAWDIR)/pdb_chain_uniprot.tsv.gz $(SCAMPIDIR)/TMs.3line
# 	# Generate up lists of uniprot IDs to map against PDB and ss from RCSB
# 	./bin/gen_uniprot_list.sh $(SCAMPIDIR)/globular.3line > $(SCAMPIDIR)/globularlist.txt
# 	# Actually generate the 3line from the list, use the -t H flag to use H for alpha helix topology
# 	./bin/make_3line_from_uniprot_list.py $(SCAMPIDIR)/globularlist.txt $(RAWDIR)/pdb_chain_uniprot.tsv.gz $(RAWDIR)/ss.txt.gz $(SCAMPIDIR)/Globular.3line -t H
# 	# Make sure all proteins contain at least one helix
# 	./bin/remove_nonTMs_from_3line.py $(SCAMPIDIR)/Globular.3line
# 	# Clean up
# 	rm -f $(SCAMPIDIR)/membranes.3line
# 	rm -f $(SCAMPIDIR)/globularlist.txt
# 	rm -f $(SCAMPIDIR)/globular.3line


$(3LINES): $(3LINESPDBTM)
	cp $(3LINESPDBTM) $(PROCDIR)/
	cp $(PDBTMDIR)/pdbtm_bridges.pickle $(PROCDIR)/

$(3LINESRED): $(3LINESPDBTMRED)
	cp $< $@

$(FASTAS) : %.fa: %.3line | $(3LINES) $(3LINESGLOB)
	sed '3~3d' $< > $@

$(CLUST) : %.clust.fa: %.fa | $(FASTAS)
	cd-hit -i $< -o $@ -c 0.4 -n 2 -T 0 -M 0 -d 0

$(3LINESCLUST) : %.clust.3line: %.clust.fa | $(CLUST)
	./bin/add_topo_from_3line.py $< $(subst .clust,,$@) > $@

$(MEMS) : %.clust.mems.pickle: %.clust.3line | $(3LINESCLUST)
	./bin/make_mems_from_3line.py $< $@

$(A3MMEMS) : $(3LINESPDBTM) $(3LINESGLOB)
	./bin/make_mems_from_a3m.py $(PROCDIR)/pdbtm.3line $(DATADIR)/pdbtm_a3m $(PROCDIR)/pdbtm_a3m.mems.pickle
	./bin/make_mems_from_a3m.py $(PROCDIR)/scop_glob.clust.3line $(DATADIR)/scop_a3m $(PROCDIR)/scop_glob_a3m.mems.pickle

$(MEMSRED) : %.mems.pickle: %.3line | $(3LINESRED)
	./bin/make_mems_from_3line.py $< $@
	# ./bin/make_mems_from_3line.py $< $(PROCDIR)/pdbtm_redundant_tolerant.mems.pickle --t True

$(MEMSGLOB) : %.mems.pickle: %.3line | $(3LINESGLOB)
	./bin/make_mems_from_3line.py $< $@

$(STATS) : $(STATDIR)/%_stats.txt : $(PROCDIR)/%.clust.mems.pickle | $(MEMS)
	./bin/statsCharges.py $< > $@
	./bin/dataset_analysis.py $(PROCDIR)/pdbtm.clust.mems.pickle $(PROCDIR)/pdbtm_bridges.pickle $(PROCDIR)/pdbtm.clust.3line
	./bin/dataset_analysis.py $(PROCDIR)/scop_glob.clust.mems.pickle $(PROCDIR)/scop_glob.clust.mems_bridges.pickle $(PROCDIR)/scop_glob.clust.3line
	./bin/dataset_analysis.py $(PROCDIR)/pdbtm_a3m.mems.pickle $(PROCDIR)/pdbtm_bridges.pickle $(PROCDIR)/pdbtm.3line -a3m True
	./bin/dataset_analysis.py $(PROCDIR)/scop_glob_a3m.mems.pickle $(PROCDIR)/scop_glob.clust.mems_bridges.pickle $(PROCDIR)/scop_glob.clust.3line -a3m True
	./bin/csv_plots.py $(PROCDIR)/pdbtm.clust.aas.csv $(PROCDIR)/pdbtm.clust.pairs.csv
	./bin/csv_plots.py $(PROCDIR)/pdbtm_a3m.aas.csv $(PROCDIR)/pdbtm_a3m.pairs.csv
	./bin/csv_plots.py $(PROCDIR)/scop_glob.clust.aas.csv $(PROCDIR)/scop_glob.clust.pairs.csv
	./bin/csv_plots.py $(PROCDIR)/scop_glob_a3m.aas.csv $(PROCDIR)/scop_glob_a3m.pairs.csv
	./bin/manual_pairs.py $(PROCDIR)/pdbtm.clust.pairs.csv $(PROCDIR)/pdbtm.clust.mems.pickle $(PROCDIR)/pdbtm_bridges.pickle > $(STATDIR)/manual_pairs_info.txt
	# ./bin/statsCharges.py $(PROCDIR)/pdbtm_redundant.mems.pickle > $(STATDIR)/pdbtm_redundant_stats.txt

# $(RCSB) : %_info.txt : %_pairs.txt | $(PROTS)
# 	./bin/rcsb_info.sh $< > $@

# $(LOGODDS) : $(PROCDIR)/Topcons_TMs.clust.mems.pickle $(PROCDIR)/pdbtm.clust.mems.pickle
# 	./bin/logodd_barchart.py $(PROCDIR)/Topcons_TMs.clust.mems.pickle
# 	./bin/logodd_barchart.py $(PROCDIR)/pdbtm.clust.mems.pickle

$(LISTS) : $(MEMSRED) $(3LINESRED)
	# ./bin/gen_potential_list.py $(PROCDIR)/pdbtm.clust.mems.pickle $(PROCDIR)/pdbtm.clust.3line -b 1 > $(STATDIR)/pdbtm_clust_1_list.txt
	#./bin/gen_potential_list.py $(PROCDIR)/pdbtm_a3m.mems_all_a3m.pickle $(PROCDIR)/pdbtm.3line -b 1 -s True -a3m True > $(STATDIR)/pdbtm_a3m_1_list.txt
	# ./bin/gen_potential_list.py $(PROCDIR)/pdbtm_redundant.mems.pickle $(PROCDIR)/pdbtm_redundant.3line -b 1 > $(STATDIR)/pdbtm_redundant_1_list.txt
	./bin/gen_potential_list.py $(PROCDIR)/pdbtm_redundant_bridges.pickle $(PROCDIR)/pdbtm_redundant.mems.pickle $(PROCDIR)/pdbtm_redundant.3line -b 1 > $(STATDIR)/pdbtm_redundant_1_list.txt
	# ./bin/gen_potential_list.py $(PROCDIR)/pdbtm_redundant.mems.pickle $(PROCDIR)/pdbtm_redundant.3line -b 1 -t False > $(STATDIR)/pdbtm_redundant_strict_1_list.txt
	./bin/gen_potential_list.py $(PROCDIR)/pdbtm_redundant_bridges.pickle $(PROCDIR)/pdbtm_redundant.mems.pickle $(PROCDIR)/pdbtm_redundant.3line -b 1 -t False > $(STATDIR)/pdbtm_redundant_strict_1_list.txt
	# ./bin/gen_potential_list.py $(PROCDIR)/scop_glob.clust.mems.pickle $(PROCDIR)/scop_glob.clust.3line -b 1 -s True > $(STATDIR)/scop_clust_1_list.txt
	#	./bin/gen_potential_list.py $(PROCDIR)/scop_glob_a3m.mems_all_a3m.pickle $(PROCDIR)/scop_glob.clust.3line -b 1 -s True -a3m True> $(STATDIR)/scop_a3m_1_list.txt
# 	./bin/gen_potential_list.py $(PROCDIR)/pdbtm_redundant.mems.pickle $(PROCDIR)/pdbtm_redundant.3line -b 1 -g 0 > $(STATDIR)/pdbtm_redundant_1_all_list.txt
# 	./bin/gen_potential_list.py $(PROCDIR)/pdbtm_redundant_tolerant.mems.pickle $(PROCDIR)/pdbtm_redundant.3line -b 1 > $(STATDIR)/pdbtm_redundant_1_tolerant_list.txt -t True
# 	./bin/gen_potential_list.py $(PROCDIR)/pdbtm_redundant_tolerant.mems.pickle $(PROCDIR)/pdbtm_redundant.3line -b 1 -g 0 > $(STATDIR)/pdbtm_redundant_1_tolerant_all_list.txt -t True
# 	./bin/gen_potential_list.py $(PROCDIR)/pdbtm_redundant.mems.pickle $(PROCDIR)/pdbtm_redundant.3line -b 2 > $(STATDIR)/pdbtm_redundant_2_list.txt
# 	./bin/gen_potential_list.py $(PROCDIR)/pdbtm_redundant.mems.pickle $(PROCDIR)/pdbtm_redundant.3line -b 2 -g 0 > $(STATDIR)/pdbtm_redundant_2_all_list.txt

$(GLOBCHARGES) : $(MEMSGLOB)
	./bin/make_charges_from_mems.py $< $@

# $(MEMCHARGES) : %.clust.charges.pickle : %.clust.mems.pickle | $(MEMS)
# 	./bin/make_charges_from_mems.py $< $@

$(A3MCHARGES) : %.charges.pickle : %.mems.pickle | $(A3MMEMS)
	./bin/make_charges_from_mems.py $< $@

$(VISIMAGES) : $(A3MCHARGES) $(GLOBCHARGES)
	./bin/visualizeAAs.py $(PROCDIR)/pdbtm_a3m.charges.pickle $(IMAGEDIR)/pdbtm_vis.svg ""
	./bin/visualizeAAs.py $(PROCDIR)/pdbtm_a3m.charges.pickle $(IMAGEDIR)/pdbtm_vis.png ""
	./bin/visualizeAAs.py $(PROCDIR)/pdbtm_a3m.charges.pickle $(IMAGEDIR)/pdbtm_vis.pdf ""
	./bin/visualizeAAs.py $(PROCDIR)/scop_glob_a3m.charges.pickle $(IMAGEDIR)/scop_glob_vis.svg ""
	./bin/visualizeAAs.py $(PROCDIR)/scop_glob_a3m.charges.pickle $(IMAGEDIR)/scop_glob_vis.png ""
	./bin/visualizeAAs.py $(PROCDIR)/scop_glob_a3m.charges.pickle $(IMAGEDIR)/scop_glob_vis.pdf ""
	./bin/visualize_pairs.py $(PROCDIR)/pdbtm.clust.mems.pickle $(PROCDIR)/pdbtm.clust.pairs.csv
	# ./bin/visualizeAAs.py $(A3MCHARGES) $(IMAGEDIR)/pdbtm_vis.svg "Charges"
	./bin/visualizeAAs_mem.py $(PROCDIR)/pdbtm_a3m.charges.pickle $(IMAGEDIR)/mem_cluster ""
	# ./bin/visualizeAAs_Compare.py $(GLOBCHARGES) $(A3MCHARGES) $(IMAGEDIR)/charges_vis.svg "pdbtm vs Globular"
	# ./bin/visualizeAAs_Compare_all.py $(GLOBCHARGES) $(A3MCHARGES) $(IMAGEDIR)/mem_vs_glob.svg "pdbtm vs Globular"
	./bin/flow_image.py $(PROCDIR)/pdbtm.clust.mems.stats.pickle "PDBTM non-redundant" $(IMAGEDIR)/pdbtm_stats.png
	./bin/flow_image.py $(PROCDIR)/pdbtm.clust.mems.stats.pickle "PDBTM non-redundant" $(IMAGEDIR)/pdbtm_stats.svg 
	# ./bin/flow_image.py $(STATDIR)/pdbtm_redundant.mems.stats.pickle "PDBTM redundant" $(IMAGEDIR)/pdbtm_stats_red.png
	# ./bin/flow_image.py $(STATDIR)/pdbtm_redundant.mems.stats.pickle "PDBTM redundant" $(IMAGEDIR)/pdbtm_stats_red.svg 

.PHONY: clean deepclean
clean:
	rm -rf $(SCAMPIDIR)
	rm -rf $(PDBTMDIR)
	rm -rf $(PROCDIR)
	rm -r $(IMAGEDIR)/*
	rm -r $(STATS)
	rm -r $(LISTS)
deepclean: clean
	rm -rf $(OPMDIR)
	rm -r $(IMAGEDIR)/*
	rm -r $(STATS)
	rm -r $(LISTS)
	rm -rf $(RAWDIR)/pdb_chain_uniprot.tsv.gz
	rm -rf $(RAWDIR)/scop*
	rm -rf $(RAWDIR)/ss.txt.gz
	rm -rf $(RAWDIR)/pdbtm_alpha_entries.xml
	rm -rf $(RAWDIR)/pdbtm_non_redundant_alpha_list.txt
