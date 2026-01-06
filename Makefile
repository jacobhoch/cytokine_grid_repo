GRID1_STEM=241120_primeseq
GRID2_STEM=240722_primeseq
TB_STEM=250201_TBinfection

ALL_STEMS = $(GRID1_STEM) $(GRID2_STEM) $(TB_STEM)

GRID1_INPUTS := $(foreach n, data/raw_counts/$(GRID1_STEM)/*.rds, $n)
GRID1_CSVS := $(foreach n, data/raw_counts/$(GRID1_STEM)/*.csv, $n)

GRID2_INPUTS := $(foreach n, data/raw_counts/$(GRID2_STEM)/*.rds, $n)
GRID2_CSVS := $(foreach n, data/raw_counts/$(GRID2_STEM)/*.csv, $n)

TB_INPUTS := $(foreach n, data/raw_counts/$(TB_STEM)/*.rds, $n)
TB_CSVS := $(foreach n, data/raw_counts/$(TB_STEM)/*.csv, $n)

### COMBINE ZUMIS OUTPUTS ###
data/raw_dds/$(GRID1_STEM).Rds: src/loadData.R $(GRID1_INPUTS) $(GRID1_CSVS)
	Rscript $< -i $(GRID1_INPUTS) -m $(GRID1_CSVS) -s $(GRID1_STEM) -d Donor Stim 

data/raw_dds/$(GRID2_STEM).Rds: src/loadData.R $(GRID2_INPUTS) $(GRID2_CSVS)
	Rscript $< -i $(GRID2_INPUTS) -m $(GRID2_CSVS) -s $(GRID2_STEM) -d Donor Stim

data/raw_dds/$(TB_STEM).Rds: src/loadData.R $(TB_INPUTS) $(TB_CSVS)
	Rscript $< -i $(TB_INPUTS) -m $(TB_CSVS) -s $(TB_STEM) -d Donor Stim

raw := $(foreach n, $(ALL_STEMS), $(addprefix data/raw_dds/, $(addprefix $n, .Rds)))
raw_dds: $(raw)

### RUN QC AND COLLAPSE TECHNICAL REPLICATES ###
data/clean_dds/241120_primeseq.Rds: src/runQC.R data/raw_dds/241120_primeseq.Rds
	Rscript $< -i $(word 2, $^) -r data/241120_primeseq_remove.txt

data/clean_dds/240722_primeseq.Rds: src/runQC.R data/raw_dds/240722_primeseq.Rds
	Rscript $< -i $(word 2, $^) -r data/240722_primeseq_remove.txt #sample I messed up

data/clean_dds/250201_TBinfection.Rds: src/runQC.R data/raw_dds/250201_TBinfection.Rds
	Rscript $< -i $(word 2, $^) -r data/250201_TBinfection_remove.txt

clean := $(foreach n, $(ALL_STEMS), $(addprefix data/clean_dds/, $(addprefix $n, .Rds)))
clean_dds: $(clean)

### RUN DESEQ2 AND GENERATE DEG LISTS ###
data/DE_results/%_list.Rds: src/runDE.R data/clean_dds/%.Rds data/comparisons/comparisons_%.txt
	Rscript $< -i $(word 2, $^) -c $(word 3, $^)

DE := $(foreach n, $(ALL_STEMS), $(addprefix data/DE_results/, $(addprefix $n, _list.Rds)))
DE_dds: $(DE)
	date +"finish time = %F %T"

### MISC ####
# PCA
fig/%_PCA.png: src/PCA.R data/DE_results/%_list.Rds
	Rscript $< -i $(word 2, $^)

PCA := $(foreach n, $(ALL_STEMS), $(addprefix fig/, $(addprefix $n, _PCA.png)))
PCA_plot: $(PCA)

# DEG grid
fig/241120_primeseq_percentDEGgrid.png: src/generateDEGgrid.R data/DE_results/241120_primeseq/*_up.txt
	Rscript $< -s "241120_primeseq"

# truth tables
fig/241120_primeseq_truthtable: src/DEGtruthtables.R data/DE_results/241120_primeseq/*_full.csv
	Rscript $< -s "241120_primeseq" -a "IFNb" -b "IL4"
