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

grid1_de := $(foreach n, $(shell cat data/comparisons/comparisons_$(GRID1_STEM).txt), $(addprefix data/DE_results/$(GRID1_STEM)/, $(addsuffix _full.csv, $n)))
grid2_de := $(foreach n, $(shell cat data/comparisons/comparisons_$(GRID2_STEM).txt), $(addprefix data/DE_results/$(GRID2_STEM)/, $(addsuffix _full.csv, $n)))
tb_de := $(foreach n, $(shell cat data/comparisons/comparisons_$(TB_STEM).txt), $(addprefix data/DE_results/$(TB_STEM)/, $(addsuffix _full.csv, $n)))
DE_tables: $(grid1_de) $(grid2_de) $(tb_de)

DE := $(foreach n, $(ALL_STEMS), $(addprefix data/DE_results/, $(addprefix $n, _list.Rds)))
DE_dds: $(DE)

### FIGURES 1 & 2 CYTO GRID ###
fig/cytokine_grid/correlation_heatmap.pdf: src/grid_correlation.R data/DE_results/$(GRID1_STEM)_list.Rds
	Rscript $< -d $(word 2, $^)

fig/cytokine_grid/grid_PCA.pdf: src/grid_PCA.R data/DE_results/$(GRID1_STEM)_list.Rds
	Rscript $< -d $(word 2, $^)

fig/cytokine_grid/indep_DEG_grid.pdf: src/generateDEGgrid.R $(grid1_de)
	Rscript $< -s $(GRID1_STEM)

fig/cytokine_grid/upset_IFNb_dep.pdf: src/upset.R $(grid1_de)
	Rscript $< -s $(GRID1_STEM) -c IFNb

fig/cytokine_grid/upset_IL4_dep.pdf: src/upset.R $(grid1_de)
	Rscript $< -s $(GRID1_STEM) -c IL4

fig1: fig/cytokine_grid/correlation_heatmap.pdf fig/cytokine_grid/grid_PCA.pdf

fig2: fig/cytokine_grid/indep_DEG_grid.pdf fig/cytokine_grid/upset_IFNb_dep.pdf fig/cytokine_grid/upset_IL4_dep.pdf


fig: fig1 fig2


### MISC ####

# truth tables
fig/241120_primeseq_truthtable: src/DEGtruthtables.R data/DE_results/241120_primeseq/*_full.csv
	Rscript $< -s "241120_primeseq" -a "IFNb" -b "IL4"
