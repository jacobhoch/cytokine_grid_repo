GRID1_STEM=241120_primeseq
GRID2_STEM=240722_primeseq

ALL_STEMS = $(GRID1_STEM) $(GRID2_STEM)

GRID1_INPUTS := $(foreach n, data/raw_counts/$(GRID1_STEM)/*.rds, $n)
GRID1_CSVS := $(foreach n, data/raw_counts/$(GRID1_STEM)/*.csv, $n)

GRID2_INPUTS := $(foreach n, data/raw_counts/$(GRID2_STEM)/*.rds, $n)
GRID2_CSVS := $(foreach n, data/raw_counts/$(GRID2_STEM)/*.csv, $n)


### COMBINE ZUMIS OUTPUTS
data/raw_dds/$(GRID1_STEM).Rds: src/loadData.R $(GRID1_INPUTS) $(GRID1_CSVS)
	Rscript $< -i $(GRID1_INPUTS) -m $(GRID1_CSVS) -s $(GRID1_STEM) -d Donor Stim 

data/raw_dds/$(GRID2_STEM).Rds: src/loadData.R data/raw_counts/$(GRID2_STEM)/*.rds
	Rscript $< -i $(GRID2_INPUTS) -m $(GRID2_CSVS) -s $(GRID2_STEM) -d Donor Stim

raw := $(foreach n, $(ALL_STEMS), $(addprefix data/raw_dds/, $(addprefix $n, .Rds)))
raw_dds: $(raw)

### RUN QC
data/clean_dds/%.Rds: src/runQC.R data/raw_dds/%.Rds
	Rscript src/runQC.R -i $(word 2, $^)

clean := $(foreach n, $(ALL_STEMS), $(addprefix data/clean_dds/, $(addprefix $n, .Rds)))
clean_dds: $(clean)

### DE
data/DE_results/%.Rds: src/runDE.R data/clean_dds/%.Rds data/comparisons/comparisons_%.txt
	Rscript src/runDE.R -i $(word 2, $^) -c $(word 3, $^) ### ADD DESIGN VARIABLE

DE := $(foreach n, $(ALL_STEMS), $(addprefix data/DE_results/, $(addprefix $n, .Rds)))
DE_dds: $(DE)
