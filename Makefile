GRID1_STEM=241120_primeseq
GRID2_STEM=240722_primeseq

ALL_STEMS := $(GRID1_STEM) $(GRID2_STEM)

### COMBINE ZUMIS OUTPUTS
data/raw_dds/$(GRID1_STEM)_Stimonly.Rds: src/loadData.R data/raw_counts/$(GRID1_STEM)/%.Rds data/raw_counts/$(GRID1_STEM)/%.csv
	Rscript $< -i $(word 2, $^) -s $(word 3, $^) -d Donor Stim

data/raw_dds/$(GRID1_STEM)_day1day2.Rds: src/loadData.R data/raw_counts/$(GRID1_STEM)/%.Rds
	Rscript $< -i $(word 2, $^) -s $(word 3, $^) -d Donor day_1*day_2

data/raw_dds/$(GRID2_STEM).Rds: src/loadData.R data/raw_counts/$(GRID2_STEM)/%.Rds
	Rscript $< -i $(word 2, $^) -s $(word 3, $^) -d Donor Stim

raw := $(foreach n, $(ALL_STEMS), $(addprefix data/raw_dds/, $(addprefix $n, .Rds)))
raw_dds: $(raw)

### RUN QC
data/clean_dds/%.Rds: src/runQC.R data/raw_dds/%.Rds
	Rscript src/runQC.R -i $(word 2, $^)

clean := $(foreach n, $(ALL_STEMS), $(addprefix data/clean_dds/, $addprefix $n, .Rds)))
clean_dds: $(clean)

### DE
data/DE_results/%.Rds: src/runDE.R data/clean_dds/%.Rds data/comparisons/comparisons_%.txt
	Rscript src/runDE.R -i $(word 2, $^) -c $(word 3, $^)

DE := $(foreach n, $(ALL_STEMS), $(addprefix data/DE_results/, $(addprefix $n, .Rds)))
DE_dds: $(DE)
