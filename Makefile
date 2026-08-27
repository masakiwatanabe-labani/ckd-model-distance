.PHONY: all refs data analysis figures clean help

PY ?= python

REFS     = 00_fetch_refs 00b_fetch_external
CORE     = 01_load 02_de 03_crossspecies 04_modules 05_stage 06_compartment 07_candidates
CROSS    = 09_dataset_matrix 10_human_validation 11_ercb 12_model_distance
TIME     = 13_time_axis 14_time_axis_full 15_within_model_time 16_entry_model_checks
EXTRA    = 18_monotonic_genes 19_distance_distributions 20_null_diagnostics 21_table1
MODULE   = 24_module_level_distance 25_module_power 28_table2
FIGURES  = 17_fig3 22_fig1 23_fig2 26_fig6 27_fig5 08_figures

STEPS = $(CORE) $(CROSS) $(TIME) $(EXTRA) $(MODULE) $(FIGURES)

help:
	@echo "make refs      - download reference and external data (network)"
	@echo "make analysis  - run the analysis, assuming data are in place"
	@echo "make all       - refs then analysis"
	@echo "make clean     - remove results/ and logs/ contents"

refs:
	@for s in $(REFS); do echo "=== $$s ==="; $(PY) src/$$s.py || exit 1; done

analysis:
	@for s in $(STEPS); do echo "=== $$s ==="; $(PY) src/$$s.py || exit 1; done

figures:
	@for s in $(FIGURES); do echo "=== $$s ==="; $(PY) src/$$s.py || exit 1; done

all: refs analysis

clean:
	rm -rf results/* logs/*.log
	mkdir -p results/figures results/supplementary
	touch logs/.gitkeep
