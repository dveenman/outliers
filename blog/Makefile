# If you are new to Makefiles: https://makefiletutorial.com

SIMULATION_VIZ := output/simulation_visuals.rda \
	output/fig_sim_base.png \
	output/ani_sim_bad_leverage.gif \
	output/ani_sim_vert_outlier.gif \
	output/fig_sim_irwls_base.png \
	output/ani_sim_irwls.gif \
	output/fig_sim_nl_100runs.png \
	output/fig_sim_nl_base.png \
	output/ani_sim_nl_irwls.gif

CFO_SAMPLE_VIZ := output/cfo_sample_visuals.rda \
	output/fig_cfo_scatter.png \
	output/ani_cfo_linear_irwls.gif \
	output/tab_cfo_rrweights_by_ind.html \
	output/fig_cfo_rrweights_by_size.png \
	output/fig_cfo_rrweights_by_sales_growth.png \
	output/fig_cfo_linear.png \
	output/fig_cfo_iacted.png

TARGETS :=  $(SIMULATION_VIZ) $(CFO_SAMPLE_VIZ)

EXTERNAL_DATA := data/external/fama_french_12_industries.csv \
	data/external/fama_french_48_industries.csv

WRDS_DATA := data/pulled/cstat_global_fund.rds \
	data/pulled/cstat_global_sec.rds

GENERATED_DATA := data/generated/sample.rds

RSCRIPT := Rscript --encoding=UTF-8

.phony: all clean very-clean dist-clean

all: $(TARGETS)

clean:
	rm -f $(TARGETS)
	rm -f $(GENERATED_DATA)

very-clean: clean
	rm -f $(WRDS_DATA)

dist-clean: very-clean
	rm config.csv

config.csv:
	@echo "To start, you need to copy _config.csv to config.csv and edit it"
	@false

$(WRDS_DATA)&: code/pull_wrds_data.R code/read_config.R config.csv
	$(RSCRIPT) code/pull_wrds_data.R

$(GENERATED_DATA): $(WRDS_DATA) $(EXTERNAL_DATA) code/prepare_sample.R
	$(RSCRIPT) code/prepare_sample.R

$(SIMULATION_VIZ)&:	code/animate_irwls.R \
		code/prepare_simulation_visuals.R
	$(RSCRIPT) code/prepare_simulation_visuals.R

$(CFO_SAMPLE_VIZ)&:	$(GENERATED_DATA) code/animate_irwls.R \
		code/prepare_cfo_sample_visuals.R
	$(RSCRIPT) code/prepare_cfo_sample_visuals.R

