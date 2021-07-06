## Code to prepare the visuals for our outlier blog post 

The blog post itself can be found [here](https://joachim-gassen.github.io/2021/07/outliers/) or [here](https://arc.eaa-online.org/blog/taking-outlier-treatment-next-level). The code in this folder has been used to prepare the visuals presented in the blog.


### Setup

0. This will need a current R/RStudio environment to run.
1. Copy the file _config.csv to config.csv in the `blog` directory. Edit it by adding your WRDS credentials. 
3. In `blog/code/prepare_sample.R`, you can decide which sample you want to base the visuals on. Currently there are two options: `World` (all non-financial firm observations with key variables including returns and fiscal years prior to 2020) and `Eurozone` (the `World` sample limited to firms headquartered in a Eurozone country and fiscal years 2005 to 2019). The `World` sample is denoted in U.S.-$ while the `Eurozone` sample is denoted in Euro. The code defaults to the `Eurozone` sample as this is what we use in the blog post.
2. Assuming that you have `make`installed, you can tun 'make all' either via the console or by identifying the 'Build All' button in the 'Build' tab (normally in the upper right quadrant of the RStudio screen). 
3. Eventually, you will be greeted with sample data in `data/generated` and 
output in `output`.


### Exploration

After running make, if you want to have a quick look at the data, you can 
run ExPanD:

```
# install.packages("ExPanDaR")
smp <- readRDS("blog/data/generated/sample.rds") %>% select(loc:ret)
ExPanDaR::ExPanD(smp, cs_id = c("gvkey", "conm"), ts_id = "fyear")
```


### Folder structure 

- `code`: This directory holds program scripts that are being called to download data from WRDS, prepare the data, run the analysis and create the output files

- `data`: A directory where data is stored. The data in the sub-folders `pulled` and `generated` will be populated by running the code.

- `output`: Will take output files but is empty. Why? Because you will create the output locally on your computer, if you want.


### Troubleshooting

If you do not see 'Build' tab this is most likely because you do not have 'make' installed on your system. 
  - For Windows: Install Rtools: https://cran.r-project.org/bin/windows/Rtools/
  - For MacOS: You need to install the Mac OS developer tools. Open a terminal and run `xcode-select --install` Follow the instructions
  - On Linux: I have never seen a Unix environment without 'make'. 

If you do not want to go through this process, you can instead source the code files by hand in the following order:

- `pull_wrds_data.R`
- `prepare_sample.R`
- `prepare_cfo_sample_visuals.R`
- `prepare_simulation_visuals.R`

If you have problems because of missing packages, try to run this code:

```
# Code to install packages to your system
install_package_if_missing <- function(pkg) {
  if (! pkg %in% installed.packages()[, "Package"]) install.packages(pkg)
}
install_package_if_missing("tidyverse")
install_package_if_missing("gganimate")
install_package_if_missing("animate")
install_package_if_missing("countrycode")
install_package_if_missing("gifski")
install_package_if_missing("ExPanDaR")
install_package_if_missing("RPostgres")
install_package_if_missing("DBI")
```

### Notes and references

This repository was built based on the ['treat' template for reproducible research](https://github.com/trr266/treat).


