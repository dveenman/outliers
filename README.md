# "Outliers and Robust Inference in Archival Accounting Research"

by Joachim Gassen and David Veenman

This repository accompanies our paper on outliers (https://ssrn.com/abstract=3880942) and provides the following:

- ado program files that can be installed and used in Stata, and  
- the code that creates the visuals of our [outlier blog post](https://joachim-gassen.github.io/2021/07/outliers/).


## Overview

- `robcluster2`: this program produces twoway cluster-robust standard errors for a robust regression estimator. It computes standard errors for robust estimators based on equation 3.11 of [Croux, Dhaene, and Hoorelbeke 2003](https://ideas.repec.org/p/ete/ceswps/ces0316.html). Twoway clustered standard errors are obtained based on (1) the intuition of [Zeileis, Koll, and Graham (2020)](https://www.jstatsoft.org/article/view/v095i01) for sandwich estimators to obtain cluster-robust standard errors and (2) the formula for twoway cluster-robust standard errors from [Thompson (2011)](https://doi.org/10.1016/j.jfineco.2010.08.016) and [Cameron, Gelbach, and Miller (2011)](https://doi.org/10.1198/jbes.2010.07136). By default, the robust estimator in `robcluster2` is the MM-estimator defined by `robreg mm`, but the program requires an explicit input for the desired estimation efficiency when regression errors are normally distributed. Other options are an M-estimator (with Huber or biweight objective function), an S-estimator, or median regression. Significance levels of reported *t*-values are based on the degrees of freedom defined by the cluster dimension with the lowest number of unique clusters (G-1).
- `roboot`: this program produces bootstrapped standard errors for MM-estimators using the fast bootstrap procedure described by [Salibian-Barrera and Zamar (2002)](https://www.jstor.org/stable/2699968#metadata_info_tab_contents). The robust estimator in `roboot` is the MM-estimator defined by `robreg mm`, but the program requires an explicit input for the desired estimation efficiency when regression errors are normally distributed. The program produces MM-estimates with bootstrapped standard errors without cluster-adjustment, one-way cluster-robust standard errors, or twoway cluster-robust standard errors. The bootstrapped standard errors are derived from the standard deviations of B bootstrap samples drawn from the sample with replacement. One-way cluster-robust bootstrapped standard errors are obtained similarly by drawing entire clusters of data from the sample with replacement. Twoway cluster-robust standard errors are obtained using the formula from [Thompson (2011)](https://doi.org/10.1016/j.jfineco.2010.08.016) and [Cameron, Gelbach, and Miller (2011)](https://doi.org/10.1198/jbes.2010.07136). The cluster-robust versions rely on the computational efficiency of drawing low-dimensional matrices per cluster as described by [MacKinnon (2022)](https://www.sciencedirect.com/science/article/abs/pii/S2452306221001404), which further speeds up the bootstrap procedure.
- `robfmb`: this program can be used to run Fama-Macbeth type robust regressions. Although the typical application is to derive estimates from the time-series of coefficient estimates derived from cross-sectional regressions, the program can be used to perform any type of sample splitting (e.g., [Conley, Gonçalves, and Hansen 2018](https://doi.org/10.1111/1475-679X.12219)). When the estimations are performed for each unit of time (e.g., by year), the program additionally allows for the calculation of a Newey-West correction of the standard errors.
- `robalance`: this program produces a table of output in which input variables are compared between the full sample and the reweighted sample of observations, based on the weights stored from a previous robust regression estimator. For continuous variables, the program compares the mean, 5th percentile value, first quartile (p25), median (p50), third quartile (p75), and 95th percentile value of the variable distribution between the full sample and the reweighted sample. Besides the difference in these distributional characteristics, the program provides z-statistics based on bootstrapped standard errors derived for the relevant characteristic. Bootstrapped standard errors can additionally be adjusted for clustering in one or two dimensions. For indicator variables, the program computes the difference in the unweighted versus the weighted mean of the variables. 
- `blog`: Folder that contains the code for preparing the blog post visuals. See its README for additional information.


## Installation and dependencies

The programs can be installed by simply saving the \*.ado files into your local Stata directory that contains additional ado programs. To identify this folder, type and execute "sysdir" in your Stata command window and go to the folder listed at "PLUS:". Make sure you place the program in the folder with the correct starting letter of the program (i.e., the folder names "r" for all programs) and the file extension is correctly specified as \*.ado.

The `robcluster2` and `roboot` programs require the latest versions of `moremata` and `robreg` to be installed in Stata:
```
ssc inst moremata, replace
ssc inst robreg, replace
```

## Program details: robcluster2

Syntax:
**robcluster2** *depvar* *indepvars*, *cluster(varlist)* [*options*]

This program can be used to estimate a robust regression with cluster-robust standard errors in the two dimensions specified in *cluster(varlist)*.

- The options available in [*options*] are:  
  - *eff(value)*: specifies the desired level of estimation efficiency of the estimator in case the regression errors follow a normal distribution. This option is required for the default MM-estimator and the M-estimator when using the *m* option.
  - *m* (optional): estimates an M-estimator (`robreg m`) instead of the default MM-estimator `robreg mm` (cannot be combined with *s* or *median* as options)  
  - *s* (optional): estimates an S-estimator (`robreg s`) instead of the default MM-estimator `robreg mm` (cannot be combined with *m* or *median* as options)  
  - *median* (optional): estimates median regression (`robreg q`) instead of the default MM-estimator `robreg mm` (cannot be combined with *m* or *s* as options)  
  - *biweight* (optional): when combined with option *m*, this option replaces the default Huber objective function in the M-estimator with the biweight objective function.
  - *tol(value)* (optional): specifies the tolerance for iterative algorithms (default is 1e-6).
  - *sopts(options)* (optional): specifies additional options to be passed to the initial S-estimator in MM-estimation.
  - *weightvar(name)* (optional): stores the robust regression weights in new variable ''name''. 

**Output**: 
The program produces standard regression output, with the number of observations used in the estimation, the degrees of freedom used to compute the critical values of the *t*-statistics (G-1, where G is the number of clusters in the dimension with lowest number of clusters), as well as the number of clusters in each dimension. Below the regression output, the program also lists information on the average weights assigned to observations in the final weighted least squares estimation, as well as the fraction of observations with weights below 0.5 and weights equal to zero. Illustration of output:

![image](https://user-images.githubusercontent.com/65561067/200134167-eba04f0c-7225-4158-a71b-b24be65ca360.png)

**Examples**: 
```
robcluster2 y x x2 x3 x4 x5 x6, eff(95) cluster(firm year)
robcluster2 y x x2 x3 x4 x5 x6, eff(70) cluster(firm year) m
robcluster2 y x x2 x3 x4 x5 x6, eff(95) cluster(firm year) m biweight
robcluster2 y x x2 x3 x4 x5 x6, cluster(firm year) s
robcluster2 y x x2 x3 x4 x5 x6, cluster(firm year) median
```

## Program details: roboot
Syntax:
**roboot** *depvar* *indepvars*, *nboot(value)* [*cluster(varlist)*] [*options*]

This program uses fast bootstrapping procedures to produce bootstrapped standard errors for MM-estimators based on B bootstrap samples defined in *nboot(value)*, with a potential adjustment for clustering up to two dimensions specified in option *cluster(varlist)*. 

- The options available in [*options*] are: 
  - *nboot(value)*: the number of bootstrap samples to be drawn.
  - *eff(value)*: specifies the desired level of estimation efficiency of the estimator in case the regression errors follow a normal distribution.
  - *tol(value)* (optional): specifies the tolerance for iterative algorithms (default is 1e-6).
  - *sopts(options)* (optional): specifies additional options to be passed to the initial S-estimator in MM-estimation.
  - *weightvar(name)* (optional): stores the robust regression weights in new variable ''name''.
  - *seed(value)* (optional): sets the seed to ensure reproducibility of the results. Not required but highly recommended.

**Output**:
The program produces standard regression output, with the number of observations used in the estimation and the number of clusters in each dimension. The bootstrapped standard errors are adjusted for a finite-sample correction similar to **robcluster2** and p-values rely on critical values from a t-distribution with G-1 degrees of freedom following the recommendations of [Cameron and Miller (2015)](http://cameron.econ.ucdavis.edu/research/Cameron_Miller_JHR_2015_February.pdf). Below the regression output, the program also lists information on the average weights assigned to observations in the final weighted least squares estimation, as well as the fraction of observations with weights below 0.5 and weights equal to zero. Illustration of output:

![image](https://user-images.githubusercontent.com/65561067/200160575-697e458b-ddc4-4d89-9b30-e10b6999704a.png)

**Examples**:
```
roboot y x x2 x3 x4 x5 x6, eff(95) cluster(firm year) nboot(1000) 
roboot y x x2 x3 x4 x5 x6, eff(70) cluster(firm) nboot(1000)
roboot y x x2 x3 x4 x5 x6, eff(50) nboot(1000) seed(1234)
``` 

## Program details: robfmb
Syntax:
**robfmb** *depvar* *indepvars*, *splitvar(varname)* [*options*]

- *varname* in *splitvar(varname)* contains the sample-splitting variable based on which the separate estimations should be performed.

- The options available in [*options*] are:  
  - *lag(value)*: specifies the number of lags for Newey-West correction. This option should be used only when the split variable is a time variable.
  - *eff(value)*: specifies the desired level of estimation efficiency of the estimator in case the regression errors follow a normal distribution. This option is required for the default MM-estimator and the M-estimator when using the *m* option.
  - *m*: estimates an M-estimator (`robreg m`) instead of the default MM-estimator `robreg mm` (cannot be combined with *s* or *median* as options)  
  - *s*: estimates an S-estimator (`robreg s`) instead of the default MM-estimator `robreg mm` (cannot be combined with *m* or *median* as options)  
  - *median*: estimates median regression (`robreg q`) instead of the default MM-estimator `robreg mm` (cannot be combined with *m* or *s* as options)  
  - *biweight*: when combined with option *m*, this option replaces the default Huber objective function in the M-estimator with the biweight objective function.
  - *tol(value)*: specifies the tolerance for iterative algorithms.

**Example and output**:

![image](https://user-images.githubusercontent.com/65561067/124244725-8d898080-db1f-11eb-852d-ae9d40b93d0c.png)

## Program details: robalance 
Syntax:
**robalance** *varlist*, *weightvar(varname)* *nboot(value)*  [*cluster(varlist)* *dummies(varlist)* *box(string)*]

For all variables included in *varlist*, this program provides a table comparing the distributional properties (mean and percentiles 5, 25, 50, 75, and 95) between the base sample including outliers and the reweighted sample after a robust regression. The program uses as input the variable specified in *weightvar(varname)*, which is obtained from a previous robust regression estimation. For indicator variables, which should be specified in option *dummies(varlist)*, only a mean comparison is provided. Standard errors to derive z-statistics for differences in the statistic of interest between the base sample and the reweighted sample are obtained using a bootstrapping procedure, where the number of bootstrap samples B is defined in *nboot(value)*. Bootstrapped standard errors can additionally be adjusted for one- or twoway clustering based on the (maximum two) variables specified in option *cluster(varlist)*. To create and export PDF box plots (using Stata's standard box plots, excluding outside values) the option *box(string)* can be added, where *string* refers to the path of the folder to which the PDFs should be exported.

**Example and output**:

![image](https://user-images.githubusercontent.com/65561067/124132541-36cd6980-da81-11eb-9f8f-947023aba06c.png)

