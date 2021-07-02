/*
This code can be used to test whether the programs cluster2mm and bootmm function properly
The dataset is simulated to mimic a panel dataset with 100 firms and 50 years and a dependence 
structure in x and the error term across both dimensions similar to Gow, Ormazabal, and Taylor (2010)
*/

clear all
set seed 1234
local firms=100
local years=50
local obs=`firms'*`years'
set obs `obs'

gen n=_n
gen firm=ceil(n/`years')
gen year=n-(firm-1)*`years'

* Induce dependence structure:
local rho=0.8
local ccor=0.5
local sigma_x=sqrt(1)
local sigma_e=sqrt(2)
local sigma_xt=`ccor'*`sigma_x'
local sigma_xf=sqrt(`sigma_x'^2-`sigma_xt'^2)
local sigma_et=`ccor'*`sigma_e'
local sigma_ef=sqrt(`sigma_e'^2-`sigma_et'^2)
	
gen xt=rnormal()*`sigma_xt' if firm==1
egen max=max(xt), by(year)
replace xt=max
drop max
gen et=rnormal()*`sigma_et' if firm==1
egen max=max(et), by(year)
replace et=max
drop max

gen vx=sqrt(1-`rho'^2)*rnormal()*`sigma_xf'
gen xf=rnormal()*`sigma_xf' if year==1
replace xf=`rho'*xf[_n-1]+vx if year>=2
gen ve=sqrt(1-`rho'^2)*rnormal()*`sigma_ef'
gen ef=rnormal()*`sigma_ef' if year==1
replace ef=`rho'*ef[_n-1]+ve if year>=2

gen x=xt+xf
gen e=et+ef

* Induce vertical outliers:
gen double r=rnormal()
sort r
replace e=e*3 if _n<=(`obs'/10)
gen y=x+e	

sum y x e
reg y x

gen x2=rnormal()
gen x3=rnormal()
gen x4=rnormal()
gen x5=rnormal()
gen x6=rnormal()

robcluster2 y x x2 x3 x4 x5 x6, eff(95) cluster(firm year)
robcluster2 y x x2 x3 x4 x5 x6, eff(70) cluster(firm year) m
robcluster2 y x x2 x3 x4 x5 x6, eff(95) cluster(firm year) m biweight
robcluster2 y x x2 x3 x4 x5 x6, eff(70) cluster(firm year) s
robcluster2 y x x2 x3 x4 x5 x6, eff(95) cluster(firm year) median

roboot y x x2 x3 x4 x5 x6, eff(95) cluster(firm year) nboot(400) 
roboot y x x2 x3 x4 x5 x6, eff(70) cluster(firm) nboot(400) m
roboot y x x2 x3 x4 x5 x6, eff(70) nboot(400) m biweight
roboot y x x2 x3 x4 x5 x6, cluster(firm) nboot(400) s 
roboot y x x2 x3 x4 x5 x6, cluster(year) nboot(400) median

