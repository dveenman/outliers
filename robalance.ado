*! version 1.0 01july2021 David Veenman

program define robalance, sortpreserve
	syntax varlist [in] [if], [cluster(varlist)] nboot(integer) weightvar(varname) [dummies(varlist)] [box(string)] 
	tokenize `varlist'
	marksample touse

	local nc: word count `cluster'
	if (`nc'>2){
	    di as text "ERROR: Maximum number of dimensions to cluster on is two"
		exit
	}
	if (`nc'>0){
		local fcluster: word 1 of `cluster'
	}
	if (`nc'>1){
		local tcluster: word 2 of `cluster'
	}
	
	foreach var of local dummies{
		qui sum `var'
		if (r(max)!=1 | r(min)!=0){
			di as text "ERROR: variable ---" "`var'" "--- is not a dummy variable"
			exit
		}
		else{
			qui sum `var' if `var'>0 & `var'<1
			if r(N)!=0{
				di as text "ERROR: variable ---" "`var'" "--- is not a dummy variable"
				exit
			}
		}
	
		if strpos("`varlist'","`var'"){
		}
		else{
			di as text "ERROR: dummy variable ---" "`var'" "--- must be included in variable list"
			exit
		}
	}
	
	di ""
	di as text "Continuous variables:", _continue
	foreach var of local varlist{
		if strpos("`dummies'","`var'"){
		}
		else{
			di "`var'", _continue
		}
	}
	di ""
	di as text "Dummy variables: " "`dummies'"

	tempvar bw n n2 intersection
	qui gen `bw'=.
	qui gsort -`touse'
	qui gen `n'=_n if `touse'
	qui gen `n2'=.
	qui egen `intersection'=group(`fcluster' `tcluster') if `touse'
	
	/* Estimation without clustering */
	if (`nc'==0){	
	local j=1
	foreach var of local varlist{
		tempvar p5_`j'
		qui gen `p5_`j''=.
		tempvar p5w_`j'
		qui gen `p5w_`j''=.
		tempvar p5d_`j'
		qui gen `p5d_`j''=.
		tempvar p25_`j'
		qui gen `p25_`j''=.
		tempvar p25w_`j'
		qui gen `p25w_`j''=.
		tempvar p25d_`j'
		qui gen `p25d_`j''=.
		tempvar p50_`j'
		qui gen `p50_`j''=.
		tempvar p50w_`j'
		qui gen `p50w_`j''=.
		tempvar p50d_`j'
		qui gen `p50d_`j''=.
		tempvar p75_`j'
		qui gen `p75_`j''=.
		tempvar p75w_`j'
		qui gen `p75w_`j''=.
		tempvar p75d_`j'
		qui gen `p75d_`j''=.
		tempvar p95_`j'
		qui gen `p95_`j''=.
		tempvar p95w_`j'
		qui gen `p95w_`j''=.
		tempvar p95d_`j'
		qui gen `p95d_`j''=.
		tempvar mean_`j'
		qui gen `mean_`j''=.
		tempvar meanw_`j'
		qui gen `meanw_`j''=.
		tempvar meand_`j'
		qui gen `meand_`j''=.
		tempvar full_p5_`j'
		qui gen `full_p5_`j''=.
		tempvar full_p5w_`j'
		qui gen `full_p5w_`j''=.
		tempvar full_p5d_`j'
		qui gen `full_p5d_`j''=.
		tempvar full_p25_`j'
		qui gen `full_p25_`j''=.
		tempvar full_p25w_`j'
		qui gen `full_p25w_`j''=.
		tempvar full_p25d_`j'
		qui gen `full_p25d_`j''=.
		tempvar full_p50_`j'
		qui gen `full_p50_`j''=.
		tempvar full_p50w_`j'
		qui gen `full_p50w_`j''=.
		tempvar full_p50d_`j'
		qui gen `full_p50d_`j''=.
		tempvar full_p75_`j'
		qui gen `full_p75_`j''=.
		tempvar full_p75w_`j'
		qui gen `full_p75w_`j''=.
		tempvar full_p75d_`j'
		qui gen `full_p75d_`j''=.
		tempvar full_p95_`j'
		qui gen `full_p95_`j''=.
		tempvar full_p95w_`j'
		qui gen `full_p95w_`j''=.
		tempvar full_p95d_`j'
		qui gen `full_p95d_`j''=.
		tempvar full_mean_`j'
		qui gen `full_mean_`j''=.
		tempvar full_meanw_`j'
		qui gen `full_meanw_`j''=.
		tempvar full_meand_`j'
		qui gen `full_meand_`j''=.
		local j=`j'+1
	}

	forvalues i=1(1)`nboot'{
		bsample if `touse', weight(`bw')
			qui expand `bw'
			qui sort `n'
			qui by `n': replace `n2'=_n
			local j=1
			foreach var of local varlist{
				qui sum `var' if `bw'>0 & `touse', d
				qui replace `p5_`j''=r(p5) if `n'==`i'
				qui replace `p25_`j''=r(p25) if `n'==`i'
				qui replace `p50_`j''=r(p50) if `n'==`i'
				qui replace `p75_`j''=r(p75) if `n'==`i'
				qui replace `p95_`j''=r(p95) if `n'==`i'
				qui replace `mean_`j''=r(mean) if `n'==`i'
				qui sum `var' [aw=`weightvar'] if `bw'>0 & `touse', d
				qui replace `p5w_`j''=r(p5) if `n'==`i'
				qui replace `p25w_`j''=r(p25) if `n'==`i'
				qui replace `p50w_`j''=r(p50) if `n'==`i'
				qui replace `p75w_`j''=r(p75) if `n'==`i'
				qui replace `p95w_`j''=r(p95) if `n'==`i'
				qui replace `meanw_`j''=r(mean) if `n'==`i'
				qui replace `p5d_`j''=`p5w_`j''-`p5_`j'' if `n'==`i'
				qui replace `p25d_`j''=`p25w_`j''-`p25_`j'' if `n'==`i'
				qui replace `p50d_`j''=`p50w_`j''-`p50_`j'' if `n'==`i'
				qui replace `p75d_`j''=`p75w_`j''-`p75_`j'' if `n'==`i'
				qui replace `p95d_`j''=`p95w_`j''-`p95_`j'' if `n'==`i'
				qui replace `meand_`j''=`meanw_`j''-`mean_`j'' if `n'==`i'
				local j=`j'+1
			}
			qui keep if `n2'==1
			qui replace `n2'=.		
		di ".", _continue
	}
	}

	/* Estimation with clustering on one dimension*/
	if (`nc'==1){	
	local j=1
	foreach var of local varlist{
		tempvar p5_`j'
		qui gen `p5_`j''=.
		tempvar p5w_`j'
		qui gen `p5w_`j''=.
		tempvar p5d_`j'
		qui gen `p5d_`j''=.
		tempvar p25_`j'
		qui gen `p25_`j''=.
		tempvar p25w_`j'
		qui gen `p25w_`j''=.
		tempvar p25d_`j'
		qui gen `p25d_`j''=.
		tempvar p50_`j'
		qui gen `p50_`j''=.
		tempvar p50w_`j'
		qui gen `p50w_`j''=.
		tempvar p50d_`j'
		qui gen `p50d_`j''=.
		tempvar p75_`j'
		qui gen `p75_`j''=.
		tempvar p75w_`j'
		qui gen `p75w_`j''=.
		tempvar p75d_`j'
		qui gen `p75d_`j''=.
		tempvar p95_`j'
		qui gen `p95_`j''=.
		tempvar p95w_`j'
		qui gen `p95w_`j''=.
		tempvar p95d_`j'
		qui gen `p95d_`j''=.
		tempvar mean_`j'
		qui gen `mean_`j''=.
		tempvar meanw_`j'
		qui gen `meanw_`j''=.
		tempvar meand_`j'
		qui gen `meand_`j''=.
		tempvar full_p5_`j'
		qui gen `full_p5_`j''=.
		tempvar full_p5w_`j'
		qui gen `full_p5w_`j''=.
		tempvar full_p5d_`j'
		qui gen `full_p5d_`j''=.
		tempvar full_p25_`j'
		qui gen `full_p25_`j''=.
		tempvar full_p25w_`j'
		qui gen `full_p25w_`j''=.
		tempvar full_p25d_`j'
		qui gen `full_p25d_`j''=.
		tempvar full_p50_`j'
		qui gen `full_p50_`j''=.
		tempvar full_p50w_`j'
		qui gen `full_p50w_`j''=.
		tempvar full_p50d_`j'
		qui gen `full_p50d_`j''=.
		tempvar full_p75_`j'
		qui gen `full_p75_`j''=.
		tempvar full_p75w_`j'
		qui gen `full_p75w_`j''=.
		tempvar full_p75d_`j'
		qui gen `full_p75d_`j''=.
		tempvar full_p95_`j'
		qui gen `full_p95_`j''=.
		tempvar full_p95w_`j'
		qui gen `full_p95w_`j''=.
		tempvar full_p95d_`j'
		qui gen `full_p95d_`j''=.
		tempvar full_mean_`j'
		qui gen `full_mean_`j''=.
		tempvar full_meanw_`j'
		qui gen `full_meanw_`j''=.
		tempvar full_meand_`j'
		qui gen `full_meand_`j''=.
		local j=`j'+1
	}

	forvalues i=1(1)`nboot'{
		bsample if `touse', cluster(`fcluster') weight(`bw')
			qui expand `bw'
			qui sort `n'
			qui by `n': replace `n2'=_n
			local j=1
			foreach var of local varlist{
				qui sum `var' if `bw'>0 & `touse', d
				qui replace `p5_`j''=r(p5) if `n'==`i'
				qui replace `p25_`j''=r(p25) if `n'==`i'
				qui replace `p50_`j''=r(p50) if `n'==`i'
				qui replace `p75_`j''=r(p75) if `n'==`i'
				qui replace `p95_`j''=r(p95) if `n'==`i'
				qui replace `mean_`j''=r(mean) if `n'==`i'
				qui sum `var' [aw=`weightvar'] if `bw'>0 & `touse', d
				qui replace `p5w_`j''=r(p5) if `n'==`i'
				qui replace `p25w_`j''=r(p25) if `n'==`i'
				qui replace `p50w_`j''=r(p50) if `n'==`i'
				qui replace `p75w_`j''=r(p75) if `n'==`i'
				qui replace `p95w_`j''=r(p95) if `n'==`i'
				qui replace `meanw_`j''=r(mean) if `n'==`i'
				qui replace `p5d_`j''=`p5w_`j''-`p5_`j'' if `n'==`i'
				qui replace `p25d_`j''=`p25w_`j''-`p25_`j'' if `n'==`i'
				qui replace `p50d_`j''=`p50w_`j''-`p50_`j'' if `n'==`i'
				qui replace `p75d_`j''=`p75w_`j''-`p75_`j'' if `n'==`i'
				qui replace `p95d_`j''=`p95w_`j''-`p95_`j'' if `n'==`i'
				qui replace `meand_`j''=`meanw_`j''-`mean_`j'' if `n'==`i'
				local j=`j'+1
			}
			qui keep if `n2'==1
			qui replace `n2'=.		
		di ".", _continue
	}
	}

	/* Estimation with clustering on two dimensions*/
	if (`nc'==2){	
	local j=1
	foreach var of local varlist{
		forvalues u=1(1)3{
			tempvar p5_`j'_`u'
			qui gen `p5_`j'_`u''=.
			tempvar p5w_`j'_`u'
			qui gen `p5w_`j'_`u''=.
			tempvar p5d_`j'_`u'
			qui gen `p5d_`j'_`u''=.
			tempvar p25_`j'_`u'
			qui gen `p25_`j'_`u''=.
			tempvar p25w_`j'_`u'
			qui gen `p25w_`j'_`u''=.
			tempvar p25d_`j'_`u'
			qui gen `p25d_`j'_`u''=.
			tempvar p50_`j'_`u'
			qui gen `p50_`j'_`u''=.
			tempvar p50w_`j'_`u'
			qui gen `p50w_`j'_`u''=.
			tempvar p50d_`j'_`u'
			qui gen `p50d_`j'_`u''=.
			tempvar p75_`j'_`u'
			qui gen `p75_`j'_`u''=.
			tempvar p75w_`j'_`u'
			qui gen `p75w_`j'_`u''=.
			tempvar p75d_`j'_`u'
			qui gen `p75d_`j'_`u''=.
			tempvar p95_`j'_`u'
			qui gen `p95_`j'_`u''=.
			tempvar p95w_`j'_`u'
			qui gen `p95w_`j'_`u''=.
			tempvar p95d_`j'_`u'
			qui gen `p95d_`j'_`u''=.
			tempvar mean_`j'_`u'
			qui gen `mean_`j'_`u''=.
			tempvar meanw_`j'_`u'
			qui gen `meanw_`j'_`u''=.
			tempvar meand_`j'_`u'
			qui gen `meand_`j'_`u''=.
		}
		tempvar full_p5_`j'
		qui gen `full_p5_`j''=.
		tempvar full_p5w_`j'
		qui gen `full_p5w_`j''=.
		tempvar full_p5d_`j'
		qui gen `full_p5d_`j''=.
		tempvar full_p25_`j'
		qui gen `full_p25_`j''=.
		tempvar full_p25w_`j'
		qui gen `full_p25w_`j''=.
		tempvar full_p25d_`j'
		qui gen `full_p25d_`j''=.
		tempvar full_p50_`j'
		qui gen `full_p50_`j''=.
		tempvar full_p50w_`j'
		qui gen `full_p50w_`j''=.
		tempvar full_p50d_`j'
		qui gen `full_p50d_`j''=.
		tempvar full_p75_`j'
		qui gen `full_p75_`j''=.
		tempvar full_p75w_`j'
		qui gen `full_p75w_`j''=.
		tempvar full_p75d_`j'
		qui gen `full_p75d_`j''=.
		tempvar full_p95_`j'
		qui gen `full_p95_`j''=.
		tempvar full_p95w_`j'
		qui gen `full_p95w_`j''=.
		tempvar full_p95d_`j'
		qui gen `full_p95d_`j''=.
		tempvar full_mean_`j'
		qui gen `full_mean_`j''=.
		tempvar full_meanw_`j'
		qui gen `full_meanw_`j''=.
		tempvar full_meand_`j'
		qui gen `full_meand_`j''=.
		local j=`j'+1
	}

	forvalues i=1(1)`nboot'{
		bsample if `touse', cluster(`fcluster') weight(`bw')
			qui expand `bw'
			qui sort `n'
			qui by `n': replace `n2'=_n
			local j=1
			foreach var of local varlist{
				qui sum `var' if `bw'>0 & `touse', d
				qui replace `p5_`j'_1'=r(p5) if `n'==`i'
				qui replace `p25_`j'_1'=r(p25) if `n'==`i'
				qui replace `p50_`j'_1'=r(p50) if `n'==`i'
				qui replace `p75_`j'_1'=r(p75) if `n'==`i'
				qui replace `p95_`j'_1'=r(p95) if `n'==`i'
				qui replace `mean_`j'_1'=r(mean) if `n'==`i'
				qui sum `var' [aw=`weightvar'] if `bw'>0 & `touse', d
				qui replace `p5w_`j'_1'=r(p5) if `n'==`i'
				qui replace `p25w_`j'_1'=r(p25) if `n'==`i'
				qui replace `p50w_`j'_1'=r(p50) if `n'==`i'
				qui replace `p75w_`j'_1'=r(p75) if `n'==`i'
				qui replace `p95w_`j'_1'=r(p95) if `n'==`i'
				qui replace `meanw_`j'_1'=r(mean) if `n'==`i'
				qui replace `p5d_`j'_1'=`p5w_`j'_1'-`p5_`j'_1' if `n'==`i'
				qui replace `p25d_`j'_1'=`p25w_`j'_1'-`p25_`j'_1' if `n'==`i'
				qui replace `p50d_`j'_1'=`p50w_`j'_1'-`p50_`j'_1' if `n'==`i'
				qui replace `p75d_`j'_1'=`p75w_`j'_1'-`p75_`j'_1' if `n'==`i'
				qui replace `p95d_`j'_1'=`p95w_`j'_1'-`p95_`j'_1' if `n'==`i'
				qui replace `meand_`j'_1'=`meanw_`j'_1'-`mean_`j'_1' if `n'==`i'
				local j=`j'+1
			}
			qui keep if `n2'==1
			qui replace `n2'=.		
		bsample if `touse', cluster(`tcluster') weight(`bw')
			qui expand `bw'
			qui sort `n'
			qui by `n': replace `n2'=_n
			local j=1
			foreach var of local varlist{
				qui sum `var' if `bw'>0 & `touse', d
				qui replace `p5_`j'_2'=r(p5) if `n'==`i'
				qui replace `p25_`j'_2'=r(p25) if `n'==`i'
				qui replace `p50_`j'_2'=r(p50) if `n'==`i'
				qui replace `p75_`j'_2'=r(p75) if `n'==`i'
				qui replace `p95_`j'_2'=r(p95) if `n'==`i'
				qui replace `mean_`j'_2'=r(mean) if `n'==`i'
				qui sum `var' [aw=`weightvar'] if `bw'>0 & `touse', d
				qui replace `p5w_`j'_2'=r(p5) if `n'==`i'
				qui replace `p25w_`j'_2'=r(p25) if `n'==`i'
				qui replace `p50w_`j'_2'=r(p50) if `n'==`i'
				qui replace `p75w_`j'_2'=r(p75) if `n'==`i'
				qui replace `p95w_`j'_2'=r(p95) if `n'==`i'
				qui replace `meanw_`j'_2'=r(mean) if `n'==`i'
				qui replace `p5d_`j'_2'=`p5w_`j'_2'-`p5_`j'_2' if `n'==`i'
				qui replace `p25d_`j'_2'=`p25w_`j'_2'-`p25_`j'_2' if `n'==`i'
				qui replace `p50d_`j'_2'=`p50w_`j'_2'-`p50_`j'_2' if `n'==`i'
				qui replace `p75d_`j'_2'=`p75w_`j'_2'-`p75_`j'_2' if `n'==`i'
				qui replace `p95d_`j'_2'=`p95w_`j'_2'-`p95_`j'_2' if `n'==`i'
				qui replace `meand_`j'_2'=`meanw_`j'_2'-`mean_`j'_2' if `n'==`i'
				local j=`j'+1
			}
			qui keep if `n2'==1
			qui replace `n2'=.		
		bsample if `touse', cluster(`intersection') weight(`bw')
			qui expand `bw'
			qui sort `n'
			qui by `n': replace `n2'=_n
			local j=1
			foreach var of local varlist{
				qui sum `var' if `bw'>0 & `touse', d
				qui replace `p5_`j'_3'=r(p5) if `n'==`i'
				qui replace `p25_`j'_3'=r(p25) if `n'==`i'
				qui replace `p50_`j'_3'=r(p50) if `n'==`i'
				qui replace `p75_`j'_3'=r(p75) if `n'==`i'
				qui replace `p95_`j'_3'=r(p95) if `n'==`i'
				qui replace `mean_`j'_3'=r(mean) if `n'==`i'
				qui sum `var' [aw=`weightvar'] if `bw'>0 & `touse', d
				qui replace `p5w_`j'_3'=r(p5) if `n'==`i'
				qui replace `p25w_`j'_3'=r(p25) if `n'==`i'
				qui replace `p50w_`j'_3'=r(p50) if `n'==`i'
				qui replace `p75w_`j'_3'=r(p75) if `n'==`i'
				qui replace `p95w_`j'_3'=r(p95) if `n'==`i'
				qui replace `meanw_`j'_3'=r(mean) if `n'==`i'
				qui replace `p5d_`j'_3'=`p5w_`j'_3'-`p5_`j'_3' if `n'==`i'
				qui replace `p25d_`j'_3'=`p25w_`j'_3'-`p25_`j'_3' if `n'==`i'
				qui replace `p50d_`j'_3'=`p50w_`j'_3'-`p50_`j'_3' if `n'==`i'
				qui replace `p75d_`j'_3'=`p75w_`j'_3'-`p75_`j'_3' if `n'==`i'
				qui replace `p95d_`j'_3'=`p95w_`j'_3'-`p95_`j'_3' if `n'==`i'
				qui replace `meand_`j'_3'=`meanw_`j'_3'-`mean_`j'_3' if `n'==`i'
				local j=`j'+1
			}
			qui keep if `n2'==1
			qui replace `n2'=.		
		di ".", _continue
	}
	}
	
	local j=1
	foreach var of local varlist{
		qui sum `var' if `touse', d
		qui replace `full_p5_`j''=r(p5) if `n'==1
		qui replace `full_p25_`j''=r(p25) if `n'==1
		qui replace `full_p50_`j''=r(p50) if `n'==1
		qui replace `full_p75_`j''=r(p75) if `n'==1
		qui replace `full_p95_`j''=r(p95) if `n'==1
		qui replace `full_mean_`j''=r(mean) if `n'==1
		qui sum `var' [aw=`weightvar'] if `touse', d
		qui replace `full_p5w_`j''=r(p5) if `n'==1
		qui replace `full_p25w_`j''=r(p25) if `n'==1
		qui replace `full_p50w_`j''=r(p50) if `n'==1
		qui replace `full_p75w_`j''=r(p75) if `n'==1
		qui replace `full_p95w_`j''=r(p95) if `n'==1
		qui replace `full_meanw_`j''=r(mean) if `n'==1
		qui replace `full_p5d_`j''=`full_p5w_`j''-`full_p5_`j'' if `n'==1
		qui replace `full_p25d_`j''=`full_p25w_`j''-`full_p25_`j'' if `n'==1
		qui replace `full_p50d_`j''=`full_p50w_`j''-`full_p50_`j'' if `n'==1
		qui replace `full_p75d_`j''=`full_p75w_`j''-`full_p75_`j'' if `n'==1		
		qui replace `full_p95d_`j''=`full_p95w_`j''-`full_p95_`j'' if `n'==1		
		qui replace `full_meand_`j''=`full_meanw_`j''-`full_mean_`j'' if `n'==1				
		local j=`j'+1
	}
	
	di ""
	di ""
	di "Balancing table of differences in distributional characteristics after robust re-weighting of sample"
	di "------------------------------------------------------------------------------------------------------"
	di as text "Statistic -->" ///
		_column(25) "mean" ///
		_column(37) "p5" ///
		_column(49) "p25" ///
		_column(61) "p50" ///
		_column(73) "p75" ///
		_column(85) "p95" 	    
	di "------------------------------------------------------------------------------------------------------"
	
	local j=1
	foreach var of local varlist{
		qui sum `full_p5d_`j'' if `n'==1
		local m_p5=r(mean)
		if (`nc'==2){
			qui sum `p5d_`j'_1'
			local v1=r(Var)
			qui sum `p5d_`j'_2'
			local v2=r(Var)
			qui sum `p5d_`j'_3'
			local v3=r(Var)
			local sd_p5=sqrt(`v1'+`v2'-`v3')
		}
		else{
			qui sum `p5d_`j''
			local sd_p5=r(sd)
		}
		if `sd_p5'==0{
			local z_p5=0
		}
		else{
			local z_p5=`m_p5'/`sd_p5'
		}
		qui sum `full_p25d_`j'' if `n'==1
		local m_p25=r(mean)
		if (`nc'==2){
			qui sum `p25d_`j'_1'
			local v1=r(Var)
			qui sum `p25d_`j'_2'
			local v2=r(Var)
			qui sum `p25d_`j'_3'
			local v3=r(Var)
			local sd_p25=sqrt(`v1'+`v2'-`v3')
		}
		else{
			qui sum `p25d_`j''
			local sd_p25=r(sd)
		}
		if `sd_p25'==0{
			local z_p25=0
		}
		else{
			local z_p25=`m_p25'/`sd_p25'
		}
		qui sum `full_p50d_`j'' if `n'==1
		local m_p50=r(mean)
		if (`nc'==2){
			qui sum `p50d_`j'_1'
			local v1=r(Var)
			qui sum `p50d_`j'_2'
			local v2=r(Var)
			qui sum `p50d_`j'_3'
			local v3=r(Var)
			local sd_p50=sqrt(`v1'+`v2'-`v3')
		}
		else{
			qui sum `p50d_`j''
			local sd_p50=r(sd)
		}
		if `sd_p50'==0{
			local z_p50=0
		}
		else{
			local z_p50=`m_p50'/`sd_p50'
		}
		qui sum `full_p75d_`j'' if `n'==1
		local m_p75=r(mean)
		if (`nc'==2){
			qui sum `p75d_`j'_1'
			local v1=r(Var)
			qui sum `p75d_`j'_2'
			local v2=r(Var)
			qui sum `p75d_`j'_3'
			local v3=r(Var)
			local sd_p75=sqrt(`v1'+`v2'-`v3')
		}
		else{
			qui sum `p75d_`j''
			local sd_p75=r(sd)
		}
		if `sd_p75'==0{
			local z_p75=0
		}
		else{
			local z_p75=`m_p75'/`sd_p75'
		}
		qui sum `full_p95d_`j'' if `n'==1
		local m_p95=r(mean)
		if (`nc'==2){
			qui sum `p95d_`j'_1'
			local v1=r(Var)
			qui sum `p95d_`j'_2'
			local v2=r(Var)
			qui sum `p95d_`j'_3'
			local v3=r(Var)
			local sd_p95=sqrt(`v1'+`v2'-`v3')
		}
		else{
			qui sum `p95d_`j''
			local sd_p95=r(sd)
		}
		if `sd_p95'==0{
			local z_p95=0
		}
		else{
			local z_p95=`m_p95'/`sd_p95'
		}
		qui sum `full_meand_`j'' if `n'==1
		local m_mean=r(mean)
		if (`nc'==2){
			qui sum `meand_`j'_1'
			local v1=r(Var)
			qui sum `meand_`j'_2'
			local v2=r(Var)
			qui sum `meand_`j'_3'
			local v3=r(Var)
			local sd_mean=sqrt(`v1'+`v2'-`v3')
		}
		else{
			qui sum `meand_`j''
			local sd_mean=r(sd)
		}
		if `sd_mean'==0{
			local z_mean=0
		}
		else{
			local z_mean=`m_mean'/`sd_mean'
		}
		
		if strpos("`dummies'","`var'"){
			di "`var'" " (dummy)" 
			di _column(5) "Base" ///
				_column(25) %5.3f `full_mean_`j'' 
			di 	_column(5) "Weighted" ///
				_column(25) %5.3f `full_meanw_`j'' 
			di _column(5) "z-stat:"  ///
				_column(25) %5.2f `z_mean'  
			di "------------------------------------------------------------------------------------------------------"
		}
		else{
			di "`var'" 
			di _column(5) "Base" ///
				_column(25) %5.3f `full_mean_`j'' ///
				_column(37) %5.3f `full_p5_`j'' ///
				_column(49) %5.3f `full_p25_`j'' ///
				_column(61) %5.3f `full_p50_`j'' ///
				_column(73) %5.3f `full_p75_`j'' ///
				_column(85) %5.3f `full_p95_`j'' 
			di _column(5) "Weighted" ///
				_column(25) %5.3f `full_meanw_`j'' ///
				_column(37) %5.3f `full_p5w_`j'' ///
				_column(49) %5.3f `full_p25w_`j'' ///
				_column(61) %5.3f `full_p50w_`j'' ///
				_column(73) %5.3f `full_p75w_`j'' ///
				_column(85) %5.3f `full_p95w_`j'' 
			di _column(5) "z-stat:"  ///
				_column(25) %5.2f `z_mean' ///
				_column(37) %5.2f `z_p5' ///
				_column(49) %5.2f `z_p25' ///
				_column(61) %5.2f `z_p50' ///
				_column(73) %5.2f `z_p75' ///
				_column(85) %5.2f `z_p95' 
			di "------------------------------------------------------------------------------------------------------"
		}
		local j=`j'+1
	}
	
	if ("`box'"!=""){
		preserve
		qui expand 2 if `touse'
		qui sort `n'
		qui by `n': gen d=_n-1 if `touse'
		qui sum d
		qui drop if d==1 & w_mm70<.5
		local j=1
		foreach var of local varlist{
			if strpos("`dummies'","`var'"){
			}
			else{
				graph hbox ``j'', over(d, gap(10) label(nolabels) axis(noline)) nooutsides intensity(0) ylabel(none) ysize(2.5) ytitle("") graphregion(color(white) margin(zero)) bgcolor(white) yscale(off) note("") medtype(cline) medline(lwidth(vthick) lcolor(maroon)) plotregion(margin(2 2 0 0))
				set printcolor asis
				graph export "`box'\box`j'.pdf", replace
				local j=`j'+1
			}
		}	
		restore
	}
	
end
