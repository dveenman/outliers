*! version 1.1 20211011 David Veenman

/* 
20211011: 1.1	Added tolerance option 
				Added option nohaus to s and mm estimator for faster execution
20210701: 1.0	First version
*/

program define robfmb, eclass sortpreserve 
	syntax varlist [in] [if], splitvar(varname) [eff(real 0)] [lag(int 0)] [m|s|median] [biweight] [tol(real 0)]
	
	if (`tol'==0){
	    local tolerance=1e-6
	}
	else {
	    local tolerance=`tol'
	}
		
	if ("`median'"=="" & "`s'"=="" & `eff'==0) {
	    if ("`m'"!=""){
			di as text "ERROR: You must specify the desired estimation efficiency in option eff() with M estimation"
		}
		else{
		    di as text "ERROR: You must specify the desired estimation efficiency in option eff() with MM estimation"
		}
		exit
	}
	if ("`m'"!="" & "`s'"!="") {
	    di as text "ERROR: Cannot specify both S and M estimation"
		exit
	}
	if ("`median'"!="" & "`s'"!="") {
	    di as text "ERROR: Cannot specify both S estimation and median regression"
		exit
	}
	if ("`median'"!="" & "`s'"!="") {
	    di as text "ERROR: Cannot specify both M estimation and median regression"
		exit
	}
	if "`s'"!="" {
	    local est="s"
	    local options="nohaus tol(`tolerance')"
	}
	if ("`m'"!="" & "`biweight'"!="") {
		local est="m"
	    local options="eff(`eff') biw tol(`tolerance')"
	}
	if ("`m'"!="" & "`biweight'"=="") {
		local est="m"
	    local options="eff(`eff') tol(`tolerance')"
	}
	if ("`m'"=="" & "`s'"=="") {
	    local est="mm"
	    local options="eff(`eff') nohaus tol(`tolerance')"	    
	}
	if "`median'"!="" {
	    local est="q"
	    local options="tol(`tolerance')"	    
	}
    local est0="robreg"
	
	tokenize `varlist'
	marksample touse
	local depv `"`1'"'
	macro shift 1
    local indepv "`*'"
	
	/* Create temporary variables: */
	tempvar b_0 nn splitvarnumeric splitcount pseudor2 
	* The following are variables that will contain the coefficients
	local j=1
	foreach var of local indepv{
		tempvar b_`j'
		qui gen `b_`j''=.
		local j=`j'+1
	}
	
	qui egen `splitvarnumeric'=group(`splitvar') if `touse'
	qui sum `splitvarnumeric' if `touse'
	scalar e_Ntot=r(N)
	local maxsplit=r(max)
	scalar e_N=r(max)
	scalar e_df_r=e_N-1
	qui gen `b_0'=.
	qui gen `nn'=_n
	
	qui reg `varlist' if `touse'
	local nparm = e(df_m)+1
	qui egen `splitcount'=count(`splitvarnumeric') if `touse', by(`splitvarnumeric')
	qui sum `splitcount' if `touse'
	if (`nparm'>r(min)) {
	    di as text "ERROR: Minimum number of observations in each split must exceed number of parameters"
		exit
	}
	
	di ""
	di as text "Number of estimations by variable `splitvar': `maxsplit'"
	
	qui gen `pseudor2'=.
	/* Estimate regression for each split (e.g., year) */
	forvalues i=1(1)`maxsplit'{
		qui capture `est0' `est' `varlist' if `touse' & `splitvarnumeric'==`i', `options'
		qui replace `b_0'=_b[_cons] if `nn'==`i'	
		local j=1
		foreach var of local indepv{
			qui replace `b_`j''=_b[``j''] if `nn'==`i'	
			local j=`j'+1
		}
		qui replace `pseudor2'=e(r2_p) if `nn'==`i'	
		di ".", _continue
	}
	
	qui preserve
	qui keep if `nn'<=`maxsplit'
	
	/* Store coefficient means and variances per variable in matrix */
	local j=1
	foreach var of local indepv{
		/* Normal standard errors */
		if (`lag'==0){
			qui reg `b_`j''
		}
		/* Newey-West standard errors with t lags (should only be used when split is based on a time variable)*/
		else{
			qui tsset `nn'
			qui newey `b_`j'', lag(`lag')
		}
		matrix b_`j'=e(b)
		matrix v_`j'=e(V)
		matrix colnames b_`j' = ``j''
		matrix rownames v_`j' = ``j''
		if (`j'==1){
			matrix coef=(b_`j')
			matrix Vmatrix=(v_`j')
		}
		else{
			matrix coef=(coef, b_`j')
			matrix Vmatrix=(Vmatrix \ v_`j')
		}
		local j=`j'+1
	}
	
	/* Add average and variance of intercept coefficients */	
	qui sum `b_0'
	matrix cons_b=r(mean)	
	matrix cons_V=r(Var)	
	matrix colnames cons_b = _cons
	matrix rownames cons_V = _cons
	matrix coef=(coef, cons_b)
	matrix Vmatrix=(Vmatrix \ cons_V)

	/* Transform variance matrix to symmetric diagonal matrix */		
	matrix Vmatrix=diag(Vmatrix)

	/* Post results in e() */
	qui sum `pseudor2'
	scalar e_r2_p=r(mean)
	ereturn clear
	tempname b V
	matrix `b' = coef
	matrix `V' = Vmatrix
	ereturn post `b' `V'
	ereturn scalar Ntot=e_Ntot
	ereturn scalar N=e_N
	ereturn scalar df_r=e_df_r
	ereturn scalar r2_p=e_r2_p
	ereturn local depvar "`depv'"

	/* Print results */
	di " "
	di " "
	if("`m'"=="" & "`s'"=="" & "`median'"==""){
		di in green "MM-estimator with `eff'% efficiency and Fama-Macbeth SEs" 
	}
	if("`m'"!="" & "`biweight'"==""){
		di in green "M-estimator (Huber) with `eff'% efficiency and Fama-Macbeth SEs" 
	}
	if("`m'"!="" & "`biweight'"!=""){
		di in green "M-estimator (biweight) with `eff'% efficiency and Fama-Macbeth SEs" 
	}
	if("`s'"!=""){
		di in green "S-estimator with Fama-Macbeth SEs" 
	}
	if("`median'"!=""){
		di in green "Median regression with Fama-Macbeth SEs" 
	}
	if (`lag'!=0){
		di "Newey-West correction for SEs based on `lag' lags"
	}
	di in green "Sample splitted for separate estimations by variable `splitvar'"
	di " "
	di _column(33) in green "Number of observations full sample =   " %7.0f in yellow e(Ntot)
	di _column(33) in green "Number of sample splits/estimations =  " %7.0f in yellow e(N)
	di _column(33) in green "Average pseudo R2 =                    " %7.4f in yellow e(r2_p)
			
	ereturn display

	scalar drop e_N e_Ntot e_df_r e_r2_p
	matrix drop coef Vmatrix cons_b cons_V
	local j=1
	foreach var of local indepv{
		matrix drop b_`j'
		matrix drop v_`j'
		local j=`j'+1
	}
	
	qui restore
end
