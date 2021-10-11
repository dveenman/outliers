*! version 1.1 20211011 David Veenman

/* 
20211011: 1.1	Dropped capture to allow for hard exit from loop
				Added tolerance option for initial estimates of coefficients
				Added option nohaus to s and mm estimator, and nor2 to all, for faster execution
				Small fix in scalar drop at bottom; "drop _all" causes all scalars to be dropped outside the program too
20210701: 1.0	First version
*/

program define roboot, eclass sortpreserve
	syntax varlist [in] [if], [cluster(varlist)] nboot(integer) [eff(real 0)] [m|s|median] [biweight] [tol(real 0)]
		
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
	    local options0="nohaus tol(`tolerance')"	    
	    local options="tolerance(1e-3) nose nor2 nohaus"
	}
	if ("`m'"!="" & "`biweight'"!="") {
		local est="m"
	    local options0="eff(`eff') biw tol(`tolerance')"	    
	    local options="tolerance(1e-3) eff(`eff') nose biw nor2"
	}
	if ("`m'"!="" & "`biweight'"=="") {
		local est="m"
	    local options0="eff(`eff') tol(`tolerance')"	    
	    local options="tolerance(1e-3) eff(`eff') nose nor2"
	}
	if ("`m'"=="" & "`s'"=="") {
	    local est="mm"
	    local options0="eff(`eff') nohaus tol(`tolerance')"	    
	    local options="tolerance(1e-3) eff(`eff') sopts(tolerance(1e-3)) nor2 nohaus"	    
	}
	if "`median'"!="" {
	    local est="q"
	    local options0="tol(`tolerance')"	    
	    local options="nor2"	    
	}
    local est0="robreg"
	
	tokenize `varlist'
	marksample touse
	local depv `"`1'"'
	macro shift 1
    local indepv "`*'"
	
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
	
	/* Obtain vector of coefficients: */
	qui capture `est0' `est' `varlist' if `touse', `options0' 
	scalar e_N=e(N)
	scalar e_r2_p=e(r2_p)
	if "`median'"=="" {
		tempvar robweight zero w50
		qui predict `robweight' if `touse', weights
		qui reg `varlist' if `touse' [aw=`robweight']
		matrix coef=e(b)
		qui gen `zero'=0 if `touse'
		qui replace `zero'=1 if `touse' & `robweight'==0
		qui gen `w50'=0 if `touse'
		qui replace `w50'=1 if `touse' & `robweight'<.5
		qui sum `robweight' if `touse' 
		local avweight=r(mean)
		qui sum `zero' if `touse'
		local fraczero=r(mean)
		qui sum `w50' if `touse'
		local fraclow=r(mean)
	}
	else{
		matrix coef=e(b)
	}	
	
	/* Create temporary variables: */
	tempvar bf_0 bt_0 bft_0 bw n unique_obs
	* The following are variables that will contain the bootstrapped coefficients
	local j=1
	foreach var of local indepv{
		tempvar bf_`j'
		qui gen `bf_`j''=.
		tempvar bt_`j'
		qui gen `bt_`j''=.
		tempvar bft_`j'
		qui gen `bft_`j''=.
		local j=`j'+1
	}
	
	qui gen `bf_0'=.
	qui gen `bt_0'=.
	qui gen `bft_0'=.
	qui gen `bw'=.
	qui gen `n'=_n
	
	/* Create variable for intersection between two clustering dimensions: */
	if (`nc'==2){	    
		tempvar intersection
		qui egen `intersection'=group(`fcluster' `tcluster') if `touse'
	}	
	/* Store number of clusters by cluster dimension */
	if (`nc'>0){
		qui reg `varlist' if `touse', robust cluster(`fcluster')
		local nfcluster=e(N_clust)
	}
	if (`nc'>1){
		qui reg `varlist' if `touse', robust cluster(`tcluster')
		local ntcluster=e(N_clust)
	}
	
	di ""
	
	/* Estimation without clustering */
	if (`nc'==0){
		/* Perform estimations for each bootstrap sample  */
		di in green "Performing bootstrap resampling:"
		forvalues i=1(1)`nboot'{
			qui preserve
			qui bsample if `touse'
			qui `est0' `est' `varlist' if `touse', `options'
			qui restore
			qui replace `bf_0'=_b[_cons] if `n'==`i'	
			local j=1
			foreach var of local indepv{
				qui replace `bf_`j''=_b[``j''] if `n'==`i'	
				local j=`j'+1
			}
			di ".", _continue
		}

		/* Store bootstrap sample variances of coeff per variable in matrix Vf */
		local j=1
		foreach var of local indepv{
			qui sum `bf_`j''
			matrix vf_`j'=r(Var)
			matrix rownames vf_`j' = ``j''
			if (`j'==1){
				matrix Vf=(vf_`j')
			}
			else{
				matrix Vf=(Vf \ vf_`j')
			}
			local j=`j'+1
		}
		/* Add variance of intercept coefficients */	
		qui sum `bf_0'
		matrix cons_f=r(Var)	
		matrix rownames cons_f = _cons
		matrix Vf=(Vf \ cons_f)
		/* Transform to symmetric diagonal matrices */		
		matrix Vc=diag(Vf)
	}

	/* Estimation with cluster-adjustment in one dimension */
	if (`nc'==1){
		/* Perform estimations for each bootstrap sample  */
		di in green "Performing bootstrap resampling:"
		forvalues i=1(1)`nboot'{
			qui preserve
			qui bsample if `touse', cluster(`fcluster') 
			qui `est0' `est' `varlist' if `touse', `options'
			qui restore
			qui replace `bf_0'=_b[_cons] if `n'==`i'	
			local j=1
			foreach var of local indepv{
				qui replace `bf_`j''=_b[``j''] if `n'==`i'	
				local j=`j'+1
			}
			di ".", _continue
		}

		/* Store bootstrap sample variances of coeff per variable in matrices Vf */
		local j=1
		foreach var of local indepv{
			qui sum `bf_`j''
			matrix vf_`j'=r(Var)
			matrix rownames vf_`j' = ``j''
			if (`j'==1){
				matrix Vf=(vf_`j')
			}
			else{
				matrix Vf=(Vf \ vf_`j')
			}
			local j=`j'+1
		}
		/* Add variance of intercept coefficients */	
		qui sum `bf_0'
		matrix cons_f=r(Var)	
		matrix rownames cons_f = _cons
		matrix Vf=(Vf \ cons_f)
		/* Transform to symmetric diagonal matrices */		
		matrix Vc=diag(Vf)
	}
	
	/* Estimation with cluster-adjustment in two dimensions */
	if (`nc'==2){	    
		/* Perform estimations for each bootstrap sample 
		(1) Cluster bootstrap samples by first dimension (e.g., firm) to derive Vf 
		(2) Cluster bootstrap samples by second dimension (e.g., time) to derive Vt
		(3) Cluster bootstrap samples by intersection of the dimensions to derive Vft */
		di in green "Performing bootstrap resampling:"
		forvalues i=1(1)`nboot'{
			/* (1) First dimension */
			qui preserve
			qui bsample if `touse', cluster(`fcluster') 
			qui `est0' `est' `varlist' if `touse', `options'
			qui restore
			qui replace `bf_0'=_b[_cons] if `n'==`i'	
			local j=1
			foreach var of local indepv{
				qui replace `bf_`j''=_b[``j''] if `n'==`i'	
				local j=`j'+1
			}
			/* (2) Second dimension */
			qui preserve
			qui bsample if `touse', cluster(`tcluster') 
			qui `est0' `est' `varlist' if `touse', `options'
			qui restore
			qui replace `bt_0'=_b[_cons] if `n'==`i'	
			local j=1
			foreach var of local indepv{
				qui replace `bt_`j''=_b[``j''] if `n'==`i'	
				local j=`j'+1
			}
			/* (3) Intersection */
			qui preserve
			qui bsample if `touse', cluster(`intersection') 
			qui `est0' `est' `varlist' if `touse', `options'
			qui restore
			qui replace `bft_0'=_b[_cons] if `n'==`i'	
			local j=1
			foreach var of local indepv{
				qui replace `bft_`j''=_b[``j''] if `n'==`i'	
				local j=`j'+1
			}
			di ".", _continue
		}

		/* Store bootstrap sample variances of coeff per variable in matrices Vf, Vt, and Vft */
		local j=1
		foreach var of local indepv{
			qui sum `bf_`j''
			matrix vf_`j'=r(Var)
			matrix rownames vf_`j' = ``j''
			qui sum `bt_`j''
			matrix vt_`j'=r(Var)
			matrix rownames vt_`j' = ``j''
			qui sum `bft_`j''
			matrix vft_`j'=r(Var)
			matrix rownames vft_`j' = ``j''
			if (`j'==1){
				matrix Vf=(vf_`j')
				matrix Vt=(vt_`j')
				matrix Vft=(vft_`j')
			}
			else{
				matrix Vf=(Vf \ vf_`j')
				matrix Vt=(Vt \ vt_`j')
				matrix Vft=(Vft \ vft_`j')
			}
			local j=`j'+1
		}
		/* Add variance of intercept coefficients */	
		qui sum `bf_0'
		matrix cons_f=r(Var)	
		matrix rownames cons_f = _cons
		qui sum `bt_0'
		matrix cons_t=r(Var)	
		matrix rownames cons_t = _cons
		qui sum `bft_0'
		matrix cons_ft=r(Var)	
		matrix rownames cons_ft = _cons	
		matrix Vf=(Vf \ cons_f)
		matrix Vt=(Vt \ cons_t)
		matrix Vft=(Vft \ cons_ft)
		/* Transform to symmetric diagonal matrices */	
		matrix Vf=diag(Vf)
		matrix Vt=diag(Vt)
		matrix Vft=diag(Vft)
		/* Create 2-dimension cluster-adjusted variance matrix */
		matrix Vc=Vf+Vt-Vft
	}
		
	/* Post results in e() */
	ereturn clear
	tempname b V
	matrix `b' = coef
	matrix `V' = Vc
	ereturn post `b' `V'
	ereturn scalar N=e_N
	ereturn scalar r2_p=e_r2_p
	ereturn local depvar "`depv'"

	/* Print results */
	di " "
	di " "
	if (`nc'==0){
		if("`m'"=="" & "`s'"=="" & "`median'"==""){
			di in green "MM-estimator with `eff'% efficiency and bootstrapped SEs" 
		}
		if("`m'"!="" & "`biweight'"==""){
			di in green "M-estimator (Huber) with `eff'% efficiency and bootstrapped SEs" 
		}
		if("`m'"!="" & "`biweight'"!=""){
			di in green "M-estimator (biweight) with `eff'% efficiency and bootstrapped SEs" 
		}
		if("`s'"!=""){
			di in green "S-estimator with bootstrapped SEs" 
		}
		if("`median'"!=""){
			di in green "Median regression with bootstrapped SEs" 
		}
	}
	if (`nc'==1){
		if("`m'"=="" & "`s'"=="" & "`median'"==""){
			di in green "MM-estimator with `eff'% efficiency and bootstrapped clustered SEs" 
		}
		if("`m'"!="" & "`biweight'"==""){
			di in green "M-estimator (Huber) with `eff'% efficiency and bootstrapped clustered SEs" 
		}
		if("`m'"!="" & "`biweight'"!=""){
			di in green "M-estimator (biweight) with `eff'% efficiency and bootstrapped clustered SEs" 
		}
		if("`s'"!=""){
			di in green "S-estimator with bootstrapped clustered SEs" 
		}
		if("`median'"!=""){
			di in green "Median regression with bootstrapped clustered SEs" 
		}
	}
	if (`nc'==2){
		if("`m'"=="" & "`s'"=="" & "`median'"==""){
			di in green "MM-estimator with `eff'% efficiency and bootstrapped 2D clustered SEs" 
		}
		if("`m'"!="" & "`biweight'"==""){
			di in green "M-estimator (Huber) with `eff'% efficiency and bootstrapped 2D clustered SEs" 
		}
		if("`m'"!="" & "`biweight'"!=""){
			di in green "M-estimator (biweight) with `eff'% efficiency and bootstrapped 2D clustered SEs" 
		}
		if("`s'"!=""){
			di in green "S-estimator with bootstrapped 2D clustered SEs" 
		}
		if("`median'"!=""){
			di in green "Median regression with bootstrapped 2D clustered SEs" 
		}
	}
	di " "
	di in green "Number of bootstrap samples = " _column(31) %5.0f in yellow `nboot' ///
		_column(56) in green "Number of obs = " %7.0f in yellow e(N)
	if (`nc'==0){
	    di _column(56) in green "Pseudo R2 =     " %7.4f in yellow e(r2_p)
	}
	if (`nc'>0){
		di in green "Number of clusters (`fcluster') = " _column(31) %5.0f in yellow `nfcluster' ///
			_column(56) in green "Pseudo R2 =     " %7.4f in yellow e(r2_p)
	}
	if (`nc'>1){
		di in green "Number of clusters (`tcluster') = " _column(31) %5.0f in yellow `ntcluster' 
	}
			
	ereturn display
	if (`nc'==0){
	    di "SE not adjusted for clustering" 
	}
	if (`nc'==1){
	    di "SE clustered by " "`fcluster'" 
	}
	if (`nc'==2){
	    di "SE clustered by " "`fcluster'" " and " "`tcluster'" 
	}	
	di " " 

	if "`median'"=="" {
		di as text "Average weight assigned to observations:               " %7.4f in yellow `avweight'
		di as text "Fraction of observations with weight < 0.50:           " %7.4f in yellow `fraclow'
		di as text "Fraction of observations with weight equal to zero:    " %7.4f in yellow `fraczero'
	}
	di " " 
	
	scalar drop e_N e_r2_p
	matrix drop _all

end
