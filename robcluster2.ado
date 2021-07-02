*! version 1.0 01july2021 David Veenman

program define robcluster2, eclass sortpreserve
	syntax varlist [in] [if], cluster(varlist) [eff(real 0)] [m|s|median] [biweight]
	
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
	    local options=""
	}
	if ("`m'"!="" & "`biweight'"!="") {
		local est="m"
	    local options="eff(`eff') biw"
	}
	if ("`m'"!="" & "`biweight'"=="") {
		local est="m"
	    local options="eff(`eff')"
	}
	if ("`m'"=="" & "`s'"=="") {
	    local est="mm"
	    local options="eff(`eff')"	    
	}
	if "`median'"!="" {
	    local est="q"
	    local options=" "	    
	}
    local est0="robreg"
	
	tokenize `varlist'
	marksample touse
	local depv `"`1'"'
	macro shift 1
    local indepv "`*'"
	
	local nc: word count `cluster'
	if (`nc'!=2){
	    di as text "ERROR: Number of dimensions to cluster on is two"
		exit
	}
	local fcluster: word 1 of `cluster'
	local tcluster: word 2 of `cluster'
		
	/* Create temporary variables: */
	tempvar bf_0 bt_0 bft_0 bw n intersection unique_obs
	qui egen `intersection'=group(`fcluster' `tcluster') if `touse'
			
	/* Estimation with cluster-adjustment in two dimensions 
	(1) Cluster bootstrap samples by first dimension to derive Vf 
	(2) Cluster bootstrap samples by second dimension to derive Vt
	(3) Cluster bootstrap samples by intersection of dimensions to derive Vft */

	/* (1) First dimension */
	qui `est0' `est' `varlist' if `touse', `options' cluster(`fcluster')
	if "`median'"=="" {
		tempvar weight zero w50
		qui predict `weight' if `touse', weights
		qui gen `zero'=0 if `touse'
		qui replace `zero'=1 if `touse' & `weight'==0
		qui gen `w50'=0 if `touse'
		qui replace `w50'=1 if `touse' & `weight'<.5
		qui sum `weight' if `touse' 
		local avweight=r(mean)
		qui sum `zero' if `touse'
		local fraczero=r(mean)
		qui sum `w50' if `touse'
		local fraclow=r(mean)
	}	
	scalar e_N=e(N)
	scalar e_df_r1=e(df_r)
	scalar e_r2_p=e(r2_p)
	local nfcluster=e(N_clust)
	local j=1
	foreach var of local indepv{
		matrix b_`j'=_b[`var']
		matrix colnames b_`j' = `var'
		matrix v_`j'=(_se[`var'])^2
		matrix colnames v_`j' = `var'
		if (`j'==1){
			matrix coef=(b_`j')
			matrix Vf=(v_`j')
		}
		else{
			matrix coef=(coef, b_`j')
			matrix Vf=(Vf, v_`j')
		}
		local j=`j'+1
	}
	matrix cons=_b[_cons]	
	matrix colnames cons = _cons
	matrix coef=(coef, cons)
	matrix consV=(_se[_cons])^2	
	matrix colnames consV = _cons
	matrix Vf=(Vf, consV)
	matrix Vf=diag(Vf)
				
	/* (2) Second dimension */
	qui `est0' `est' `varlist' if `touse', `options' cluster(`tcluster')
	scalar e_df_r2=e(df_r)
	local ntcluster=e(N_clust)
	local j=1
	foreach var of local indepv{
		matrix v_`j'=(_se[`var'])^2
		matrix colnames v_`j' = `var'
		if (`j'==1){
			matrix Vt=(v_`j')
		}
		else{
			matrix Vt=(Vt, v_`j')
		}
		local j=`j'+1
	}
	matrix consV=(_se[_cons])^2	
	matrix colnames consV = _cons
	matrix Vt=(Vt, consV)
	matrix Vt=diag(Vt)

	/* (3) Intersection */
	qui `est0' `est' `varlist' if `touse', `options' cluster(`intersection')
	local j=1
	foreach var of local indepv{
		matrix v_`j'=(_se[`var'])^2
		matrix colnames v_`j' = `var'
		if (`j'==1){
			matrix Vft=(v_`j')
		}
		else{
			matrix Vft=(Vft, v_`j')
		}
		local j=`j'+1
	}
	matrix consV=(_se[_cons])^2	
	matrix colnames consV = _cons
	matrix Vft=(Vft, consV)
	matrix Vft=diag(Vft)

	/* Create 2-dimension cluster-adjusted variance matrix */
	matrix Vc=Vf+Vt-Vft
		
	/* Post resuls in e() */
	ereturn clear
	tempname b V
	matrix `b'=coef
	matrix `V'=Vc
	ereturn post `b' `V'
	ereturn scalar N=e_N
	ereturn scalar r2_p=e_r2_p
	ereturn local depvar "`depv'"
	if e_df_r1<e_df_r2{
		scalar e_df_r=e_df_r1
	}
	else{
	    scalar e_df_r=e_df_r2
	}
	ereturn scalar df_r=e_df_r

	/* Print results */
	di " "
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
	di " "
	di in green "Number of clusters (`fcluster') = " _column(31) %5.0f in yellow `nfcluster' ///
		_column (56) in green "Number of obs = " %7.0f in yellow e(N)
	di in green "Number of clusters (`tcluster') = " _column(31) %5.0f in yellow `ntcluster' ///
		_column (56) in green "Degr. freedom = " %7.0f in yellow e(df_r)
	di 	_column (56) in green "Pseudo R2 =     " %7.4f in yellow e(r2_p)
			
	ereturn display
    di "SE clustered by " "`fcluster'" " and " "`tcluster'" 
	di " " 

	if "`median'"=="" {
		di as text "Average weight assigned to observations:               " %7.4f in yellow `avweight'
		di as text "Fraction of observations with weight < 0.50:           " %7.4f in yellow `fraclow'
		di as text "Fraction of observations with weight equal to zero:    " %7.4f in yellow `fraczero'
	}
	di " " 

	scalar drop _all
	matrix drop _all
end
