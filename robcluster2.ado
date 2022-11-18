*! version 2.0.2 20221118 David Veenman

/* 
20221118: 2.0.2     Fixed small inconsistency with S-estimator
20221109: 2.0.1     Some housekeeping: changed scalars to locals
20221106: 2.0.0     Major update by replacing separate estimations with cluster-robust variance matrices for more stable estimates and computationally efficient estimation 
                    Added option to store weights in new variable
                    Scalars with weight info now stored in ereturn
20211228: 1.1.1     Allow for factor variables 
                    Added S options to be passed to MM estimator
20211011: 1.1.0     Added tolerance option
                    Improved small-sample correction for 2way variance matrix to be consistent with reghdfe
                    Added option nohaus to s and mm estimator for faster execution
                    Small fix in scalar drop at bottom; previous "drop _all" caused all scalars to be dropped outside the program too
20210701: 1.0.0     First version
*/

program define robcluster2, eclass sortpreserve
	syntax varlist(numeric fv) [in] [if], cluster(varlist) [eff(real 0)] [m|s|median] [biweight] [tol(real 0)] [sopts(str)] [weightvar(str)]
	
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
	    local options="eff(`eff') nohaus tol(`tolerance') sopts(`sopts')"	    
	}
	if "`median'"!="" {
	    local est="q"
	    local options="tol(`tolerance')"	    
	}
    local est0="robreg"
	
	if ("`est'"=="mm" | "`est'"=="s"){
		local biweight="biweight"
	}
	
	tokenize `varlist'
	marksample touse
	local depv `"`1'"'
	// Ensure dv is not a factor variable:
	_fv_check_depvar `depv'
	macro shift 1
    local indepv "`*'"
	// Check and expand factor variable list:
	fvexpand `indepv'
	local indepv `r(varlist)'
	local fvcheck `r(fvops)'
	
	local nc: word count `cluster'
	if (`nc'!=2){
	    di as text "ERROR: Number of dimensions to cluster on is two"
		exit
	}
	local clusterdim1: word 1 of `cluster'
	local clusterdim2: word 2 of `cluster'
		
	// Create temporary variables: 
	tempvar bf_0 bt_0 bft_0 bw n intersection unique_obs res w
	qui egen `intersection'=group(`clusterdim1' `clusterdim2') if `touse'

	if ("`median'"=="") {
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		// Step 0: Obtain coefficients, scale, and residuals from robust estimation:
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		di ""
		di as text "STEP 1: Obtaining robust regression estimates, residuals, and scale estimate..."
		di ""
		capture `est0' `est' `varlist' if `touse', `options' nose

		if _rc>0 {
			di "Estimation failed - trying again... try 2"
			capture noisily `est0' `est' `varlist' if `touse', `options' nose
		}
		if _rc>0 {
			di "Estimation failed - trying again... try 3"
			capture noisily `est0' `est' `varlist' if `touse', `options' nose
		}
		if _rc>0 {
			di "Estimation failed - trying again... try 4"
			capture noisily `est0' `est' `varlist' if `touse', `options' nose
		}
		if _rc>0 {
			di "Estimation failed - trying again... final try"
			capture noisily `est0' `est' `varlist' if `touse', `options' nose
		}
		
		local e_N=e(N)
		local e_r2_p=e(r2_p)
		predict double `res', res
		predict double `w', weights
		scalar scale=e(scale)
		scalar krob=e(k)
		local indepvnames "`indepv' _cons"
		if ("`est'"=="mm" | "`est'"=="s") {
			local j=1
			foreach var of local indepvnames{
				matrix b_`j'=_b[`var']
				matrix colnames b_`j' = `var'
				if (`j'==1){
					matrix beta=(b_`j')
				}
				else{
					matrix beta=(beta, b_`j')
				}
				local j=`j'+1
			}
		}		
		else{
			matrix beta=e(b)
		}
		if ("`weightvar'"!="") {
			capture drop `weightvar'
			predict `weightvar', w 
		}
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
		
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		// Step 1: Getting cluster-robust VCE for first dimension 
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		di as text "STEP 2: Obtaining cluster-robust variance estimators..."
		di "STEP 2A: First clustering dimension..."
		tempvar clusterid
		egen `clusterid'=group(`clusterdim1') if `touse'
		sort `clusterid' 
		local res "`res'"	
		local cvar "`clusterid'"	
		mata: _vce_cluster()
		local nclusterdim1=mata_nclusters
		local e_df_r1=mata_nclusters-1
		matrix colnames Vclust=`indepvnames'
		matrix rownames Vclust=`indepvnames'	
		matrix V1=Vclust
		
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		// Step 2: Getting cluster-robust VCE for second dimension 
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		di "STEP 2B: Second clustering dimension..."
		tempvar clusterid
		egen `clusterid'=group(`clusterdim2') if `touse'
		sort `clusterid' 
		local cvar "`clusterid'"	
		mata: _vce_cluster()
		local nclusterdim2=mata_nclusters
		local e_df_r2=mata_nclusters-1
		matrix colnames Vclust=`indepvnames'
		matrix rownames Vclust=`indepvnames'	
		matrix V2=Vclust
		
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		// Step 3: Getting cluster-robust VCE for second dimension 
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		di "STEP 2C: Intersection of the two clustering dimensions..."
		tempvar clusterid
		egen `clusterid'=group(`intersection') if `touse'
		sort `clusterid' 
		local cvar "`clusterid'"	
		mata: _vce_cluster()
		matrix colnames Vclust=`indepvnames'
		matrix rownames Vclust=`indepvnames'	
		matrix V3=Vclust
	}
	
	else{
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		// Steps 1-3 for median regression 
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		qui `est0' `est' `varlist' if `touse', `options' cluster(`clusterdim1')
		local e_N=e(N)
		local e_r2_p=e(r2_p)
		matrix beta=e(b)
		matrix V1=e(V)/((e(N_clust)/(e(N_clust)-1))*((e(N)-1)/(e(N)-e(rank))))
		local nclusterdim1=e(N_clust)
		local e_df_r1=e(N_clust)-1
		qui `est0' `est' `varlist' if `touse', `options' cluster(`clusterdim2')
		matrix V2=e(V)/((e(N_clust)/(e(N_clust)-1))*((e(N)-1)/(e(N)-e(rank))))
		local nclusterdim2=e(N_clust)
		local e_df_r2=e(N_clust)-1
		qui `est0' `est' `varlist' if `touse', `options' cluster(`intersection')
		matrix V3=e(V)/((e(N_clust)/(e(N_clust)-1))*((e(N)-1)/(e(N)-e(rank))))
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	// Step 4: Combining VCEs and applying small sample correction: 
	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	matrix Vc=V1+V2-V3

	local N=`e_N'
	if "`fvcheck'"=="true"{
		local K=rowsof(Vc)-1
	}
	else{
		local K=rowsof(Vc) 
	}
	local factor1=(`nclusterdim1'/(`nclusterdim1'-1))*((`N'-1)/(`N'-`K'))
	local factor2=(`nclusterdim2'/(`nclusterdim2'-1))*((`N'-1)/(`N'-`K'))
				
	if `nclusterdim1'<`nclusterdim2'{
		local factormin=`factor1'
	}
	else{
		local factormin=`factor2'
	}
	matrix Vc=`factormin'*Vc
	
	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	// Step 5: Post resuls in e() and print results:
	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	ereturn clear
	tempname b V
	matrix `b'=beta
	matrix `V'=Vc
	
	ereturn post `b' `V'
	ereturn scalar N=`e_N'
	ereturn scalar r2_p=`e_r2_p'
	ereturn local depvar "`depv'"
	if `e_df_r1'<`e_df_r2'{
		local e_df_r=`e_df_r1'
	}
	else{
	    local e_df_r=`e_df_r2'
	}
	ereturn scalar df_r=`e_df_r'

	di " "
	di "---------------------------------------------------------"
	if("`m'"=="" & "`s'"=="" & "`median'"==""){
		di in green "MM-estimator with `eff'% efficiency and 2D clustered SEs" 
	}
	if("`m'"!="" & "`biweight'"==""){
		di in green "M-estimator (Huber) with `eff'% efficiency and 2D clustered SEs" 
	}
	if("`m'"!="" & "`biweight'"!=""){
		di in green "M-estimator (biweight) with `eff'% efficiency and 2D clustered SEs" 
	}
	if("`s'"!=""){
		di in green "S-estimator with 2D clustered SEs" 
	}
	if("`median'"!=""){
		di in green "Median regression with 2D clustered SEs" 
	}
	di " "
	di in green "Number of clusters (`clusterdim1') = " _column(31) %5.0f in yellow `nclusterdim1' ///
		_column (56) in green "Number of obs = " %7.0f in yellow e(N)
	di in green "Number of clusters (`clusterdim2') = " _column(31) %5.0f in yellow `nclusterdim2' ///
		_column (56) in green "Degr. freedom = " %7.0f in yellow e(df_r)
	di 	_column (56) in green "Pseudo R2 =     " %7.4f in yellow e(r2_p)
			
	ereturn display
    di "SE clustered by " "`clusterdim1'" " and " "`clusterdim2'" 
	di " " 

	if "`median'"=="" {
		ereturn scalar avweight=`avweight'
		ereturn scalar fraclow=`fraclow'
		ereturn scalar fraczero=`fraczero'
		di as text "Average weight assigned to observations:               " %7.4f in yellow `avweight'
		di as text "Fraction of observations with weight < 0.50:           " %7.4f in yellow `fraclow'
		di as text "Fraction of observations with weight equal to zero:    " %7.4f in yellow `fraczero'
	}
	di " " 

	matrix drop _all
end

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
// Mata program
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

mata:
	mata clear
	void _vce_cluster() {
		st_view(y=., ., st_local("depv"), st_local("touse"))
		st_view(X=., ., tokens(st_local("indepv")), st_local("touse"))
		st_view(w=., ., tokens(st_local("w")), st_local("touse"))
		st_view(r=., ., tokens(st_local("res")), st_local("touse"))
		st_view(cvar=., ., tokens(st_local("cvar")), st_local("touse"))
		st_view(cvarn=., ., tokens(st_local("cvarn")), st_local("touse"))
		scale=st_numscalar("scale")		
		krob=st_numscalar("krob")
		bi=st_local("biweight")
		X=(X,J(rows(X),1,1))
		
		n1=nonmissing(y)
		n2=nonmissing(rowsum(X, 1))
		if (n1<n2) n=n1 
		else n=n2
		k=cols(X)

		XWXinv=invsym(quadcross(X,w,X))
		z=r:/scale
		if (bi=="biweight"){
			psi=mm_biweight_psi(z,krob)
			phi=mm_biweight_phi(z,krob)			
		}
		else{
			psi=mm_huber_psi(z,krob)
			phi=mm_huber_phi(z,krob)			
		}
				
		XphiXinv=invsym(quadcross(X,phi,X))
        info=panelsetup(cvar, 1)
        nc=rows(info)
        M=J(k,k,0)
        dfr=nc-1
        for(i=1; i<=nc; i++) {
            xi=panelsubmatrix(X,i,info)
            psii=panelsubmatrix(psi,i,info)
			mi=xi:*psii
			M=M+xi'*(psii*psii')*xi
        }
		
        Vclust=makesymmetric(scale^2*XphiXinv*M*XphiXinv)
		
		st_matrix("Vclust",Vclust)
		st_numscalar("mata_nclusters",nc)
	}
end
	