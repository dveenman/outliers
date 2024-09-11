*! version 2.1.0 20240911 David Veenman

/* 
20240911: 2.1.0     Adjusted syntax for consistency with robreg (e.g., "robcluster2 m y x" instead of "robcluster2 y x, m")
                    Added 3-way clustering option
                    Added option for sampling weights (pw)
                    Added noconstant option                    
                    Added Cameron/Gelbach/Miller (2011) adjustment from Gu and Yoo (2019, DOI: 10.1177/1536867X19893637)                    
                    Added multicollinearity check and omit collinear variables with _rmcoll					
                    Changed sopts() option to m() for consistency between S- and MM-estimators
                    Small Mata efficiency improvements					
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
	syntax [anything] [in] [if] [pw], cluster(varlist) [eff(real 0)] [biweight] [noconstant] [tol(real 0)] [m(str)] [weightvar(str)] 

	marksample touse
	
	if "`anything'"=="" {
		if ("`e(cmd)'"!="robreg" & "`e(cmd)'"!="robcluster2") {
			di as err "You should estimate using robreg first before using the postestimation version of the program"
			exit
		}
		if (`eff'!=0 | "`biweight'"!="" | "`constant'"!="" | `tol'!=0 | "`m'"!="" | "`weightvar'"!="") {
			di as err "Only option cluster() is allowed with the postestimation version of the program"
			exit
		}
		local postestimation="postestimation"
		local constant "`e(noconstant)'"
		local eff=round(`e(efficiency)')
		local depv "`e(depvar)'"
		local indepv "`e(indepvars)'"
		local subcmd="`e(subcmd)'"
		if ("`subcmd'"=="mm" | "`subcmd'"=="s"){
			local biweight="biweight"
		}
		local wexp "`e(wexp)'"
		local avweight "`e(avweight)'"
		local fraclow "`e(fraclow)'"
		local fraczero "`e(fraczero)'"
	}
	
	else {
		
		local postest=""
		
		tokenize `anything'
		local subcmd `"`1'"'

		local cmdlist "q m s mm"
		if !`:list subcmd in cmdlist' {
			di as err `"Invalid subcommand: `subcmd'"'
			exit 
		}
		
		macro shift 1
		local depv `"`1'"'
		local varlist `"`*'"'
		
		// Ensure dv is not a factor variable:
		_fv_check_depvar `depv'
		macro shift 1
		local indepv "`*'"
		
		// Check for collinearity and omit one variable if needed:
		_rmcoll `indepv'
		local indepv "`r(varlist)'"
		if ("`m'"!="") {
			qui _rmcoll `m'
			local m "`r(varlist)'"
		}

		// Check and expand factor variable list:
		fvexpand `indepv'
		local indepv `r(varlist)'
		local fvcheck `r(fvops)'
		
		if (`tol'!=0){
			local tolerance="tol(`tol')"
		}
		else {
			local tolerance=""
		}	
		if ("`subcmd'"=="m" & `eff'==0) {
			di as text "ERROR: You must specify the desired estimation efficiency in option eff() with M estimation"
			exit
		}
		if ("`subcmd'"=="mm" & `eff'==0) {
			di as text "ERROR: You must specify the desired estimation efficiency in option eff() with MM estimation"
			exit
		}
		if ("`subcmd'"!="mm" & "`subcmd'"!="s" & "`m'"!="") {
			di as text "ERROR: You can only combine option m() with MM or S estimation"
			exit
		}
		if ("`subcmd'"=="s" & `eff'!=0) {
			di as text "ERROR: Option eff() not allowed with S estimation"
			exit
		}
		if ("`subcmd'"=="q" & `eff'!=0) {
			di as text "ERROR: Option eff() not allowed with median regression"
			exit
		}
		if ("`subcmd'"=="s") {
			if "`m'"!="" {
				local options="nohaus `tolerance' `constant' m(`m')"
			}
			else {
				local options="nohaus `tolerance' `constant'"
			}
		}
		if ("`subcmd'"=="mm") {
			if "`m'"!="" {
				local options="eff(`eff') nohaus `tolerance' `constant' sopts(m(`m'))"
			}
			else {
				local options="eff(`eff') nohaus `tolerance' `constant'"
			}
		}
		if ("`subcmd'"=="m" & "`biweight'"!="") {
			local options="eff(`eff') biw `tolerance' `constant'"
		}
		if ("`subcmd'"=="m" & "`biweight'"=="") {
			local options="eff(`eff') `tolerance' `constant'"
		}
		if ("`subcmd'"=="q") {
			local options="`tolerance' `constant'"	    
		}	
		if ("`subcmd'"=="mm" | "`subcmd'"=="s"){
			local biweight="biweight"
		}
		
	}
	
	local nc: word count `cluster'
	if (`nc'!=2 & `nc'!=3){
	    di as text "ERROR: Number of dimensions to cluster on is two or three"
		exit
	}
	local clusterdim1: word 1 of `cluster'
	local clusterdim2: word 2 of `cluster'
	if `nc'==3 {
		local clusterdim3: word 3 of `cluster'
	}
	
	// Create temporary variables: 
	tempvar bf_0 bt_0 bft_0 bw n clus1 clus2 clus12 unique_obs res wvar
	qui egen `clus1'=group(`clusterdim1') if `touse'
	qui egen `clus2'=group(`clusterdim2') if `touse'
	qui egen `clus12'=group(`clusterdim1' `clusterdim2') if `touse'
	
	if `nc'==3 {
		tempvar clus3 clus13 clus23 clus123
		qui egen `clus3'=group(`clusterdim3') if `touse'
		qui egen `clus13'=group(`clusterdim1' `clusterdim3') if `touse'
		qui egen `clus23'=group(`clusterdim2' `clusterdim3') if `touse'
		qui egen `clus123'=group(`clusterdim1' `clusterdim2' `clusterdim3') if `touse'
	}
	
	if "`postestimation'"=="" {
		if "`weight'"!="" {
			qui gen double `wvar' `exp' if `touse'
			local woption [`weight'=`wvar']
		}
		else {
			qui gen double `wvar'=1 if `touse'
		}
	}
	else {
		if "`e(wexp)'"!="" {
			qui gen double `wvar' `e(wexp)' if `touse'
		}
		else {
			qui gen double `wvar'=1 if `touse'
		}
	}
	
	if ("`subcmd'"!="q") {
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		// Step 0: Obtain coefficients, scale, and residuals from robust estimation:
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		if "`postestimation'"=="" {
			di ""
			di as text "STEP 1: Obtaining robust regression estimates, residuals, and scale estimate..."
			di ""
			capture robreg `subcmd' `depv' `indepv' if `touse' `woption', `options' nose

			if _rc>0 {
				di "Estimation failed - trying again... try 2"
				capture noisily robreg `subcmd' `depv' `indepv' if `touse' `woption', `options' nose
			}
			if _rc>0 {
				di "Estimation failed - trying again... try 3"
				capture noisily robreg `subcmd' `depv' `indepv' if `touse' `woption', `options' nose
			}
			if _rc>0 {
				di "Estimation failed - trying again... try 4"
				capture noisily robreg `subcmd' `depv' `indepv' if `touse' `woption', `options' nose
			}
			if _rc>0 {
				di "Estimation failed - trying again... final try"
				capture noisily robreg `subcmd' `depv' `indepv' if `touse' `woption', `options' nose
			}			
		}
		
		local e_N=e(N)
		local e_r2_p=e(r2_p)
		scalar scale=e(scale)
		scalar krob=e(k)
				
		if "`e(cmd)'"=="robreg" {
			predict double `res', res
		}
			
		if ("`postestimation'"!="" & "`e(cmd)'"=="robreg") {
			capture drop _robcluster2_res
			qui predict double _robcluster2_res if `touse', res
		}
		
		if ("`postestimation'"!="" & "`e(cmd)'"=="robcluster2") {
			qui gen double `res'=_robcluster2_res if `touse'
		}	
		
		if "`postestimation'"=="" {
			if "`constant'"=="" {
				local indepvnames "`indepv' _cons"
			}
			else {
				local indepvnames "`indepv'"
			}
		}
		else {
			if "`e(noconstant)'"=="" {
				local indepvnames "`e(indepvars)' _cons"
			}
			else {
				local indepvnames "`e(indepvars)'"
			}
		}
		
		if ("`subcmd'"=="mm" | "`subcmd'"=="s") {
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
		
		if "`e(cmd)'"=="robreg" {
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
		
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		// Step 1: Getting cluster-robust VCE for first dimension 
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		if "`postestimation'"=="" {
			di as text "STEP 2: Obtaining cluster-robust variance estimator"
		}
		else {
			di as text "Obtaining cluster-robust variance estimator"
		}
		sort `clus1' 
		local res "`res'"	
		local cvar "`clus1'"	
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
		sort `clus2' 
		local cvar "`clus2'"	
		mata: _vce_cluster()
		local nclusterdim2=mata_nclusters
		local e_df_r2=mata_nclusters-1
		matrix colnames Vclust=`indepvnames'
		matrix rownames Vclust=`indepvnames'	
		matrix V2=Vclust
		
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		// Step 3: Getting cluster-robust VCE for interaction between dimension 1 and 2
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		sort `clus12' 
		local cvar "`clus12'"	
		mata: _vce_cluster()
		matrix colnames Vclust=`indepvnames'
		matrix rownames Vclust=`indepvnames'	
		matrix V12=Vclust
				
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		// Step 3*: Getting additional VCEs for 3-way clustering
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		if `nc'==3 {
			// Dimension 3:
			sort `clus3' 
			local cvar "`clus3'"	
			mata: _vce_cluster()
			local nclusterdim3=mata_nclusters
			local e_df_r3=mata_nclusters-1
			matrix colnames Vclust=`indepvnames'
			matrix rownames Vclust=`indepvnames'	
			matrix V3=Vclust
			// Dimension 1/3:
			sort `clus13' 
			local cvar "`clus13'"	
			mata: _vce_cluster()
			matrix colnames Vclust=`indepvnames'
			matrix rownames Vclust=`indepvnames'	
			matrix V13=Vclust
			// Dimension 2/3:
			sort `clus23' 
			local cvar "`clus23'"	
			mata: _vce_cluster()
			matrix colnames Vclust=`indepvnames'
			matrix rownames Vclust=`indepvnames'	
			matrix V23=Vclust
			// Dimension 1/2/3:
			sort `clus123' 
			local cvar "`clus123'"	
			mata: _vce_cluster()
			matrix colnames Vclust=`indepvnames'
			matrix rownames Vclust=`indepvnames'	
			matrix V123=Vclust
		}
		
	}
	
	else{
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		// Steps for median regression 
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		qui robreg q `depv' `indepv' if `touse' `woption', `options' cluster(`clus1')
		local e_N=e(N)
		local e_r2_p=e(r2_p)
		matrix beta=e(b)
		matrix V1=e(V)/((e(N_clust)/(e(N_clust)-1))*((e(N)-1)/(e(N)-e(rank))))
		local nclusterdim1=e(N_clust)
		local e_df_r1=e(N_clust)-1
		qui robreg q `depv' `indepv' if `touse' `woption', `options' cluster(`clus2')
		matrix V2=e(V)/((e(N_clust)/(e(N_clust)-1))*((e(N)-1)/(e(N)-e(rank))))
		local nclusterdim2=e(N_clust)
		local e_df_r2=e(N_clust)-1
		qui robreg q `depv' `indepv' if `touse' `woption', `options' cluster(`clus12')
		matrix V12=e(V)/((e(N_clust)/(e(N_clust)-1))*((e(N)-1)/(e(N)-e(rank))))
		
		if `nc'==3 {
			qui robreg q `depv' `indepv' if `touse' `woption', `options' cluster(`clus3')
			matrix V3=e(V)/((e(N_clust)/(e(N_clust)-1))*((e(N)-1)/(e(N)-e(rank))))
			local nclusterdim3=e(N_clust)
			local e_df_r3=e(N_clust)-1
			qui robreg q `depv' `indepv' if `touse' `woption', `options' cluster(`clus13')
			matrix V13=e(V)/((e(N_clust)/(e(N_clust)-1))*((e(N)-1)/(e(N)-e(rank))))
			qui robreg q `depv' `indepv' if `touse' `woption', `options' cluster(`clus23')
			matrix V23=e(V)/((e(N_clust)/(e(N_clust)-1))*((e(N)-1)/(e(N)-e(rank))))
			qui robreg q `depv' `indepv' if `touse' `woption', `options' cluster(`clus123')
			matrix V123=e(V)/((e(N_clust)/(e(N_clust)-1))*((e(N)-1)/(e(N)-e(rank))))
		}
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	// Step 4: Combining VCEs and applying small sample correction: 
	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	if `nc'==2 {
		matrix Vc=V1+V2-V12
	}
	else {
		matrix Vc=V1+V2+V3-V12-V13-V23+V123
	}
	
	mata: _check_vce()
	if (negative==1) {
		matrix colnames Vcadj=`indepvnames'
		matrix rownames Vcadj=`indepvnames'	
		matrix Vc=Vcadj
		di ""
		di "Note: adjustment from Cameron, Gelbach, and Miller (2011) applied to non-positive semi-definite VCE"
	}
	
	local N=`e_N'
	if "`fvcheck'"=="true"{
		local K=rowsof(Vc)-1
	}
	else{
		local K=rowsof(Vc) 
	}
	local factor1=(`nclusterdim1'/(`nclusterdim1'-1))*((`N'-1)/(`N'-`K'))
	local factor2=(`nclusterdim2'/(`nclusterdim2'-1))*((`N'-1)/(`N'-`K'))
	if `nc'==3 {
		local factor3=(`nclusterdim3'/(`nclusterdim3'-1))*((`N'-1)/(`N'-`K'))
	}
	
	local factormin=`factor1'
	if `factor2'>`factor1' {
		local factormin=`factor2'
	}
	if `nc'==3 {
		if `factor3'>`factor2' {
			local factormin=`factor3'
		}
	}
		
	matrix Vc=`factormin'*Vc
	
	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	// Step 5: Post results in e() and print results:
	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	ereturn clear
	tempname b V
	matrix `b'=beta
	matrix `V'=Vc
	
	ereturn post `b' `V'
	
	ereturn scalar N=`e_N'
	ereturn scalar r2_p=`e_r2_p'
	if ("`subcmd'"!="q") {
		ereturn scalar k=krob
		ereturn scalar scale=scale
		ereturn scalar avweight=`avweight'
		ereturn scalar fraclow=`fraclow'
		ereturn scalar fraczero=`fraczero'
	}
	
	ereturn local depvar "`depv'"
	ereturn local indepvars "`indepv'"
	ereturn local cmd "robcluster2"
	ereturn local subcmd "`subcmd'"
	ereturn local noconstant "`constant'"
	ereturn local wexp "`wexp'"
	if ("`subcmd'"!="q") {
		ereturn local efficiency "`eff'"		
	}
	
	local e_df_r=`e_df_r1'
	if `e_df_r2'<`e_df_r1' {
		local e_df_r=`e_df_r2'
	}
	if `nc'==3 {
		if `e_df_r3'<`e_df_r2' {
			local e_df_r=`e_df_r3'
		}			
	}
	ereturn scalar df_r=`e_df_r'

	di " "
	di " "
	di "---------------------------------------------------------"
	if("`subcmd'"=="mm"){
		if `nc'==2 {
			di in green "MM-estimator with `eff'% efficiency and 2D clustered SEs" 
		}
		else {
			di in green "MM-estimator with `eff'% efficiency and 3D clustered SEs" 
		}
	}
	if("`subcmd'"=="m" & "`biweight'"==""){
		if `nc'==2 {
			di in green "M-estimator (Huber) with `eff'% efficiency and 2D clustered SEs" 
		}
		else {
			di in green "M-estimator (Huber) with `eff'% efficiency and 3D clustered SEs" 			
		}
	}
	if("`subcmd'"=="m" & "`biweight'"!=""){
		if `nc'==2 {
			di in green "M-estimator (biweight) with `eff'% efficiency and 2D clustered SEs" 
		}
		else {
			di in green "M-estimator (biweight) with `eff'% efficiency and 3D clustered SEs" 			
		}
	}
	if("`subcmd'"=="s"){
		if `nc'==2 {
			di in green "S-estimator with 2D clustered SEs" 
		}
		else {
			di in green "S-estimator with 3D clustered SEs" 			
		}
	}
	if("`subcmd'"=="q"){
		if `nc'==2 {
			di in green "Median regression with 2D clustered SEs" 
		}
		else {
			di in green "Median regression with 3D clustered SEs" 			
		}
	}
	di " "
	di in green "Number of clusters (`clusterdim1') = " _column(31) %5.0f in yellow `nclusterdim1' ///
		_column (56) in green "Number of obs = " %7.0f in yellow e(N)
	di in green "Number of clusters (`clusterdim2') = " _column(31) %5.0f in yellow `nclusterdim2' ///
		_column (56) in green "Degr. freedom = " %7.0f in yellow e(df_r)
	if `nc'==2 {
		di _column (56) in green "Pseudo R2 =     " %7.4f in yellow e(r2_p)
	}
	else {
		di in green "Number of clusters (`clusterdim3') = " _column(31) %5.0f in yellow `nclusterdim3' ///
		_column (56) in green "Pseudo R2 =     " %7.4f in yellow e(r2_p)
	}	
			
	ereturn display
	if `nc'==2 {
		di "SE clustered by " "`clusterdim1'" " and " "`clusterdim2'" 
	}
	else {
		di "SE clustered by " "`clusterdim1'" ", " "`clusterdim2', and " "`clusterdim3'" 
	}
    
	di " " 

	di as text "Average weight assigned to observations:               " %7.4f in yellow `avweight'
	di as text "Fraction of observations with weight < 0.50:           " %7.4f in yellow `fraclow'
	di as text "Fraction of observations with weight equal to zero:    " %7.4f in yellow `fraczero'

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
		st_view(w=., ., tokens(st_local("wvar")), st_local("touse"))
		st_view(r=., ., tokens(st_local("res")), st_local("touse"))
		st_view(cvar=., ., tokens(st_local("cvar")), st_local("touse"))
		scale=st_numscalar("scale")		
		krob=st_numscalar("krob")
		bi=st_local("biweight")
		noc=st_local("constant")
		n=st_local("e_N")
		
		if (noc!="noconstant"){
			X=(X,J(rows(X),1,1))
		}
				
		n1=nonmissing(y)
		n2=nonmissing(rowsum(X, 1))
		if (n1<n2) n=n1 
		else n=n2
		k=cols(X)

		z=r:/scale
		if (bi=="biweight"){
			psi=w:*mm_biweight_psi(z,krob)
			phi=w:*mm_biweight_phi(z,krob)			
		}
		else{
			psi=w:*mm_huber_psi(z,krob)
			phi=w:*mm_huber_phi(z,krob)			
		}
	
		XphiXinv=invsym(quadcross(X,phi,X))

		info=panelsetup(cvar, 1)
        nc=rows(info)
        M=J(k,k,0)
        dfr=nc-1
		
		// Loop over clusters if G<N:
		if (nc<n){
			for(i=1; i<=nc; i++) {
				xi=panelsubmatrix(X,i,info)
				psii=panelsubmatrix(psi,i,info)
				M=M+(xi'*psii)*(psii'*xi) 
			}
		}
		// Else use heteroskedasticity-robust version:
		else{
			psi2=psi:*psi
			M=quadcross(X,psi2,X)			
		}
		
		Vclust=makesymmetric(scale^2*XphiXinv*M*XphiXinv)
		
		st_matrix("Vclust",Vclust)
		st_numscalar("mata_nclusters",nc)
		
		printf(".")
	}
	
	void _check_vce() {
		real scalar neg
		real matrix diag, EVEC, eval
		
		neg=0
		diag=diagonal(st_matrix("Vc"))
		for (i=1; i<=rows(diag); i++) { 
			if (diag[i]<=0) {
				neg=1
			}
		}
		st_numscalar("negative", neg)
		if (neg==1) {
			// Cameron, Gelbach, and Miller (2011) adjustment from Gu and Yoo (2019, DOI: 10.1177/1536867X19893637) 
			symeigensystem(st_matrix("Vc"), EVEC = ., eval = .)
			eval = eval :* (eval :> 0)
			st_matrix("Vcadj", EVEC*diag(eval)*EVEC')
		}
	}
	
end
	
	
	