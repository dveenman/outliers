*! version 2.0.1 20221109 David Veenman

/* 
20221109: 2.0.1     Some housekeeping: changed scalars to locals
20221106: 2.0.0     Complete new version based on fast bootstrap for MM estimators following Salibian-Barrera and Zamar (2002) [SZ2002]
                    Option for cluster-bootstrap up to two dimensions by adjusting bootstrap in SZ2002 to pairs-cluster bootstrap
                    Including computationally efficient procedure from MacKinnon (2022) by drawing lower-dimensional matrices in cluster bootstrap
20220106: 1.1.1     Added finite-sample correction and significance testing based on t(G-1)
20211011: 1.1.0     Dropped capture to allow for hard exit from loop
                    Added tolerance option for initial estimates of coefficients
                    Added option nohaus to s and mm estimator, and nor2 to all, for faster execution
                    Small fix in scalar drop at bottom; "drop _all" causes all scalars to be dropped outside the program too
20210701: 1.0.0     First version
*/

program define roboot, eclass sortpreserve
	syntax varlist(numeric fv) [in] [if], [cluster(varlist)] nboot(integer) eff(real) [tol(real 0)] [sopts(str)] [weightvar(str)] [seed(integer 0)]
	
	if (`seed'==0){
		di ""
		di as text "Careful: make sure to use seed() option or the 'set seed' command in Stata to obtain reproducable results!"
	}
	else{
		set seed `seed'
	}
	
	if (`tol'==0){
	    local tolerance=1e-6
	}
	else {
	    local tolerance=`tol'
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
	
	local ncdim: word count `cluster'
	if (`ncdim'>2){
	    di as text "ERROR: Maximum number of dimensions to cluster on is two"
		exit
	}
	if (`ncdim'>0){
		local clusterdim1: word 1 of `cluster'
	}
	if (`ncdim'>1){
		local clusterdim2: word 2 of `cluster'
	}
	
	/* Create variable for intersection between two clustering dimensions: */
	if (`ncdim'==2){	    
		tempvar intersection
		qui egen `intersection'=group(`clusterdim1' `clusterdim2') if `touse'
	}	
	
	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	// Step 0: Obtain coefficients, scale, and residuals from robust estimation:
	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	di ""
	di as text "STEP 1: Obtaining robust regression estimates, residuals, and scale estimate..."
	di ""
	capture robreg mm `varlist' if `touse', eff(`eff') nohaus tol(`tolerance') nose sopts(`sopts')
	
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
	if (`ncdim'==0){
		local e_df_r=e(df_r)
	}
	
	if "`fvcheck'"=="true"{
		local param=(e(df_m)+2)
	}
	else{
		local param=(e(df_m)+1)
	}
	matrix bmm0=(e(b)[1,1..`param'])'
	
	matrix bmm=J(`param',1,0)
	forvalues i=1(1)`param'{
		matrix bmm[`i', 1]=bmm0[`i', 1]
	}
	matrix bs=(e(b)[1,(`param'+1)..(2*`param')])'
	
	scalar scale=e(scale)
	scalar kmm=e(k)
	scalar ks=e(kS)

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
	
	local indepvnames "`indepv' _cons"
	if (`ncdim'>=1){	
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		// Step 1: Getting cluster-robust SEs for first dimension  
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		di as text "STEP 2: Obtaining cluster-robust bootstrapped standard errors..."
		if (`ncdim'==2){
			di "STEP 2A: First clustering dimension..."
		}
		tempvar clusterid
		egen `clusterid'=group(`clusterdim1') if `touse'
		sort `clusterid' 
		qui by `clusterid': gen _cvarn_temp=_n if `touse'
		qui replace _cvarn_temp=. if _cvarn_temp>1		
		local cvar "`clusterid'"	
		local cvarn "`clusterid' _cvarn_temp"
		mata: _vce_cluster(`nboot')
		local nclusterdim1=mata_nclusters
		local e_df_r1=mata_nclusters-1
		matrix colnames Vclust=`indepvnames'
		matrix rownames Vclust=`indepvnames'	
		matrix V1=Vclust
		drop _cvarn_temp
	}
	
	if (`ncdim'==2){		
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		// Step 2: Getting cluster-robust SEs for second dimension 
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		di "STEP 2B: Second clustering dimension..."
		tempvar clusterid
		egen `clusterid'=group(`clusterdim2') if `touse'
		sort `clusterid' 
		qui by `clusterid': gen _cvarn_temp=_n if `touse'
		qui replace _cvarn_temp=. if _cvarn_temp>1
		local cvar "`clusterid'"	
		local cvarn "`clusterid' _cvarn_temp"
		mata: _vce_cluster(`nboot')
		local nclusterdim2=mata_nclusters
		local e_df_r2=mata_nclusters-1
		matrix colnames Vclust=`indepvnames'
		matrix rownames Vclust=`indepvnames'	
		matrix V2=Vclust
		drop _cvarn_temp
			
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		// Step 3: Getting cluster-robust SEs for second dimension 
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		di "STEP 2C: Intersection of the two clustering dimensions..."
		tempvar clusterid
		egen `clusterid'=group(`intersection') if `touse'
		sort `clusterid' 
		qui by `clusterid': gen _cvarn_temp=_n if `touse'
		qui replace _cvarn_temp=. if _cvarn_temp>1
		local cvar "`clusterid'"	
		local cvarn "`clusterid' _cvarn_temp"
		mata: _vce_cluster_inter(`nboot')
		matrix colnames Vclust=`indepvnames'
		matrix rownames Vclust=`indepvnames'	
		matrix V3=Vclust
		drop _cvarn_temp
			
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		// Step 4: Combining VCEs and applying small sample correction: 
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		matrix Vc=V1+V2-V3
		matrix colnames Vc=`indepvnames'
		matrix rownames Vc=`indepvnames'
			
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
		if `e_df_r1'<`e_df_r2'{
			local e_df_r=`e_df_r1'
		}
		else{
			local e_df_r=`e_df_r2'
		}
	}
	else{
		if (`ncdim'==1) {
			matrix Vc=V1
			local indepvnames "`indepv' _cons"
			matrix colnames Vc=`indepvnames'
			matrix rownames Vc=`indepvnames'
				
			local N=`e_N'
			if "`fvcheck'"=="true"{
				local K=rowsof(Vc)-1
			}
			else{
				local K=rowsof(Vc) 
			}

			local factor=(`nclusterdim1'/(`nclusterdim1'-1))*((`N'-1)/(`N'-`K'))
			matrix Vc=`factor'*Vc
			local e_df_r=`e_df_r1'			
		}		
		if (`ncdim'==0) {
			di "STEP 2: Estimating bootstrapped standard errors..."
			mata: _vce_nocl(`nboot')
			matrix colnames Vnoclust=`indepvnames'
			matrix rownames Vnoclust=`indepvnames'	
			matrix Vc=Vnoclust
		}
	}
		
	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	// Step 5: Post resuls in e() and print results:
	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	ereturn clear
	tempname b V

	matrix rownames bmm=`indepvnames'
	matrix colnames Vc=`indepvnames'
	matrix rownames Vc=`indepvnames'	
	matrix `b'=bmm'
	matrix `V'=Vc
	
	ereturn post `b' `V'
	
	ereturn scalar N=`e_N'
	ereturn scalar r2_p=`e_r2_p'
	ereturn local depvar "`depv'"
	ereturn scalar df_r=`e_df_r'
		
	di " "
	di " "
	if (`ncdim'==0){
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
	if (`ncdim'==1){
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
	if (`ncdim'==2){
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
	if (`ncdim'==0){
	    di _column(56) in green "Pseudo R2 =     " %7.4f in yellow e(r2_p)
	}
	if (`ncdim'>0){
		di in green "Number of clusters (`clusterdim1') = " _column(31) %5.0f in yellow `nclusterdim1' ///
			_column(56) in green "Pseudo R2 =     " %7.4f in yellow e(r2_p)
	}
	if (`ncdim'>1){
		di in green "Number of clusters (`clusterdim2') = " _column(31) %5.0f in yellow `nclusterdim2' 
	}
			
	ereturn display
	if (`ncdim'==0){
	    di "SE not adjusted for clustering" 
	}
	if (`ncdim'==1){
	    di "SE clustered by " "`clusterdim1'" 
	}
	if (`ncdim'==2){
	    di "SE clustered by " "`clusterdim1'" " and " "`clusterdim2'" 
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
	void _vce_cluster(real scalar B) {
		st_view(y=., ., st_local("depv"), st_local("touse"))
		st_view(X=., ., tokens(st_local("indepv")), st_local("touse"))
		st_view(cvar=., ., tokens(st_local("cvar")), st_local("touse"))
		st_view(cvarn=., ., tokens(st_local("cvarn")), st_local("touse"))
		bmm=st_matrix("bmm")
		bs=st_matrix("bs")
		s=st_numscalar("scale")		
		kmm=st_numscalar("kmm")		
		ks=st_numscalar("ks")
		X=(X,J(rows(X),1,1))

		n1=nonmissing(y)
		n2=nonmissing(rowsum(X, 1))
		if (n1<n2) n=n1 
		else n=n2
		k=cols(X)
		
		// Create residual vectors for MM and S estimator:
		rmm=y-X*bmm
		rs=y-X*bs
		zmm=rmm:/s
		zs=rs:/s
		
		// Compute the constant needed for minimization of the M-scale (K=b in eq 2.4 in Salibian-Barrera and Zamar (2002)):
		K=.5*mm_biweight_rho(ks,ks)
		
		// Compute weights w and v to be bootstrapped (no need to recalculate each bootstrap):
		// Note: w should be scaled by standardized residual, not raw residual in eq 3.1:
		w=mm_biweight_psi(zmm, kmm):/rmm
		v=(s/(n*K))*((mm_biweight_rho(zs, ks)):/rs)
		
		// Compute correction factors per eq 3.6-3.8 in Salibian-Barrera and Zamar (2002):
		phi1=mm_biweight_phi(zmm, kmm)
		M=s*invsym(cross(X,phi1,X))*cross(X,w,X)
		psi0=mm_biweight_psi(zs, ks)
		aa=(psi0:*rs)/s
		a=(1/(n*K))*colsum(aa)
		d=(1/a)*invsym(cross(X,phi1,X))*cross(X,phi1,rmm)

		// Cluster setup:
        info=panelsetup(cvar, 1)
        nc=rows(info)

		// Obtain the relevant information from the clusters and store in lower-dimensional matrices:
		nb=0
		XXg=J(0,k,0)
		Xyg=J(0,1,0)
		vrg=J(0,1,0)
		for(i=1; i<=nc; i++) {
			// Cluster-level data needed for coefficients:
			xg=panelsubmatrix(X,i,info)
			ng=nonmissing(rowsum(xg, 1))
			wg=panelsubmatrix(w,i,info)
			XXg=(XXg \ cross(xg,wg,xg))
			nb=nb+ng
			yg=panelsubmatrix(y,i,info)
			Xyg=(Xyg \ cross(xg,wg,yg))
			// Cluster-level data needed for scale estimate:
			vg=panelsubmatrix(v,i,info)
			rsg=panelsubmatrix(rs,i,info)
			vrg=(vrg \ colsum(vg:*rsg))
		}
		
		// Perform bootstrapping by drawing the lower-dimensional matrices for calculation of coefficients and scale:
		beta=J(0,k,0)
		bresults=J(0,k,.)
		for(b=1; b<=B; b++) {
			XXb=J(k,k,0)
			Xyb=J(k,1,0)
			vrb=J(1,1,0)
			submatx=J(k,k,0)
			submaty=J(n,1,0)
			bsample_c=mm_sample(nc,nc)
			for(i=1; i<=nc; i++) {
				cluster=bsample_c[i]
				pos0=1+k*(cluster-1)
				pos=(pos0::pos0+k-1)
				submatx=XXg[pos,.]
				XXb=XXb+submatx
				submaty=Xyg[pos,.]
				Xyb=Xyb+submaty
				vrb=vrb+vrg[i]
			}
			betabmm=cross(invsym(XXb),Xyb)
			sb=vrb
			
			// Apply correction and store results:
			bbcenter=(M*(betabmm-bmm)+d*(sb-s))'
			bresults=(bresults \ bbcenter)
		}
				
        Vclust=variance(bresults)
		st_matrix("Vclust",Vclust)
		st_numscalar("mata_nclusters",nc)
	}
	
	void _vce_cluster_inter(real scalar B) {
		st_view(y=., ., st_local("depv"), st_local("touse"))
		st_view(X=., ., tokens(st_local("indepv")), st_local("touse"))
		st_view(cvar=., ., tokens(st_local("cvar")), st_local("touse"))
		st_view(cvarn=., ., tokens(st_local("cvarn")), st_local("touse"))
		bmm=st_matrix("bmm")
		bs=st_matrix("bs")
		s=st_numscalar("scale")		
		kmm=st_numscalar("kmm")		
		ks=st_numscalar("ks")
		X=(X,J(rows(X),1,1))

		n1=nonmissing(y)
		n2=nonmissing(rowsum(X, 1))
		if (n1<n2) n=n1 
		else n=n2
		k=cols(X)

		// Create residual vectors for MM and S estimator:
		rmm=y-X*bmm
		rs=y-X*bs
		zmm=rmm:/s
		zs=rs:/s
		
		// Compute the constant needed for minimization of the M-scale (K=b in eq 2.4 in Salibian-Barrera and Zamar (2002)):
		K=.5*mm_biweight_rho(ks,ks)
		
		// Compute weights w and v to be bootstrapped (no need to recalculate each bootstrap):
		// Per eq 3.1 in Salibian-Barrera and Zamar (2002):
		w=mm_biweight_psi(zmm, kmm):/rmm
		v=(s/(n*K))*((mm_biweight_rho(zs, ks)):/rs)
		
		// Compute correction factors per eq 3.6-3.8 in Salibian-Barrera and Zamar (2002):
		phi1=mm_biweight_phi(zmm, kmm)
		M=s*invsym(cross(X,phi1,X))*cross(X,w,X)
		psi0=mm_biweight_psi(zs, ks)
		aa=(psi0:*rs)/s
		a=(1/(n*K))*colsum(aa)
		d=(1/a)*invsym(cross(X,phi1,X))*cross(X,phi1,rmm)

		// Cluster setup:
        info=panelsetup(cvar, 1)
        nc=rows(info)
		mm_panels(cvar, Cinfo=.)
				
		bresults=J(0,k,.)
		for(b=1; b<=B; b++) {
			// Create permutation vector and draw the relevant variables for bootstrap sample:
			bsample_n=mm_sample(nc,.,Cinfo)
				
			yb=y[bsample_n]
			Xb=X[bsample_n,.]		
			wb=w[bsample_n]
			vb=v[bsample_n]
			rbs=rs[bsample_n]
			
			// Compute bootstrap-sample betas and scale:
			betabmm=invsym(cross(Xb,wb,Xb))*cross(Xb,wb,yb)
			sb=colsum(vb:*rbs)

			// Apply correction and store results:
			bbcenter=(M*(betabmm-bmm)+d*(sb-s))'
			bresults=(bresults \ bbcenter)
		}
				
        Vclust=variance(bresults)
		st_matrix("Vclust",Vclust)
	}

	void _vce_nocl(real scalar B) {
		st_view(y=., ., st_local("depv"), st_local("touse"))
		st_view(X=., ., tokens(st_local("indepv")), st_local("touse"))
		bmm=st_matrix("bmm")
		bs=st_matrix("bs")
		s=st_numscalar("scale")		
		kmm=st_numscalar("kmm")		
		ks=st_numscalar("ks")
		X=(X,J(rows(X),1,1))

		n1=nonmissing(y)
		n2=nonmissing(rowsum(X, 1))
		if (n1<n2) n=n1 
		else n=n2
		k=cols(X)

		// Create residual vectors for MM and S estimator:
		rmm=y-X*bmm
		rs=y-X*bs
		zmm=rmm:/s
		zs=rs:/s
		
		// Compute the constant needed for minimization of the M-scale (K=b in eq 2.4 in Salibian-Barrera and Zamar (2002)):
		K=.5*mm_biweight_rho(ks,ks)
		
		// Compute weights w and v to be bootstrapped (no need to recalculate each bootstrap):
		// Per eq 3.1 in Salibian-Barrera and Zamar (2002):
		w=mm_biweight_psi(zmm, kmm):/rmm
		v=(s/(n*K))*((mm_biweight_rho(zs, ks)):/rs)
		
		// Compute correction factors per eq 3.6-3.8 in Salibian-Barrera and Zamar (2002):
		phi1=mm_biweight_phi(zmm, kmm)
		M=s*invsym(cross(X,phi1,X))*cross(X,w,X)
		psi0=mm_biweight_psi(zs, ks)
		aa=(psi0:*rs)/s
		a=(1/(n*K))*colsum(aa)
		d=(1/a)*invsym(cross(X,phi1,X))*cross(X,phi1,rmm)

		bresults=J(0,k,.)
		for(b=1; b<=B; b++) {
			// Create permutation vector and draw the relevant variables for bootstrap sample:
			bsample_n=mm_sample(n,n)
				
			yb=y[bsample_n]
			Xb=X[bsample_n,.]		
			wb=w[bsample_n]
			vb=v[bsample_n]
			rbs=rs[bsample_n]
			
			// Compute bootstrap-sample betas and scale:
			betabmm=invsym(cross(Xb,wb,Xb))*cross(Xb,wb,yb)
			sb=colsum(vb:*rbs)

			// Apply correction and store results:
			bbcenter=(M*(betabmm-bmm)+d*(sb-s))'
			bresults=(bresults \ bbcenter)
		}
				
        Vnoclust=variance(bresults)
		st_matrix("Vnoclust",Vnoclust)
	}
	
end
	