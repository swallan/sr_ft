ptukey_function(q,nmeans,df,nranges=1)
# J.H.Maindonald, DSIR Physical Sciences, Applied Mathematics Group, 1990
# Calls the Fortran code that is described in
# Copenhaver, Margaret DiPonzio and Holladn, Burt S: Computation of
#   the Distribution of the Maximum Studentized Range Statistic with
#   Application to Multiple Significance Testing of Simple Effects.
#   Journal of Statistical Computation and Simulation 30, pp.1-15, 1988.
{
	nval_max(c(length(q),length(nmeans),length(df),length(nranges)))
	if(nval>1){
	if(length(q)==1)q_rep(q,nval)
	
	if(length(nmeans)==1)nmeans_rep(nmeans,nval)
		if(length(df)==1)df_rep(df,nval)
	if(length(nranges)==1)nranges_rep(nranges,nval)
	}
	pval_rep(0,nval)
	for(i in seq(along=q)){
	
    ir_c(0,0);p_0
    junk_.Fortran("qprob2",as.double(q[i]),as.double(nranges[i]),
	as.double(nmeans[i]),as.double(df[i]),as.integer(ir),p=as.double(p))
    pval[i]_junk$p
}
	pval
	}

qtukey_function(p,nmeans,df,nranges=1)
# J.H.Maindonald, DSIR Physical Sciences, Applied Mathematics Group, 1990
# Calls the Fortran code that is described in
# Copenhaver, Margaret DiPonzio and Holladn, Burt S: Computation of
#   the Distribution of the Maximum Studentized Range Statistic with
#   Application to Multiple Significance Testing of Simple Effects.
#   Journal of Statistical Computation and Simulation 30, pp.1-15, 1988.
{
	nval_max(c(length(p),length(nmeans),length(df),length(nranges)))
	if(nval>1){
	if(length(p)==1)p_rep(p,nval)
	
	if(length(nmeans)==1)nmeans_rep(nmeans,nval)
		if(length(df)==1)df_rep(df,nval)
	if(length(nranges)==1)nranges_rep(nranges,nval)
	}
	qval_rep(0,nval)
	for(i in seq(along=p)){
    ir_c(0,0,0); q_0
    junk_.Fortran("cv2",as.double(p[i]),as.double(nranges[i]),
	as.double(nmeans[i]),as.double(df[i]),as.integer(ir),q=as.double(q))
    qval[i]_junk$q
	}
	qval
}

