#***********************************************************************
# Changes NA's to single precision missing value code
.na.to.snglcode_function(x,mvcode){
      x[is.na(x)]_as.double(mvcode)
  x}
#***********************************************************************
# Changes missing value code to NA
.code.to.na_function(x,mvcode){
      x[x==mvcode]_NA
  x}
#***********************************************************************
#  Perform preliminary manipulations on matrix of continuous data.  
#  Rows are sorted by missing data pattern.
prelim.norm_function(x){
# get dimensions of x
  if(is.vector(x)) x_matrix(x,length(x),1)
  n_nrow(x); p_ncol(x); storage.mode(x)_"double"
# find missingness patterns
  r_1*is.na(x)
  nmis_as.integer(apply(r,2,sum))
  names(nmis)_dimnames(x)[[2]]
# index the missing data patterns
  mdp_as.integer((r%*%(2^((1:ncol(x))-1)))+1)
# do row sort
  ro_order(mdp)
  x_matrix(x[ro,],n,p)
  mdp_mdp[ro]
  r_matrix(r[ro,],n,p)
  ro_order(ro)
# compress missing data patterns
  mdpst_as.integer(seq(along=mdp)[!duplicated(mdp)])
  mdp_unique(mdp); npatt_length(mdpst)
# create r-matrix for display purposes
  r_1-r; r_matrix(r[mdpst,],npatt,p)
  if(npatt==1) tmp_format(n)
  if(npatt>1)  tmp_format(c(mdpst[2:npatt],n+1)-mdpst)
  dimnames(r)_list(tmp,dimnames(x)[[2]])
  storage.mode(r)_"integer"
# center and scale the columns of x
  if(sum(is.na(x))<length(x)){
    mvcode_as.double(max(x[!is.na(x)])+1000)
    x_.na.to.snglcode(x,mvcode)
    tmp_.Fortran("ctrsc",x,n,p,numeric(p),numeric(p),mvcode)
    x_tmp[[1]]; xbar_tmp[[4]]; sdv_tmp[[5]]
    x_.code.to.na(x,mvcode)}
  if(sum(is.na(x))==length(x)){
    xbar_rep(0,p); sdv_rep(1,p)}  
# form matrix of packed storage indices
  d_as.integer((2+3*p+p^2)/2)
  psi_.Fortran("mkpsi",p,matrix(as.integer(0),p+1,p+1))[[2]]
# other bookkeeping quantities
  if(npatt>1) nmdp_as.integer(c(mdpst[-1],n+1)-mdpst)
  if(npatt==1) nmdp_n
  sj_.Fortran("sjn",p,npatt,r,integer(p))[[4]]
  nmon_.Fortran("nmons",p,npatt,r,nmdp,sj,integer(p))[[6]]
  last_.Fortran("lasts",p,npatt,sj,integer(npatt))[[4]]
  tmp_.Fortran("layers",p,sj,integer(p),integer(1))
  layer_tmp[[3]]; nlayer_tmp[[4]]
# return list
  list(x=x,n=n,p=p,r=r,nmis=nmis,ro=ro,mdpst=mdpst,
    nmdp=nmdp,npatt=npatt,xbar=xbar,sdv=sdv,d=d,psi=psi,sj=sj,
    nmon=nmon,last=last,layer=layer,nlayer=nlayer)}
#***********************************************************************
# Retrieves means and covariances from theta. If corr=F, returns
# a list containing a vector of means and a covariance matrix. If
# corr=T, returns a list containing a vector of means, a vector of
# standard deviations, and a correlation matrix.
getparam.norm_function(s,theta,corr=F){
  mu_theta[s$psi[1,2:(s$p+1)]]*s$sdv + s$xbar
  names(mu)_dimnames(s$x)[[2]]
  sigma_theta[s$psi[2:(s$p+1),2:(s$p+1)]]
  sigma_matrix(sigma,s$p,s$p)
  tmp_matrix(s$sdv,s$p,s$p)
  sigma_sigma*tmp*t(tmp)
  dimnames(sigma)_list(names(mu),names(mu))
  if(corr){
    sdv_sqrt(diag(sigma)); names(sdv)_names(mu)
    tmp_matrix(sdv,s$p,s$p)
    r_sigma/(tmp*t(tmp)); dimnames(r)_list(names(mu),names(mu))
    result_list(mu=mu,sdv=sdv,r=r)}
  else result_list(mu=mu,sigma=sigma)
  result}
#***********************************************************************
# Makes a theta vector out of a list of specified parameters.
makeparam.norm_function(s,thetalist){
  result_numeric(s$d); result[1]_-1
  xbar_s$xbar;sdv_s$sdv
  mu_(thetalist[[1]]-xbar)/sdv
  result[2:(s$p+1)]_mu
  if(length(thetalist)==3){
    tmp_matrix(thetalist[[2]],s$p,s$p)
    sigma_thetalist[[3]]*tmp*t(tmp)}
  else sigma_thetalist[[2]]
  tmp_matrix(sdv,s$p,s$p)
  sigma_sigma/(tmp*t(tmp))
  tmp_as.vector(s$psi[2:(s$p+1),2:(s$p+1)])
  result[tmp]_as.vector(sigma)
  result}
#***********************************************************************
# Finds posterior mode of theta under the multivariate
# normal model. If no prior is specified, finds the mle.
em.norm_function(s,start,showits=T,maxits=1000,criterion=.0001,
     prior){
  s$x_.na.to.snglcode(s$x,999)
  if(missing(start)){
    start_.Fortran("stvaln",s$d,numeric(s$d),s$p,s$psi)[[2]]}
  if(missing(prior)){
    mle_as.integer(1)
    tau_numeric(1); m_numeric(1); mu0_numeric(s$p);
    lambdainv_matrix(0,s$p,s$p)}
  if(!(missing(prior))){
    mle_as.integer(0)
    tau_as.numeric(prior[[1]]); m_as.numeric(prior[[2]])
    mu0_as.numeric(prior[[3]])
    lambdainv_as.numeric(prior[[4]])}
  tmp_as.integer(numeric(s$p))
  tobs_.Fortran("tobsn",s$d,numeric(s$d),s$p,s$psi,s$n,s$x,s$npatt,
    s$r,s$mdpst,s$nmdp,tmp)[[2]]
# iterate to mle
  it_0; converged_F
  if(showits) cat(paste("Iterations of EM:","\n"))
  while((!converged)&(it<maxits)){
  old_start
  start_.Fortran("emn",s$d,old,start,tobs,s$p,s$psi,s$n,
    s$x,s$npatt,s$r,s$mdpst,s$nmdp,tmp,tmp,numeric(s$p),
    mle,tau,m,mu0,lambdainv)[[3]]
# print iteration number
  it_it+1; if(showits) cat(paste(format(it),"...",sep=""))
  converged_max(abs(old-start))<=criterion}
  if(showits)cat("\n")
  start}       
#***********************************************************************
# Calculates log observed-data posterior at theta
logpost.norm_function(s,theta,prior){
  s$x_.na.to.snglcode(s$x,999)
  l1_.Fortran("lobsn",s$d,theta,numeric(s$d),s$p,s$psi,s$n,s$x,
    s$npatt,s$r,s$mdpst,s$nmdp,as.integer(numeric(s$p)),numeric(s$p),
    0)[[14]]
  if(!(missing(prior))){
  tau_as.numeric(prior[[1]]); m_as.numeric(prior[[2]])
  mu0_as.numeric(prior[[3]]); lambdainv_as.numeric(prior[[4]])}
  if(missing(prior)){
    tau_as.numeric(0); m_as.numeric(-1); mu0_numeric(s$p)
    lambdainv_matrix(0,s$p,s$p)}
  l2_.Fortran("lprin",s$d,theta,s$p,s$psi,numeric(s$p),tau,m,mu0,
    lambdainv,0)[[10]]
  l1+l2}
#***********************************************************************
# Calculates observed-data loglikelihood at theta
loglik.norm_function(s,theta){
  s$x_.na.to.snglcode(s$x,999)
  .Fortran("lobsn",s$d,theta,numeric(s$d),s$p,s$psi,s$n,s$x,s$npatt,
    s$r,s$mdpst,s$nmdp,as.integer(numeric(s$p)),numeric(s$p),0)[[14]]}
#***********************************************************************
# Simulate a value of theta from a normal-inverted Wishart
# distribution
ninvwish_function(s,params){
    tau_as.numeric(params[[1]]); m_as.numeric(params[[2]])
    mu0_params[[3]]; lambdainv_params[[4]]
    tmpi_as.integer(numeric(s$p)); tmpr_as.double((numeric(s$p)))
    pri_numeric(s$d); pri[2:(s$p+1)]_mu0
    tmp_as.vector(s$psi[2:(s$p+1),2:(s$p+1)])
    pri[tmp]_as.vector(lambdainv)
    .Fortran("ninvwn",s$d,pri,tau,m,s$p,s$psi,
       numeric((s$p)^2),tmpr,numeric(s$d),tmpi)[[2]]}
#***********************************************************************
# Data augmentation for the multivariate normal.
# Produces a new draw of theta from its posterior distribution. 
da.norm_function(s,start,prior,steps=1,showits=F,return.ymis=F){
  s$x_.na.to.snglcode(s$x,999)
  tmpi_as.integer(numeric(s$p)); tmpr_as.double((numeric(s$p)))
  tobs_.Fortran("tobsn",s$d,numeric(s$d),s$p,s$psi,s$n,s$x,s$npatt,
    s$r,s$mdpst,s$nmdp,tmpi)[[2]]
  if(missing(prior)){
    tau_0; m_-1; mu0_numeric(s$p);lambdainv_matrix(0,s$p,s$p)}
  if(!(missing(prior))){
    tau_as.numeric(prior[[1]]); m_as.numeric(prior[[2]])
    mu0_prior[[3]]; lambdainv_prior[[4]]}
  pri_numeric(s$d); pri[2:(s$p+1)]_mu0
  tmp_as.vector(s$psi[2:(s$p+1),2:(s$p+1)])
  pri[tmp]_as.vector(lambdainv)
  if(showits) cat(paste("Steps of Data Augmentation:","\n"))
  for(i in 1:steps){
    if(showits) cat(paste(format(i),"...",sep=""))
    tmp_.Fortran("is1n",s$d,start,start,tobs,s$p,s$psi,s$n,s$x,
      s$npatt,s$r,s$mdpst,s$nmdp,tmpi,tmpi,tmpr,start)
    start_tmp[[3]]
    start_.Fortran("ps1n",s$d,start,m,tau,pri,s$p,s$psi,s$n,
      matrix(0,s$p,s$p),tmpr,start,tmpi)[[5]]}
  if(showits)cat("\n")
  if(return.ymis){
    ymis_tmp[[8]]*matrix(s$sdv,s$n,s$p,T)+matrix(s$xbar,s$n,s$p,T)
    ymis_ymis[s$ro,];ymis_ymis[s$x[s$ro,]==999]
    start_list(parameter=start,ymis=ymis)}
  start}
#***********************************************************************
# Generates a single imputed dataset under theta
imp.norm_function(s,theta,x){
  s$x_.na.to.snglcode(s$x,999)
  tmpi_as.integer(numeric(s$p)); tmpr_as.double((numeric(s$p)))
  tobs_.Fortran("tobsn",s$d,numeric(s$d),s$p,s$psi,s$n,s$x,s$npatt,
    s$r,s$mdpst,s$nmdp,tmpi)[[2]]
  s$x_.Fortran("is1n",s$d,theta,theta,tobs,s$p,s$psi,s$n,s$x,
      s$npatt,s$r,s$mdpst,s$nmdp,tmpi,tmpi,tmpr,theta)[[8]]
  s$x_s$x*matrix(s$sdv,s$n,s$p,T)+matrix(s$xbar,s$n,s$p,T)
  s$x_s$x[s$ro,]
  if(!missing(x))x[is.na(x)]_s$x[is.na(x)]
  else{x_s$x; storage.mode(x)_"double"}
  x}
#***********************************************************************
# Monotone data augmentation for the multivariate normal.
# Produces a new draw of theta from its posterior distribution. 
mda.norm_function(s,theta,steps=1,showits=F){
  s$x_.na.to.snglcode(s$x,999)
  tobs_.Fortran("tobsmn",s$p,s$psi,s$n,s$x,s$npatt,s$r,s$mdpst,s$nmdp,
    s$last,integer(s$p),s$sj,s$layer,s$nlayer,s$d,
    matrix(0,s$nlayer,s$d))[[15]]
  if(showits) cat(paste("Steps of Monotone Data Augmentation:",
    "\n")) 
  for(i in 1:steps){
    if(showits) cat(paste(format(i),"...",sep=""))
    s$x_.Fortran("is2n",s$d,theta,s$p,s$psi,s$n,s$x,s$npatt,
      s$r,s$mdpst,s$nmdp,s$sj,s$last,integer(s$p),integer(s$p),
      double(s$p),theta)[[6]]
    theta_.Fortran("ps2n",s$p,s$psi,s$n,s$x,s$npatt,s$r,s$mdpst,
      s$nmdp,integer(s$p),integer(s$p),s$nmon,s$sj,s$nlayer,s$d,
      tobs,numeric(s$d),numeric(s$d),numeric(s$p+1),numeric(s$d))[[19]]}
  if(showits)cat("\n")
  theta}
#***********************************************************************
# Initializes random number generator seed. Argument should be a 
# positive integer
rngseed_function(seed){
  seed_as.integer(seed)
  if(seed<=0)stop("Seed must be a positive integer")
  tmp_.Fortran("rngs",seed)
  invisible()}
#***********************************************************************
# multiple imputation inference
mi.inference_function(est,std.err,confidence=.95){
  qstar_est[[1]]
  for(i in 2:length(est)){qstar_cbind(qstar,est[[i]])}
  qbar_apply(qstar,1,mean)
  u_std.err[[1]]
  for(i in 2:length(std.err)){u_cbind(u,std.err[[i]])}
  dimnames(u)[[1]]_dimnames(qstar)[[1]]
  u_u^2
  ubar_apply(u,1,mean)
  bm_apply(qstar,1,var)
  m_dim(qstar)[2]
  tm_ubar+((1+(1/m))*bm)
  rem_(1+(1/m))*bm/ubar
  nu_(m-1)*(1+(1/rem))**2
  alpha_1-(1-confidence)/2
  low_qbar-qt(alpha,nu)*sqrt(tm)
  up_qbar+qt(alpha,nu)*sqrt(tm)
  pval_2*(1-pt(abs(qbar/sqrt(tm)),nu))
  fminf_(rem+2/(nu+3))/(rem+1)
  result_list(est=qbar,std.err=sqrt(tm),df=nu,signif=pval,lower=low,
  upper=up,r=rem,fminf=fminf)
  result}
#***********************************************************************

# Added function:

".First.lib" <-
function(lib, pkg) {
   library.dynam("norm", pkg, lib) }
