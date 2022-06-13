library(MASS)
library(glmnet)
llog <- function (z, eps=1/50)
{
  z[is.na(z)] <- eps
  ans <- z
  lo <- (z < eps)
  ans[lo] <- log(eps) - 1.5 + 2 * z[lo]/eps - 0.5 * (z[lo]/eps)^2
  ans[!lo] <- log(z[!lo])
  ans
}

llogp <- function (z, eps=1/50)
{
  z[is.na(z)] <- eps
  ans <- z
  lo <- (z < eps)
  ans[lo] <- 2/eps - z[lo]/eps^2
  ans[!lo] <- 1/z[!lo]
  ans
}

llogpp <- function (z, eps)
{
  ans <- z
  lo <- (z < eps)
  ans[lo] <- -1/eps^2
  ans[!lo] <- -1/z[!lo]^2
  ans
}
logelr.ee <- function (mx, lam)
{
  mx <- as.matrix(mx)
  n <- nrow(mx)
  arg <- 1 + mx %*% lam
  return(-sum(llog(arg, 1/n)))
}

m=function(theta,xy) {
  x=xy[,-1]
  y=xy[,1]
  x*drop(y-x%*%theta)
}
scad=function(mu, scadlambda, Dim=length(mu), scada=3.7){
  mu=abs(mu)
  tmp=mu<scadlambda
  s=rep(0,Dim)
  s[tmp]=scadlambda*mu[tmp]
  tmp=(mu>scadlambda)&(mu<=scada*scadlambda)
  s[tmp]=-(mu[tmp]*mu[tmp]-2*scada*scadlambda*mu[tmp]+
             scadlambda*scadlambda)/2/(scada-1)
  tmp=mu>=scada*scadlambda
  s[tmp]=(scada+1)*scadlambda*scadlambda/2
  s
}

scad.grad=function(mu, scadlambda, Dim=length(mu), scada=3.7){
  mu=abs(mu)
  t2=ifelse((mu<=scadlambda),1,0)
  t3=(scada*scadlambda-mu)
  t3=ifelse(t3>=0,t3,0)
  t3=t3/((scada-1)*scadlambda)
  t3=ifelse(mu>scadlambda,t3,0)
  ss=scadlambda*(t2+t3)
  ss[mu==0]=0
  ss
}


el.ee.tang <- function (mx, lam, maxit=25, gradtol=1e-7, svdtol=1e-9,
                        itertrace=F, TINY=sqrt(.Machine$double.xmin))
{
  mx <- as.matrix(mx)
  n <- nrow(mx)
  r <- ncol(mx)
  
  #if (n <= r)
  #stop("Need more observations than length(mu) in el.test().")
  z <- mx
  scale <- mean(abs(z)) + TINY
  z <- z/scale
  if (!missing(lam)) {
    lam <- as.vector(lam)
    lam <- lam * scale
    if (logelr.ee(z, lam) > 0)
      lam <- rep(0, r)
  }
  if (missing(lam))
    lam <- rep(0, r)
  if (svdtol < TINY)
    svdtol <- TINY
  if (gradtol < TINY)
    gradtol <- TINY
  nwts <- c(3^-c(0:3), rep(0, 12))
  gwts <- 2^(-c(0:(length(nwts) - 1)))
  gwts <- (gwts^2 - nwts^2)^0.5
  gwts[12:16] <- gwts[12:16] * 10^-c(1:5)
  nits <- 0
  gsize <- gradtol + 1
  while (nits < maxit && gsize > gradtol) {
    arg <- 1 + z %*% lam
    wts1 <- as.vector(llogp(arg, 1/n))
    wts2 <- as.vector(-llogpp(arg, 1/n))^0.5
    grad <- as.matrix(-z * wts1)
    grad <- as.vector(rowsum(grad, rep(1, nrow(grad))))
    gsize <- mean(abs(grad))
    hess <- z * wts2
    svdh <- svd(hess)
    if (min(svdh$d) < max(svdh$d) * svdtol)
      svdh$d <- svdh$d + max(svdh$d) * svdtol
    nstep <- svdh$v %*% (t(svdh$u)/svdh$d)
    nstep <- as.vector(nstep %*% matrix(wts1/wts2, n, 1))
    gstep <- -grad
    if (sum(nstep^2) < sum(gstep^2))
      gstep <- gstep * (sum(nstep^2)^0.5/sum(gstep^2)^0.5)
    ologelr <- -sum(llog(arg, 1/n))
    ninner <- 0
    for (i in 1:length(nwts)) {
      nlogelr <- logelr.ee(z, lam + nwts[i] * nstep + gwts[i] * gstep)
      if (nlogelr < ologelr) {
        lam <- lam + nwts[i] * nstep + gwts[i] * gstep
        ninner <- i
        break
      }
    }
    nits <- nits + 1
    #if (ninner == 0)
    if (ninner == 0 | abs(nlogelr-ologelr)<1e-4)
      nits <- maxit
    if (itertrace)
      #print(c(lam, nlogelr, gsize, ninner))
      print(c(ologelr, nlogelr, ninner, sum(wts1*z),sum(1/wts1)))
  }
  list(`-2LLR` = -2 * nlogelr, Pval = 1 - pchisq(-2 * nlogelr, df = r),
       lambda = lam/scale, grad = grad * scale, hess = t(hess) %*%
         hess * scale^2, wts = wts1, nits = nits)
}

el.ee.penalized.svd <- function (mx, lam, nu, maxit=25, gradtol=1e-7, svdtol=1e-9, thres.v=0.05,
                                 itertrace=F, TINY=sqrt(.Machine$double.xmin))
{
  mx <- as.matrix(mx)
  n <- nrow(mx)
  r <- ncol(mx)
  
  #if (n <= r)
  #stop("Need more observations than length(mu) in el.test().")
  z <- mx
  scale <- mean(abs(z)) + TINY
  z <- z/scale
  if (!missing(lam)) {
    lam <- as.vector(lam)
    lam <- lam * scale
    if (logelr.ee(z, lam) > 0)
      lam <- rep(0, r)
  }
  if (missing(lam))
    lam <- rep(0, r)
  if (svdtol < TINY)
    svdtol <- TINY
  if (gradtol < TINY)
    gradtol <- TINY
  nwts <- c(3^-c(0:3), rep(0, 12))
  gwts <- 2^(-c(0:(length(nwts) - 1)))
  gwts <- (gwts^2 - nwts^2)^0.5
  gwts[12:16] <- gwts[12:16] * 10^-c(1:5)
  nits <- 0
  gsize <- gradtol + 1
  while (nits < maxit && gsize > gradtol) {
    arg <- 1 + z %*% lam
    wts1 <- as.vector(llogp(arg, 1/n))
    wts2 <- as.vector(-llogpp(arg, 1/n))^0.5
    grad <- as.matrix(-z * wts1)
    grad <- as.vector(rowsum(grad, rep(1, nrow(grad))))
    gsize <- mean(abs(grad))
    hess <- z * wts2
    svdh <- svd(hess)
    if (min(svdh$d) < max(svdh$d) * svdtol)
      svdh$d <- svdh$d + max(svdh$d) * svdtol
    
    # truncation
    ind=which(abs(svdh$d)<=nu)
    if (length(ind)>0){
      nstep <- svdh$v[,-ind] %*% (t(svdh$u[,-ind])/svdh$d[-ind])
    }else{
      nstep <- svdh$v %*% (t(svdh$u)/svdh$d)
    }
    
    nstep <- as.vector(nstep %*% matrix(wts1/wts2, n, 1))
    gstep <- -grad
    if (sum(nstep^2) < sum(gstep^2))
      gstep <- gstep * (sum(nstep^2)^0.5/sum(gstep^2)^0.5)
    ologelr <- -sum(llog(arg, 1/n))
    ninner <- 0
    for (i in 1:length(nwts)) {
      nlogelr <- logelr.ee(z, lam + nwts[i] * nstep + gwts[i] * gstep)
      if (nlogelr <= ologelr) {
        lam <- lam + nwts[i] * nstep + gwts[i] * gstep
        zero=abs(lam)<thres.v
        lam[zero]=0
        
        ninner <- i
        break
      }
    }
    zero=abs(lam)<thres.v
    lam[zero]=0
    
    nits <- nits + 1
    if (ninner == 0 | abs(nlogelr-ologelr)<1e-4)
      nits <- maxit
    if (itertrace)
      #print(c(lam, nlogelr, gsize, ninner))
      print(c(ologelr, nlogelr, ninner, sum(wts1*z),sum(1/wts1)))
  }
  zero=abs(lam)<thres.v
  lam[zero]=0
  
  list(`-2LLR` = -2 * nlogelr, Pval = 1 - pchisq(-2 * nlogelr, df = r),
       lambda = lam/scale, grad = grad * scale, hess = t(hess) %*%
         hess * scale^2, wts = wts1, nits = nits)
}



elhood.ee.lqa.2penalty=function(m,x,init.theta,scadlambda,nu,thres.v,maxiter=30){
  n=nrow(x)
  p=length(init.theta)
  mx=m(init.theta,x)
  r=ncol(mx)
  
  init.lambda=rep(0, r)
  
  elhood=function(theta){
    
    theta.all=rep(0,p)
    theta.all[!zero]=theta
    mx=m(theta.all,x)
    lambda=el.ee.penalized.svd(mx,init.lambda,nu)$lambda
    z=mx
    arg=1+z%*%lambda
    nlogelr=sum(llog(arg,1/n))
    
    tmp=scad.grad(init.theta, scadlambda)/abs(init.theta)  
    vf=nlogelr+sum(tmp*theta*theta)*n/2
    
    return(vf)
  }
  
  counter=1
  init.theta.all=init.theta
  zero=abs(init.theta)<thres.v
  if(sum(zero)==p) zero[2]=0
  init.theta=init.theta[!zero]
  init.theta.all[zero]=0
  
  repeat {
    zz=nlm(elhood, init.theta, print.level=0,iterlim=maxiter,
           check.analyticals=F)
    dff=zz$estimate-init.theta
    init.theta=zz$estimate
    
    init.theta.all[!zero]=init.theta
    
    zero1=abs(init.theta)<thres.v
    init.theta=init.theta[!zero1]
    
    zero=abs(init.theta.all)<thres.v
    init.theta.all[zero]=0
    
    mx=m(init.theta.all,x)
    lambda=el.ee.penalized.svd(mx,init.lambda,nu)$lambda
    init.lambda=lambda
    
    if(sum(zero)>p-1 | max(abs(dff))<1e-4)
      break
    
    counter=counter+1
    if (counter >= 20) {
      warning("iteration fails to converge")
      break
    }
  }
  #minithetam value
  est=init.theta.all
  mx=m(est,x)
  arg=1+mx%*%lambda
  val=2*sum(llog(arg,1/n))
  
  list(estimate=est, lambda=lambda, val=val)
}





bic.ee.el.2penalty=function(m,x,theta0,tauv,nuv,thres.v,method){
  nn=nrow(x)
  p=length(theta0)
  r=ncol(m(theta0,x))
  n1=length(tauv)
  n2=length(nuv)
  if(method=="newton"){n2=1}
  val=matrix(rep(0,n1*n2),nrow=n1)
  est=array(0,dim=c(n1,p,n2)) #n1*p*n2
  lambda=array(0,dim=c(n1,r,n2))
  theta.init=theta0
  for(i in 1:n1) {
    for(j in 1:n2) {
      tryCatch({
        cat(paste("i=",i,"j=",j,"tau=",tauv[i],"nu=",nuv[j],sep="\t","\n"))
        #if(method=="newton"){
        #  fit=elhood.ee.lqa.newton(m,x,theta.init,m_1dev,m_2dev,tauv[i],thres.v)
        #}else if(method=="2penalty"){
        #  fit=elhood.ee.lqa.2penalty(m,x,theta.init,tauv[i],nuv[j],thres.v)
        #}else if(method=="2penalty.tang"){
        fit=elhood.ee.lqa.2penalty(m,x,theta.init,tauv[i],nuv[j],thres.v)
        #}else if(method=="2penalty.newton"){
        #  fit=elhood.ee.lqa.2penalty.newton(m,x,theta.init,m_1dev,m_2dev,tauv[i],nuv[j],thres.v)
        #}else if(method=="2penalty.cd"){
        #  fit=elhood.ee.lqa.2penalty.cd(m,x,theta.init,tauv[i],nuv[j],thres.v)
        #}
        val[i,j]=fit$val
        est[i,,j]=fit$estimate
        lambda[i,,j]=fit$lambda
        #theta.init=est[i,,j] #warm start
      },error=function(e){
        cat("ERROR :",conditionMessage(e), "\n")
        val[i,j]=999999
        est[i,j]=0
        lambda[i,j]=0
      })
      cat(paste(est[i,,j],sep="\t"),"\n")
      cat(paste(lambda[i,,j],sep="\t"),"\n")
    }
  }
  dof=apply(abs(est)>0,1,sum)
  bic=val+dof*log(nn)
  bicc=val+dof*log(nn)*max(1,log(log(ncol(x))))
  ebic=val+dof*(log(nn)+2*log(ncol(x)))
  min.ind=which(bic==min(bic),arr.ind=T)
  if(nrow(min.ind)>1){min.ind=min.ind[1,]}
  est.bic=est[min.ind[1],,min.ind[2]]
  lambda.bic=lambda[min.ind[1],,min.ind[2]]
  min.ind=which(bicc==min(bicc),arr.ind=T)
  if(nrow(min.ind)>1){min.ind=min.ind[1,]}
  est.bicc=est[min.ind[1],,min.ind[2]]
  lambda.bicc=lambda[min.ind[1],,min.ind[2]]
  min.ind=which(ebic==min(ebic),arr.ind=T)
  if(nrow(min.ind)>1){min.ind=min.ind[1,]}
  est.ebic=est[min.ind[1],,min.ind[2]]
  lambda.ebic=lambda[min.ind[1],,min.ind[2]]
  list(est=est,lambda=lambda,bic=bic,bicc=bicc,ebic=ebic,dof=dof,
       est.bic=est.bic,est.bicc=est.bicc,est.ebic=est.ebic,
       lambda.bic=lambda.bic,lambda.bicc=lambda.bicc,lambda.ebic=lambda.ebic,
       l1=min(bic),l2=min(bicc),l3=min(ebic),vv=val)
}