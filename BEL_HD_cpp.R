setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
require("Rcpp")
sourceCpp("Act_MH.cpp")
sourceCpp("grad.cpp")
require("picasso")
require("statmod")
require("mvtnorm")
require("MASS")

# Pseudo-logarithm function
llog <- function (z, eps=1/50)
{
  z[is.na(z)] <- eps
  ans <- z
  lo <- (z < eps)
  ans[lo] <- log(eps) - 1.5 + 2 * z[lo]/eps - 0.5 * (z[lo]/eps)^2
  ans[!lo] <- log(z[!lo])
  ans
}
# Gradient of pseudo-logarithm function, which avoids the problem of a vanishing denominator
llogp <- function (z, eps=1/50)
{
  z[is.na(z)] <- eps
  ans <- z
  lo <- (z < eps)
  ans[lo] <- 2/eps - z[lo]/eps^2
  ans[!lo] <- 1/z[!lo]
  ans
}

#'@title Bayesian variable selection based on empirical likelihood for ultra-high dimensional data
#'@description A great amount of literature has shown that development of variable selection methods 
#'can enable efficient and interpretable analysis of high dimensional data that are ubiquitous in 
#'recent years. However, much less attention has been paid to Bayesian empirical likelihood for 
#'variable selection involving ultrahigh-dimensional data where the number of covariates p is 
#'(much) large than the sample size n. In the semi-parametric domain, under the ultra-high 
#'dimensional setting, we propose a Bayesian variable selection method based on empirical likelihood, 
#'which requires only estimating equations and no distributional assumptions.
#'@param X n*p covariate matrix, with each column representing a covariate and standardized 
#'to have mean 0 and standard deviation 1
#'@param Y the outcome vector of length n, centered to have mean 0
#'@param Oiter Total number of iterations for outerloop. Default is 300.
#'@param Iiter Total number of iterations for innerloop. Default is 100.
#'@param Oburnin Number of burnin cycles for outerloop. Default is 150.
#'@param Iburnin Number of burnin cycles for innerloop. Default is 20.
#'@param thres_s pre-specified thresholds, which can be set to be a small positive value 
#'and then adjusted by the trial and error approach or by cross-validation.
#'@param Quick logical variable. "TRUE" would simplify the algorithm, which is normally 
#'applied to ultra-high dimensional problem. Default is TRUE. 
#'@details Motivated by Chang et al. (2018) on doubly penalized empirical likelihood (EL), 
#'we introduce priors to regularize both regression parameters and Lagrange multipliers 
#'associated with the estimating equations. We further develop an efficient Markov chain 
#'Monte Carlo sampling algorithm based on the active set idea, first proposed by 
#'Friedman et al. (2007). The proposed method not only inherits merits from the Bayesian 
#'paradigm and shares the attractiveness of EL approaches, but also has superior performance 
#'in both the prediction and variable selection, as shown in our numerical studies.
#'@return The function returns a list of elements including: a list of MCMC samples and the 
#'final estimate for regression coefficients. 
#'@export
#'@import MASS
#'@import statmod
#'@import mvtnorm
#'@import picasso
#'@import Rcpp
BEL_HD = function(X,Y,Oiter = 300,Iiter = 100, Oburnin = 150,Iburnin = 20,thres_s=0.01,Quick=TRUE){
  ptm <- proc.time()
  p = ncol(X)
  n = length(Y)
  # Initialization
  Beta = matrix(0,Oiter,p)
  Lambda = matrix(0,Oiter,p)
  index = matrix(0,Oiter,p)
  gam1= numeric(Oiter)
  gam2= numeric(Oiter)
  
  gam_1 = rexp(1,rate=0.5)
  gam_2 = rexp(1,rate=0.5)
  gam1[1] = gam_1
  gam2[1] = gam_2
  
  beta = picasso(X, Y, nlambda=1000, type.gaussian="naive",lambda.min.ratio = 0.01)$beta[,1000]+0.0001
  
  if(p<3200 & Quick==FALSE){
    TT = X*as.vector(Y-X%*%beta)
    
    U_l = function(q){
      return(sum(llog(1+TT%*%q)))}
    
    grad_U_l = function(q){
      return(colSums(TT*as.vector(llogp(1+TT%*%q))))
    }
    
    lambda = optim(rep(0.001,p), U_l,grad_U_l, control=list(fnscale=-1),method = "BFGS")$par
    
  }else{
    lambda = rep(0.00001,p)
  }
  
  Beta[1,] = beta
  Lambda[1,] = lambda
  index_act = which(abs(beta)>thres_s)
  index_act = sort(index_act)
  index[1,] = (abs(beta)>thres_s)
  p_new = length(index_act)
  Tau_sum = sum(1/statmod::rinvgauss(rep(1, p_new),sqrt(gam_1^2/beta[index_act]^2),gam_1^2))# squared
  V_sum = sum(1/statmod::rinvgauss(rep(1, p_new),sqrt(gam_2^2/lambda[index_act]^2),gam_2^2))
  for(iter in 2:Oiter){
    if(iter%%10==0) print(iter)
    # first check the length of active set
    if (length(index_act) == 0){
      print("No active set")
      break
    }else if (length(index_act) == 1){
      cat("The only active set is:", index_act)
      break
    }
    # Step 1: update gamma
    gam1[iter] = sqrt(rgamma(1, shape=p_new+.1, rate = Tau_sum/2+.1))
    gam2[iter] = sqrt(rgamma(1, shape=p_new+.1, rate = V_sum/2+.1))
    
    
    # step 2: inner loop to update beta and lambda
    res_new = Inner_loop(X,Y,index_act,Beta[iter-1,],Lambda[iter-1,],gam1[iter],gam2[iter],Iiter=Iiter,Iburnin=Iburnin)
    
    print(res_new$acp_b)
    
    if (res_new$acp_b<0.02 & thres_s<1.0){
      thres_s = thres_s+0.01
    }
    
    Beta[iter,] = res_new$beta
    Lambda[iter,]=res_new$lambda
    
    # step 3: update active set
    if(p<3200 & Quick==FALSE){
      grad = def_grad(X,Y,Lambda[iter,],Beta[iter,])
      potential_ind_grad = which(abs(grad)>3*max(abs(grad[index_act])))
      index_act = unique(c(which(abs(Beta[iter,])>thres_s),
                           which(abs(grad)>quantile(abs(grad[potential_ind_grad]),0.9))))
    }else{
      index_act = which(abs(Beta[iter,])>thres_s)
    }
    index[iter,] = (abs(Beta[iter,])>thres_s)
    
    p_new=length(index_act)
    Tau_sum = res_new$Tau_sum
    V_sum = res_new$V_sum
    
  }
  # Output
  reg = list()
  if (length(index_act) == 0 ){
    reg$beta = rep(0,p)
  }else if (length(index_act) == 1){
    beta_final = rep(0,p)
    A_s= t(X[,index_act])%*%X[,index_act]/n
    B_s = (t(X[,index_act])%*%Y)/n
    beta_final[index_act] = ginv(A_s)%*%(B_s)
    reg$beta = beta_final
  }else{
    burn = Oburnin
    beta_final = numeric(p)
    beta_final[index_act] = apply(Beta[-(1:burn),index_act],2,mean)
    beta_final[-index_act]=0
    reg$beta = beta_final
    reg$Beta = Beta
    reg$Lambda = Lambda
    reg$index=index
  }
  print(proc.time() - ptm)
  return(reg)
}

# Inner loop function
Inner_loop = function(X,Y,index_act,beta_int,lambda_int,gamma_1,gamma_2,Iiter=100,Iburnin=20){
  X_new = X[,index_act]
  p_new = length(index_act)
  n = length(Y)
  # Initialization for inner loop
  Beta_I = matrix(0,Iiter,p_new)
  Lambda_I = matrix(0,Iiter,p_new)
  Tau_I = matrix(0,Iiter,p_new)
  V_I = matrix(0,Iiter,p_new)
  gam1=numeric(Iiter)
  gam2=numeric(Iiter)
  accept=rep(0,Iiter)
  accept[1]=1
  A_I = t(X_new)%*%X_new/n
  B_I = (t(X_new)%*%Y)/n
  beta_ols_I = ginv(A_I)%*%(B_I)
  Beta_I[1,] = beta_ols_I
  
  U = X_new*as.vector(Y-X_new%*%beta_ols_I)
  U_l = function(q){
    return(sum(llog(1+U%*%q)))}
  grad_U_l = function(q){
    return(colSums(U*as.vector(llogp(1+U%*%q))))
  }
  lambda = optim(rep(0.0001,p_new), U_l,grad_U_l, control=list(fnscale=-1),method = "BFGS")$par
  
  Lambda_I[1,] = lambda
  Eps_I = diag(as.vector(Y-X_new%*%beta_ols_I)^2)
  Tn_I = t(X_new)%*%Eps_I%*%X_new/n
  
  for (iter_i in 2:Iiter) {
    # step 1: sample tau and nu
    Tau_I[iter_i,] = 1/statmod::rinvgauss(rep(1, p_new),sqrt(gamma_1^2/Beta_I[iter_i-1,]^2),gamma_1^2) # squared
    V_I[iter_i,] = 1/statmod::rinvgauss(rep(1, p_new),sqrt(gamma_2^2/Lambda_I[iter_i-1,]^2),gamma_2^2) # squared
    # step 2: sample beta and lambda
    Beta_I[iter_i,] = Beta_I[iter_i-1,] 
    Lambda_I[iter_i,] = Lambda_I[iter_i-1,] 
    
    Sigma_lp = ginv(n*t(A_I)%*%ginv(Tn_I)%*%A_I + diag(1/Tau_I[iter_i,]))
    mu_lp = t(t(B_I)%*%ginv(Tn_I)%*%A_I%*%Sigma_lp)*n
    U_old = X_new*as.vector(Y-X_new%*%Beta_I[iter_i-1,])
    
    lambda_cpp_int = Lambda_I[iter_i-1,]
    beta_cpp_int = Beta_I[iter_i-1,]
    tau_cpp_int = Tau_I[iter_i,]
    v_cpp_int = V_I[iter_i,]
    out = dp_cpp(U_old,X_new,Sigma_lp,lambda_cpp_int,beta_cpp_int,tau_cpp_int,v_cpp_int,Y,mu_lp,p_new)
    if(out$accept==1){
      Beta_I[iter_i,] = out$beta
      accept[iter_i] = 1
      Lambda_I[iter_i,] = out$lambda
    }
  }
  beta_final = beta_int
  beta_final[index_act] = colMeans(Beta_I[Iburnin:Iiter,])
  lambda_final = lambda_int
  lambda_final[index_act] = colMeans(Lambda_I[Iburnin:Iiter,])
  Tau_sum = sum(colMeans(Tau_I[Iburnin:Iiter,]))
  V_sum = sum(colMeans(V_I[Iburnin:Iiter,]))
  # Output
  res = list(beta=beta_final,lambda=lambda_final,
             acp_b = sum(accept)/Iiter,Tau_sum=Tau_sum,
             V_sum=V_sum)
}

# Penalty function
prior_new = function(beta,tau){
  pi_beta = - sum(beta^2/tau/2)
  return(sum(pi_beta))
}
# Log likelihood function
lik_new <- function(U,lambda){
  a = llog(1/(1+U%*%lambda))
  if(sum(is.na(a))<3){return(sum(a,na.rm=TRUE))}else{return(sum(a))}
}
# Penalized log likelihood function
penlik_new <- function(beta,tau,U,lambda){
  lik_new(U,lambda) + prior_new(beta,tau)
}

