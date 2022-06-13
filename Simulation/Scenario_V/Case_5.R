require("nleqslv")
require("rmutil")
require("statmod")
require("bayesreg")
require("mvtnorm")
require("monomvn")
require("ncvreg")
require("glmnet")
require("coda")
require("BayesS5")
require("picasso")
require("SIS")
library("MASS")


res_1 = matrix(NA,28,100)
for (ii in 1:100) {
  set.seed(ii)
  n=100
  p=3200
  sde=1 
  sigma = matrix(rep(0.5,p*p),p)+diag(0.5,p)
  x =rmvnorm(n=n, mean=rep(0,p), sigma=sigma)
  b = matrix(0,p,1)
  
  b[1:20] =  rnorm(20,5,3)
  error_term = rt(n,3)
  y = x %*% b + (rbinom(n,1,.5)*2-1)*3 + error_term/sd(error_term)*sde
  
  
  X = x
  mu = apply(X,2,mean)
  X = t(t(X) - mu)
  Y = y
  
  #SCAD
  ptm_SCAD <- proc.time()
  res_SCAD = cv.ncvreg(X, Y,penalty="SCAD")
  beta_SCAD <- res_SCAD$fit$beta[,res_SCAD$min]
  time_SCAD=proc.time() - ptm_SCAD
  MSE1 = sum((beta_SCAD[-1]-b)^2)/p
  cor_1= sum(abs(beta_SCAD[c(2:21)])>0)
  wro_1 = sum(abs(beta_SCAD[-c(1:21)])>0)
  
  # lasso
  ptm_la <- proc.time()
  cv <- lars::cv.lars(X, Y, plot.it = FALSE, mode = "step")
  idx <- which.max(cv$cv - cv$cv.error <= min(cv$cv))
  beta_la=coef(lars::lars(X, Y))[idx,]
  time_la=proc.time() - ptm_la
  
  MSE2 = sum((beta_la-b)^2)/p
  cor_2= sum(abs(beta_la[c(1:20)])>0)
  wro_2 = sum(abs(beta_la[-c(1:20)])>0)
  
  # EL_DP
  A = t(X)%*%X/n
  B = (t(X)%*%Y)/n
  beta_ols = picasso(X, Y, nlambda=1000, type.gaussian="naive",lambda.min.ratio = 0.01)$beta[,1000]
  
  tauv2=exp(seq(-3,-1,by=1))
  nuv2=exp(seq(-3,-2,by=1))
  thres.v=0.05;
  xy=cbind(Y,X)
  ptm_ELDP <- proc.time()
  res3 =bic.ee.el.2penalty(m,xy,beta_ols,tauv2,nuv2,thres.v,method="2penalty.tang")
  time_ELDP=proc.time() - ptm_ELDP
  MSE3=sum((res3$est.bic-b)^2)/p 
  cor_3= sum(abs(res3$est.bic[c(1:20)])>0)
  wro_3 = sum(abs(res3$est.bic[-c(1:20)])>0)
  
  #SIS
  ptm_SIS <- proc.time()
  res4 = SIS(X, Y, family='gaussian', penalty = 'lasso', tune='bic',perm=TRUE,
             greedy=TRUE)
  time_SIS=proc.time() - ptm_SIS
  b_est_sis2 = rep(0,p)
  b_est_sis2[res4$ix]= res4$coef.est[-1]
  
  MSE4 = sum((b_est_sis2 - b)^2)/p
  cor_4= sum(abs(b_est_sis2[c(1:20)])>0)
  wro_4 = sum(abs(b_est_sis2[-c(1:20)])>0)
  
  # BEL_HD
  ptm_HD <- proc.time()
  res5 = BEL_HD(X,Y,thres_s =0.01,Quick=FALSE)
  time_HD=proc.time() - ptm_HD
  MSE5 = sum((res5$beta - b)^2)/p
  cor_5= sum(abs(res5$beta[c(1:20)])>0)
  wro_5 = sum(abs(res5$beta[-c(1:20)])>0)
  
  #PICASSO
  ptm_pic <- proc.time()
  res6 = picasso(X, Y, nlambda=1000, type.gaussian="naive",lambda.min.ratio = 0.01)$beta[,1000]
  time_pic=proc.time() - ptm_pic
  MSE6 = sum((res6 - b)^2)/p
  cor_6= sum(abs(res6[c(1:20)])>0)
  wro_6 = sum(abs(res6[-c(1:20)])>0)
  
  #Bayes_S5
  ptm_S5 <- proc.time()
  fit_default = S5(X,Y)
  res_default = result(fit_default)
  est.LS = result_est_LS(res_default,X,Y) # Averged over the Least 
  time_S5=proc.time() - ptm_S5
  MSE7 = sum((est.LS$beta.MAP-b)^2)/p
  cor_7= sum(abs(est.LS$beta.MAP[c(1:20)])>0)
  wro_7 = sum(abs(est.LS$beta.MAP[-c(1:20)])>0)
  
  res_1[,ii]=c(MSE1,MSE2,MSE3,MSE4,MSE5,MSE6,MSE7,
               cor_1,cor_2,cor_3,cor_4,cor_5,cor_6,cor_7,
               wro_1,wro_2,wro_3,wro_4,wro_5,wro_6,wro_7,
               time_SCAD[3],time_la[3],time_ELDP[3],time_SIS[3],
               time_HD[3],time_pic[3],time_S5[3])
}
rownames(res_1) = c("SCAD","LASSO","EL_DP","SIS","BEL_HD","PICASSO","BayesS5",
                    "SCAD_cor","LASSO_cor","EL_DP_cor","SIS_cor","BEL_HD_cor",
                    "PICASSO_cor","BayesS5_cor",
                    "SCAD_wro","LASSO_wro","EL_DP_wro","SIS_wro","BEL_HD_wro",
                    "PICASSO_wro","BayesS5_wro",
                    "SCAD_time","LASSO_time","EL_DP_time","SIS_time","BEL_HD_time",
                    "PICASSO_time","BayesS5_time") 
write.csv(res_1, file = "Result_100_3200.csv")



