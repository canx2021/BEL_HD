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


getcorx <- function(n,rho,x1) {
  # desired correlation = cos(angle)
  theta <- acos(rho)             # corresponding angle
  x2    <- rnorm(n, 0, 1)      # new random data
  X     <- cbind(x1, x2)         # matrix
  Xctr  <- scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)
  Id   <- diag(n)                               # identity matrix
  Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
  P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
  x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
  Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
  Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1
  x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]
  return(x)
}

y = UScrime[,16]
n=47
p_rand=4*n-15*2
p=4*n
corr=0.3
res_1 = matrix(NA,18,100)
# Run 100 replicates
set.seed(874135)
for (ii in 1:100) {
  X_rand_1 = getcorx(n,corr,x=UScrime[,1])
  X_rand_2 = getcorx(n,corr,x=UScrime[,2])
  X_rand_3 = getcorx(n,corr,x=UScrime[,3])
  X_rand_4 = getcorx(n,corr,x=UScrime[,4])
  X_rand_5 = getcorx(n,corr,x=UScrime[,5])
  X_rand_6 = getcorx(n,corr,x=UScrime[,6])
  X_rand_7 = getcorx(n,corr,x=UScrime[,7])
  X_rand_8 = getcorx(n,corr,x=UScrime[,8])
  X_rand_9 = getcorx(n,corr,x=UScrime[,9])
  X_rand_10 = getcorx(n,corr,x=UScrime[,10])
  X_rand_11 = getcorx(n,corr,x=UScrime[,11])
  X_rand_12 = getcorx(n,corr,x=UScrime[,12])
  X_rand_13 = getcorx(n,corr,x=UScrime[,13])
  X_rand_14 = getcorx(n,corr,x=UScrime[,14])
  X_rand_15 = getcorx(n,corr,x=UScrime[,15])
  sigma = matrix(rep(0.1,p_rand*p_rand),p_rand)+diag(0.9,p_rand)
  X_rand_16 = rmvnorm(n=n, mean=rep(0,p_rand), sigma=sigma)
  x=cbind(UScrime[,-16],X_rand_1,X_rand_2,X_rand_3,X_rand_4,
          X_rand_5,X_rand_6,X_rand_7,X_rand_8,X_rand_9,X_rand_10,
          X_rand_11,X_rand_12,X_rand_13,X_rand_14,
          X_rand_15,X_rand_16)
  X = scale(x)
  Y = log(y)-mean(log(y))
  N_test = 5
  N_total = length(y)
  N_train = N_total - N_test
  train_ind = sample(1:N_total,N_train)
  X_train = X[train_ind,]
  X_test = X[-train_ind,]
  # SCAD
  res_SCAD = cv.ncvreg(X_train, Y[train_ind],penalty="SCAD")
  beta_SCAD <- res_SCAD$fit$beta[,res_SCAD$min]
  cor_1= sum(abs(beta_SCAD[c(2:16)])>0)
  wro_1 = sum(abs(beta_SCAD[-c(1:16)])>0)
  mspe_1 = mean((Y[-train_ind] - X_test%*%beta_SCAD[-1])^2)
  # lasso
  cv <- lars::cv.lars(X_train, Y[train_ind], plot.it = FALSE,mode = "step")
  idx <- which.max(cv$cv - cv$cv.error <= min(cv$cv))
  beta_la=coef(lars::lars(X_train, Y[train_ind]))[idx,]
  mspe_2 = mean((Y[-train_ind] - X_test%*%beta_la)^2)
  cor_2= sum(abs(beta_la[c(1:15)])>0)
  wro_2 = sum(abs(beta_la[-c(1:15)])>0)
  #SIS
  res4 = SIS(X_train, Y[train_ind], family='gaussian', penalty = 'lasso', tune='bic',perm=TRUE)
  b_est_sis2 = rep(0,p)
  b_est_sis2[res4$ix]= res4$coef.est[-1]
  mspe_4 = mean((Y[-train_ind] - X_test%*%b_est_sis2)^2)
  cor_4= sum(abs(b_est_sis2[c(1:15)])>0)
  wro_4 = sum(abs(b_est_sis2[-c(1:15)])>0)
  # BEL_HD
  res5 = BEL_HD(X_train, Y[train_ind],thres_s=0.08,Quick=FALSE)
  cor_5= sum(abs(res5$beta[c(1:15)])>0)
  wro_5 = sum(abs(res5$beta[-c(1:15)])>0)
  mspe_5 = mean((Y[-train_ind] - X_test%*%res5$beta)^2)
  # PICASSO
  res6 = picasso(X_train, Y[train_ind], nlambda=1000, type.gaussian="naive")$beta[,1000]
  mspe_6 = mean((Y[-train_ind] - X_test%*%res6)^2)
  cor_6= sum(abs(res6[c(1:15)])>0)
  wro_6 = sum(abs(res6[-c(1:15)])>0)
  # Bayes_S5
  fit_default = S5(X_train, Y[train_ind])
  res_default = result(fit_default)
  est.LS = result_est_LS(res_default,X_train, Y[train_ind]) 
  mspe_7 = mean((Y[-train_ind] - X_test%*%est.LS$beta.MAP)^2)
  cor_7= sum(abs(est.LS$beta.MAP[c(1:15)])>0)
  wro_7 = sum(abs(est.LS$beta.MAP[-c(1:15)])>0)
  
  res_1[,ii]=c(cor_1,cor_2,cor_4,cor_5,cor_6,cor_7,wro_1,wro_2,wro_4,wro_5,wro_6,wro_7,mspe_1,mspe_2,mspe_4,mspe_5,mspe_6,mspe_7)}

write.csv(res_1, file = "Result_crime_cor03.csv")



