library(GEOquery)
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


data=getGEO("GSE3330") 

datExpr1 = exprs(data[[1]])

datExpr2 = exprs(data[[2]])

x = rbind(datExpr1,datExpr2)

y = x[296,]
X = t(x[-296,])

X = scale(X)
Y = y -mean(y)
p=dim(x)[1]
N_total = length(y)
N_train = 55
N_test = N_total-N_train

p=dim(X)[2]
res_1 = matrix(NA,5,100)
set.seed(1020)

for (ii in 1:100) {
  
  train_ind = sample(1:N_total,N_train)
  X_train = X[train_ind,]
  X_test = X[-train_ind,]
  
  res_SCAD = cv.ncvreg(X_train, Y[train_ind],penalty="SCAD")
  beta_SCAD <- res_SCAD$fit$beta[-1,res_SCAD$min]
  
  mspe_1 = mean((Y[-train_ind] - X_test%*%beta_SCAD)^2)
 
  res4 = SIS(X_train, Y[train_ind], family='gaussian', penalty = 'lasso', tune='bic',perm=TRUE,
             greedy=TRUE)
  
  b_est_sis2 = rep(0,p)
  b_est_sis2[res4$ix]= res4$coef.est[-1]
  mspe_4 = mean((Y[-train_ind] - X_test%*%b_est_sis2)^2)
  
  res6 = picasso(X_train, Y[train_ind], nlambda=1000, type.gaussian="naive")$beta[,1000]
  
  mspe_6 = mean((Y[-train_ind] - X_test%*%res6)^2)
  
  res5 = BEL_HD(X_train, Y[train_ind],Oiter =1000,Oburnin = 500,thres_s=0.01)
  
  mspe_5 = mean((Y[-train_ind] - X_test%*%res5$beta)^2)
  
  
  fit_default = S5(X_train, Y[train_ind])
  res_default = result(fit_default)
  est.LS = result_est_LS(res_default,X_train, Y[train_ind]) 
  mspe_7 = mean((Y[-train_ind] - X_test%*%est.LS$beta.MAP)^2)
  
  res_1[,ii]=c(mspe_1,mspe_4,mspe_5,mspe_6,mspe_7)}

write.csv(res_1, file = "Gene_all.csv")


