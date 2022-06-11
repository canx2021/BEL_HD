# Example: gene expression data
library(GEOquery)
data=getGEO("GSE3330") 
datExpr1 = exprs(data[[1]])
datExpr2 = exprs(data[[2]])
x = rbind(datExpr1,datExpr2)
y = x[296,]
X = t(x[-296])
X = scale(X)
Y = y -mean(y)

N_total = length(y)
N_train = 55
N_test = N_total-N_train
p=dim(X)[2]
res_1 = numeric(100)
set.seed(1020)

for (ii in 1:100) {
  
  train_ind = sample(1:N_total,N_train)
  X_train = X[train_ind,]
  X_test = X[-train_ind,]
  
  res5 = BEL_HD(X_train, Y[train_ind],Oiter =1000,Oburnin = 500,thres_s=0.1)
  
  mspe_5 = mean((Y[-train_ind] - X_test%*%res5$beta)^2)
  
  res_1[ii]=mspe_5}

mean(res_1)

