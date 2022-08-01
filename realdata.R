
rm(list = ls())
library("tidyverse") #the package to delect NA
#library('nloptr')

setwd("/vhome/shidapeng/ARIC")
ARIC = read.csv("ARIC_pheno.csv")



ARIC = (ARIC %>% drop_na(gender, v1age01, bmi01, hypert05))
exposure = c("bmi01")
outcome = c("hypert05")

signal = read.table("1e-4.txt")

sample.id.old = as.vector(unlist(read.table("ID2.txt"))) #All the individuals' ID
#sample.id.new = sample.id.old[sample.id.old %in% ARIC[,2]]
sample.id.new = sample.id.old[sample.id.old %in% ARIC[,2]]

valid.sample.number = length(sample.id.new)

computation.id = sample.id.new

#regress on the covariate of "age" and "gender"
ARIC.gender = as.character(ARIC[ ,"gender"])
head(ARIC.gender)
ARIC.gender[ARIC.gender == "M"] = c("1")
ARIC.gender[ARIC.gender == "F"] = c("0")

ARIC.gender = as.numeric(ARIC.gender)

ARIC.age = ARIC[ ,"v1age01"]
ARIC.exposure = ARIC[, exposure]

ARIC.resi = data.frame(y = ARIC.exposure,  gender = ARIC.gender, age = ARIC.age)
coeff = lm(y ~ -1 + gender + age, data = ARIC.resi)$coefficients
#compute the residual
exposure.resi = ARIC.exposure - as.matrix(ARIC.resi[,2:3])%*%coeff
ARIC[, exposure] = exposure.resi


genotype.orig = NULL
for(i in 1:22)
{
  path = paste0("chr", i, ".DS.txt")
  genotype.orig = cbind(genotype.orig, t(read.table(file = path)))
}

geno.name = genotype.orig[1,]                              #name of snps
row = dim(genotype.orig)[1]                             #number of snps
genotype = data.frame(id = c(sample.id.old), geno.name = genotype.orig[2:row,])                                  # integrate samples and snps 
names(genotype) = c("id", geno.name)# modify its name
genotype = genotype[genotype[,1] %in% sample.id.new,]


genotype.computation = genotype[genotype[,1] %in% computation.id,]


ARIC.computation = ARIC[ARIC[,2]%in% computation.id, ]

pheno_x.computation = ARIC.computation[, colnames(ARIC.computation) %in% exposure ]


pheno_x.inveranktrans.computation = qnorm((rank(pheno_x.computation, na.last="keep")-0.5)/sum(!is.na(pheno_x.computation)))
pheno_y.computation = ARIC.computation[, colnames(ARIC.computation) %in% outcome ]  #phenotype y

  
tmp1 = as.data.frame(lapply(genotype.computation[,2:dim(genotype.computation)[2]], as.numeric))
genotype.center.computation = cbind(genotype.computation[,1], tmp1 - apply(tmp1, 2, mean)) #centering gene matrix
names(genotype.center.computation) = names(genotype)

valid.geno.name = names(genotype[2:length(names(genotype))])
IV.name = intersect(valid.geno.name, signal[,1])









ARIC = (ARIC %>% drop_na(sbpa22, hypert05))
ARIC.male = ARIC
signal = read.table("dbp.1e-3.txt")

sample.id.old = as.vector(unlist(read.table("ID2.txt"))) #All the individuals' ID
#sample.id.new = sample.id.old[sample.id.old %in% ARIC[,2]]
sample.id.new.male = sample.id.old[sample.id.old %in% ARIC.male[,2]]

valid.sample.number = length(sample.id.new.male)

computation.id = sample.id.new.male

genotype.orig = NULL
for(i in 1:22)
{
  path = paste0("chr", i, ".DS.txt")
  genotype.orig = cbind(genotype.orig, t(read.table(file = path)))
}


geno.name = genotype.orig[1,]                              #name of snps
row = dim(genotype.orig)[1]                             #number of snps
genotype = data.frame(id = c(sample.id.old), geno.name = genotype.orig[2:row,])                                  # integrate samples and snps 
names(genotype) = c("id", geno.name)# modify its name
genotype = genotype[genotype[,1] %in% sample.id.new.male,]


genotype.computation = genotype[genotype[,1] %in% computation.id,]


ARIC.computation = ARIC[ARIC[,2]%in% computation.id, ]

pheno_x.computation = ARIC.computation[, colnames(ARIC.computation) %in% c("sbpa22")]


pheno_x.inveranktrans.computation = qnorm((rank(pheno_x.computation, na.last="keep")-0.5)/sum(!is.na(pheno_x.computation)))
pheno_y.computation = ARIC.computation[, colnames(ARIC.computation) %in% c("hypert05")]  #phenotype y

  
tmp1 = as.data.frame(lapply(genotype.computation[,2:dim(genotype.computation)[2]], as.numeric))
genotype.center.computation = cbind(genotype.computation[,1], tmp1 - apply(tmp1, 2, mean)) #centering gene matrix
names(genotype.center.computation) = names(genotype)

valid.geno.name = names(genotype[2:length(names(genotype))])
IV.name = intersect(valid.geno.name, signal[,1])




P_u = function(r)
{
  exp(-0.5*r^2)/sqrt(2*pi)
}




#EM2 pre-training for suitable initial values
grad2_train = function(para1, para2, train.data)
{
  
  train.length = dim(train.data)[1]
  X = train.data[, "X"]
  Y = train.data[, "Y"]
  Z = train.data[, Zk]
  beta1.new = beta1
  
  gradbeta0 = 0
  gradgamma = 0
  gradlogsigma1 = 0
  gradlogsigma2 = 0
  gradrho = 0
  
  beta0 = para1[1]
  gamma = para1[2]
  logsigma1 = para1[3]
  logsigma2 = para1[4]
  rho = para1[5]
  
  beta0.new = para2[1]
  gamma.new = para2[2]
  logsigma1.new = para2[3]
  logsigma2.new = para2[4]
  rho.new = para2[5]
  
  
  for (i in 1:train.length)
  {
    
    K = exp(logsigma1)/(exp(logsigma2)+1e-8)*rho*(X[i] - gamma*Z[i])
    J = function(r) beta0 + X[i]*beta1 + K + sqrt(1 - rho^2)*exp(logsigma1)*r
    P_yx = function(r) exp((Y[i] - 1)*J(r) - log(1 + exp(-J(r))) - 0.5*log(2*pi) - r^2/2)
    tmp = function(r) P_yx(r)
    I = integrate(tmp, -Inf, Inf)$value
    
    
    K.new = exp(logsigma1.new - logsigma2.new)*rho.new*(X[i] - gamma.new*Z[i])
    J.new = function(r)  beta0.new + X[i]*beta1.new + K.new + sqrt(1 - rho.new^2)*exp(logsigma1.new)*r
    
    #compute derivitives
    tmp1 = function(r) tmp(r)/(1 + exp(-J.new(r)))
    gradbeta0 = gradbeta0 + Y[i] - integrate(tmp1, -Inf, Inf)$value/I
    
    tmp2 = function(r) - exp(logsigma1.new - logsigma2.new)*rho.new*Z[i]*tmp(r)/(1 + exp(-J.new(r)))
    gradgamma = gradgamma - exp(logsigma1.new - logsigma2.new)*rho.new*Z[i]*Y[i] +(X[i] - gamma.new*Z[i])/exp(2*logsigma2.new)*Z[i] - integrate(tmp2, -Inf, Inf)$value/I 
    
    tmp3 = function(r) (Y[i] - 1/(1 + exp(-J.new(r))))*(rho.new*exp(logsigma1.new - logsigma2.new)*(X[i] - gamma.new*Z[i]) + exp(logsigma1.new)*sqrt(1-rho.new^2)*r)*tmp(r)
    gradlogsigma1 = gradlogsigma1 + integrate(tmp3, -Inf, Inf)$value/I
    
    
    tmp4 = function(r) -rho.new*exp(logsigma1.new - logsigma2.new)*(X[i] - gamma.new*Z[i])*tmp(r)/(1 + exp(-J.new(r)))
    gradlogsigma2 = gradlogsigma2 - Y[i]*rho.new*exp(logsigma1.new - logsigma2.new)*(X[i] - gamma.new*Z[i])- 1 + (X[i] - gamma.new*Z[i])^2*exp(-2*logsigma2.new) - integrate(tmp4, -Inf, Inf)$value/I
    
    tmp5 = function(r) (Y[i] - 1/(1 + exp(-J.new(r))))*(exp(logsigma1.new - logsigma2.new)*(X[i] - gamma.new*Z[i]) - rho.new*exp(logsigma1.new)*r/sqrt(1-rho.new^2))*tmp(r)
    gradrho = gradrho + integrate(tmp5, -Inf, Inf)$value/I
  }
  grad = -c(gradbeta0, gradgamma ,gradlogsigma1, gradlogsigma2 ,gradrho )
  return(grad)
}

conlike2_train = function(para1, para2, train.data)
{
  
  train.length = dim(train.data)[1]
  X = train.data[, "X"]
  Y = train.data[, "Y"]
  Z = train.data[, Zk]
  
  
  beta01 = para1[1]
  gamma1 = para1[2]
  logsigma11 = para1[3]
  logsigma21 = para1[4]
  rho1 = para1[5]
  
  beta11 = beta1
  beta12 = beta1
  
  flink = 0
  
  for (i in 1:train.length)
  {
    
    K1 = exp(logsigma11)/(exp(logsigma21)+1e-8)*rho1*(X[i] - gamma1*Z[i])
    P_yx = function(r)
    {    
      exp((Y[i]-1)*(beta01 + X[i]*beta11 + K1 + sqrt(abs(1 - rho1^2))*exp(logsigma11)*r) - 
            log(1 + exp(-(beta01 + X[i]*beta11 + K1 + sqrt(abs(1 - rho1^2))*exp(logsigma11)*r))) - 0.5*log(2*pi) - r^2/2)
    }
    
    Jointlike = function(r) P_yx(r)
    MargInteg2 =  integrate(Jointlike, -Inf, Inf)$value
    #MargInteg2 =  distrExIntegrate(Jointlike, -10, 10,  subdivisions = 200, abs.tol = rel.tol)
    
    beta02 = para2[1]
    gamma2 = para2[2]
    logsigma12 = para2[3]
    logsigma22 = para2[4]
    rho2 = para2[5]
    
    func = function(r)
    {
      K2 = exp(logsigma12)/(exp(logsigma22)+1e-8)*rho2*(X[i] - gamma2*Z[i])
      (Y[i]-1)*(beta02 + X[i]*beta12 + K2 + sqrt(abs(1 - rho2^2))*exp(logsigma12)*r) - 
        log(1 + exp(-(beta02 + X[i]*beta12 + K2 + sqrt(abs(1 - rho2^2))*exp(logsigma12)*r)))  - 0.5*log(2*pi) - r^2/2
    } 
    
    opt1 = function(r)  P_yx(r)*func(r)
    #MargInteg1 = integrate(opt1, -10, 10)$value
    MargInteg1 = distrExIntegrate(opt1, -20, 20, subdivisions = 200, abs.tol = rel.tol)
    flink =  flink - MargInteg1/MargInteg2 + 0.5*(X[i] - gamma2*Z[i])^2/exp(2*logsigma22) + logsigma22
  }  
  return(flink)
}


EM_algorithm2_train = function(beta00, gamma0, sigma10, sigma20, rho0, train.data, maxit = 100)
{
  
  par.old <<- c(beta00, gamma0, sigma10, sigma20, rho0)
  
  PostExp = function(para)
  {
    conlike2_train(par.old, para, train.data) 
  }
  grad.iter = function(para)
  {
    grad2_train(par.old, para, train.data)
  }
  

   par.new = optim(par = par.old, fn = PostExp, grad.iter, method =  "L-BFGS-B",
                  lower = c(par.old[1] - 10, par.old[2] - 2, -Inf, -Inf, -0.9), 
                    upper = c(par.old[1] + 10, par.old[2] + 2, Inf,  Inf, 0.9))$par

  iter = 1
  print(sqrt(sum((par.new - par.old)^2)) )
  while (sqrt(sum((par.new - par.old)^2)) > 1e-4 & iter < maxit) 
  { 

    par.old <<-  par.new
    PostExp = function(para)
    {
      conlike2_train(par.old, para, train.data) 
    }
    grad.iter = function(para)
    {
      grad2_train(par.old, para, train.data)
    }
    

   par.new = optim(par = par.old, fn = PostExp, grad.iter, method =  "L-BFGS-B",
                  lower = c(par.old[1] - 10, par.old[2] - 2, -Inf, -Inf, -0.9), 
                    upper = c(par.old[1] + 10, par.old[2] + 2, Inf,  Inf, 0.9))$par
    iter = iter + 1   
  }
  
  if (sqrt(sum((par.new - par.old)^2)) <= 1e-4)
  {conver2 = 1} else
  {conver2 = 0}
  return(c(par.new[1:5], conver2, sqrt(sum((par.new - par.old)^2))))
  
}




#EM algorithm for all  parameters

conlike1_train = function(para1, para2, train.data)
{ 
  train.length = dim(train.data)[1]
  X = train.data[, "X"]
  Y = train.data[, "Y"]
  Z = train.data[, Zk]
  
  beta01 = para1[1]
  beta11 = para1[2]
  gamma1 = para1[3]
  logsigma11 = para1[4]
  logsigma21 = para1[5]
  rho1 = para1[6]
  
  flink = 0
  
  for (i in 1:train.length)
  {
    
    K1 = exp(logsigma11)/(exp(logsigma21)+1e-8)*rho1*(X[i] - gamma1*Z[i])
    P_yx = function(r)
    {    
      exp((Y[i]-1)*(beta01 + X[i]*beta11 + K1 + sqrt(abs(1 - rho1^2))*exp(logsigma11)*r) - 
            log(1 + exp(-(beta01 + X[i]*beta11 + K1 + sqrt(abs(1 - rho1^2))*exp(logsigma11)*r)))- 0.5*log(2*pi) - r^2/2)
    }
    
    Jointlike = function(r) {P_yx(r)}
    MargInteg2 =  integrate(Jointlike, -Inf, Inf)$value
    #MargInteg2 =  distrExIntegrate(Jointlike, -10, 10,  subdivisions = 200, abs.tol = rel.tol)
    
    beta02 = para2[1]
    beta12 = para2[2]
    gamma2 = para2[3]
    logsigma12 = para2[4]
    logsigma22 = para2[5]
    rho2 = para2[6]
    
    func = function(r)
    {
      K2 = exp(logsigma12)/(exp(logsigma22)+1e-8)*rho2*(X[i] - gamma2*Z[i])
      (Y[i]-1)*(beta02 + X[i]*beta12 + K2 + sqrt(abs(1 - rho2^2))*exp(logsigma12)*r) - 
        log(1 + exp(-(beta02 + X[i]*beta12 + K2 + sqrt(abs(1 - rho2^2))*exp(logsigma12)*r)))- 0.5*log(2*pi) - r^2/2
    } 
    
    opt1 = function(r)  P_yx(r)*func(r)
    #MargInteg1 = integrate(opt1, -10, 10)$value
    MargInteg1 = distrExIntegrate(opt1, -10, 10, subdivisions = 100, abs.tol = rel.tol)
    flink =  flink - MargInteg1/MargInteg2 + 0.5*(X[i] - gamma2*Z[i])^2/exp(2*logsigma22) + logsigma22
  }  
  return(flink)
}

grad1_train = function(para1, para2, train.data)
{
  
  train.length = dim(train.data)[1]
  X = train.data[, "X"]
  Y = train.data[, "Y"]
  Z = train.data[, Zk]
  
  gradbeta0 = 0
  gradbeta1 = 0
  gradgamma = 0
  gradlogsigma1 = 0
  gradlogsigma2 = 0
  gradrho = 0
  
  beta0 = para1[1]
  beta1 = para1[2]
  gamma = para1[3]
  logsigma1 = para1[4]
  logsigma2 = para1[5]
  rho = para1[6]
  
  beta0.new = para2[1]
  beta1.new = para2[2]
  gamma.new = para2[3]
  logsigma1.new = para2[4]
  logsigma2.new = para2[5]
  rho.new = para2[6]
  
  
  for (i in 1:train.length)
  {
    
    K = exp(logsigma1)/(exp(logsigma2)+1e-8)*rho*(X[i] - gamma*Z[i])
    J = function(r) beta0 + X[i]*beta1 + K + sqrt(1 - rho^2)*exp(logsigma1)*r
    P_yx = function(r) exp((Y[i] - 1)*J(r) - log(1 + exp(-J(r)))  - 0.5*log(2*pi) - r^2/2 )
    tmp = function(r) P_yx(r)
    I = integrate(tmp, -Inf, Inf)$value
    
    
    K.new = exp(logsigma1.new - logsigma2.new)*rho.new*(X[i] - gamma.new*Z[i])
    J.new = function(r)  beta0.new + X[i]*beta1.new + K.new + sqrt(1 - rho.new^2)*exp(logsigma1.new)*r
    
    #compute derivitives
    tmp1 = function(r) tmp(r)/(1 + exp(-J.new(r)))
    gradbeta0 = gradbeta0 + Y[i] - integrate(tmp1, -Inf, Inf)$value/I
    
    tmp2 = function(r) X[i]*tmp(r)/(1 + exp(-J.new(r)))
    gradbeta1 = gradbeta1 + Y[i]*X[i] - integrate(tmp2, -Inf, Inf)$value/I
    
    tmp3 = function(r) - exp(logsigma1.new - logsigma2.new)*rho.new*Z[i]*tmp(r)/(1 + exp(-J.new(r)))
    gradgamma = gradgamma - exp(logsigma1.new - logsigma2.new)*rho.new*Z[i]*Y[i] +(X[i] - gamma.new*Z[i])/exp(2*logsigma2.new)*Z[i] - integrate(tmp3, -Inf, Inf)$value/I 
    
    tmp4 = function(r) (Y[i] - 1/(1 + exp(-J.new(r))))*(rho.new*exp(logsigma1.new - logsigma2.new)*(X[i] - gamma.new*Z[i]) + exp(logsigma1.new)*sqrt(1-rho.new^2)*r)*tmp(r)
    gradlogsigma1 = gradlogsigma1 + integrate(tmp4, -Inf, Inf)$value/I
    
    
    tmp5 = function(r) -rho.new*exp(logsigma1.new - logsigma2.new)*(X[i] - gamma.new*Z[i])*tmp(r)/(1 + exp(-J.new(r)))
    gradlogsigma2 = gradlogsigma2 - Y[i]*rho.new*exp(logsigma1.new - logsigma2.new)*(X[i] - gamma.new*Z[i])- 1 + (X[i] - gamma.new*Z[i])^2*exp(-2*logsigma2.new) - integrate(tmp5, -Inf, Inf)$value/I
    
    tmp6 = function(r) (Y[i] - 1/(1 + exp(-J.new(r))))*(exp(logsigma1.new - logsigma2.new)*(X[i] - gamma.new*Z[i]) - rho.new*exp(logsigma1.new)*r/sqrt(1-rho.new^2))*tmp(r)
    gradrho = gradrho + integrate(tmp6, -Inf, Inf)$value/I
  }
  grad = -c(gradbeta0 ,gradbeta1, gradgamma ,gradlogsigma1, gradlogsigma2 ,gradrho )
  return(grad)
}

EM_algorithm1_train = function(beta00, beta10, gamma0, sigma10, sigma20, rho0, train.data, maxit = 100)
{
  
  par.old <<- c(beta00, beta10, gamma0, sigma10, sigma20, rho0)
  PostExp = function(para)
  {
    postexp = conlike1_train(par.old, para, train.data) 
    postexp
  }
  grad.iter = function(para)
  {
    grad1_train(par.old, para, train.data)
  }
  
   par.new = optim(par = par.old, fn = PostExp, grad.iter , method =  "L-BFGS-B",
                  lower = c(par.old[1] - 10, par.old[2] - 10, par.old[3] - 2, -Inf, -Inf, -0.9), 
                  upper = c(par.old[1] + 10, par.old[2] + 10, par.old[3] + 2, Inf,  Inf, 0.9))$par 

  iter = 1
  while (sqrt(sum((par.new - par.old)^2)) > 1e-4 & iter < maxit) 
  { 
    par.old <<-  par.new
    PostExp = function(para)
    {
      postexp = conlike1_train(par.old, para, train.data) 
      postexp
    }
    grad.iter = function(para)
    {
      grad1_train(par.old, para, train.data)
    }
    
    par.new = optim(par = par.old, fn = PostExp, grad.iter , method =  "L-BFGS-B",
                  lower = c(par.old[1] - 10, par.old[2] - 10, par.old[3] - 2, -Inf, -Inf, -0.9), 
                  upper = c(par.old[1] + 10, par.old[2] + 10, par.old[3] + 2, Inf,  Inf, 0.9))$par 
    iter = iter + 1
  }
  
  if (sqrt(sum((par.new - par.old)^2)) <= 1e-4)
  {conver1 = 1} else
  {conver1 = 0}
  return(c(par.new[1:6], conver1, sqrt(sum((par.new - par.old)^2))))
  
}


#EM algorithm for fixed beta1 = beta1.true.


#compute the likelihood value


like.value = function(para)
{ 
  beta0 = para[1]
  beta1 = para[2]
  gamma = para[3]
  sigma1 = para[4]
  sigma2 = para[5]
  rho = para[6]
  
  f = 0
  for ( i in 1:n )
  { 
    likelihood = function(r)
    {
      K = sigma1/sigma2*rho*(X[i] - gamma*Z[i])
      exp(Y[i]*(beta0 + X[i]*beta1 + K + sqrt(1 - rho^2)*sigma1*r) - 
            log(1 + exp(beta0 + X[i]*beta1 + K + sqrt(1 - rho^2)*sigma1*r)) - 0.5*log(2*pi) - r^2/2)
      
    } 
    f = f + log(integrate(likelihood, -Inf, Inf)$value) - 0.5*log(2*pi*sigma2^2) - 0.5*(X[i] - gamma*Z[i])^2/sigma2^2
  }
  return(f)
}




library(mvtnorm)
library(foreach)
library(doParallel)
library(distrEx)

#find out the significant snp
setwd("/vhome/shidapeng/BinaryMRshi/sbpreal")          #set print out path

out = NULL
n = length(pheno_x.inveranktrans.computation )
X = pheno_x.inveranktrans.computation 
Y = pheno_y.computation
  
print(IV.name)

#traing data
test.matrix = data.frame(computation.id, pheno_y.computation, pheno_x.inveranktrans.computation, genotype.center.computation[, IV.name])
rownames(test.matrix) = NULL
colnames(test.matrix) = c("id" ,"Y", "X", as.character(IV.name))


  for (k in 6:length(IV.name))
  {
    Zk = as.character(IV.name[k])   #the location of the kth significant snp
    Z = genotype.center.computation[, Zk]
    x_z = data.frame(x = X, z = Z)
    xz.coeff = lm(x ~ -1 + z, data = x_z)$coefficients
    Xhat = Z*xz.coeff
    
    y_z =data.frame(y = Y, z = Z)
    yz.coeff = glm(y ~ z, family = binomial(), data = y_z)$coefficients[2]
    
    
    #directed estimator
    y_x = data.frame(y = Y, x = X)
    direc = coef(summary(glm(y~x, family = binomial(), data = y_x)))[2,]
    beta1.direc = direc[1]
    beta1.direc.sd = direc[2]
    beta1.direc.pval = direc[4]
    
    
    #two stage estimator
    y_xhat = data.frame(y = Y, xhat = Xhat)
    twostage = coef(summary(glm(y ~ xhat, family = binomial(), data = y_xhat)))[2,]
    beta1.twostage = twostage[1]
    beta1.twostage.sd = twostage[2]
    beta1.twostage.pval = twostage[4]
    
    
    #adjusted IV
    Vhat = X - Xhat
    y_xhat_vhat = data.frame(y = Y, xhat = Xhat, vhat = Vhat)
    beta1.adjIV = glm(y ~ xhat + vhat, family = binomial(), data = y_xhat_vhat)$coefficients[2]
    adjIV = coef(summary(glm(y ~ xhat + vhat, family = binomial(), data = y_xhat_vhat)))[2,]
    beta1.adjIV = adjIV[1]
    beta1.adjIV.sd = adjIV[2]
    beta1.adjIV.pval = adjIV[4]
    
    beta1 = 0   # Null Hypothesis -- H_0: \beta1 = 0
    gamma.Mle = sum(X*Z)/sum(Z^2)
    sigma2.Mle = sqrt(1/n*(sum(X^2) - 1/sum(Z^2)*(sum(X*Z))^2))
    print(c(gamma.Mle, sigma2.Mle))
    
    #pre-training for local optimization
    train.number = 200  #number of pre-training
    bach.number = 1000 #number of training sample  for every pre-train
    max.train.step = 100 # max iterative steps for training data
    total.maxstep = 200  # max` iterative steps for complete data
    
    
    #parallel computing for pre-training
    pre.likelihood = NULL  #report every pre-training outcome
    
    cl <- makeCluster(100)
    registerDoParallel(cl)
    pre.likelihood = foreach ( l = 1:train.number, .combine = rbind, .packages =  c('distrEx','mvtnorm', 'nloptr') ) %dopar%
    {
      #pre-training for local optimization
       #report every pre-training outcome
      
      set.seed(100+l)
      train.l = rep(FALSE, n)
      train.l[sample(n, bach.number)] = TRUE
      train.data.l = test.matrix[train.l, ]    
            
      pre.likelihood = NULL

      tryCatch({
        mle2.train = EM_algorithm2_train(runif(1, -5, 5),  gamma.Mle, runif(1, log(1),log(5)), log(sigma2.Mle),  runif(1, -0.3, 0.3), train.data.l, max.train.step)
        Partpara.mle = c(mle2.train[1], beta1, mle2.train[2], exp(mle2.train[3:4]), mle2.train[5])
        
        mle1.train = EM_algorithm1_train(mle2.train[1], runif(1, -5, 5), mle2.train[2], mle2.train[3], mle2.train[4], mle2.train[5], train.data.l, max.train.step)
         #mle1.train = EM_algorithm1_train(mle2.train[1], beta1, mle2.train[2], mle2.train[3], mle2.train[4], mle2.train[5], train.data.l, max.train.step)
        Allpara.mle = c(mle1.train[1:3], exp(mle1.train[4:5]), mle1.train[6])  
        
        pre.likelihood = c(like.value(Partpara.mle), like.value(Allpara.mle), Partpara.mle, Allpara.mle)

      }, warning = function(w){
        # 这里是出现warning状态时，应该怎么做，可以用print打印出来，可以执行其它命令
      }, error = function(e){ 
        # 这里时出现Error状态时，应该怎么做，可以用print打印出来，也可以执行其它命令
      },finally = {
        # 这是运行正常时，应该怎么做，可以用print打印出来，也可以执行其它命令
      })
    
     return(pre.likelihood)
    }
    
    stopImplicitCluster()
    stopCluster(cl)
    
    warmup.Partpara = pre.likelihood[order(-pre.likelihood[,1]), 3:8]
    warmup.Allpara = pre.likelihood[order(-pre.likelihood[,2]), 9:dim(pre.likelihood)[2]]
    
    #print(warmup.Partpara[1,])
    #print(warmup.Allpara[1,])
    #warmup.Partpara = pre.likelihood[which(pre.likelihood[,1] == max(pre.likelihood[,1])), 3:8]
    #warmup.Allpara = pre.likelihood[which(pre.likelihood[,2] == max(pre.likelihood[,2])), 9:dim(pre.likelihood)[2]]

    total.data = data.frame( Y, X, Z)
    names(total.data) = c("Y", "X", Zk)
    
   LRT = 0
    z = 1
    while (LRT == 0 && z < 5) { 

      tryCatch({
        
        mle2 = EM_algorithm2_train(warmup.Partpara[z,1], warmup.Partpara[z,3], log(warmup.Partpara[z,4]), log(warmup.Partpara[z,5]), warmup.Partpara[z,6],  total.data, total.maxstep)
        opt.Partpara = c(mle2[1], beta1, mle2[2], exp(mle2[3:4]), mle2[5])
        
        mle1 = EM_algorithm1_train(warmup.Allpara[z,1], warmup.Allpara[z,2], warmup.Allpara[z,3], log(warmup.Allpara[z,4]), log(warmup.Allpara[z,5]), warmup.Allpara[z,6], total.data, total.maxstep)
        opt.Allpara = c(mle1[1:3], exp(mle1[4:5]), mle1[6])
        
        LRT = -2*(like.value(opt.Partpara) - like.value(opt.Allpara))
        
        if(LRT < 0)
        {
          mle1 = EM_algorithm1_train(mle2[1], beta1, mle2[2], mle2[3], mle2[4], mle2[5], total.data, total.maxstep)
          opt.Allpara = c(mle1[1:3], exp(mle1[4:5]), mle1[6])
          LRT = -2*(like.value(opt.Partpara) - like.value(opt.Allpara))
        }
          
        }, warning = function(w){
        
      LRT = 0
      z = z + 1
        
        # 这里是出现warning状态时，应该怎么做，可以用print打印出来，可以执行其它命令
      }, error = function(e){ 

      LRT = 0
      z = z + 1

        # 这里时出现Error状态时，应该怎么做，可以用print打印出来，也可以执行其它命令
      },finally = {
       z  = z  + 1
        # 这是运行正常时，应该怎么做，可以用print打印出来，也可以执行其它命令
      })
}
   
    
    
    beta1.MLE = opt.Allpara[2]
    beta1.MLE.sd = sqrt(abs((beta1.MLE - beta1)^2/LRT))
    beta1.LRT.pval = 1 - pchisq(LRT, df = 1, ncp=0, lower.tail = TRUE, log.p = FALSE)
    #summary = matrix(c(Partpara1.mle,  Allpara.mle), 2,6, byrow =  TRUE)
    
    
    beta1.estimate = as.numeric(c(beta1.direc, beta1.twostage, beta1.adjIV, beta1.MLE))
    beta1.sd = as.numeric(c(beta1.direc.sd, beta1.twostage.sd, beta1.adjIV.sd, beta1.MLE.sd))
    beta1.pval = as.numeric(c(beta1.direc.pval, beta1.twostage.pval, beta1.adjIV.pval, beta1.LRT.pval))
    
    
    write.table(data.frame("SNP.location" = IV.name[k], "valid.sample" = 
                             valid.sample.number,  "EM.mle" = t(opt.Allpara),  "estimate" = t(beta1.estimate), "pvalue" = t(beta1.pval), "LRT" = LRT), file = "sbpa22.to.hypert05.RealData.txt", append = TRUE, sep = ",", col.names = NA, qmethod = "double")
  
  }
  
  
 





