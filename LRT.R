
rm(list = ls())

#In the full parameter space to compute the Mle
conlike1 = function(para1,para2)
{
  beta01 = para1[1]
  beta11 = para1[2]
  gamma1 = para1[3]
  logsigma11 = para1[4]
  logsigma21 = para1[5]
  rho1 = para1[6]
  
  flink = 0
  
  for (i in 1:n)
  {
    
    K1 = exp(logsigma11 - logsigma21)*rho1*(X[i] - gamma1*Z[i])
    P_yx = function(r)
    {    
      exp((Y[i]-1)*(beta01 + X[i]*beta11 + K1 + sqrt(abs(1 - rho1^2))*exp(logsigma11)*r) - 
            log(1 + exp(-(beta01 + X[i]*beta11 + K1 + sqrt(abs(1 - rho1^2))*exp(logsigma11)*r))))
    }
    
    Jointlike = function(r) {P_yx(r)*P_u(r)}
    MargInteg2 =  integrate(Jointlike, -Inf, Inf)$value
    
      
    beta02 = para2[1]
    beta12 = para2[2]
    gamma2 = para2[3]
    logsigma12 = para2[4]
    logsigma22 = para2[5]
    rho2 = para2[6]
    
    func = function(r)
    {
      K2 = exp(logsigma12 - logsigma22)*rho2*(X[i] - gamma2*Z[i])
      (Y[i]-1)*(beta02 + X[i]*beta12 + K2 + sqrt(abs(1 - rho2^2))*exp(logsigma12)*r) - 
      log(1 + exp(-(beta02 + X[i]*beta12 + K2 + sqrt(abs(1 - rho2^2))*exp(logsigma12)*r))) - 0.5*(X[i] - gamma2*Z[i])^2/exp(2*logsigma22) - logsigma22
    } 
    
    opt1 = function(r)  P_yx(r)*P_u(r)*func(r)
    #MargInteg1 = integrate(opt1, -10, 10)$value
    MargInteg1 = distrExIntegrate(opt1, -10, 10, subdivisions = 100, abs.tol = rel.tol)
    flink =  flink - MargInteg1/MargInteg2
  }  
  return(flink)
}


P_u = function(r)
{
  exp(-0.5*r^2)/sqrt(2*pi)
}


#All parameter gradient
grad1 = function(para1, para2)
{
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
  
  
  for (i in 1:n)
  {
    
    K = exp(logsigma1 - logsigma2)*rho*(X[i] - gamma*Z[i])
    J = function(r) beta0 + X[i]*beta1 + K + sqrt(1 - rho^2)*exp(logsigma1)*r
    P_yx = function(r) exp((Y[i] - 1)*J(r) - log(1 + exp(-J(r))))
    tmp = function(r) P_yx(r)*P_u(r)
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

#fix beta1 gradient
grad2 = function(para1, para2)
{
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
  
  
  for (i in 1:n)
  {
    
    K = exp(logsigma1 - logsigma2)*rho*(X[i] - gamma*Z[i])
    J = function(r) beta0 + X[i]*beta1 + K + sqrt(1 - rho^2)*exp(logsigma1)*r
    P_yx = function(r) exp((Y[i] - 1)*J(r) - log(1 + exp(-J(r))))
    tmp = function(r) P_yx(r)*P_u(r)
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


conlike2 = function(para1,para2)
{
  beta01 = para1[1]
  gamma1 = para1[2]
  logsigma11 = para1[3]
  logsigma21 = para1[4]
  rho1 = para1[5]
  
  beta11 = beta1
  beta12 = beta1
  
  flink = 0
  
  for (i in 1:n)
  {
    
    K1 = exp(logsigma11 - logsigma21)*rho1*(X[i] - gamma1*Z[i])
    P_yx = function(r)
    {    
      exp((Y[i]-1)*(beta01 + X[i]*beta11 + K1 + sqrt(abs(1 - rho1^2))*exp(logsigma11)*r) - 
            log(1 + exp(-(beta01 + X[i]*beta11 + K1 + sqrt(abs(1 - rho1^2))*exp(logsigma11)*r))))
    }
    
    Jointlike = function(r) {P_yx(r)*P_u(r)}
    MargInteg2 =  integrate(Jointlike, -Inf, Inf)$value
    
      
    beta02 = para2[1]
    gamma2 = para2[2]
    logsigma12 = para2[3]
    logsigma22 = para2[4]
    rho2 = para2[5]
    
    func = function(r)
    {
      K2 = exp(logsigma12 - logsigma22)*rho2*(X[i] - gamma2*Z[i])
      (Y[i]-1)*(beta02 + X[i]*beta12 + K2 + sqrt(abs(1 - rho2^2))*exp(logsigma12)*r) - 
      log(1 + exp(-(beta02 + X[i]*beta12 + K2 + sqrt(abs(1 - rho2^2))*exp(logsigma12)*r))) - 0.5*(X[i] - gamma2*Z[i])^2/exp(2*logsigma22) - logsigma22
    } 
    
    opt1 = function(r)  P_yx(r)*P_u(r)*func(r)
    #MargInteg1 = integrate(opt1, -10, 10)$value
    MargInteg1 = distrExIntegrate(opt1, -20, 20, subdivisions = 200, abs.tol = rel.tol)
    flink =  flink - MargInteg1/MargInteg2
  }  
  return(flink)
}




#EM algorithm for all  parameters

EM_algorithm1 = function(beta00, beta10, gamma0, sigma10, sigma20, rho0, maxit = 20)
{
  
  par.old <<- c(beta00, beta10, gamma0, sigma10, sigma20, rho0)
  PostExp = function(para)
  {
    postexp = conlike1(par.old, para) 
    postexp
  }
  grad.iter = function(para)
  {
    grad1(par.old, para)
  }
  
  par.new = optim(par = par.old, fn = PostExp, grad.iter, method =  "L-BFGS-B",
                  lower = c(beta0 - 5, beta1 -5, gamma - 1, -Inf, -Inf, -0.99), 
                  upper = c(beta0 + 5, beta1 + 5, gamma + 1, Inf,  Inf, 0.99))$par
  iter = 1
  while (sqrt(sum((par.new - par.old)^2)) > 1e-4 & iter < maxit) 
  { 
    par.old <<-  par.new
    PostExp = function(para)
    {
      postexp = conlike1(par.old, para) 
      postexp
    }
    grad.iter = function(para)
    {
      grad1(par.old, para)
    }
    par.new = optim(par = par.old, fn = PostExp, grad.iter , method =  "L-BFGS-B",
                    lower = c(beta0 - 5, beta1 -5, gamma - 1, -Inf,- Inf, -0.99), 
                    upper = c(beta0 + 5, beta1 + 5, gamma + 1, Inf,  Inf, 0.99))$par 
    iter = iter + 1

  }
  
  if (sqrt(sum((par.new - par.old)^2)) <= 1e-4)
  {conver1 = 1} else
  {conver1 = 0}
  return(c(par.new[1:6], conver1, sqrt(sum((par.new - par.old)^2))))

}


#EM algorithm for fixed beta1 = beta1.true.

EM_algorithm2 = function(beta00, gamma0, sigma10, sigma20, rho0, maxit = 100)
{
  
  par.old <<- c(beta00, gamma0, sigma10, sigma20, rho0)
  PostExp = function(para)
  {
    postexp = conlike2(par.old, para) 
    postexp
  }
  grad.iter = function(para)
  {
    grad2(par.old, para)
  }
  
  par.new = optim(par = par.old, fn = PostExp, grad.iter, method =  "L-BFGS-B",
                  lower = c(beta0 - 5, gamma - 1, -Inf, -Inf, -0.99), 
                  upper = c(beta0 + 5, gamma + 1, Inf,  Inf, 0.99))$par
  iter = 1
  while (sqrt(sum((par.new - par.old)^2)) > 1e-4 & iter < maxit) 
  { 
    par.old <<-  par.new
    PostExp = function(para)
    {
      postexp = conlike2(par.old, para) 
      postexp
    }
    grad.iter = function(para)
    {
      grad2(par.old, para)
    }
    par.new = optim(par = par.old, fn = PostExp, grad.iter , method =  "L-BFGS-B",
                    lower = c(beta0 - 5, gamma - 1, -Inf,- Inf, -0.99), 
                    upper = c(beta0 + 5, gamma + 1, Inf,  Inf, 0.99))$par 
    iter = iter + 1

  }
  
  if (sqrt(sum((par.new - par.old)^2)) <= 1e-4)
  {conver2 = 1} else
  {conver2 = 0}
  return(c(par.new[1:5], conver2, sqrt(sum((par.new - par.old)^2))))

}


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
      exp(Y[i]*(beta0 + X[i]*beta1 + r) - log(1 + exp(beta0 + X[i]*beta1 + r)) - 
            log(sigma2) - 0.5*log(1 - rho^2)  -
            0.5*log(2*pi) -
            (X[i] - gamma*Z[i] - rho*sigma2/sigma1*r)^2/(2*sigma2^2*(1 - rho^2)) - 
            0.5*log(2*pi) -log(sigma1) - 0.5*r^2/sigma1^2)
      
    } 
    f = f + log(integrate(likelihood, -Inf, Inf)$value)
  }
    return(f)
}




library(mvtnorm)
library(foreach)
library(doParallel)
library(distrEx)



#set.seed(123)
n = 1000
maf = 0.3
beta0 = 2
beta1 =  0
gamma = 1


sigma1 = 4
sigma2 = 4
rho = 0.8

Z = rbinom(n, 2, maf)
out = matrix(rep(0,12000),1000,12)


cl <- makeCluster(100)
registerDoParallel(cl)

out = foreach ( k = 1:1000, .combine = rbind, .packages =  c('distrEx','mvtnorm') ) %dopar%
{
  #set.seed(k)

  C = rmvnorm(n, c(0, 0), sigma = matrix(c(sigma1^2, rho*sigma1*sigma2, rho*sigma1*sigma2, sigma2^2), 2, 2))
  U = C[, 1]
  V = C[, 2]
  X = Z*gamma  +  V
  logitp = beta0 + X*beta1 + U
  p = 1/(1 + exp(-logitp))
  Y=rbinom(n, 1, p)
  
  x_z = data.frame(x = X, z = Z)
  xz.coeff = lm(x ~ -1 + z, data = x_z)$coefficients
  Xhat = Z*xz.coeff
  
  y_z =data.frame(y = Y, z = Z)
  yz.coeff = glm(y ~ z, family = binomial(), data = y_z)$coefficients[2]
  
  
  #directed estimator
  y_x = data.frame(y = Y, x = X)
  direc = coef(summary(glm(y~x, family = binomial(), data = y_x)))[2,]
  beta1.direc = direc[1]
  beta1.direc.pval = direc[4]
  
  
  #two stage estimator
  y_xhat = data.frame(y = Y, xhat = Xhat)
  twostage = coef(summary(glm(y ~ xhat, family = binomial(), data = y_xhat)))[2,]
  beta1.twostage = twostage[1]
  beta1.twostage.pval = twostage[4]
  
  
  #adjusted IV
  Vhat = X - Xhat
  y_xhat_vhat = data.frame(y = Y, xhat = Xhat, vhat = Vhat)
  beta1.adjIV = glm(y ~ xhat + vhat, family = binomial(), data = y_xhat_vhat)$coefficients[2]
  adjIV = coef(summary(glm(y ~ xhat + vhat, family = binomial(), data = y_xhat_vhat)))[2,]
  beta1.adjIV = adjIV[1]
  beta1.adjIV.pval = adjIV[4]
  
  mle2 = EM_algorithm2(runif(1, beta0 - 0.5, beta0+ 0.5),  runif(1, gamma - 0.2, gamma + 0.2), log(sigma1+ 0.3), runif(1, log(sigma2 - 0.1), log(sigma2 + 0.2)), runif(1, rho-0.1, rho), 100)
  Partpara.mle = c(mle2[1:2], exp(mle2[3:4]), mle2[5])
  Partpara1.mle = c(Partpara.mle[1], beta1, Partpara.mle[2:5])
    res2 = mle2[6:7]
  
  #Allpara.Mle = EM_algorithm(Propara1.Mle[1], Propara1.Mle[2], Propara1.Mle[3], Propara1.Mle[4], Propara1.Mle[5], Propara1.Mle[6],  50)
  mle1 = EM_algorithm1(runif(1, beta0 - 0.5, beta0+ 0.5), beta1, runif(1, gamma - 0.2, gamma + 0.2), log(sigma1+ 0.3), runif(1, log(sigma2 - 0.1), log(sigma2 + 0.2)), runif(1, rho-0.1, rho), 100)
  Allpara.mle = c(mle1[1:3], exp(mle1[4:5]), mle1[6])
  beta1.MLE = Allpara.mle[2]
  res1 = mle1[7:8]

  LRT = -2*(like.value(Partpara1.mle) - like.value(Allpara.mle))
  beta1.LRT.pval = 1 - pchisq(LRT, df = 1, ncp=0, lower.tail = TRUE, log.p = FALSE)

  
  beta1.estimate = as.numeric(c(beta1.direc, beta1.twostage, beta1.adjIV, beta1.MLE))
  beta1.pval = as.numeric(c(beta1.direc.pval, beta1.twostage.pval, beta1.adjIV.pval, beta1.LRT.pval))
  
  write.table(list("Times " = k, "estimate" = beta1.estimate, "pval" = beta1.pval", file = “LRT.txt", append = TRUE, sep = ",", col.names = NA, qmethod = "double")
  return(c(beta1.estimate, beta1.sd, beta1.pval))
}

stopImplicitCluster()
stopCluster(cl)










