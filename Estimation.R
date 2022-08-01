rm(list = ls())

conlike1 = function(para1,para2)
{
  beta01 = para1[1]
  beta11 = para1[2]
  gamma1 = para1[3]
  logsigma11 = para1[4]
  logsigma21 = para1[5]
  rho1 = para1[6]
  alpha1 = para1[7]
  
  
  flink = 0
  
  for (i in 1:n)
  {
    
    K1 = exp(logsigma11 - logsigma21)*rho1*(X[i] - gamma1*Z[i])
    P_yx = function(r)
    {    
      exp((Y[i]-1)*(beta01 + X[i]*beta11 + K1 + sqrt(abs(1 - rho1^2))*exp(logsigma11)*alpha1*r) - 
            log(1 + exp(-(beta01 + X[i]*beta11 + K1 + sqrt(abs(1 - rho1^2))*exp(logsigma11)*alpha1*r))))
    }
    
    Jointlike = function(r) {P_yx(r)*P_u(r)}
    MargInteg2 =  integrate(Jointlike, -Inf, Inf)$value
    
      
    beta02 = para2[1]
    beta12 = para2[2]
    gamma2 = para2[3]
    logsigma12 = para2[4]
    logsigma22 = para2[5]
    rho2 = para2[6]
    alpha2 = para2[7]
    
    func = function(r)
    {
      K2 = exp(logsigma12 - logsigma22)*rho2*(X[i] - gamma2*Z[i])
      (Y[i]-1)*(beta02 + X[i]*beta12 + K2 + sqrt(abs(1 - rho2^2))*exp(logsigma12)*alpha2*r) - 
      log(1 + exp(-(beta02 + X[i]*beta12 + K2 + sqrt(abs(1 - rho2^2))*exp(logsigma12)*alpha2*r))) - 0.5*(X[i] - gamma2*Z[i])^2/exp(2*logsigma22) - logsigma22
    } 
    
    opt1 = function(r)  P_yx(r)*P_u(r)*func(r)
    #MargInteg1 = integrate(opt1, -10, 10)$value
    MargInteg1 = distrExIntegrate(opt1, -20, 20, subdivisions = 200, abs.tol = rel.tol)
    flink =  flink - MargInteg1/MargInteg2
  }  
  return(flink)
}


P_u = function(r)
{
  exp(-0.5*r^2)/sqrt(2*pi)
}


grad = function(para1, para2)
{
  gradbeta0 = 0
  gradbeta1 = 0
  gradgamma = 0
  gradlogsigma1 = 0
  gradlogsigma2 = 0
  gradrho = 0
  gradalpha = 0
  
  beta0 = para1[1]
  beta1 = para1[2]
  gamma = para1[3]
  logsigma1 = para1[4]
  logsigma2 = para1[5]
  rho = para1[6]
  alpha = para1[7]
  
  beta0.new = para2[1]
  beta1.new = para2[2]
  gamma.new = para2[3]
  logsigma1.new = para2[4]
  logsigma2.new = para2[5]
  rho.new = para2[6]
  alpha.new = para2[7]
  
  for (i in 1:n)
  {
    
    K = exp(logsigma1 - logsigma2)*rho*(X[i] - gamma*Z[i])
    J = function(r) beta0 + X[i]*beta1 + K + sqrt(1 - rho^2)*exp(logsigma1)*alpha*r
    P_yx = function(r) exp((Y[i] - 1)*J(r) - log(1 + exp(-J(r))))
    tmp = function(r) P_yx(r)*P_u(r)
    I = integrate(tmp, -Inf, Inf)$value
    
    
    K.new = exp(logsigma1.new - logsigma2.new)*rho.new*(X[i] - gamma.new*Z[i])
    J.new = function(r)  beta0.new + X[i]*beta1.new + K.new + sqrt(1 - rho.new^2)*exp(logsigma1.new)*alpha.new*r
    
    #compute derivitives
    tmp1 = function(r) tmp(r)/(1 + exp(-J.new(r)))
    gradbeta0 = gradbeta0 + Y[i] - integrate(tmp1, -Inf, Inf)$value/I
    
    tmp2 = function(r) X[i]*tmp(r)/(1 + exp(-J.new(r)))
    gradbeta1 = gradbeta1 + Y[i]*X[i] - integrate(tmp2, -Inf, Inf)$value/I
    
    tmp3 = function(r) - exp(logsigma1.new - logsigma2.new)*rho.new*Z[i]*tmp(r)/(1 + exp(-J.new(r)))
    gradgamma = gradgamma - exp(logsigma1.new - logsigma2.new)*rho.new*Z[i]*Y[i] +(X[i] - gamma.new*Z[i])/exp(2*logsigma2.new)*Z[i] - integrate(tmp3, -Inf, Inf)$value/I 
    
    tmp4 = function(r) (Y[i] - 1/(1 + exp(-J.new(r))))*(rho.new*exp(logsigma1.new - logsigma2.new)*(X[i] - gamma.new*Z[i]) + exp(logsigma1.new)*sqrt(1-rho.new^2)*alpha.new*r)*tmp(r)
    gradlogsigma1 = gradlogsigma1 + integrate(tmp4, -Inf, Inf)$value/I
    
    
    tmp5 = function(r) -rho.new*exp(logsigma1.new - logsigma2.new)*(X[i] - gamma.new*Z[i])*tmp(r)/(1 + exp(-J.new(r)))
    gradlogsigma2 = gradlogsigma2 - Y[i]*rho.new*exp(logsigma1.new - logsigma2.new)*(X[i] - gamma.new*Z[i])- 1 + (X[i] - gamma.new*Z[i])^2*exp(-2*logsigma2.new) - integrate(tmp5, -Inf, Inf)$value/I
    
    tmp6 = function(r) (Y[i] - 1/(1 + exp(-J.new(r))))*(exp(logsigma1.new - logsigma2.new)*(X[i] - gamma.new*Z[i]) - rho.new*exp(logsigma1.new)*alpha.new*r/sqrt(1-rho.new^2))*tmp(r)
    gradrho = gradrho + integrate(tmp6, -Inf, Inf)$value/I
    
    tmp7 = function(r) (Y[i] - 1/(1 + exp(-J.new(r))))*( exp(logsigma1.new)*sqrt(1-rho.new^2)*r)*tmp(r)
    gradalpha = gradalpha + integrate(tmp7, -Inf, Inf)$value/I
    
  }
  grad = -c(gradbeta0 ,gradbeta1, gradgamma ,gradlogsigma1, gradlogsigma2 ,gradrho, gradalpha)
  return(grad)
}



EM_algorithm = function(beta00, beta10, gamma0, sigma10, sigma20, rho0, maxit = 100)
{
  alpha0 = 1
  par.old <<- c(beta00, beta10, gamma0, sigma10, sigma20, rho0, alpha0)
  PostExp = function(para)
  {
    postexp = conlike1(par.old, para) 
    postexp
  }
  grad.iter = function(para)
  {
    grad(par.old, para)
  }
  par.new = optim(par = par.old, fn = PostExp, grad.iter, method =  "L-BFGS-B",
                  lower = c(beta0 - 5, beta1 -5, gamma - 1, -Inf, -Inf, -0.99, 0.5), 
                  upper = c(beta0 + 5, beta1 + 5, gamma + 1, Inf,  Inf, 0.99, Inf))$par
  iter = 1
  while (sqrt(sum((par.new - par.old)^2)) > 1e-4 & iter < maxit) 
  { 
    sigma10 = par.new[4]*par.new[7]
    alpha0 = 1
    par.old <<- c(par.new[1:3], sigma10, par.new[5:6], alpha0)
    PostExp = function(para)
    {
      postexp = conlike1(par.old, para) 
      postexp
    }
    grad.iter = function(para)
    {
      grad(par.old, para)
    }
    par.new = optim(par = par.old, fn = PostExp, grad.iter , method =  "L-BFGS-B",
                    lower = c(beta0 - 5, beta1 -5, gamma - 1, -Inf,- Inf, -0.99, 0.5), 
                    upper = c(beta0 + 5, beta1 + 5, gamma + 1, Inf,  Inf, 0.99, Inf))$par 
    iter = iter + 1
  }

  if (sqrt(sum((par.new - par.old)^2)) <= 1e-4)
  {conver = 1} else
  {conver = 0}
  #write.table(conver , file = "conver2.txt", append = TRUE, sep = ",", col.names = NA, qmethod = "double")
  return(c(par.new[1:6], conver, sqrt(sum((par.new - par.old)^2))))
}





library(mvtnorm)
library(foreach)
library(doParallel)
#library(pacman)
#library(Rsolnp)
library(distrEx)


#set.seed(12)
n = 1000
maf = 0.3
beta0 = 2
beta1 =  -1
gamma = 1


sigma1 = 1
sigma2 = 2
rho = 0.8

Z = rbinom(n, 2, maf)


cl <- makeCluster(100)
registerDoParallel(cl)

Paraout = matrix(rep(0,4000),1000,4)
Paraout = foreach ( k = 1:1000, .combine = rbind, .packages =  c('distrEx','mvtnorm') ) %dopar%
  {
     set.seed(k)
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

    #direct type estimator
    y_x = data.frame(y = Y, x = X)
    beta1.direc = glm(y~x, family = binomial(), data = y_x)$coefficient[2]

    #two stage estimator
    y_xhat = data.frame(y = Y, xhat = Xhat)
    beta1.twostage = glm(y ~ xhat, family = binomial(), data = y_xhat)$coefficients[2]

    #adjusted IV
    Vhat = X - Xhat
    y_xhat_vhat = data.frame(y = Y, xhat = Xhat, vhat = Vhat)
    beta1.adjIV = glm(y ~ xhat + vhat, family = binomial(), data = y_xhat_vhat)$coefficients[2]

    mle = EM_algorithm(beta0 - 1, beta1 - 0.5, gamma - 0.2, log(sigma1 - 0.5), runif(1, log(sigma2 - 0.3), log(sigma2+ 0.1)), rho -0.2, 100)   
    mle1 = c(mle[1:3], exp(mle[4:5]), mle[6])
    beta1.estimate = as.numeric(c(beta1.direc,beta1.twostage,beta1.adjIV, mle1[2]))

    write.table(data.frame(k, beta1.direc,beta1.twostage,beta1.adjIV, mle1[2],mle[7], mle[8]) , file = "Mle.beta-1rho0.8gam1sig1equ1sig2equ2.txt", append = TRUE, sep = ",", col.names = F, qmethod = "double")
    return(beta1.estimate)
  }

stopImplicitCluster()
stopCluster(cl)




















