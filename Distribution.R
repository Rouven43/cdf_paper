dnormexp <- function(x, lambda = 1, sigma_v = 1, log = FALSE, N=10000){
  
  out <-lambda*exp(lambda*x+(sigma_v^2*lambda^2)/2)*pnorm(-x/sigma_v-lambda*sigma_v)
  
  return(out)
}

pnormexp <- function(x, lambda = 1, sigma_v = 1, meth="num", log = FALSE, N=10000, E=NULL){
  out<-matrix(NA,nrow=length(x),ncol=1)
  if(any(abs(x)==Inf)){
    out[x == Inf,] <- 1
    out[x == -Inf,] <- 0
    xx<-x[abs(x)!=Inf,]
  } else {
    xx<-x
  }
  if(meth == "num"){
    value <- unlist(
      lapply(1:length(xx),
             function(i) tryCatch(pracma::quadinf(f=dnormexp, xa = -Inf, xb= xx[i], lambda=lambda,  sigma_v = sigma_v)$Q,
                                  error=function(e) NA )))
  }
  
  if(meth=="ana"){
    value <-1/2*(erf(xx/(sqrt(2)*sigma_v))+exp(1/2*lambda^2*sigma_v^2+lambda*xx)*erfc((lambda*sigma_v+xx/sigma_v)/sqrt(2)))+1/2
  }
  
  if(meth=="ana2"){
    a<--lambda*sigma_v
    b<-1
    y<--(xx+lambda*sigma_v^2)/(sigma_v) 
    value <- -lambda*sigma_v*exp((sigma_v^2*lambda^2)/2+lambda^2*sigma_v^2)*(1/a*exp(a*y)*pnorm(b*y)-1/a*exp(a^2/(2*b^2))*pnorm(b*y-a/b))
  }
  
  if(meth=="emg"){
    value <-1-pemg(-xx,mu=0,sigma=sigma_v,lambda=lambda, log = log)
  }
  
  if(meth=="empir"){
    if(is.null(E)){
      U<-rtruncnorm(N, a = 0, b = Inf, mean = mu, sd = sigma_u)
      print(paste0(mu," + ", sigma_u, " + ", sigma_v))
      V<-rnorm(N, 0, sigma_v)
      E<-V-U
    }
    value<-ecdf(E)(xx)
  }
  
  if(any(abs(x)==Inf)){
    out[abs(x)!=Inf,]<-value
  } else {
    out<-value 
  }
  
  out[out>1]<-1
  out[out<0]<-0
  return(out)
}

dnormtnorm <- function(x, mu = 1, sigma_u = 1, sigma_v = 1, log = FALSE, N=10000){
  
  out <- 1 / (sqrt(sigma_v^2 + sigma_u^2) * pnorm(mu / sigma_u)) * 
    dnorm((x + mu) / (sqrt(sigma_v^2 + sigma_u^2))) * pnorm((mu * sigma_v^2 - x * 
        sigma_u^2) / (sqrt(sigma_v^2 + sigma_u^2) * sigma_v * sigma_u))
  
  return(out)
}


gfun <- function(nomin,denomin){
 
  if(sign(nomin) / sign(denomin) < 0){
    out <- - 1 / 2
  } else {
    out <- 0
  }

  return(out)
}


pnormtnorm <- function(x, mu = 1, sigma_u = 1, sigma_v = 1, meth="ana", N=10000, E=NULL, log = FALSE){
  # if(x == Inf){
    out<-matrix(NA,nrow=length(x),ncol=1)
    if(any(abs(x)==Inf)){
      out[x == Inf,] <- 1
      out[x == -Inf,] <- 0
      xx<-x[abs(x)!=Inf,]
    } else {
      xx<-x
    }
    a <- mu * sqrt(sigma_v^2 + sigma_u^2) / (sigma_v * sigma_u)
    b <- - sigma_u / sigma_v
    y <- (xx + mu) / sqrt(sigma_v^2 + sigma_u^2)
    
    if(meth == "num"){
      value <- unlist(
        lapply(1:length(xx),
          function(i) tryCatch(pracma::quadinf(f=dnormtnorm, xa = -Inf, xb= xx[i], mu = mu, sigma_u = sigma_u, sigma_v = sigma_v)$Q,
            error=function(e) NA )))
    }
    
  if(meth=="ana"){
    value <- unlist(
        lapply(1:length(xx),
            function(i) tryCatch( 1 / pnorm(mu / sigma_u) * (1 / 2 * pnorm(y[i]) + 1 / 2 * pnorm(a / sqrt(1 + b^2)) + gfun(a, y[i] * sqrt(1 + b^2)) -# ifelse(sign(a/(y[i] * sqrt(1 + b^2)))<0,-1/2,0)- #gfun(a, y[i] * sqrt(1 + b^2)) -
              T.Owen(y[i], (a + b * y[i]) / y[i]) - T.Owen(a / sqrt(1 + b^2), (a * b + y[i] * (1 + b^2)) / a)),
              error=function(e) NA )))
    }
  
    if(meth=="bvn"){
      value <- unlist(
        lapply(1:length(xx),
               function(i) tryCatch( 1 / pnorm(mu / sigma_u) * pbivnorm(x = a / sqrt(1 + b^2), y = y[i], rho = -b / sqrt(1 + b^2), recycle = TRUE), 
                error=function(e) NA )))
    } 
    
    if(meth=="empir"){
      if(is.null(E)){
        U<-rtruncnorm(N, a = 0, b = Inf, mean = mu, sd = sigma_u)
        V<-rnorm(N, 0, sigma_v)
        E<-V-U
      }
      value<-ecdf(E)(xx)
    }

    if(any(abs(x)==Inf)){
      out[abs(x)!=Inf,]<-value
    } else {
     out<-value 
    }
    
    out[out>1]<-1
    out[out<0]<-0
    
  return(out)
} 

sim_cdf<-function(mu = c(1), sigma_u = c(1), sigma_v = c(1), lambda = c(1), u="tn", N = 1000, quan=c(0.01,0.05,0.1,0.5,0.9,0.95,0.99), MC_reps=1, benchmark){
  
  if(u!="exp"){
    E<-rnorm(N, 0, sd=sigma_v)-rtruncnorm(N, a = 0, b = Inf, mean = mu, sd = sigma_u)
  } else {
    E<-rnorm(N, 0, sd=sigma_v)-rexp(N, rate=lambda)
  }
  
  eval_grid<-quantile(E,quan)
  
  if(benchmark==T){
    if(u!="exp"){
      num_cdf  <- microbenchmark(pnormtnorm(x = eval_grid, mu = mu, sigma_u = sigma_u, sigma_v = sigma_v, meth="num",  N=N, E=E, log = FALSE))
      
      ana_cdf <- microbenchmark(pnormtnorm(x = eval_grid, mu = mu, sigma_u = sigma_u, sigma_v = sigma_v, meth="ana",  N=N, E=E, log = FALSE))
      
      bvn_cdf  <- microbenchmark(pnormtnorm(x = eval_grid, mu = mu, sigma_u = sigma_u, sigma_v = sigma_v, meth="bvn",  N=N, E=E, log = FALSE))
      
      empir_cdf <- microbenchmark(pnormtnorm(x = eval_grid, mu = mu, sigma_u = sigma_u, sigma_v = sigma_v, meth="empir",  N=N, E=E, log = FALSE))
    } else {
      num_cdf  <- microbenchmark(pnormexp(x = eval_grid, lambda = lambda, sigma_v = sigma_v, meth="num",  N=N, E=E, log = FALSE))
      
      ana_cdf <- microbenchmark(pnormexp(x = eval_grid, lambda = lambda, sigma_v = sigma_v, meth="ana",  N=N, E=E, log = FALSE))
      
      emg_cdf  <- microbenchmark(pnormexp(x = eval_grid, lambda = lambda, sigma_v = sigma_v, meth="emg",  N=N, E=E, log = FALSE))
      
      empir_cdf <- microbenchmark(pnormexp(x = eval_grid, lambda = lambda, sigma_v = sigma_v, meth="empir",  N=N, E=E, log = FALSE))
    }
    
  } else {
    if(u!="exp"){
      num_cdf  <- pnormtnorm(x = eval_grid, mu = mu, sigma_u = sigma_u, sigma_v = sigma_v, meth="num",  N=N, E=E, log = FALSE)
      
      ana_cdf <- pnormtnorm(x = eval_grid, mu = mu, sigma_u = sigma_u, sigma_v = sigma_v, meth="ana",  N=N, E=E, log = FALSE)
      
      bvn_cdf  <- pnormtnorm(x = eval_grid, mu = mu, sigma_u = sigma_u, sigma_v = sigma_v, meth="bvn",  N=N, E=E, log = FALSE)
      
      empir_cdf <- pnormtnorm(x = eval_grid, mu = mu, sigma_u = sigma_u, sigma_v = sigma_v, meth="empir",  N=N, E=E, log = FALSE)
    } else {
      num_cdf  <- pnormexp(x = eval_grid, lambda = lambda, sigma_v = sigma_v, meth="num",  N=N, E=E, log = FALSE)
      
      ana_cdf <- pnormexp(x = eval_grid, lambda = lambda, sigma_v = sigma_v, meth="ana",  N=N, E=E, log = FALSE)
      
      emg_cdf  <- pnormexp(x = eval_grid, lambda = lambda, sigma_v = sigma_v, meth="emg",  N=N, E=E, log = FALSE)
      
      empir_cdf <- pnormexp(x = eval_grid, lambda = lambda, sigma_v = sigma_v, meth="empir",  N=N, E=E, log = FALSE)
    }
  }
  if(u!="exp"){
    out<-list(num_cdf,ana_cdf,bvn_cdf,empir_cdf,eval_grid)
  } else {
    out<-list(num_cdf,ana_cdf,emg_cdf,empir_cdf,eval_grid)
  }
  return(out)
} 
