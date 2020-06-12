## Server Wrapper Script
if(F){
  rm(list=ls())
  setwd("A:/ownCloud/SFA_CET/R_Code")
}

##Set System Up ----
##Get Process Starting Time ----
Sys.time()
timestart <- Sys.time()
sessionInfo()
sys <- Sys.info()

#load packages
library('parallel')
library("truncnorm")
library("pracma")
library("sn")
library("pbivnorm")
library("emg")
library("xtable")

#load code
source("Distribution_CET_20200505.R")

## Set Seed
set.seed(1337)

## Makes the functions run on server without change
if (sys[1] == "Windows") {
  numCores <- 1 # Number of cores used
  MC_reps <- 1 # Number of MC reps
} else {
  numCores <- 200 # Number of cores used
  MC_reps <- 1 # Number of MC reps
}


#Set parameters
mu =  c(-8, -4,-2,-1, 1, 2, 4, 8)
sigma_u = c(0.25,0.5,1,2,4)
sigma_v = c(0.25,0.5,1,2,4)
lambda = c(0.25,0.5,1,2,4, 8)
N = 1e8
quan=c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99)
u<-"exp"
u<-"tn"
benchmark = F

if(u!="exp"){
  par_grid_tn<-cbind(rep(rep(mu, each=length(sigma_u)*length(sigma_v)),each=MC_reps),
                     rep(rep(rep(sigma_u, each=length(sigma_v)),length(mu)),each=MC_reps),
                     rep(rep(sigma_v,length(sigma_u)*length(mu)),each=MC_reps),
                     rep(1:MC_reps, length(mu)*length(sigma_u)*length(sigma_v)))
} else {
  par_grid_exp<-cbind(rep(rep(lambda, each=length(sigma_v)),each=MC_reps), 
                      rep(rep(sigma_v,length(lambda)),each=MC_reps),
                      rep(1:MC_reps, length(lambda)*length(sigma_v)))
}



if(u!="exp"){
  results_tn<-mclapply(1:nrow(par_grid_tn), function(i)   sim_cdf(mu = par_grid_tn[i,1], sigma_u = par_grid_tn[i,2], sigma_v = par_grid_tn[i,3], lambda = c(1), u="tn", N = N, quan=quan, MC_reps=1, benchmark=benchmark), mc.cores = numCores)
  
  # save(results_tn, file = paste0("results_tn.Rdata")) #save results
  
  # load("results_tn.Rdata")
  # 
  if(benchmark){
    err_list<-list()
    rel<-F
    for(j in 1:4){
      err<-matrix(0, ncol=length(quan),nrow=nrow(par_grid_tn))
      for(i in 1:nrow(par_grid_tn)){
        # if(rel){
        #   err[i,]<-(results_tn[[i]][[3]]-quan)/quan
        #   err[i,][quan>0.5]<-((results_tn[[i]][[3]]-quan)/(1-quan))[quan>0.5]
        # } else {
          err[i,]<-results_tn[[i]][[j]]-quan#abs(results_tn[[i]][[j]]-quan)
        # }
          # err[i,]<-out
      }
      err_list[[j]]<-err
    }
  } else {
    err_list<-list()
    rel<-F
    # for(j in 1:4){
      err<-matrix(0, ncol=1,nrow=nrow(par_grid_tn))
      for(i in 1:nrow(par_grid_tn)){
        # if(rel){
        #   err[i,]<-(results_tn[[i]][[3]]-quan)/quan
        #   err[i,][quan>0.5]<-((results_tn[[i]][[3]]-quan)/(1-quan))[quan>0.5]
        # } else {
          err[i,]<-mean(results_tn[[i]][[1]]$time/results_tn[[i]][[3]]$time)#abs(results_tn[[i]][[j]]-quan)
        # }
          # err[i,]<-out
      }
      # err_list[[j]]<-err
  }
} else {
  results_exp<-mclapply(1:nrow(par_grid_exp), function(i) sim_cdf(mu = 1, sigma_u = 1, sigma_v = par_grid_exp[i,2], lambda = par_grid_exp[i,1], u="exp", N = N, quan=quan, MC_reps=1, benchmark=benchmark), mc.cores = numCores)
  
  # save(results_exp, file = paste0("results_exp_benchmark.Rdata")) #save results
  
  # load("results_exp.Rdata")
  # 
  if(benchmark){
    err_list<-list()
    rel<-F
    for(j in 1:4){
      err<-matrix(0, ncol=length(quan),nrow=nrow(par_grid_exp))
      for(i in 1:nrow(par_grid_exp)){
        # if(rel){
        #   err[i,]<-(results_tn[[i]][[3]]-quan)/quan
        #   err[i,][quan>0.5]<-((results_tn[[i]][[3]]-quan)/(1-quan))[quan>0.5]
        # } else {
        err[i,]<-results_exp[[i]][[j]]-quan#abs(results_tn[[i]][[j]]-quan)
        # }
        # err[i,]<-out
      }
      err_list[[j]]<-err
    }
  } else {
    err_list<-list()
    rel<-F
    # for(j in 1:4){
    err<-matrix(0, ncol=1,nrow=nrow(par_grid_exp))
    for(i in 1:nrow(par_grid_exp)){
      # if(rel){
      #   err[i,]<-(results_tn[[i]][[3]]-quan)/quan
      #   err[i,][quan>0.5]<-((results_tn[[i]][[3]]-quan)/(1-quan))[quan>0.5]
      # } else {
      err[i,]<-mean(results_exp[[i]][[1]]$time/results_exp[[i]][[3]]$time)#abs(results_tn[[i]][[j]]-quan)
      # }
      # err[i,]<-out
    }
    # err_list[[j]]<-err
  }
}

if(F){
  par(mfrow=c(1,2))
  boxplot(err_list[[2]],
          ylab=expression(paste("F"[epsilon],"(","Q"^"*","(p))","-p")),
          xlab="p",xaxt = 'n',
          main="Owen CDF")
  axis(side=1,at=1:9,labels=quan)
  
  boxplot(err_list[[3]], 
          ylab=expression(paste("F"[epsilon],"(","Q"^"*","(p))","-p")),
          xlab="p",xaxt = 'n',
          main="BvN CDF")
  axis(side=1,at=1:9,labels=quan)
  
  
  par(mfrow=c(1,2))
  boxplot(err_list[[2]],
          ylab=expression(paste("F"[epsilon],"(","Q"^"*","(p))","-p")),
          xlab="p",xaxt = 'n',
          main="Erf CDF")
  axis(side=1,at=1:9,labels=quan)
  boxplot(err_list[[3]], 
          ylab=expression(paste("F"[epsilon],"(","Q"^"*","(p))","-p")),
          xlab="p",xaxt = 'n',
          main="EMG CDF")
  axis(side=1,at=1:9,labels=quan)
  
  
  dat<-data.frame(cbind(par_grid_tn[,-4]),err)
  dat<-data.frame(cbind(par_grid_exp[,-3]),err)
  
  colnames(dat)<-c("mu","sigma_u","sigma_v",paste0("p=",quan))
  colnames(dat)<-c("lambda","sigma_v",paste0("p=",quan))
  
  print(xtable(dat, digits=c(0,0,1,2,rep(6,length(quan))), type = "latex"),include.rownames=F)
  
}