#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
dirname <- args[1]
gc()
dir <- "/home/projects/feariel/mishash/Development/ThreeFold/" 
dirname <- "addGammaSum"
.libPaths("/home/projects/feariel/mishash/R/x86_64-pc-linux-gnu-library/4.4")
if (1==0)
{
  system("rm -r /home/projects/feariel/mishash/R/x86_64-pc-linux-gnu-library/4.4/armaInv/")
  setwd("/home/projects/feariel/mishash/Development/ThreeFold/")
  system("rm -r /home/projects/feariel/mishash/Development/ThreeFold/armaInv/")
  RcppArmadillo::RcppArmadillo.package.skeleton(name = "armaInv")
  # copy the cpp file to /home/projects/feariel/mishash/Development/ThreeFold/armaInv/src/  !!!!!!
  system("cp /home/projects/feariel/mishash/Development/ThreeFold/src/armaInv.cpp /home/projects/feariel/mishash/Development/ThreeFold/armaInv/src/")
  Rcpp::compileAttributes("/home/projects/feariel/mishash/Development/ThreeFold/armaInv/")
  install.packages("/home/projects/feariel/mishash/Development/ThreeFold/armaInv/", repos = NULL, type = "source",force=TRUE, INSTALL_opts = '--no-lock')
  library(armaInv, lib.loc = "/home/projects/feariel/mishash/R/x86_64-pc-linux-gnu-library/4.4")
}

prefactor <- function(rhoM)
{
  n <- nrow(rhoM)
  if (n==2) {return(sum(rhoM))}
  else
  {
    A <- 0
    for (j in 1:ncol(rhoM))
    {
      foo <- prefactor(rhoM[-j,-j])
      for (i in 1:nrow(rhoM))
      {
        A <- A + rhoM[i,j]*foo
      }
    }
    return(A)
  }
}


## Intro ----
# Simulation <- "Simulation/" # for simuation
Simulation <- "" # for empirical
print(dirname)
set.seed(Sys.time())
rTh <- 1000
rB <- 10^seq(log10(rTh-0.5),6,0.05) 
rV <- 0.5*(rB[-length(rB)]+rB[-1])
dr <- rV*0.01
drm2 <- matrix(1/dr^2,nrow=1)
ncores <- as.numeric(Sys.getenv("LSB_DJOB_NUMPROC")) 
print(paste0(" found ", ncores, " cores"))
IndR <- which(rV<10^40.5)
rVsim <- 10^seq(0,7,0.05)
rhoD <- gammaD <- data.frame()
print("here0")
# system("export OPENBLAS_NUM_THREADS=1") # do before starting R !
library(circlize)
library(foreach) 
print("here0.6")
# library(matlib)
library(iterators)
library(doParallel)
require(intervals)
library(sets)
require(data.table)
library(matrixStats)
library(plotrix)
# library(rhdf5)
library(ggplot2)
library(MASS)
library(gtools)
library(expm)  
library(GA)
library(optimParallel)
# library(RcppArmadillo)
# library(Rcpp)
library(resample)
library(mvtnorm)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(purrr)
library(RColorBrewer)
library(stringr)
library(parallel)
library(rray)

.libPaths("/home/projects/feariel/mishash/R/x86_64-pc-linux-gnu-library/4.4")


if (is.na(ncores)) {ncores <- 2}
print(paste0("Starting cluster with ", ncores, " cores"))
cl <- parallel::makeCluster(ncores, type = "PSOCK")
print("generated the cluster")
registerDoParallel(cl)
setDefaultCluster(cl=cl)
clusterEvalQ(cl, library(foreach));
clusterEvalQ(cl, library(rray));

print("Done starting cluster")



# auxillary function----

normalizeMLD1 <- function(MLD1,MLD,IndR,getL1=FALSE)
{
  L1 <- 0
  count <- 0
  for (n in 2:length(MLD1))
  {
    for (s in 1:nrow(MLD[[n]])) 
    {
      if (any(MLD[[n]][s,]!=0)) 
      {
        count <- count+1
        Ind <- which(MLD[[n]][s,IndR]>0 & MLD1[[n]][s,IndR]>0)
        L1 <- L1 + mean(log10(MLD[[n]][s,IndR][Ind]/MLD1[[n]][s,IndR][Ind]))
      }
    }
  }
  L1 <- 10^(L1/count); #L1 <- 1; # here I remove L1 dependence
  for (n in 2:length(MLD1))
  {
    MLD1[[n]] <- L1*MLD1[[n]]
  }
  if (!getL1) {return(MLD1)}
  else {return(list(MLD1,L1))}
}

distMLDmoments <- function(MLD1,MLD,IndR)
{
  FigureOfMerit <- 0
  for (n in 2:min(c(length(MLD1),length(MLD))))
  {
    for (s in 1:nrow(MLD[[n]]))
    {
      for (m in seq(0,n-1,1))
      {
        mom <- sum(diff(rB)[IndR]*rV[IndR]^m*(MLD[[n]][s,]))
        mom1 <- sum(diff(rB)[IndR]*rV[IndR]^m*(MLD1[[n]][s,]))
        FigureOfMerit <- FigureOfMerit + (mom-mom1)^2/mom^2/n
      }
    }
  }
  return(FigureOfMerit)
}

distMLD <- function(MLD1,MLD,IndR)
{
  FigureOfMerit <- 0
  for (n in 2:8)
  {
    # Ind <- which(MLD1[[n]] > 0 & MLD[[n]] > 0)
    for (s in 1:nrow(MLD[[n]]))
    {    
      a <- MLD1[[n]][s,IndR]
      b <- MLD[[n]][s,IndR]
      Ind <- which(b>0)
      FigureOfMerit <- FigureOfMerit + mean(log10(a[Ind]/b[Ind])^2,trim=0.2)/nrow(MLD[[n]])
    }
  }
  return(FigureOfMerit/(length(MLD1)-1))
}

# function to plot results----
source("/home/projects/feariel/mishash/Development/ThreeFold/src/plotResults.R")
# plotAlpha()
# plotResults()
# fjft
# loop over different folds ----
print(paste0("FOUND ", ncores, " CORES!"))
species=c("Cronobacter","Raoultella","Proteus","Serratia","Citrobacter","Vibrio","Enterobacter","Salmonella","Klebsiella","Escherichia")

print(species)
nSp <- length(species)
nG <- 0
spList <- as.list(species)

dir <- "/home/projects/feariel/mishash/Development/ThreeFold/" # Afalina
print("Plotting MLDs")
system(paste0("mkdir ",dir,"data/processed/",dirname))

# get maps for indices ----
ijk_from_set <- list()
set_from_ijk <- list()
for (n in 1:(nSp-nG))
{
  ijk_from_set[[n]] <- combinations(n = nSp-nG, r = n, v = 1:(nSp-nG), repeats.allowed = FALSE,set=TRUE)
  set_from_ijk[[n]] <- 1:nrow(ijk_from_set[[n]])
  names(set_from_ijk[[n]]) <- sapply(1:nrow(ijk_from_set[[n]]),function(i){paste(ijk_from_set[[n]][i,],collapse='_')})
}

sList <- s_i_insteadof_jList <- iList <- jList <- list()
for (n in 1:(nSp-nG))
{
  sList[[n]] <- s_i_insteadof_jList[[n]] <- iList[[n]] <- jList[[n]] <- 1;
  for (s in 1:length(set_from_ijk[[n]]))
  {
    ijk <- ijk_from_set[[n]][s,]
    not_ijk <- setdiff(1:(nSp-nG),ijk)
    for (i in not_ijk)
    {
      for (j in ijk)
      {
        i_insteadof_j <- ijk; i_insteadof_j[which(i_insteadof_j==j)] <- i;  s_i_insteadof_j <- set_from_ijk[[n]][paste(sort(i_insteadof_j),collapse='_')]
        sList[[n]] <- c(sList[[n]],s)
        s_i_insteadof_jList[[n]] <- c(s_i_insteadof_jList[[n]],s_i_insteadof_j)
        iList[[n]] <- c(iList[[n]],i)
        jList[[n]] <- c(jList[[n]],j)
      }
    }
  }
  sList[[n]] <- sList[[n]][-1]; s_i_insteadof_jList[[n]] <- s_i_insteadof_jList[[n]][-1]; iList[[n]] <- iList[[n]][-1]; jList[[n]] <- jList[[n]][-1]
}



# get empirical MLD from full one  ####
MLD <- readRDS(file = paste0("/home/projects/feariel/mishash/Development/ThreeFold/data/processed/global/MLDplasmids.rds"))

# initial values ----
L10 <- 1e4
rhoM0 <- matrix(0,nSp,nSp)
colnames(rhoM0) <- rownames(rhoM0) <- species
for (s in 1:((nSp-1)*nSp/2))
{
  i <- ijk_from_set[[2]][s,1]
  j <- ijk_from_set[[2]][s,2]
  rhoM0[i,j] <- rhoM0[j,i] <- sum(MLD[[2]][s,])/sum(rV^(-3))/L10/2
}
rhoM <- rhoM0

gammaD0 <- read.table("/home/projects/feariel/mishash/Development/ThreeFold/data/processed/remove4backup/gammaD.tsv", header = TRUE); 
L20 <- median(gammaD0$L2)
gamma0 <- rep(10,nSp)
names(gamma0) <- species


# gamma opt----

LLgamma <- function(par,returnMLD=FALSE,sign=1)
{
  L2 <- L20*10^par[1]; par <- par[-1]
  
  gamma <- gamma0*10^par
  names(gamma) <- species
  MLD2 <- MLD
  for (n in 5:length(MLD))
  {
    for (s in 1:nrow(MLD[[n]]))
    {
      ijk <- ijk_from_set[[n]][s,]
      MLD2[[n]][s,] <- prod(gamma[ijk])*n*(n+1)/rV^(n+2)*L2/sum(gamma[ijk])
    }
  }

  
  foo <- distMLD(MLD2,MLD,IndR)
  return(sign*ifelse(!is.na(foo),foo,1e5))
}
dat <- data.frame()
par <- c(rep(0, 1+nSp))
par <- par + 1e-5*(runif(par)-0.5)
lower <- par-2
upper <- par+2

print("adding variables to cluster")
clusterExport(cl, list("LLgamma","set_from_ijk","ijk_from_set","rV","MLD","IndR","nSp","nG","rTh","rB","getMLD","dr","drm2","L20","gamma0","distMLD","species"));
print("starting GA optimizing")
library(GA)
# res <- ga(
#   type = "real-valued",
#   fitness = LLgamma,
#   sign = -1, # GA maximizes by default, so negate the function
#   lower = lower,
#   upper = upper,
#   popSize = ncores,
#   maxiter = 200,
#   run = 500,
#   parallel = cl,
#   monitor = TRUE,
#   optim=FALSE,
#   optimArgs = list(method = "L-BFGS-B", poptim=1/50, control = list(maxit = 20))
# )
# par <- res@solution[1,]
# res <- optimParallel(par=par, fn=LLgamma,method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent")[4],
#                      lower=lower, upper=upper,
#                      control = list(REPORT=1,maxit=300,trace=1,pgtol=1e-9,factr=1e-9,factr=1e-9,pgtol=1e-9,fnscale=1), hessian = FALSE,
#                      parallel = list(cl=cl),returnMLD=FALSE)
# par <- res$par
L2 <- L20*10^par[1]; par <- par[-1]
L20 <- L2
gamma <- gamma0*10^par

gamma <- readRDS(file = paste0(dir,"data/processed/",dirname,"/gamma.rds"))
L1 <- readRDS(file = paste0(dir,"data/processed/",dirname,"/L1.rds"))
L2 <- readRDS(file = paste0(dir,"data/processed/",dirname,"/L2.rds"))


gamma0 <- gamma
names(gamma) <- species
MLD2 <- MLDanalytic <- MLD
for (n in 2:length(MLD))
{
  for (s in 1:nrow(MLD[[n]]))
  {
    ijk <- ijk_from_set[[n]][s,]
    MLD2[[n]][s,] <- prod(gamma[ijk])*n*(n+1)/rV^(n+2)*L2/sum(gamma[ijk])
  }
}
L10 <- 1e4
rhoM0 <- matrix(0,nSp,nSp)
colnames(rhoM0) <- rownames(rhoM0) <- species
for (s in 1:((nSp-1)*nSp/2))
{
  i <- ijk_from_set[[2]][s,1]
  j <- ijk_from_set[[2]][s,2]
  Ind <- which(MLD[[2]][s,]>0)
  rhoM0[i,j] <- rhoM0[j,i] <- mean((MLD[[2]][s,Ind]-MLD2[[2]][s,Ind])*rV[Ind]^3,trim=0.2)/L10/2

  if (rhoM0[i,j]<0)
  {
    rhoM0[i,j] <- rhoM0[j,i] <- mean(MLD[[2]][s,]*rV^3,trim=0.2)/L10/2
  }
}
rhoM <- rhoM0

# full opt ----
clusterExport(cl, list("rhoM0","set_from_ijk","ijk_from_set","rV","MLD","IndR","nSp","nG","rhoM","rTh","rB","getMLD","dr","drm2","L10","L20","gamma0","distMLD","species","prefactor"));
LL <- function(par,returnMLD=FALSE,sign=1)
{
  rhoM[upper.tri(rhoM, diag=FALSE)] <- rhoM0[upper.tri(rhoM, diag=FALSE)]*10^par[1:(nSp*(nSp-1)/2)]; par <- par[-(1:(nSp*(nSp-1)/2))]
  rhoM[lower.tri(rhoM, diag=FALSE)] <- 0; rhoM <- rhoM + t(rhoM) #symmetric
  diag(rhoM) <- 0
  colnames(rhoM) <- rownames(rhoM) <- species
  
  L1 <- L10*10^par[1]; par <- par[-1]
  L2 <- L20*10^par[1]; par <- par[-1]
  
  gamma <- gamma0*10^(par)
  names(gamma) <- species
  MLD1 <- MLD
  if (1==1) {MLD1 <- getMLD(rhoM,N=5,IndR); for (n in 6:10) {MLD1[[n]] <- MLD[[n]]*0}}
  else
  {
    for (n in 2:length(MLD)){MLD1[[n]] <- MLD[[n]]*0}
    for (s in 1:nrow(MLD[[2]]))
    {
      ijk <- ijk_from_set[[2]][s,]
      MLD1[[2]][s,] <- prefactor(rhoM[ijk,ijk])/rV^3
    }
    for (s in 1:nrow(MLD[[3]]))
    {
      ijk <- ijk_from_set[[3]][s,]
      MLD1[[3]][s,] <- prefactor(rhoM[ijk,ijk])/rV^4
    }
    for (s in 1:nrow(MLD[[4]]))
    {
      ijk <- ijk_from_set[[4]][s,]
      MLD1[[4]][s,] <- prefactor(rhoM[ijk,ijk])/rV^5/2
    }
    for (s in 1:nrow(MLD[[5]]))
    {
      ijk <- ijk_from_set[[5]][s,]
      MLD1[[5]][s,] <- prefactor(rhoM[ijk,ijk])/rV^6/6
    }
    for (s in 1:nrow(MLD[[6]]))
    {
      ijk <- ijk_from_set[[6]][s,]
      MLD1[[6]][s,] <- prefactor(rhoM[ijk,ijk])/rV^7/24
    }
  }
  
  
  MLD2 <- MLDanalytic <- MLD
  
  for (n in 2:length(MLD))
  {
    MLD1[[n]] <- MLD1[[n]]*L1
    for (s in 1:nrow(MLD[[n]]))
    {
      ijk <- ijk_from_set[[n]][s,]
      MLD2[[n]][s,] <- prod(gamma[ijk])*n*(n+1)/rV^(n+2)*L2/sum(gamma[ijk])
    }
    MLDanalytic[[n]] <- MLD1[[n]] + MLD2[[n]]
  }
  
  foo <- distMLD(MLDanalytic,MLD,IndR)
  if (!returnMLD)
  {return(sign*ifelse(!is.na(foo),foo,1e5))}
  else
  {return(list(value=foo,MLDanalytic=MLDanalytic,MLD1=MLD1,MLD2=MLD2,rhoM,gamma,L1,L2))}
}
clusterExport(cl, list("LL"));
# optimization ----
par <- c(rep(0, (nSp-1)*nSp/2+2+nSp))
par <- par + 1e-15*(runif(par)-0.5)
lower <- par-1.5
upper <- par+1.5

print("starting full BFGS optimization")
res <- optimParallel(par=par, fn=LL,method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent")[4],
                     lower=lower, upper=upper,
                     control = list(REPORT=1,maxit=200,trace=1,pgtol=1e-6,factr=1e-6,factr=1e-6,pgtol=1e-6,fnscale=1), hessian = FALSE,
                     parallel = list(cl=cl),returnMLD=FALSE)
par <- res$par

listfoo <- LL(par,returnMLD=TRUE)
MLDanalytic <- listfoo[[2]]; MLD1 <- listfoo[[3]]; MLD2 <- listfoo[[4]]; rhoM <- listfoo[[5]]; gamma <- listfoo[[6]]; L1 <- listfoo[[7]]; L2 <- listfoo[[8]];
saveRDS(MLDanalytic, file = paste0(dir,"data/processed/",dirname,"/MLDanalytic.rds"))    
saveRDS(MLD1, file = paste0(dir,"data/processed/",dirname,"/MLDdirect.rds"))   
saveRDS(MLD2, file = paste0(dir,"data/processed/",dirname,"/MLDghost.rds"))
saveRDS(rhoM, file = paste0(dir,"data/processed/",dirname,"/rhoM.rds"))   
saveRDS(ijk_from_set, file = paste0(dir,"data/processed/",dirname,"/ijk_from_set.rds"))  
saveRDS(set_from_ijk, file = paste0(dir,"data/processed/",dirname,"/set_from_ijk.rds"))  
saveRDS(MLD, file = paste0(dir,"data/processed/",dirname,"/MLD.rds"))
saveRDS(gamma, file = paste0(dir,"data/processed/",dirname,"/gamma.rds"))
saveRDS(L1, file = paste0(dir,"data/processed/",dirname,"/L1.rds"))
saveRDS(L2, file = paste0(dir,"data/processed/",dirname,"/L2.rds"))
# plotResults(dirname=dirname,name=1,rhoM=rhoM,nSp=nSp,nG=nG,MLDsim=MLDanalytic,MLDanalytic=MLDanalytic,MLDdirect=MLD1,MLDghost=MLD2,ijk_from_set=ijk_from_set,set_from_ijk=set_from_ijk,MLD=MLD)
plotResults()

for (i in 1:(nrow(rhoM)-1))
{
  for (j in (i+1):ncol(rhoM))
  {
    if (nrow(rhoD)>0)
    {
      rhoD <- rbind(rhoD,data.frame(species1=species[i],species2=species[j],pair=paste0(sort(c(species[i],species[j]))[1],"_vs_",sort(c(species[i],species[j]))[2]),rate=rhoM[species[i],species[j]],L1=L1))
    }
    else
    {
      rhoD <- data.frame(species1=species[i],species2=species[j],pair=paste0(sort(c(species[i],species[j]))[1],"_vs_",sort(c(species[i],species[j]))[2]),rate=rhoM[species[i],species[j]],L1=L1)
    }
  }
}
write.table(rhoD,paste0(dir,"data/processed/",dirname,"/rhoD.tsv"),row.names=FALSE);

for (i in 1:nSp)
{
  if (nrow(gammaD)>0)
  {
    gammaD <- rbind(gammaD,data.frame(species1=species[i],rate=gamma[species[i]],L2=L2))
  }
  else
  {
    gammaD <- data.frame(species1=species[i],rate=gamma[species[i]],L2=L2)
  }
}
write.table(gammaD,paste0(dir,"data/processed/",dirname,"/gammaD.tsv"),row.names=FALSE);


pdf(paste0(dir,"data/processed/",dirname,"/rates.pdf"),width=15,height=10)
p <- ggplot(rhoD, aes(x = as.factor(pair), y = log10(rate),fill=as.factor(pair))) + 
  geom_violin(scale = "width") +
  geom_jitter(shape=16, position=position_jitter(0.2),cex=0.5) +
  # scale_y_continuous(breaks=seq(0,12,1), expand = c(0, 0),limits=c(2,12)) +
  scale_x_discrete() +
  theme(panel.grid.major = element_line(color = "grey",linewidth = 0.25,linetype = 2),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none")
print(p)
p <- ggplot(gammaD, aes(x = as.factor(species1), y = log10(rate),fill=as.factor(species1))) + 
  geom_violin(scale = "width") +
  geom_jitter(shape=16, position=position_jitter(0.2),cex=0.5) +
  # scale_y_continuous(breaks=seq(0,12,1), expand = c(0, 0),limits=c(2,12)) +
  scale_x_discrete() +
  theme(panel.grid.major = element_line(color = "grey",linewidth = 0.25,linetype = 2),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none")
print(p)
L1V <- unique(rhoD$L1)
L2V <- unique(gammaD$L2)
dat <- data.frame(value=c(L1V,L2V),L=c(rep("L1",length(L1V)),rep("L2",length(L2V))))
p <- ggplot(dat, aes(x = as.factor(L), y = log10(value),fill=as.factor(L))) + 
  geom_violin(scale = "width") +
  geom_jitter(shape=16, position=position_jitter(0.2),cex=0.5) +
  # scale_y_continuous(breaks=seq(0,12,1), expand = c(0, 0),limits=c(2,12)) +
  scale_x_discrete() +
  theme(panel.grid.major = element_line(color = "grey",linewidth = 0.25,linetype = 2),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none")
print(p)

dev.off()
gc()





