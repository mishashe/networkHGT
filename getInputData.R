#!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# dirname <- args[1]
gc()
dir <- "/home/projects/feariel/mishash/Development/ThreeFold/" # Afalina
dirname <- "PlasmidScope"
system(paste0("mkdir ",dir,"data/processed/",dirname))

.libPaths("/home/projects/feariel/mishash/R/x86_64-pc-linux-gnu-library/4.4")
library(plotrix)
library(gtools)
library(foreach)
library(parallel)
library(doParallel)
library(MASS) # to access Animals data sets
library(scales)
library(ggplot2)
library(gridExtra)
library(Biostrings)

ncores <- 1
print(paste0("Starting cluster with ", ncores, " cores"))
cl <- parallel::makeCluster(ncores, type = "PSOCK")
# cl <- makeMPIcluster(ncores, outfile='', type='MPI')
print("generated the cluster")
registerDoParallel(cl)
setDefaultCluster(cl=cl)
clusterEvalQ(cl, library(foreach));


rTh <- 100
rB <- 10^seq(log10(rTh-0.5),6,0.05) 
rV <- 0.5*(rB[-length(rB)]+rB[-1])
nG <- 0

species=c("Cronobacter","Raoultella","Proteus","Serratia","Citrobacter","Vibrio","Enterobacter","Salmonella","Klebsiella","Escherichia")

remove <- c(1,6); species=species[-remove]
nSp <- length(species)

# get number of genomes for each genus----
nList <- foreach (i = 1:nSp, .inorder=TRUE) %dopar%
{
  as.numeric(system(paste0("grep -c '>' /home/projects/feariel/mishash/Development/ThreeFold/data/external/PlasmidScope/",species[i], ".fa"),intern = TRUE))
}
saveRDS(nList, file = paste0(dir,"data/processed/",dirname,"/nList.rds"))
nList <- readRDS(file = paste0(dir,"data/processed/",dirname,"/nList.rds"))    
print(nList)

# get mappings from s to ijk and back----
ijk_from_set <- list()
for (n in 1:(nSp-nG))
{
  ijk_from_set[[n]] <- combinations(n = nSp-nG, r = n, v = 1:(nSp-nG), repeats.allowed = FALSE,set=TRUE)
}
ijk_from_set <- list()
set_from_ijk <- list()
for (n in 1:(nSp-nG))
{
  ijk_from_set[[n]] <- combinations(n = (nSp-nG), r = n, v = 1:(nSp-nG), repeats.allowed = FALSE,set=TRUE)
  set_from_ijk[[n]] <- 1:nrow(ijk_from_set[[n]])
  names(set_from_ijk[[n]]) <- sapply(1:nrow(ijk_from_set[[n]]),function(i){paste(ijk_from_set[[n]][i,],collapse='_')})
}



# get MLD ----
MLD <- MLDplasmids <- MLDchromosomes  <- list()
fold <- 2
MLD[[fold]] <- matrix(0,choose(nSp-nG,fold),length(rV));
MLDplasmids[[fold]] <- MLDchromosomes[[fold]] <- MLD[[fold]]
wT <- rT <- c()
for (i in 1:(nSp-1-nG))
{
  for (j in (i+1):(nSp-nG))
  {
    Comparison <- paste0(species[j],"_vs_",species[i])
    file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
    comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
    colnames(comp) <- c("r","n")
    r <- comp$r
    w <- comp$n
    w <- w[r>=rTh]
    r <- r[r>=rTh]
    p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
    MLD[[fold]][set_from_ijk[[fold]][paste(c(i,j),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]
    wT <- c(wT,w)
    rT <- c(rT,r)
    
    file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
    comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
    colnames(comp) <- c("r","n")
    r <- comp$r
    w <- comp$n
    w <- w[r>=rTh]
    r <- r[r>=rTh]
    p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
    MLDplasmids[[fold]][set_from_ijk[[fold]][paste(c(i,j),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]
    
    file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
    comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
    colnames(comp) <- c("r","n")
    r <- comp$r
    w <- comp$n
    w <- w[r>=rTh]
    r <- r[r>=rTh]
    p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
    MLDchromosomes[[fold]][set_from_ijk[[fold]][paste(c(i,j),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]
  }
}
print(1+sum(wT)/sum(wT*log(rT/(rTh+0.5))))
m <- colMeans(MLD[[fold]]); m2 <- m
mPlasmids <- colMeans(MLDplasmids[[fold]])
mChromosomes <- colMeans(MLDchromosomes[[fold]])

Ind <- 1:max(which(m>1e-1*m[1]/(rV/rTh)^(fold+1)));
A <- exp(mean(log(m[Ind]))+mean(log(rV[Ind]^(fold+1))))
mTh1 <- A/rV^(fold+1)
A <- exp(mean(log(m[Ind]))+mean(log(rV[Ind]^(fold+2))))
mTh2 <- A/rV^(fold+2)
Ind <- 1:max(which(m>1e-3*mTh1));
plot2 <- ggplot(data=data.frame(m=m,rV=rV)[Ind,], aes(x = rV, y = m)) + geom_point(color = "black") +
  scale_x_log10(breaks = 10^seq(3,6,1),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +   
  # geom_point(data=data.frame(m=mPlasmids,rV=rV,mTh=mTh1)[Ind,],aes(x = rV, y = m), color = "red",shape=1) +
  # geom_point(data=data.frame(m=m,rV=rV,mTh=mTh1)[Ind,],aes(x = rV, y = m), color = "blue", shape=3) +
  # geom_line(data=data.frame(m=m,rV=rV,mTh=mPlasmids[10]*(rV[10]/rV)^(fold+1))[Ind,],aes(x = rV, y = mTh), color = "grey") +
  geom_line(data=data.frame(m=m,rV=rV,mTh=mPlasmids[10]*(rV[10]/rV)^(fold+1))[Ind,],aes(x = rV, y = mTh), color = "blue") +
  # geom_line(data=data.frame(m=m,rV=rV,mTh=mNonPlasmids[10]*(rV[10]/rV)^(fold+3))[Ind,],aes(x = rV, y = mTh), color = "grey") +
  theme_bw()
# 3-fold
if (nSp-nG>=3)
{
  fold <- 3
  MLD[[fold]] <- matrix(0,choose(nSp-nG,fold),length(rV));
  MLDplasmids[[fold]] <- MLDchromosomes[[fold]] <- MLD[[fold]]
  for (i in 1:(nSp-nG-2))
  {
    for (j in (i+1):(nSp-nG-1))
    {
      for (k in (j+1):(nSp-nG))
      {
        
        Comparison <- paste0(species[k],"_vs_",species[j],"_vs_",species[i])
        file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
        comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
        colnames(comp) <- c("r","n")
        r <- comp$r
        w <- comp$n
        w <- w[r>=rTh]
        r <- r[r>=rTh]
        p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
        MLD[[fold]][set_from_ijk[[fold]][paste(c(i,j,k),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]
        
        
        file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
        comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
        colnames(comp) <- c("r","n")
        r <- comp$r
        w <- comp$n
        w <- w[r>=rTh]
        r <- r[r>=rTh]
        p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
        MLDplasmids[[fold]][set_from_ijk[[fold]][paste(c(i,j,k),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]
        
        file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
        if (file.exists(file) & file.size(file) != 0L)
        {
          comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
          colnames(comp) <- c("r","n")
          r <- comp$r
          w <- comp$n
          w <- w[r>=rTh]
          r <- r[r>=rTh]
          p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
          MLDchromosomes[[fold]][set_from_ijk[[fold]][paste(c(i,j,k),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]
        }
        else
        {
          MLDchromosomes[[fold]][set_from_ijk[[fold]][paste(c(i,j,k),collapse='_')],]  <- rV*0
        }
      }
    }
  }
}
m <- colMeans(MLD[[fold]]); m3 <- m
mPlasmids <- colMeans(MLDplasmids[[fold]])
mChromosomes <- colMeans(MLDchromosomes[[fold]])

Ind <- 1:max(which(m>1e-1*m[1]/(rV/rTh)^(fold+2)));
A <- exp(mean(log(m[Ind]))+mean(log(rV[Ind]^(fold+1))))
mTh1 <- A/rV^(fold+1)
A <- exp(mean(log(m[Ind]))+mean(log(rV[Ind]^(fold+2))))
mTh2 <- A/rV^(fold+2)
Ind <- 1:max(which(m>1e-3*mTh2));
plot3 <- ggplot(data=data.frame(m=m,rV=rV)[Ind,], aes(x = rV, y = m)) + geom_point(color = "black") +
  scale_x_log10(breaks = 10^seq(3,6,1),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +   
  # geom_point(data=data.frame(m=m,rV=rV,mTh=mTh1)[Ind,],aes(x = rV, y = m), color = "red",shape=1) +
  # geom_point(data=data.frame(m=m,rV=rV,mTh=mTh1)[Ind,],aes(x = rV, y = m), color = "blue", shape=3) +
  # geom_line(data=data.frame(m=m,rV=rV,mTh=mPlasmids[10]*(rV[10]/rV)^(fold+2))[Ind,],aes(x = rV, y = mTh), color = "grey") +
  geom_line(data=data.frame(m=m,rV=rV,mTh=m[5]*(rV[5]/rV)^(fold+1))[Ind,],aes(x = rV, y = mTh), color = "blue") +
  geom_line(data=data.frame(m=m,rV=rV,mTh=m[5]*(rV[5]/rV)^(fold+2))[Ind,],aes(x = rV, y = mTh), color = "red") +
  theme_bw()
# 4-fold
if ((nSp-nG)>=4)
{
  fold <- 4
  MLD[[fold]] <- matrix(0,choose(nSp-nG,fold),length(rV));
  MLDplasmids[[fold]] <- MLDchromosomes[[fold]] <- MLD[[fold]]
  for (i in 1:(nSp-nG-3))
  {
    for (j in (i+1):(nSp-nG-2))
    {
      for (k in (j+1):(nSp-nG-1))
      {
        for (l in (k+1):(nSp-nG))
        {
          Comparison <- paste0(species[l],"_vs_",species[k],"_vs_",species[j],"_vs_",species[i])
          file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
          comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
          colnames(comp) <- c("r","n")
          r <- comp$r
          w <- comp$n
          w <- w[r>=rTh]
          r <- r[r>=rTh]
          p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
          MLD[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]/nList[[l]]
          
          
          file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
          comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric"));
          colnames(comp) <- c("r","n")
          r <- comp$r
          w <- comp$n
          w <- w[r>=rTh]
          r <- r[r>=rTh]
          p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
          MLDplasmids[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]/nList[[l]]
          
          file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
          if (file.exists(file) & file.size(file) != 0L)
          {
            comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
            colnames(comp) <- c("r","n")
            r <- comp$r
            w <- comp$n
            w <- w[r>=rTh]
            r <- r[r>=rTh]
            p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
            MLDchromosomes[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]/nList[[l]]
          }
          else
          {
            MLDchromosomes[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l),collapse='_')],]  <- rV*0
          }
        }
      }
    }
  }
}
m <- colMeans(MLD[[fold]]); m4 <-m
mPlasmids <- colMeans(MLDplasmids[[fold]])
mChromosomes <- colMeans(MLDchromosomes[[fold]])
Ind <- 1:max(which(m>1e-1*m[1]/(rV/rTh)^(fold+2)));
A <- exp(mean(log(m[Ind]))+mean(log(rV[Ind]^(fold+1))))
mTh1 <- A/rV^(fold+1)
A <- exp(mean(log(m[Ind]))+mean(log(rV[Ind]^(fold+2))))
mTh2 <- A/rV^(fold+2)
Ind <- 1:max(which(m>1e-3*mTh2));
plot4 <- ggplot(data=data.frame(m=m,rV=rV)[Ind,], aes(x = rV, y = m)) + geom_point(color = "black") +
  scale_x_log10(breaks = 10^seq(3,6,1),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +   
  # geom_point(data=data.frame(m=mPlasmids,rV=rV,mTh=mTh1)[Ind,],aes(x = rV, y = m), color = "red",shape=1) +
  # geom_point(data=data.frame(m=mChromosomes,rV=rV,mTh=mTh1)[Ind,],aes(x = rV, y = m), color = "blue", shape=3) +
  # geom_line(data=data.frame(m=m,rV=rV,mTh=mPlasmids[10]*(rV[10]/rV)^(fold+2))[Ind,],aes(x = rV, y = mTh), color = "grey") +
  geom_line(data=data.frame(m=m,rV=rV,mTh=m[5]*(rV[5]/rV)^(fold+2))[Ind,],aes(x = rV, y = mTh), color = "red") +
  # geom_line(data=data.frame(m=m,rV=rV,mTh=mNonPlasmids[10]*(rV[10]/rV)^(fold+3))[Ind,],aes(x = rV, y = mTh), color = "grey") +
  theme_bw()
# 5-fold
if ((nSp-nG)>=5)
{
  fold <- 5 
  MLD[[fold]] <- matrix(0,choose(nSp-nG,fold),length(rV));
  MLDplasmids[[fold]] <- MLDchromosomes[[fold]] <- MLD[[fold]]
  for (i in 1:(nSp-nG-4))
  {
    for (j in (i+1):(nSp-nG-3))
    {
      for (k in (j+1):(nSp-nG-2))
      {
        for (l in (k+1):(nSp-nG-1))
        {
          for (mm in (l+1):(nSp-nG))
          {
            Comparison <- paste0(species[mm],"_vs_",species[l],"_vs_",species[k],"_vs_",species[j],"_vs_",species[i])
            file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
            comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
            colnames(comp) <- c("r","n")
            r <- comp$r
            w <- comp$n
            w <- w[r>=rTh]
            r <- r[r>=rTh]
            p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
            MLD[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]/nList[[l]]/nList[[mm]]
            
            file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
            comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric"));
            colnames(comp) <- c("r","n")
            r <- comp$r
            w <- comp$n
            w <- w[r>=rTh]
            r <- r[r>=rTh]
            p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
            MLDplasmids[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]/nList[[l]]/nList[[mm]]
            
            file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
            if (file.exists(file) & file.size(file) != 0L)
            {
              comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
              colnames(comp) <- c("r","n")
              r <- comp$r
              w <- comp$n
              w <- w[r>=rTh]
              r <- r[r>=rTh]
              p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
              MLDchromosomes[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]/nList[[l]]/nList[[mm]]
            }
            else
            {
              MLDchromosomes[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm),collapse='_')],]  <- rV*0
            }
          }
        }
      }
    }
  }
}
m <- colMeans(MLD[[fold]]); m5 <- m
mPlasmids <- colMeans(MLDplasmids[[fold]])
mChromosomes <- colMeans(MLDchromosomes[[fold]])
Ind <- 1:max(which(m>1e-1*m[1]/(rV/rTh)^(fold+2)));
A <- exp(mean(log(m[Ind]))+mean(log(rV[Ind]^(fold+1))))
mTh1 <- A/rV^(fold+1)
A <- exp(mean(log(m[Ind]))+mean(log(rV[Ind]^(fold+2))))
mTh2 <- A/rV^(fold+2)
Ind <- 1:max(which(m>1e-3*mTh2));
plot5 <- ggplot(data=data.frame(m=m,rV=rV)[Ind,], aes(x = rV, y = m)) + geom_point(color = "black") +
  scale_x_log10(breaks = 10^seq(3,6,1),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  # geom_point(data=data.frame(m=mPlasmids,rV=rV,mTh=mTh1)[Ind,],aes(x = rV, y = m), color = "red",shape=1) +
  # geom_point(data=data.frame(m=mChromosomes,rV=rV,mTh=mTh1)[Ind,],aes(x = rV, y = m), color = "blue", shape=3) +
  # geom_line(data=data.frame(m=m,rV=rV,mTh=mPlasmids[10]*(rV[10]/rV)^(fold+2))[Ind,],aes(x = rV, y = mTh), color = "grey") +
  geom_line(data=data.frame(m=m,rV=rV,mTh=m[5]*(rV[5]/rV)^(fold+2))[Ind,],aes(x = rV, y = mTh), color = "red") +
  # geom_line(data=data.frame(m=m,rV=rV,mTh=mNonPlasmids[10]*(rV[10]/rV)^(fold+3))[Ind,],aes(x = rV, y = mTh), color = "grey") +
  theme_bw()
# # 6-fold
if ((nSp-nG)>=6)
{
  fold <- 6
  MLD[[fold]] <- matrix(0,choose(nSp-nG,fold),length(rV));
  MLDplasmids[[fold]] <- MLDchromosomes[[fold]] <- MLD[[fold]]
  for (i in 1:(nSp-nG-5))
  {
    for (j in (i+1):(nSp-nG-4))
    {
      for (k in (j+1):(nSp-nG-3))
      {
        for (l in (k+1):(nSp-nG-2))
        {
          for (mm in (l+1):(nSp-nG-1))
          {
            for (pp in (mm+1):(nSp-nG))
            {
              Comparison <- paste0(species[pp],"_vs_",species[mm],"_vs_",species[l],"_vs_",species[k],"_vs_",species[j],"_vs_",species[i])
              file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
              comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
              colnames(comp) <- c("r","n")
              r <- comp$r
              w <- comp$n
              w <- w[r>=rTh]
              r <- r[r>=rTh]
              p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
              MLD[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm,pp),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]/nList[[l]]/nList[[mm]]/nList[[pp]]
              
              file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
              comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
              colnames(comp) <- c("r","n")
              r <- comp$r
              w <- comp$n
              w <- w[r>=rTh]
              r <- r[r>=rTh]
              p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
              MLDplasmids[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm,pp),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]/nList[[l]]/nList[[mm]]/nList[[pp]]
            
              file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
              if (file.exists(file) & file.size(file) != 0L)
              {
                comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
                colnames(comp) <- c("r","n")
                r <- comp$r
                w <- comp$n
                w <- w[r>=rTh]
                r <- r[r>=rTh]
                p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
                MLDchromosomes[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm,pp),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]/nList[[l]]/nList[[mm]]/nList[[pp]]
              }
              else
              {
                MLDchromosomes[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm,pp),collapse='_')],]  <- rV*0
              }
              
            }
          }
        }
      }
    }
  }
}
m <- colMeans(MLD[[fold]]); m6 <- m
mPlasmids <- colMeans(MLDplasmids[[fold]])
mChromosomes <- colMeans(MLDchromosomes[[fold]])
Ind <- 1:max(which(m>1e-1*m[1]/(rV/rTh)^(fold+2)));
A <- exp(mean(log(m[Ind]))+mean(log(rV[Ind]^(fold+1))))
mTh1 <- A/rV^(fold+1)
A <- exp(mean(log(m[Ind]))+mean(log(rV[Ind]^(fold+2))))
mTh2 <- A/rV^(fold+2)
Ind <- 1:max(which(m>1e-3*mTh2));
plot6 <- ggplot(data=data.frame(m=m,rV=rV)[Ind,], aes(x = rV, y = m)) + geom_point(color = "black") +
  scale_x_log10(breaks = 10^seq(3,6,1),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +   
  # geom_point(data=data.frame(m=mPlasmids,rV=rV,mTh=mTh1)[Ind,],aes(x = rV, y = m), color = "red",shape=1) +
  # geom_point(data=data.frame(m=mChromosomes,rV=rV,mTh=mTh1)[Ind,],aes(x = rV, y = m), color = "blue", shape=3) +
  # geom_line(data=data.frame(m=m,rV=rV,mTh=mPlasmids[10]*(rV[10]/rV)^(fold+2))[Ind,],aes(x = rV, y = mTh), color = "grey") +
  geom_line(data=data.frame(m=m,rV=rV,mTh=m[5]*(rV[5]/rV)^(fold+2))[Ind,],aes(x = rV, y = mTh), color = "red") +
  # geom_line(data=data.frame(m=m,rV=rV,mTh=mNonPlasmids[10]*(rV[10]/rV)^(fold+3))[Ind,],aes(x = rV, y = mTh), color = "grey") +
  theme_bw()
# 7-fold
if ((nSp-nG)>=7)
{
  fold <- 7
  MLD[[fold]] <- matrix(0,choose(nSp-nG,fold),length(rV));
  MLDplasmids[[fold]] <- MLDchromosomes[[fold]] <- MLD[[fold]]
  wT <- rT <- c()
  for (i in 1:(nSp-nG-6))
  {
    for (j in (i+1):(nSp-nG-5))
    {
      for (k in (j+1):(nSp-nG-4))
      {
        for (l in (k+1):(nSp-nG-3))
        {
          for (mm in (l+1):(nSp-nG-2))
          {
            for (pp in (mm+1):(nSp-nG-1))
            {
              for (q in (pp+1):(nSp-nG))
              {
                Comparison <- paste0(species[q],"_vs_",species[pp],"_vs_",species[mm],"_vs_",species[l],"_vs_",species[k],"_vs_",species[j],"_vs_",species[i])
                file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
                comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
                colnames(comp) <- c("r","n")
                r <- comp$r
                w <- comp$n
                w <- w[r>=rTh]
                r <- r[r>=rTh]
                p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
                MLD[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm,pp,q),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]/nList[[l]]/nList[[mm]]/nList[[pp]]/nList[[q]]
                wT <- c(wT,w)
                rT <- c(rT,r)
                
                file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
                comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
                colnames(comp) <- c("r","n")
                r <- comp$r
                w <- comp$n
                w <- w[r>=rTh]
                r <- r[r>=rTh]
                p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
                MLDplasmids[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm,pp,q),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]/nList[[l]]/nList[[mm]]/nList[[pp]]/nList[[q]]
                
                file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
                if (file.exists(file) & file.size(file) != 0L)
                {
                  comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
                  colnames(comp) <- c("r","n")
                  r <- comp$r
                  w <- comp$n
                  w <- w[r>=rTh]
                  r <- r[r>=rTh]
                  p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
                  MLDchromosomes[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm,pp,q),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]/nList[[l]]/nList[[mm]]/nList[[pp]]/nList[[q]]
                }
                else
                {
                  MLDchromosomes[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm,pp,q),collapse='_')],]  <- rV*0
                }
                
              }
            }
          }
        }
      }
    }
  }
}
print(1+sum(wT)/sum(wT*log(rT/(rTh+0.5))))
m <- colMeans(MLD[[fold]]); m7 <- m
mPlasmids <- colMeans(MLDplasmids[[fold]])
mChromosomes <- colMeans(MLDchromosomes[[fold]])
Ind <- 1:max(which(m>1e-1*m[1]/(rV/rTh)^(fold+2)));
A <- exp(mean(log(m[Ind]))+mean(log(rV[Ind]^(fold+1))))
mTh1 <- A/rV^(fold+1)
A <- exp(mean(log(m[Ind]))+mean(log(rV[Ind]^(fold+2))))
mTh2 <- A/rV^(fold+2)
Ind <- 1:max(which(m>1e-3*mTh2));
plot7 <- ggplot(data=data.frame(m=m,rV=rV)[Ind,], aes(x = rV, y = m)) + geom_point(color = "black") +
  scale_x_log10(breaks = 10^seq(3,6,1),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +   
  # geom_point(data=data.frame(m=mPlasmids,rV=rV,mTh=mTh1)[Ind,],aes(x = rV, y = m), color = "red",shape=1) +
  # geom_point(data=data.frame(m=mChromosomes,rV=rV,mTh=mTh1)[Ind,],aes(x = rV, y = m), color = "blue", shape=3) +
  # geom_line(data=data.frame(m=m,rV=rV,mTh=mPlasmids[10]*(rV[10]/rV)^(fold+2))[Ind,],aes(x = rV, y = mTh), color = "grey") +
  geom_line(data=data.frame(m=m,rV=rV,mTh=m[5]*(rV[5]/rV)^(fold+2))[Ind,],aes(x = rV, y = mTh), color = "red") +
  # geom_line(data=data.frame(m=m,rV=rV,mTh=mNonPlasmids[10]*(rV[10]/rV)^(fold+3))[Ind,],aes(x = rV, y = mTh), color = "grey") +
  theme_bw()
# 8-fold
if ((nSp-nG)>=8)
{
  fold <- 8
  MLD[[fold]] <- matrix(0,choose(nSp-nG,fold),length(rV));
  MLDplasmids[[fold]] <- MLDchromosomes[[fold]] <- MLD[[fold]]
  wT <- rT <- c()
  for (i in 1:(nSp-nG-7))
  {
    for (j in (i+1):(nSp-nG-6))
    {
      for (k in (j+1):(nSp-nG-5))
      {
        for (l in (k+1):(nSp-nG-4))
        {
          for (mm in (l+1):(nSp-nG-3))
          {
            for (pp in (mm+1):(nSp-nG-2))
            {
              for (q in (pp+1):(nSp-nG-1))
              {
                for (rr in (q+1):(nSp-nG-0))
                {
                  Comparison <- paste0(species[rr],"_vs_",species[q],"_vs_",species[pp],"_vs_",species[mm],"_vs_",species[l],"_vs_",species[k],"_vs_",species[j],"_vs_",species[i])
                  file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
                  comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
                  colnames(comp) <- c("r","n")
                  r <- comp$r
                  w <- comp$n
                  w <- w[r>=rTh]
                  r <- r[r>=rTh]
                  p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
                  MLD[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm,pp,q,rr),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]/nList[[l]]/nList[[mm]]/nList[[pp]]/nList[[q]]/nList[[rr]]
                  wT <- c(wT,w)
                  rT <- c(rT,r)
                  
                  file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
                  comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
                  colnames(comp) <- c("r","n")
                  r <- comp$r
                  w <- comp$n
                  w <- w[r>=rTh]
                  r <- r[r>=rTh]
                  p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
                  MLDplasmids[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm,pp,q,rr),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]/nList[[l]]/nList[[mm]]/nList[[pp]]/nList[[q]]/nList[[rr]]
                  
                  file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
                  if (file.exists(file) & file.size(file) != 0L)
                  {
                    comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
                    colnames(comp) <- c("r","n")
                    r <- comp$r
                    w <- comp$n
                    w <- w[r>=rTh]
                    r <- r[r>=rTh]
                    p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
                    MLDchromosomes[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm,pp,q,rr),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]/nList[[l]]/nList[[mm]]/nList[[pp]]/nList[[q]]/nList[[rr]]
                  }
                  else
                  {
                    MLDchromosomes[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm,pp,q,rr),collapse='_')],]  <- rV*0
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
print(1+sum(wT)/sum(wT*log(rT/(rTh+0.5))))
m <- colMeans(MLD[[fold]]); m8 <- m
mPlasmids <- colMeans(MLDplasmids[[fold]])
mChromosomes <- colMeans(MLDchromosomes[[fold]])
Ind <- 1:max(which(m>1e-1*m[1]/(rV/rTh)^(fold+2)));
A <- exp(mean(log(m[Ind]))+mean(log(rV[Ind]^(fold+1))))
mTh1 <- A/rV^(fold+1)
A <- exp(mean(log(m[Ind]))+mean(log(rV[Ind]^(fold+2))))
mTh2 <- A/rV^(fold+2)
Ind <- 1:max(which(m>1e-3*mTh2));
plot8 <- ggplot(data=data.frame(m=m,rV=rV)[Ind,], aes(x = rV, y = m)) + geom_point(color = "black") +
  scale_x_log10(breaks = 10^seq(3,6,1),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +   
  # geom_point(data=data.frame(m=mPlasmids,rV=rV,mTh=mTh1)[Ind,],aes(x = rV, y = m), color = "red",shape=1) +
  # geom_point(data=data.frame(m=mChromosomes,rV=rV,mTh=mTh1)[Ind,],aes(x = rV, y = m), color = "blue", shape=3) +
  # geom_line(data=data.frame(m=m,rV=rV,mTh=mPlasmids[10]*(rV[10]/rV)^(fold+2))[Ind,],aes(x = rV, y = mTh), color = "grey") +
  geom_line(data=data.frame(m=m,rV=rV,mTh=m[5]*(rV[5]/rV)^(fold+2))[Ind,],aes(x = rV, y = mTh), color = "red") +
  # geom_line(data=data.frame(m=m,rV=rV,mTh=mNonPlasmids[10]*(rV[10]/rV)^(fold+3))[Ind,],aes(x = rV, y = mTh), color = "grey") +
  theme_bw()

# # 9-fold
if ((nSp-nG)>=9)
{
  fold <- 9
  MLD[[fold]] <- matrix(0,choose(nSp-nG,fold),length(rV));
  MLDplasmids[[fold]] <- MLDchromosomes[[fold]] <- MLD[[fold]]
  wT <- rT <- c()
  for (i in 1:(nSp-nG-8))
  {
    for (j in (i+1):(nSp-nG-7))
    {
      for (k in (j+1):(nSp-nG-6))
      {
        for (l in (k+1):(nSp-nG-5))
        {
          for (mm in (l+1):(nSp-nG-4))
          {
            for (pp in (mm+1):(nSp-nG-3))
            {
              for (q in (pp+1):(nSp-nG-2))
              {
                for (rr in (q+1):(nSp-nG-1))
                {
                  for (ss in (rr+1):(nSp-nG-0))
                  {
                    Comparison <- paste0(species[ss],"_vs_",species[rr],"_vs_",species[q],"_vs_",species[pp],"_vs_",species[mm],"_vs_",species[l],"_vs_",species[k],"_vs_",species[j],"_vs_",species[i])
                    file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
                    comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
                    colnames(comp) <- c("r","n")
                    r <- comp$r
                    w <- comp$n
                    w <- w[r>=rTh]
                    r <- r[r>=rTh]
                    p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
                    MLD[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm,pp,q,rr,ss),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]/nList[[l]]/nList[[mm]]/nList[[pp]]/nList[[q]]/nList[[rr]]/nList[[ss]]
                    wT <- c(wT,w)
                    rT <- c(rT,r)
                    
                    file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
                    comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
                    colnames(comp) <- c("r","n")
                    r <- comp$r
                    w <- comp$n
                    w <- w[r>=rTh]
                    r <- r[r>=rTh]
                    p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
                    MLDplasmids[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm,pp,q,rr,ss),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]/nList[[l]]/nList[[mm]]/nList[[pp]]/nList[[q]]/nList[[rr]]/nList[[ss]]
                  
                    file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
                    if (file.exists(file) & file.size(file) != 0L)
                    {
                      comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
                      colnames(comp) <- c("r","n")
                      r <- comp$r
                      w <- comp$n
                      w <- w[r>=rTh]
                      r <- r[r>=rTh]
                      p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
                      MLDchromosomes[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm,pp,q,rr,ss),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]/nList[[l]]/nList[[mm]]/nList[[pp]]/nList[[q]]/nList[[rr]]/nList[[ss]]
                    }
                    else
                    {
                      MLDchromosomes[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm,pp,q,rr,ss),collapse='_')],]  <- rV*0
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
print(1+sum(wT)/sum(wT*log(rT/(rTh+0.5))))
m <- colMeans(MLD[[fold]])
mPlasmids <- colMeans(MLDplasmids[[fold]])
mChromosomes <- colMeans(MLDchromosomes[[fold]])
Ind <- 1:max(which(m>1e-1*m[1]/(rV/rTh)^(fold+2)));
A <- exp(mean(log(m[Ind]))+mean(log(rV[Ind]^(fold+1))))
mTh1 <- A/rV^(fold+1)
A <- exp(mean(log(m[Ind]))+mean(log(rV[Ind]^(fold+2))))
mTh2 <- A/rV^(fold+2)
Ind <- 1:max(which(m>1e-3*mTh2));
LL <- function(par)
{
  return(mean((log10(m[Ind])-log10(10^par[1]*exp(-10^par[2]*rV[Ind])+10^par[3]*exp(-10^par[4]*rV[Ind])+10^par[5]*exp(-10^par[6]*rV[Ind])))^2))
}
par <- c(-0.596515,  -1.910987,  -8.618588,  -2.499669, -22.012788,  -3.669951)
print(LL(par))
# res <- optim(par=par, fn=LL,method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent")[1],
#              # lower=par-2, upper=par+2,
#              control = list(REPORT=1,maxit=2000,trace=1,pgtol=1e-8,factr=1e-86,factr=1e-8,pgtol=1e-8,fnscale=1), hessian = FALSE)
# par <- res$par
print(par)
mExp <- 10^par[1]*exp(-10^par[2]*rV)+10^par[3]*exp(-10^par[4]*rV)+10^par[5]*exp(-10^par[6]*rV)

plot9 <- ggplot(data=data.frame(m=m,rV=rV)[Ind,], aes(x = rV, y = m)) + geom_point(color = "black") +
  scale_x_log10(breaks = 10^seq(3,6,1),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +   
  geom_point(data=data.frame(m=mPlasmids,rV=rV,mTh=mTh1)[Ind,],aes(x = rV, y = m), color = "red",shape=1) +
  geom_point(data=data.frame(m=mChromosomes,rV=rV,mTh=mTh1)[Ind,],aes(x = rV, y = m), color = "blue", shape=3) +
  # geom_line(data=data.frame(m=m,rV=rV,mTh=mPlasmids[10]*(rV[10]/rV)^(fold+2))[Ind,],aes(x = rV, y = mTh), color = "grey") +
  geom_line(data=data.frame(m=m,rV=rV,mTh=m[10]*(rV[10]/rV)^(fold+2))[Ind,],aes(x = rV, y = mTh), color = "grey") +
  # geom_line(data=data.frame(m=m,rV=rV,mTh=mNonPlasmids[10]*(rV[10]/rV)^(fold+3))[Ind,],aes(x = rV, y = mTh), color = "grey") +
  theme_bw()
# # 10-fold
if ((nSp-nG)>=10)
{
  fold <- 10
  MLD[[fold]] <- matrix(0,choose(nSp-nG,fold),length(rV));
  MLDplasmids[[fold]] <- MLDchromosomes[[fold]] <- MLD[[fold]]
  wT <- rT <- c()
  for (i in 1:(nSp-nG-9))
  {
    for (j in (i+1):(nSp-nG-8))
    {
      for (k in (j+1):(nSp-nG-7))
      {
        for (l in (k+1):(nSp-nG-6))
        {
          for (mm in (l+1):(nSp-nG-5))
          {
            for (pp in (mm+1):(nSp-nG-4))
            {
              for (q in (pp+1):(nSp-nG-3))
              {
                for (rr in (q+1):(nSp-nG-2))
                {
                  for (ss in (rr+1):(nSp-nG-1))
                  {
                    for (t in (ss+1):(nSp-nG-0))
                    {
                      Comparison <- paste0(species[t],"_vs_",species[ss],"_vs_",species[rr],"_vs_",species[q],"_vs_",species[pp],"_vs_",species[mm],"_vs_",species[l],"_vs_",species[k],"_vs_",species[j],"_vs_",species[i])
                      file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
                      comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
                      colnames(comp) <- c("r","n")
                      r <- comp$r
                      w <- comp$n
                      w <- w[r>=rTh]
                      r <- r[r>=rTh]
                      p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
                      MLD[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm,pp,q,rr,ss,t),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]/nList[[l]]/nList[[mm]]/nList[[pp]]/nList[[q]]/nList[[rr]]/nList[[ss]]/nList[[t]]
                      wT <- c(wT,w)
                      rT <- c(rT,r)
                      
                      file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
                      comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
                      colnames(comp) <- c("r","n")
                      r <- comp$r
                      w <- comp$n
                      w <- w[r>=rTh]
                      r <- r[r>=rTh]
                      p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
                      MLDplasmids[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm,pp,q,rr,ss,t),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]/nList[[l]]/nList[[mm]]/nList[[pp]]/nList[[q]]/nList[[rr]]/nList[[ss]]/nList[[t]]
                    
                      file <- paste0(dir,"data/processed/Comparisons/PlasmidScope/",Comparison,"/",Comparison,".tsv")
                      if (file.exists(file) & file.size(file) != 0L)
                      {
                        comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
                        colnames(comp) <- c("r","n")
                        r <- comp$r
                        w <- comp$n
                        w <- w[r>=rTh]
                        r <- r[r>=rTh]
                        p <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
                        MLDchromosomes[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm,pp,q,rr,ss,t),collapse='_')],]  <- p$counts/diff(rB)/nList[[i]]/nList[[j]]/nList[[k]]/nList[[l]]/nList[[mm]]/nList[[pp]]/nList[[q]]/nList[[rr]]/nList[[ss]]/nList[[t]]
                      }
                      else
                      {
                        MLDchromosomes[[fold]][set_from_ijk[[fold]][paste(c(i,j,k,l,mm,pp,q,rr,ss,t),collapse='_')],]  <- rV*0
                      }
                      
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
print(1+sum(wT)/sum(wT*log(rT/(rTh+0.5))))
m <- colMeans(MLD[[fold]])
mPlasmids <- colMeans(MLDplasmids[[fold]])
mChromosomes <- colMeans(MLDchromosomes[[fold]])
Ind <- 1:max(which(m>1e-1*m[1]/(rV/rTh)^(fold+2)));
A <- exp(mean(log(m[Ind]))+mean(log(rV[Ind]^(fold+1))))
mTh1 <- A/rV^(fold+1)
A <- exp(mean(log(m[Ind]))+mean(log(rV[Ind]^(fold+2))))
mTh2 <- A/rV^(fold+2)
Ind <- 1:max(which(m>1e-4*mTh2));
LL <- function(par)
{
  return(mean((log10(m[Ind])-log10(10^par[1]*exp(-10^par[2]*rV[Ind])+10^par[3]*exp(-10^par[4]*rV[Ind])))^2))
}
par <- c(-6,-3,-13,-4)
# print(LL(par))
# res <- optim(par=par, fn=LL,method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent")[1],
#                      # lower=par-2, upper=par+2,
#                      control = list(REPORT=1,maxit=2000,trace=1,pgtol=1e-8,factr=1e-8,factr=1e-8,pgtol=1e-8,fnscale=1), hessian = FALSE)
# par <- res$par
mExp <- 10^par[1]*exp(-10^par[2]*rV)+10^par[3]*exp(-10^par[4]*rV)
print(par)


plot10 <- ggplot(data=data.frame(m=m,rV=rV)[Ind,], aes(x = rV, y = m)) + geom_point(color = "black") +
  scale_x_log10(breaks = 10^seq(3,6,1),limits=c(1e3,1e4),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +   
  geom_point(data=data.frame(m=mPlasmids,rV=rV,mTh=mTh1)[Ind,],aes(x = rV, y = m), color = "red",shape=1) +
  geom_point(data=data.frame(m=mChromosomes,rV=rV,mTh=mTh1)[Ind,],aes(x = rV, y = m), color = "blue", shape=3) +
  # geom_line(data=data.frame(m=m,rV=rV,mTh=mPlasmids[10]*(rV[10]/rV)^(fold+2))[Ind,],aes(x = rV, y = mTh), color = "grey") +
  geom_line(data=data.frame(m=m,rV=rV,mTh=m[5]*(rV[5]/rV)^(fold+2))[Ind,],aes(x = rV, y = mTh), color = "grey") +
  # geom_line(data=data.frame(m=m,rV=rV,mTh=mNonPlasmids[10]*(rV[10]/rV)^(fold+3))[Ind,],aes(x = rV, y = mTh), color = "grey") +
  # geom_line(data=data.frame(m=m,rV=rV,mTh=mExp)[Ind,],aes(x = rV, y = mTh), color = "green") +
  theme_bw()



pdf(paste0(dir,"data/processed/",dirname,"/MLDmean.pdf"),width=10,height=5)
grid.arrange(plot2, plot3, plot4, plot5, plot6, plot7, plot8, nrow = 2, ncol = 4)
dev.off()

plotall <- ggplot(data=data.frame(m=m2,rV=rV), aes(x = rV, y = m)) + geom_point(color = "black") +
  scale_x_log10(breaks = 10^seq(3,6,1),limits=c(1e3,8e4),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),limits=c(1e-21,1e-1),
                labels = trans_format("log10", math_format(10^.x))) + 
    geom_line(data=data.frame(m=m,rV=rV,mTh=m2[5]*(rV[5]/rV)^(2+1)),aes(x = rV, y = mTh), color = "black",linetype = "dashed") +
  geom_line(data=data.frame(m=m,rV=rV,mTh=m3[5]*(rV[5]/rV)^(3+1)),aes(x = rV, y = mTh), color = "red",linetype = "dashed") +
  geom_line(data=data.frame(m=m,rV=rV,mTh=m3[5]*(rV[5]/rV)^(3+2)),aes(x = rV, y = mTh), color = "red") +
  geom_line(data=data.frame(m=m,rV=rV,mTh=m4[5]*(rV[5]/rV)^(4+2)),aes(x = rV, y = mTh), color = "blue") +
  geom_line(data=data.frame(m=m,rV=rV,mTh=m5[5]*(rV[5]/rV)^(5+2)),aes(x = rV, y = mTh), color = "green") +
  geom_line(data=data.frame(m=m,rV=rV,mTh=m6[5]*(rV[5]/rV)^(6+2)),aes(x = rV, y = mTh), color = "purple") +
  geom_line(data=data.frame(m=m,rV=rV,mTh=m7[5]*(rV[5]/rV)^(7+2)),aes(x = rV, y = mTh), color = "orange") +
  geom_line(data=data.frame(m=m,rV=rV,mTh=m8[5]*(rV[5]/rV)^(8+2)),aes(x = rV, y = mTh), color = "brown") +
  geom_point(data=data.frame(m=m2,rV=rV,mTh=mTh1),aes(x = rV, y = m), color = "black",shape=16) +
  geom_point(data=data.frame(m=m3,rV=rV,mTh=mTh1),aes(x = rV, y = m), color = "red",shape=1) +
  geom_point(data=data.frame(m=m4,rV=rV,mTh=mTh1),aes(x = rV, y = m), color = "blue",shape=17) +
  geom_point(data=data.frame(m=m5,rV=rV,mTh=mTh1),aes(x = rV, y = m), color = "green",shape=2) +
  geom_point(data=data.frame(m=m6,rV=rV,mTh=mTh1),aes(x = rV, y = m), color = "purple",shape=15) +
  geom_point(data=data.frame(m=m7,rV=rV,mTh=mTh1),aes(x = rV, y = m), color = "orange",shape=0) +
  geom_point(data=data.frame(m=m8,rV=rV,mTh=mTh1),aes(x = rV, y = m), color = "brown",shape=3) +
  theme_bw() +  
  theme(
    axis.text.x = element_text(size = 12),  
    axis.text.y = element_text(size = 12) 
  )

pdf(paste0(dir,"data/processed/",dirname,"/MLDmean_all.pdf"),width=5,height=5)
grid.arrange(plotall, nrow = 1, ncol = 1)
dev.off()

saveRDS(MLD, file = paste0(dir,"data/processed/",dirname,"/MLD.rds"))    
# saveRDS(MLDplasmids, file = paste0(dir,"data/processed/",dirname,"/MLDplasmids.rds"))    
# saveRDS(MLDchromosomes, file = paste0(dir,"data/processed/",dirname,"/MLDchromosomes.rds"))    


