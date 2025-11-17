#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
sp1 <- args[1]
sp2 <- args[2]
sp3 <- args[3]

# sp1="Escherichia_coli"
# sp2="Klebsiella_pneumoniae"
# sp3="Salmonella_enterica"

library(stringr)
library(foreach)
library(iterators)
library(doParallel)
require(intervals)
library(sets)
require(data.table)
library(intervalaverage)
library(dplyr)
library(plotrix)
library(iterators)
library(rhdf5)
library(seqinr);
library(questionr)
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp(paste0("/home/projects/feariel/mishash/Development/ThreeFold/src/intersect_weighted.cpp"))

#Simulation <- "PlasmidScope/" # for simuation
Simulation <- "" # for empirical


registerDoParallel(5)
combine_tables <- function(table1,table2)
{
  names12 <- unique(c(names(table1),names(table2)))
  table1 <- table1[names12]
  setdiff21 <- setdiff(names(table2),names(table1))
  if (length(setdiff21)>0)
  {
    names(table1)[is.na(names(table1))] <- setdiff21
    table1[setdiff21] <- 0
  }
  table2 <- table2[names12]
  setdiff12 <- setdiff(names(table1),names(table2))
  if (length(setdiff12)>0)
  {
    names(table2)[is.na(names(table2))] <- setdiff12
    table2[setdiff12] <- 0
  }
  return(table1[names12] + table2[names12])
}





# dir <- "/home/misha/Documents/Development/ThreeFold/" # Bechet
dir <- "/home/projects/feariel/mishash/Development/ThreeFold/"

system(paste0("mkdir -p ",dir,"plots/"))

# sp1 <- "Escherichia_coli"
# sp2 <- "Klebsiella_pneumoniae"
# sp3 <- "Serratia"

system(paste0("mkdir -p ", dir,"data/processed/Comparisons/",Simulation,sp3,"_vs_",sp2,"_vs_",sp1,"/"))
group1 <- h5ls(paste0(dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".h5"))$group[2]
group2 <- h5ls(paste0(dir,"data/processed/Comparisons/",Simulation,sp3,"_vs_",sp1,"/",sp3,"_vs_",sp1,".h5"))$group[2]

headers1 <- h5ls(paste0(dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".h5"))$name; headers <- unique(headers1)
headers2 <- h5ls(paste0(dir,"data/processed/Comparisons/",Simulation,sp3,"_vs_",sp1,"/",sp3,"_vs_",sp1,".h5"))$name; headers <- unique(intersect(headers,headers2))
headers <- headers[!headers %in% c("comp1","comp2","comp")]
rTable <- foreach(header=headers, .inorder=FALSE, .combine=combine_tables) %dopar%
{
  intervals1 <- h5read(paste0(dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".h5"), paste0(group1,"/",header))
  intervals2 <- h5read(paste0(dir,"data/processed/Comparisons/",Simulation,sp3,"_vs_",sp1,"/",sp3,"_vs_",sp1,".h5"), paste0(group2,"/",header))
  intersects <- as.data.table(intersect_weighted(as.matrix(intervals1), as.matrix(intervals2)))
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  wtd.table(as.integer(intersects$V2-intersects$V1+1),w=as.numeric(intersects$w))
}
write.table(data.frame(r=as.integer(names(rTable)),m=as.numeric(rTable)),paste0(dir,"data/processed/Comparisons/",Simulation,sp3,"_vs_",sp2,"_vs_",sp1,"/",sp3,"_vs_",sp2,"_vs_",sp1,".tsv"),append=FALSE,row.names=FALSE,col.names=FALSE,sep = "\t",quote=FALSE)




# 
# 
# library(stringr)
# library(foreach)
# library(iterators)
# library(doParallel)
# require(intervals)
# library(sets)
# require(data.table)
# library(intervalaverage)
# library(dplyr)
# library(plotrix)
# library(iterators)
# library(rhdf5)
# library(seqinr);
# library(questionr)
# 
# dir <- "/home/msheinman/Development/ThreeFold/" # Afalina
# 
#   # sp1 <- "Escherichia_coli"
#   # sp2 <- "Escherichia_albertii"
#   # sp3 <- "Escherichia_fergusonii"
# n1 <- length(Sys.glob(file.path(paste0(dir,"data/external/",sp1), "*.fna")))
# n2 <- length(Sys.glob(file.path(paste0(dir,"data/external/",sp2), "*.fna")))
# n3 <- length(Sys.glob(file.path(paste0(dir,"data/external/",sp3), "*.fna")))
# 
# 
# print("Plotting MLDs")
# rB <- 10^seq(log10(10)-0.1,6,0.05) 
# rV <- 0.5*(rB[-length(rB)]+rB[-1])
# 
# 
# file12 <- paste0(dir,"data/processed/Comparisons/",sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,"_table.tsv")
# file13 <- paste0(dir,"data/processed/Comparisons/",sp3,"_vs_",sp1,"/",sp3,"_vs_",sp1,"_table.tsv")
# file23 <- paste0(dir,"data/processed/Comparisons/",sp3,"_vs_",sp2,"/",sp3,"_vs_",sp2,"_table.tsv")
# file123 <- paste0(dir,"data/processed/Comparisons/",sp3,"_vs_",sp2,"_vs_",sp1,"/",sp3,"_vs_",sp2,"_vs_",sp1,".tsv")
# 
# 
# print("import comparison 12")
# comp12 <- read.table(file12, header = FALSE,colClasses=c("integer","integer")); 
# colnames(comp12) <- c("r","n")
# r <- comp12$r
# w <- comp12$n
# w <- w[r>min(rB)]
# r <- r[r>min(rB)]
# p12 <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
# print("import comparison 13")
# comp13 <- read.table(file13, header = FALSE,colClasses=c("integer","integer")); 
# colnames(comp13) <- c("r","n")
# r <- comp13$r
# w <- comp13$n
# w <- w[r>min(rB)]
# r <- r[r>min(rB)]
# p13 <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
# print("import comparison 13")
# comp23 <- read.table(file23, header = FALSE,colClasses=c("integer","integer")); 
# colnames(comp23) <- c("r","n")
# r <- comp23$r
# w <- comp23$n
# w <- w[r>min(rB)]
# r <- r[r>min(rB)]
# p23 <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
# print("import comparison 123")
# comp123 <- read.table(file123, header = FALSE,colClasses=c("integer","integer")); 
# colnames(comp123) <- c("r","n")
# r <- comp123$r
# w <- comp123$n
# w <- w[r>min(rB)]
# r <- r[r>min(rB)]
# p123 <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
# 
# pdf(paste0(dir,"data/processed/Comparisons/",sp3,"_vs_",sp2,"_vs_",sp1,"/",sp3,"_vs_",sp2,"_vs_",sp1,"_MLD.pdf"),width=5,height=5)
# m12 <- p12$counts/diff(rB)/n1/n2
# plot(log10(rV),log10(m12),pch=1,ylim=c(-15,-0),xlim=log10(c(250,2e5)),cex=0.5,col="blue")
# mT <- 2*3e6*rV^(-3)
# # mT <- mT/mT[which.min(abs(rV-2000))]*m12[which.min(abs(rV-2000))]
# lines(log10(rV),log10(mT),col="blue")
# m13 <- p13$counts/diff(rB)/n1/n3
# points(log10(rV),log10(m13),pch=3,cex=0.5,col="red")
# mT <- 2*3e6*rV^(-3)
# # mT <- mT/mT[which.min(abs(rV-2000))]*m13[which.min(abs(rV-2000))]
# lines(log10(rV),log10(mT),col="red")
# m23 <- p23$counts/diff(rB)/n2/n3
# points(log10(rV),log10(m23),pch=3,cex=0.5,col="green")
# mT <- 2*3e6*rV^(-3)
# # mT <- mT/mT[which.min(abs(rV-2000))]*m13[which.min(abs(rV-2000))]
# lines(log10(rV),log10(mT),col="green")
# 
# m123 <- p123$counts/diff(rB)/n1/n2/n3
# points(log10(p123$mids),log10(m123),pch=2,cex=0.5,col="black")
# mT <- 6*3e6*36*rV^(-4)
# # mT <- mT/mT[which.min(abs(rV-2000))]*m123[which.min(abs(rV-2000))]
# lines(log10(rV),log10(mT),col="black")
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #################################################################################################################################
# 
# 
# 
# 
# 
# 
# 
# 
# 
