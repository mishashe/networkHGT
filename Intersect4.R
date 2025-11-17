#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
sp1 <- args[1]
sp2 <- args[2]
sp3 <- args[3]
sp4 <- args[4]

# sp1="Escherichia_coli"
# sp2="Klebsiella_pneumoniae"
# sp3="Salmonella_enterica"
# sp4="Serratia"



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


system(paste0("mkdir -p ",dir,"data/processed/Comparisons/",Simulation,sp4,"_vs_",sp3,"_vs_",sp2,"_vs_",sp1,"/"))
group1 <- h5ls(paste0(dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".h5"))$group[2]
group2 <- h5ls(paste0(dir,"data/processed/Comparisons/",Simulation,sp3,"_vs_",sp1,"/",sp3,"_vs_",sp1,".h5"))$group[2]
group3 <- h5ls(paste0(dir,"data/processed/Comparisons/",Simulation,sp4,"_vs_",sp1,"/",sp4,"_vs_",sp1,".h5"))$group[2]

headers1 <- h5ls(paste0(dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".h5"))$name; headers <- unique(headers1)
headers2 <- h5ls(paste0(dir,"data/processed/Comparisons/",Simulation,sp3,"_vs_",sp1,"/",sp3,"_vs_",sp1,".h5"))$name; headers <- unique(intersect(headers,headers2))
headers3 <- h5ls(paste0(dir,"data/processed/Comparisons/",Simulation,sp4,"_vs_",sp1,"/",sp4,"_vs_",sp1,".h5"))$name; headers <- unique(intersect(headers,headers3))

headers <- unique(headers[!headers %in% c("comp1","comp2","comp")])
rTable <- foreach(header=headers, .inorder=FALSE, .combine=combine_tables) %dopar%
{
  intervals1 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".h5"), paste0(group1,"/",header))))
  intervals2 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/Comparisons/",Simulation,sp3,"_vs_",sp1,"/",sp3,"_vs_",sp1,".h5"), paste0(group2,"/",header))))
  intervals3 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/Comparisons/",Simulation,sp4,"_vs_",sp1,"/",sp4,"_vs_",sp1,".h5"), paste0(group3,"/",header))))
  intervals1 <- (intervals1 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  intervals2 <- (intervals2 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  intervals3 <- (intervals3 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  
  intersects <- as.data.table(intersect_weighted(as.matrix(intervals1), as.matrix(intervals2)))
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  intersects <- as.data.table(intersect_weighted(as.matrix(intersects), as.matrix(intervals3)))
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  wtd.table(as.integer(intersects$V2-intersects$V1+1),w=as.numeric(intersects$w))
}
write.table(data.frame(r=as.integer(names(rTable)),m=as.numeric(rTable)),paste0(dir,"data/processed/Comparisons/",Simulation,sp4,"_vs_",sp3,"_vs_",sp2,"_vs_",sp1,"/",sp4,"_vs_",sp3,"_vs_",sp2,"_vs_",sp1,".tsv"),append=FALSE,row.names=FALSE,col.names=FALSE,sep = "\t",quote=FALSE)

# 
# 
# sp1 <- "Escherichia_coli"
# sp2 <- "Vibrio"
# sp3 <- "Enterobacter_hormaechei"
# sp4 <- "Serratia"
# 
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
# rTh <- 1000
# dir <- "/scratch/ws1/msheinman-msheinman/ThreeFold/" # Afalina
# 
# n1 <- length(Sys.glob(file.path(paste0(dir,"data/external/",sp1), "*.fna")))
# n2 <- length(Sys.glob(file.path(paste0(dir,"data/external/",sp2), "*.fna")))
# n3 <- length(Sys.glob(file.path(paste0(dir,"data/external/",sp3), "*.fna")))
# n4 <- length(Sys.glob(file.path(paste0(dir,"data/external/",sp4), "*.fna")))
# 
# 
# print("Plotting MLDs")
# rB <- 10^seq(log10(rTh),6,0.05) 
# rV <- 0.5*(rB[-length(rB)]+rB[-1])
# 
# 
# file12 <- paste0(dir,"data/processed/Comparisons/",sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".tsv")
# file13 <- paste0(dir,"data/processed/Comparisons/",sp3,"_vs_",sp1,"/",sp3,"_vs_",sp1,".tsv")
# file14 <- paste0(dir,"data/processed/Comparisons/",sp4,"_vs_",sp1,"/",sp4,"_vs_",sp1,".tsv")
# file23 <- paste0(dir,"data/processed/Comparisons/",sp3,"_vs_",sp2,"/",sp3,"_vs_",sp2,".tsv")
# file24 <- paste0(dir,"data/processed/Comparisons/",sp4,"_vs_",sp2,"/",sp4,"_vs_",sp2,".tsv")
# file34 <- paste0(dir,"data/processed/Comparisons/",sp4,"_vs_",sp3,"/",sp4,"_vs_",sp3,".tsv")
# 
# file123 <- paste0(dir,"data/processed/Comparisons/",sp3,"_vs_",sp2,"_vs_",sp1,"/",sp3,"_vs_",sp2,"_vs_",sp1,".tsv")
# file124 <- paste0(dir,"data/processed/Comparisons/",sp4,"_vs_",sp2,"_vs_",sp1,"/",sp4,"_vs_",sp2,"_vs_",sp1,".tsv")
# file134 <- paste0(dir,"data/processed/Comparisons/",sp4,"_vs_",sp3,"_vs_",sp1,"/",sp4,"_vs_",sp3,"_vs_",sp1,".tsv")
# file234 <- paste0(dir,"data/processed/Comparisons/",sp4,"_vs_",sp3,"_vs_",sp2,"/",sp4,"_vs_",sp3,"_vs_",sp2,".tsv")
# 
# file1234 <- paste0(dir,"data/processed/Comparisons/",sp4,"_vs_",sp3,"_vs_",sp2,"_vs_",sp1,"/",sp4,"_vs_",sp3,"_vs_",sp2,"_vs_",sp1,".tsv")
# 
# {
#   legends <- c()
#   print("import comparison 12")
#   comp12 <- read.table(file12, header = FALSE,colClasses=c("integer","numeric")); 
#   colnames(comp12) <- c("r","n")
#   r <- comp12$r
#   w <- comp12$n
#   w <- w[r>=rTh]/n1/n2
#   r <- r[r>=rTh]
#   a12 <- 1 + sum(w)/sum(w*log(r/rTh))
#   sa12 <- (a12-1)/sum(w)^0.5
#   legends <- c(legends,paste0(12,"_",round(a12,2),"+/-",round(sa12,3)))
#   p12 <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
#   print("import comparison 13")
#   comp13 <- read.table(file13, header = FALSE,colClasses=c("integer","numeric")); 
#   colnames(comp13) <- c("r","n")
#   r <- comp13$r
#   w <- comp13$n
#   w <- w[r>=rTh]/n1/n3
#   r <- r[r>=rTh]
#   a13 <- 1 + sum(w)/sum(w*log(r/rTh))
#   sa13 <- (a13-1)/sum(w)^0.5
#   legends <- c(legends,paste0(13,"_",round(a13,2),"+/-",round(sa13,3)))
#   p13 <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
#   print("import comparison 14")
#   comp14 <- read.table(file14, header = FALSE,colClasses=c("integer","numeric")); 
#   colnames(comp14) <- c("r","n")
#   r <- comp14$r
#   w <- comp14$n
#   w <- w[r>=rTh]
#   r <- r[r>=rTh]
#   a14 <- 1 + sum(w)/sum(w*log(r/rTh))
#   sa14 <- (a14-1)/sum(w)^0.5
#   legends <- c(legends,paste0(14,"_",round(a14,2),"+/-",round(sa14,3)))
#   p14 <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
#   print("import comparison 23")
#   comp23 <- read.table(file23, header = FALSE,colClasses=c("integer","numeric")); 
#   colnames(comp23) <- c("r","n")
#   r <- comp23$r
#   w <- comp23$n
#   w <- w[r>=rTh]
#   r <- r[r>=rTh]
#   a23 <- 1 + sum(w)/sum(w*log(r/rTh))
#   sa23 <- (a23-1)/sum(w)^0.5
#   legends <- c(legends,paste0(23,"_",round(a23,2),"+/-",round(sa23,3)))
#   p23 <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
#   print("import comparison 24")
#   comp24 <- read.table(file24, header = FALSE,colClasses=c("integer","numeric")); 
#   colnames(comp24) <- c("r","n")
#   r <- comp24$r
#   w <- comp24$n
#   w <- w[r>=rTh]
#   r <- r[r>=rTh]
#   a24 <- 1 + sum(w)/sum(w*log(r/rTh))
#   sa24 <- (a24-1)/sum(w)^0.5
#   legends <- c(legends,paste0(24,"_",round(a24,2),"+/-",round(sa24,3)))
#   p24 <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
#   print("import comparison 34")
#   comp34 <- read.table(file34, header = FALSE,colClasses=c("integer","numeric")); 
#   colnames(comp34) <- c("r","n")
#   r <- comp34$r
#   w <- comp34$n
#   w <- w[r>=rTh]
#   r <- r[r>=rTh]
#   a34 <- 1 + sum(w)/sum(w*log(r/rTh))
#   sa34 <- (a34-1)/sum(w)^0.5
#   legends <- c(legends,paste0(34,"_",round(a34,2),"+/-",round(sa34,3)))
#   p34 <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
# }
# {
#   print("import comparison 23")
#   comp23 <- read.table(file23, header = FALSE,colClasses=c("integer","numeric")); 
#   colnames(comp23) <- c("r","n")
#   r <- comp23$r
#   w <- comp23$n
#   w <- w[r>=rTh]
#   r <- r[r>=rTh]
#   a23 <- 1 + sum(w)/sum(w*log(r/rTh))
#   sa23 <- (a23-1)/sum(w)^0.5
#   legends <- c(legends,paste0(23,"_",round(a23,2),"+/-",round(sa23,3)))
#   p23 <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
#   print("import comparison 123")
#   comp123 <- read.table(file123, header = FALSE,colClasses=c("integer","numeric")); 
#   colnames(comp123) <- c("r","n")
#   r <- comp123$r
#   w <- comp123$n
#   w <- w[r>=rTh]
#   r <- r[r>=rTh]
#   a123 <- 1 + sum(w)/sum(w*log(r/rTh))
#   sa123 <- (a123-1)/sum(w)^0.5
#   legends <- c(legends,paste0(123,"_",round(a123,2),"+/-",round(sa123,3)))
#   p123 <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
#   print("import comparison 124")
#   comp124 <- read.table(file124, header = FALSE,colClasses=c("integer","numeric")); 
#   colnames(comp124) <- c("r","n")
#   r <- comp124$r
#   w <- comp124$n
#   w <- w[r>=rTh]
#   r <- r[r>=rTh]
#   a124 <- 1 + sum(w)/sum(w*log(r/rTh))
#   legends <- c(legends,paste0(124,"_",round(a124,2)))
#   p124 <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
#   print("import comparison 134")
#   comp134 <- read.table(file134, header = FALSE,colClasses=c("integer","numeric")); 
#   colnames(comp134) <- c("r","n")
#   r <- comp134$r
#   w <- comp134$n
#   w <- w[r>=rTh]
#   r <- r[r>=rTh]
#   a134 <- 1 + sum(w)/sum(w*log(r/rTh))
#   legends <- c(legends,paste0(134,"_",round(a134,2)))
#   p134 <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
#   print("import comparison 234")
#   comp234 <- read.table(file234, header = FALSE,colClasses=c("integer","numeric")); 
#   colnames(comp234) <- c("r","n")
#   r <- comp234$r
#   w <- comp234$n
#   w <- w[r>=rTh]
#   r <- r[r>=rTh]
#   a234 <- 1 + sum(w)/sum(w*log(r/rTh))
#   legends <- c(legends,paste0(234,"_",round(a234,2)))
#   p234 <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
# }
# {
#   print("import comparison 1234")
#   comp1234 <- read.table(file1234, header = FALSE,colClasses=c("integer","numeric")); 
#   colnames(comp1234) <- c("r","n")
#   r <- comp1234$r
#   w <- comp1234$n
#   w <- w[r>=rTh]
#   r <- r[r>=rTh]
#   a1234 <- 1 + sum(w)/sum(w*log(r/rTh))
#   legends <- c(legends,paste0(1234,"_",round(a1234,2)))
#   p1234 <- weighted.hist(x=r,breaks=rB,w=w,plot=FALSE)
# }
# 
# 
# {
#   pdf(paste0(dir,"data/processed/Comparisons/",sp4,"_vs_",sp3,"_vs_",sp2,"_vs_",sp1,"/",sp4,"_vs_",sp3,"_vs_",sp2,"_vs_",sp1,"_MLD.pdf"),width=5,height=5)
#   m12 <- p12$counts/diff(rB)/n1/n2; m12 <- m12/m12[1]*1
#   plot(log10(rV),log10(m12),pch=2,ylim=c(-15,-0),xlim=log10(c(rTh-100,2e5)),cex=0.5,col="blue")
#   m13 <- p13$counts/diff(rB)/n2/n3; m13 <- m13/m13[1]*1
#   points(log10(rV),log10(m13),pch=2,cex=0.5,col="blue")
#   m14 <- p14$counts/diff(rB)/n2/n3; m14 <- m14/m14[1]*1
#   points(log10(rV),log10(m14),pch=2,cex=0.5,col="green")
#   m23 <- p23$counts/diff(rB)/n2/n3; m23 <- m23/m23[1]*1
#   points(log10(rV),log10(m23),pch=2,cex=0.5,col="orange")
#   m24 <- p24$counts/diff(rB)/n2/n3; m24 <- m24/m24[1]*1
#   points(log10(rV),log10(m24),pch=2,cex=0.5,col="black")
#   m34 <- p34$counts/diff(rB)/n2/n3; m34 <- m34/m34[1]*1
#   points(log10(rV),log10(m34),pch=2,cex=0.5,col="purple")
#   mT <- exp(mean(log(m12[rV<1e4])))*(rV/exp(mean(log(rV[rV<1e4]))))^(-3)
#   lines(log10(rV),log10(mT),col="black")
#   
#   m123 <- p123$counts/diff(rB)/n1/n2/n3
#   points(log10(p123$mids),log10(m123),pch=3,cex=0.5,col="blue")
#   m124 <- p124$counts/diff(rB)/n1/n2/n4
#   points(log10(p124$mids),log10(m124),pch=3,cex=0.5,col="red")
#   m134 <- p134$counts/diff(rB)/n1/n3/n4
#   points(log10(p134$mids),log10(m134),pch=3,cex=0.5,col="green")
#   m234 <- p234$counts/diff(rB)/n2/n3/n4
#   mT <- exp(mean(log(m234[m234>0])))*(rV/exp(mean(log(rV[m234>0]))))^(-4)
#   lines(log10(rV),log10(mT),col="black")
#   points(log10(p234$mids),log10(m234),pch=3,cex=0.5,col="yellow")
#   lines(log10(rV),log10(mT),col="grey")
#   
#   
#   m1234 <- p1234$counts/diff(rB)/n1/n2/n3/n4
#   points(log10(p1234$mids),log10(m1234),pch=1,cex=0.5,col="black")
#   mT <- exp(mean(log(m1234[m1234>0])))*(rV/exp(mean(log(rV[m1234>0]))))^(-5)
#   lines(log10(rV),log10(mT),col="black")
#   legend(3, -6, legends,cex=0.5)
#   dev.off()
# }
# 








#################################################################################################################################









