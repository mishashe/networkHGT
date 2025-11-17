#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
sp1 <- args[1]
sp2 <- args[2]
sp3 <- args[3]
sp4 <- args[4]
sp5 <- args[5]
sp6 <- args[6]

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

registerDoParallel(2)
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


system(paste0("mkdir -p ",dir,"data/processed/Comparisons/",Simulation,sp6,"_vs_",sp5,"_vs_",sp4,"_vs_",sp3,"_vs_",sp2,"_vs_",sp1,"/"))
group1 <- h5ls(paste0(dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".h5"))$group[2]
group2 <- h5ls(paste0(dir,"data/processed/Comparisons/",Simulation,sp3,"_vs_",sp1,"/",sp3,"_vs_",sp1,".h5"))$group[2]
group3 <- h5ls(paste0(dir,"data/processed/Comparisons/",Simulation,sp4,"_vs_",sp1,"/",sp4,"_vs_",sp1,".h5"))$group[2]
group4 <- h5ls(paste0(dir,"data/processed/Comparisons/",Simulation,sp5,"_vs_",sp1,"/",sp5,"_vs_",sp1,".h5"))$group[2]
group5 <- h5ls(paste0(dir,"data/processed/Comparisons/",Simulation,sp6,"_vs_",sp1,"/",sp6,"_vs_",sp1,".h5"))$group[2]


headers1 <- h5ls(paste0(dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".h5"))$name; headers <- unique(headers1)
headers2 <- h5ls(paste0(dir,"data/processed/Comparisons/",Simulation,sp3,"_vs_",sp1,"/",sp3,"_vs_",sp1,".h5"))$name; headers <- unique(intersect(headers,headers2))
headers3 <- h5ls(paste0(dir,"data/processed/Comparisons/",Simulation,sp4,"_vs_",sp1,"/",sp4,"_vs_",sp1,".h5"))$name; headers <- unique(intersect(headers,headers3))
headers4 <- h5ls(paste0(dir,"data/processed/Comparisons/",Simulation,sp5,"_vs_",sp1,"/",sp5,"_vs_",sp1,".h5"))$name; headers <- unique(intersect(headers,headers4))
headers5 <- h5ls(paste0(dir,"data/processed/Comparisons/",Simulation,sp6,"_vs_",sp1,"/",sp6,"_vs_",sp1,".h5"))$name; headers <- unique(intersect(headers,headers5))

headers <- unique(headers[!headers %in% c("comp1","comp2","comp")])
rTable <- foreach(header=headers, .inorder=FALSE, .combine=combine_tables) %dopar%
{
  intervals1 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".h5"), paste0(group1,"/",header))))
  intervals2 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/Comparisons/",Simulation,sp3,"_vs_",sp1,"/",sp3,"_vs_",sp1,".h5"), paste0(group2,"/",header))))
  intervals3 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/Comparisons/",Simulation,sp4,"_vs_",sp1,"/",sp4,"_vs_",sp1,".h5"), paste0(group3,"/",header))))
  intervals4 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/Comparisons/",Simulation,sp5,"_vs_",sp1,"/",sp5,"_vs_",sp1,".h5"), paste0(group4,"/",header))))
  intervals5 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/Comparisons/",Simulation,sp6,"_vs_",sp1,"/",sp6,"_vs_",sp1,".h5"), paste0(group5,"/",header))))
  
  
  intervals1 <- (intervals1 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  intervals2 <- (intervals2 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  intervals3 <- (intervals3 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  intervals4 <- (intervals4 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  intervals5 <- (intervals5 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  
  
  intersects <- as.data.table(intersect_weighted(as.matrix(intervals1), as.matrix(intervals2)))
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  intersects <- as.data.table(intersect_weighted(as.matrix(intersects), as.matrix(intervals3)))
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  intersects <- as.data.table(intersect_weighted(as.matrix(intersects), as.matrix(intervals4)))
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  intersects <- as.data.table(intersect_weighted(as.matrix(intersects), as.matrix(intervals5)))
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  wtd.table(as.integer(intersects$V2-intersects$V1+1),w=as.numeric(intersects$w))
}
write.table(data.frame(r=as.integer(names(rTable)),m=as.numeric(rTable)),paste0(dir,"data/processed/Comparisons/",Simulation,sp6,"_vs_",sp5,"_vs_",sp4,"_vs_",sp3,"_vs_",sp2,"_vs_",sp1,"/",sp6,"_vs_",sp5,"_vs_",sp4,"_vs_",sp3,"_vs_",sp2,"_vs_",sp1,".tsv"),append=FALSE,row.names=FALSE,col.names=FALSE,sep = "\t",quote=FALSE)







