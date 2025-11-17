#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
sp1 <- args[1]
sp2 <- args[2]
rTh <- args[3]
.libPaths("/home/projects/feariel/mishash/R/x86_64-pc-linux-gnu-library/4.4")
Sys.setenv(R_LIBS_USER = "/home/projects/feariel/mishash/R/x86_64-pc-linux-gnu-library/4.4")
#sp1 <- "Cronobacter"; sp2 <- "Raoultella"

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
# library(seqinr);
library(questionr)
library(tidyr)
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp(paste0("/home/projects/feariel/mishash/Development/ThreeFold/src/intersect_weighted.cpp"))


# dir <- "/home/misha/Documents/Development/ThreeFold/"
dir <- "/home/projects/feariel/mishash/Development/ThreeFold/"
Simulation <- "" # for empirical
#Simulation <- "PlasmidScope/" # for plasmids
# Simulation <- "" # for empirical

if (1==0)
{
  # unzip
  system(paste0("gunzip ", dir,"data/external/",Simulation,sp1,"/*.fa.gz"),wait=TRUE)
  system(paste0("gunzip ", dir,"data/external/",Simulation,sp2,"/*.fa.gz"),wait=TRUE)
  
  # remove plasmids headers, concatenate with "n" break
  if (!file.exists(paste0(dir,"data/external/",Simulation,sp1,".fa"))) {system(paste0("sed -i -r '1 ! s/^(>[^ ]+) .*/n/' ",dir,"data/external/",Simulation,sp1,"/*.fa"),wait=TRUE)}
  if (!file.exists(paste0(dir,"data/external/",Simulation,sp2,".fa"))) {system(paste0("sed -i -r '1 ! s/^(>[^ ]+) .*/n/' ",dir,"data/external/",Simulation,sp2,"/*.fa"),wait=TRUE)}
  
  # concatenate all files to one
  if (!file.exists(paste0(dir,"data/external/",Simulation,sp1,".fa"))) {system(paste0("cat ",dir,"data/external/",Simulation,sp1,"/*.fa > ",dir,"data/external/",Simulation,sp1,".fa"),wait=TRUE)}
  if (!file.exists(paste0(dir,"data/external/",Simulation,sp2,".fa"))) {system(paste0("cat ",dir,"data/external/",Simulation,sp2,"/*.fa > ",dir,"data/external/",Simulation,sp2,".fa"),wait=TRUE)}
}
gc()

# find matches between two files
# soft <- "/home/misha/Documents/Development/Software/bfmem/bfmem "
# soft <- "/home/misha/Documents/Development/Software/copmem2-main/copmem2 "
soft <- "/home/projects/feariel/mishash/Software/bfmem/bfmem"

system(paste0("mkdir -p ",dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/"),wait=TRUE)
Command <- paste0(soft," -s b -t 1 -l ", rTh," -o ",dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".tsv -r ",dir,"data/external/",Simulation,sp1,".fna"," -q ",dir,"data/external/",Simulation,sp2,".fna")
system(Command,wait=TRUE)

print("remove query headers with > (do only once and the comment out)")
system(paste0("sed -i '/^>/d' ",paste0(dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".tsv")),wait=TRUE)

comp <- read.table(paste0(dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".tsv"), header = FALSE,  colClasses=c("character", "integer","integer","integer")); colnames(comp) <- c("header","start","foo","length")
comp <- comp[,c(1,2,4)]
comp <- (comp %>% group_by_at(vars(header,start,length)) %>% summarise(weight = n(),, .groups="keep"))
write.table(comp,paste0(dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,"_weighted.tsv"),append=FALSE,row.names=FALSE,col.names=TRUE,sep = "\t",quote=FALSE)

rTable <- wtd.table(as.integer(comp$length),w=as.numeric(comp$weight))
write.table(data.frame(r=as.integer(names(rTable)),m=as.integer(rTable)),paste0(dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".tsv"),append=FALSE,row.names=FALSE,col.names=FALSE,sep = "\t",quote=FALSE)

comp <- data.table(start=as.integer(comp$start),end=as.integer(comp$start+comp$length-1),weight=as.numeric(comp$weight),header=comp$header)
comp <- split(comp[,1:3] , f = comp$header )
system(paste0("rm -f ",dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".h5"),wait=TRUE)
h5save(comp, file=paste0(dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".h5"))



# ONLY TO SAVE INTERSECTS
group1 <- h5ls(paste0(dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".h5"))$group[2]
headers1 <- h5ls(paste0(dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".h5"))$name; headers <- unique(headers1)
headers <- unique(headers[!headers %in% c("comp1","comp2","comp")])

registerDoParallel(1)
intersects <- foreach(header=headers, .inorder=FALSE, .combine=rbind) %dopar%
{
  intervals1 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".h5"), paste0(group1,"/",header))))
  # intervals1 <- (intervals1 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  
  # intersects <- as.data.table(intersect_weighted(as.matrix(intervals1), as.matrix(intervals1)))
  intersects <- intervals1; colnames(intersects) <- c("V1","V2","V3")
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  intersects$header <- header
  intersects$length <- as.integer(intersects$V2-intersects$V1+1)
  intersects
}
intersects <- data.frame(header=intersects$header,start=intersects$V1,stop=intersects$V2,w=intersects$w,length=intersects$length)
intersects <- intersects[intersects$length>=rTh,]
intersects <- intersects[with(intersects, order(length, -w)),]
write.table(intersects,paste0(dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,"_intersects.bed"),append=FALSE,row.names=FALSE,col.names=FALSE,sep = "\t",quote=FALSE)

stop_here

seqs <- read.fasta(paste0(dir,"/data/external/",Simulation,sp1,".fa"))
locs <- read.delim(paste0(dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,"_intersects.bed"), header=F, sep="\t") # Tabs? Put sep="\t"
outf <- file(paste0(dir,"data/processed/Comparisons/",Simulation,sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,"_intersects.fa"), 'w')

doitall <- function(x) {
  seq_id <- x[[1]]; start <- x[[2]]; stop <- x[[3]]; seq <- seqs[[seq_id]]
  seqv   <- seq[start:stop]
  header <- paste(sep="", ">", attr(seq, "name"), "_", start, "_", stop, "\n")
  cat(file=outf, header)
  cat(file=outf, toupper(paste(sep="", collapse="", seqv)), "\n")
}
apply(locs, 1, doitall)
close(outf)













# 
# 
# 
# getHeaders <- function(files)
# {
#   headers <- foreach(i=1:length(files), .inorder=FALSE, .combine=c) %do%
#     {
#       headers <- readLines(con = files[i], n = 1, ok = TRUE, warn = TRUE, encoding = "unknown", skipNul = FALSE)
#       headers <- strsplit(headers," ")[[1]][1]
#       headers <- str_replace(headers,">","")
#       headers
#     }
#   return(headers)
# }
# files1 <- Sys.glob(file.path(paste0(dir,"data/external/",sp1), "*.fa"))
# headers1 <- getHeaders(files1)
# files2 <- Sys.glob(file.path(paste0(dir,"data/external/",sp2), "*.fa"))
# headers2 <- getHeaders(files2)
# 
# 
# grid <- expand.grid(headers1,headers2); colnames(grid) <- c("header1","header2")
# 
# lines <- readLines(paste0(dir,"data/processed/",sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".tsv"), n = -1)
# Ind1 <- which(substr(lines,1,1)==">")
# 
# dat <- as.data.frame(str_split_fixed(lines,"\t",4))
# colnames(dat) <- c("header1","pos1","pos2","r")
# dat$header1 <- substr(dat$header1,2,nchar(dat$header1))
# 
# dat$header2 <- NA
# dat$header2[Ind1] <- str_split_fixed(lines[Ind1]," ",3)[,2]
# dat <- dat %>% fill("header2", .direction = "down")
# dat <- dat[-Ind1,]
# write.table(dat,paste0(dir,"data/processed/",sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,"_positions.tsv"),append=FALSE,row.names=FALSE,col.names=TRUE,sep = "\t",quote=FALSE)
# 
# # plot coverage
# headersToPlot <- sample(dat$header1,10)
# pdf(paste0(dir,"data/processed/",sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,"_coverage_log10_",sp1,".pdf"),width=20,height=5)
# for (headerToPlot in headersToPlot)
# {
#   print(headerToPlot)
#   datToPlot <- dat[dat$header1==headerToPlot,]
#   datToPlot$pos1 <- as.integer(datToPlot$pos1)
#   datToPlot$r <- as.integer(datToPlot$r)
#   starts <- sort(unique(c(1,datToPlot$pos1,datToPlot$pos1+datToPlot$r)))
#   starts <- c(starts,max(max(starts)+1000,1e7))
#   maxcounts <- length(unique(dat$header2))
#   plot(0,0,xlim=c(1,max(starts)),ylim=c(0,maxcounts*1.5*0+4))
#   for (i in 1:(length(starts)-1))
#   {
#     counts <- sum(datToPlot$pos1<=starts[i] & (datToPlot$pos1+datToPlot$r)>=starts[i+1])
#     lines(c(starts[i],starts[i+1]-1),log10(1+c(counts,counts)))
#   }
# }
# dev.off()
# headersToPlot <- sample(dat$header2,10)
# pdf(paste0(dir,"data/processed/",sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,"_coverage_log10_",sp2,".pdf"),width=20,height=5)
# for (headerToPlot in headersToPlot)
# {
#   print(headerToPlot)
#   datToPlot <- dat[dat$header2==headerToPlot,]
#   datToPlot$pos2 <- as.integer(datToPlot$pos2)
#   datToPlot$r <- as.integer(datToPlot$r)
#   starts <- sort(unique(c(1,datToPlot$pos2,datToPlot$pos2+datToPlot$r)))
#   starts <- c(starts,max(max(starts)+1000,1e7))
#   maxcounts <- length(unique(dat$header1))
#   plot(0,0,xlim=c(1,max(starts)),ylim=c(0,log10(maxcounts*1.5)))
#   for (i in 1:(length(starts)-1))
#   {
#     counts <- sum(datToPlot$pos2<=starts[i] & (datToPlot$pos2+datToPlot$r)>=starts[i+1])
#     lines(c(starts[i],starts[i+1]-1),log10(1+c(counts,counts)))
#   }
# }
# dev.off()
# 
# dat_agg <- aggregate(as.integer(dat[,'r']), by=dat[,c("header1","header2")], sum);colnames(dat_agg)[3] <- "r" 
# grid$r <- 0
# dat_agg <- rbind(dat_agg,grid); colnames(dat_agg) <- c("header1","header2","r")
# dat_agg <- aggregate(dat_agg[,'r'], by=dat_agg[,c("header1","header2")], sum); colnames(dat_agg)[3] <- "r" 
# write.table(dat_agg,paste0(dir,"data/processed/",sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,"_agg.tsv"),append=FALSE,row.names=FALSE,col.names=TRUE,sep = "\t",quote=FALSE)
# 
# 
# 
# 
# 
# print("remove query headers with > (do only once and the comment out)")
# system(paste0("sed -i '/^>/d' ",paste0(dir,"data/processed/",sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".tsv")),wait=TRUE)
# 
# comp <- read.table(paste0(dir,"data/processed/",sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".tsv"), header = FALSE,  colClasses=c("character", "integer","integer","integer")); colnames(comp) <- c("header","start","foo","length")
# rTable <- table(comp$length)
# write.table(data.frame(r=as.integer(names(rTable)),m=as.integer(rTable)),paste0(dir,"data/processed/",sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,"_table.tsv"),append=FALSE,row.names=FALSE,col.names=FALSE,sep = "\t",quote=FALSE)
# comp <- data.table(start1=as.integer(comp$start),end1=as.integer(comp$start+comp$length-1),header=comp$header)
# comp <- split(comp[,1:2] , f = comp$header )
# system(paste0("rm -f ",dir,"data/processed/",sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".h5"),wait=TRUE)
# h5save(comp, file=paste0(dir,"data/processed/",sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".h5"))
# 
# 
