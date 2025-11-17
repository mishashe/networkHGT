# plot alpha----
plotAlpha <- function()
{
  What <- ""
  estimateAlpha <- function(r,w,fold)
  {
    # return(1+sum(w)/sum(w*log(r/(rTh+0.5))))
    # Ind <- which(w>(min(w)*100))
    # z <- lm(log(w[Ind]) ~ log(r[Ind]))
    # return(-z$coefficients[2])
    
    
    a <- 0
    for (k in 1:(fold-1))
    {
      Mk_w <- sum(w * r^k) / sum(w)
      ak <- 1 + (k * Mk_w) / (Mk_w - rTh^k)
      a <- a + ak/(fold-1)
    }
    a <- (a + 1+sum(w)/sum(w*log(r/(rTh+0.5))))/2

    return(a)
  }
  species=c("Cronobacter","Raoultella","Proteus","Serratia","Citrobacter","Vibrio","Enterobacter","Salmonella","Klebsiella","Escherichia")
  remove <- c(1,6); species=species[-remove]
  
  nSp <- length(species)
  nG <- 0
  alphadat <- data.frame()
  # 2-fold
  fold <- 2
  for (i in 1:(nSp-1-nG))
  {
    for (j in (i+1):(nSp-nG))
    {
      Comparison <- paste0(species[j],"_vs_",species[i])
      file <- paste0(dir,"data/processed/Comparisons/",What,"/",Comparison,".tsv")
      comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
      colnames(comp) <- c("r","n")
      r <- comp$r
      w <- comp$n
      w <- w[r>=rTh]
      r <- r[r>=rTh]
      alpha <- estimateAlpha(r,w,fold)
      alphadat <- rbind(alphadat,data.frame(n=fold,alpha=alpha, pair=Comparison))
    }
  }
  # 3-fold
  if (nSp-nG>=3)
  {
    fold <- 3
    for (i in 1:(nSp-nG-2))
    {
      for (j in (i+1):(nSp-nG-1))
      {
        for (k in (j+1):(nSp-nG))
        {
          
          Comparison <- paste0(species[k],"_vs_",species[j],"_vs_",species[i])
          file <- paste0(dir,"data/processed/Comparisons/",What,"/",Comparison,".tsv")
          comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
          colnames(comp) <- c("r","n")
          r <- comp$r
          w <- comp$n
          w <- w[r>=rTh]
          r <- r[r>=rTh]
          alpha <- estimateAlpha(r,w,fold)
          alphadat <- rbind(alphadat,data.frame(n=fold,alpha=alpha, pair=Comparison))
        }
      }
    }
  }
  # 4-fold
  if ((nSp-nG)>=4)
  {
    fold <- 4
    for (i in 1:(nSp-nG-3))
    {
      for (j in (i+1):(nSp-nG-2))
      {
        for (k in (j+1):(nSp-nG-1))
        {
          for (l in (k+1):(nSp-nG))
          {
            Comparison <- paste0(species[l],"_vs_",species[k],"_vs_",species[j],"_vs_",species[i])
            file <- paste0(dir,"data/processed/Comparisons/",What,"/",Comparison,".tsv")
            comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
            colnames(comp) <- c("r","n")
            r <- comp$r
            w <- comp$n
            w <- w[r>=rTh]
            r <- r[r>=rTh]
            alpha <- estimateAlpha(r,w,fold)
            alphadat <- rbind(alphadat,data.frame(n=fold,alpha=alpha, pair=Comparison))
          }
        }
      }
    }
  }
  # 5-fold
  if ((nSp-nG)>=5)
  {
    fold <- 5
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
              file <- paste0(dir,"data/processed/Comparisons/",What,"/",Comparison,".tsv")
              comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
              colnames(comp) <- c("r","n")
              r <- comp$r
              w <- comp$n
              w <- w[r>=rTh]
              r <- r[r>=rTh]
              alpha <- estimateAlpha(r,w,fold)
              alphadat <- rbind(alphadat,data.frame(n=fold,alpha=alpha, pair=Comparison))
            }
          }
        }
      }
    }
  }
  # # 6-fold
  if ((nSp-nG)>=6)
  {
    fold <- 6
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
                file <- paste0(dir,"data/processed/Comparisons/",What,"/",Comparison,".tsv")
                comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
                colnames(comp) <- c("r","n")
                r <- comp$r
                w <- comp$n
                w <- w[r>=rTh]
                r <- r[r>=rTh]
                alpha <- estimateAlpha(r,w,fold)
                alphadat <- rbind(alphadat,data.frame(n=fold,alpha=alpha, pair=Comparison))
              }
            }
          }
        }
      }
    }
  }
  # 7-fold
  if ((nSp-nG)>=7)
  {
    fold <- 7
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
                  file <- paste0(dir,"data/processed/Comparisons/",What,"/",Comparison,".tsv")
                  comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
                  colnames(comp) <- c("r","n")
                  r <- comp$r
                  w <- comp$n
                  w <- w[r>=rTh]
                  r <- r[r>=rTh]
                  alpha <- estimateAlpha(r,w,fold)
                  alphadat <- rbind(alphadat,data.frame(n=fold,alpha=alpha, pair=Comparison))
                }
              }
            }
          }
        }
      }
    }
  }
  # 8-fold
  if ((nSp-nG)>=8)
  {
    fold <- 8
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
                    file <- paste0(dir,"data/processed/Comparisons/",What,"/",Comparison,".tsv")
                    comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
                    colnames(comp) <- c("r","n")
                    r <- comp$r
                    w <- comp$n
                    w <- w[r>=rTh]
                    r <- r[r>=rTh]
                    alpha <- estimateAlpha(r,w,fold)
                    alphadat <- rbind(alphadat,data.frame(n=fold,alpha=alpha, pair=Comparison))
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  # # 9-fold
  if ((nSp-nG)>=9)
  {
    fold <- 9
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
                      file <- paste0(dir,"data/processed/Comparisons/",What,"/",Comparison,".tsv")
                      comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
                      colnames(comp) <- c("r","n")
                      r <- comp$r
                      w <- comp$n
                      w <- w[r>=rTh]
                      r <- r[r>=rTh]
                      alpha <- estimateAlpha(r,w,fold)
                      alphadat <- rbind(alphadat,data.frame(n=fold,alpha=alpha, pair=Comparison))
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
  # # 10-fold
  if ((nSp-nG)>=10)
  {
    fold <- 10
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
                        file <- paste0(dir,"data/processed/Comparisons/",What,"/",Comparison,".tsv")
                        comp <- read.table(file, header = FALSE,colClasses=c("integer","numeric")); 
                        colnames(comp) <- c("r","n")
                        r <- comp$r
                        w <- comp$n
                        w <- w[r>=rTh]
                        r <- r[r>=rTh]
                        alpha <- estimateAlpha(r,w,fold)
                        alphadat <- rbind(alphadat,data.frame(n=fold,alpha=alpha, pair=Comparison))
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
  
  pdf(paste0(dir,"data/processed/",dirname,"/alphaEmpirical.pdf"),width=5,height=5)
  p <- ggplot(alphadat[alphadat$n<18,], aes(x = as.factor(n), y = alpha)) + 
    geom_violin(outliers = FALSE,fill="grey", color = "black") +
    geom_jitter(aes(color = as.factor(n), shape = as.factor(n)), size = 1, position=position_jitter(0.2),cex=0.5) +
    scale_color_manual(values = c("2" = "black", "3" = "red", "4" = "blue", "5" = "green", "6" = "purple", "7" = "orange", "8" = "brown")) +
    scale_shape_manual(values = c("2" = 16, "3" = 1, "4" = 17, "5" = 2, "6" = 15, "7" = 0, "8" = 3)) +
    scale_y_continuous(breaks=seq(0,12,1), expand = c(0, 0),limits=c(2,10.25)) +
    scale_x_discrete() + 
    theme_bw() +
    theme(panel.grid.major = element_line(color = "grey",linewidth = 0.25,linetype = 3)) +
    geom_abline(intercept = 2, slope = 1, color="grey",  linetype="dashed", linewidth=0.5) +
    geom_abline(intercept = 3, slope = 1, color="grey",  linetype="solid", linewidth=0.5)
  print(p)
  dev.off()
  write.table(alphadat,paste0(dir,"data/processed/",dirname,"/alphaEmpirical.tsv"), row.names = FALSE,col.names = TRUE)
  alphadat <- read.table(paste0(dir,"data/processed/",dirname,"/alphaEmpirical.tsv"),header=TRUE)

  pdf(paste0(dir,"data/processed/",dirname,"/alphaHist.pdf"),width=1.5,height=1.5)
  for (n in 2:10)
  {
    p <- ggplot(alphadat[alphadat$n==n,], aes(x=alpha)) +
      geom_histogram(bins = 10,color="black", fill="white") + 
      geom_vline(aes(xintercept=median(alpha)),color="blue", linetype="dashed", linewidth=1) +
      xlim((median(alphadat[alphadat$n==n,"alpha"]))-1.5, (median(alphadat[alphadat$n==n,"alpha"]))+1.5) +
      theme_classic()
    
    print(p)
  }
  dev.off()
  
}

# function to plot results----
# plotResults <- function(dirname,name,rhoM,nSp,nG,MLDsim,MLDanalytic,MLDdirect,MLDghost,ijk_from_set=ijk_from_set,set_from_ijk=set_from_ijk,MLD)
# to run do:
# module load R
# R
# source("/home/projects/feariel/mishash/Development/ThreeFold/src/plotResults.R"); plotResults()

plotResults <- function()
{
  library(ggplot2)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(circlize)
  library(stringr)
  rTh <- 1000
  rB <- 10^seq(log10(rTh-0.5),6,0.05) 
  rV <- 0.5*(rB[-length(rB)]+rB[-1])
  
  dir <- "/home/projects/feariel/mishash/Development/ThreeFold/"
  dirname <- "addGammaSum"
  species=c("Cronobacter","Raoultella","Proteus","Serratia","Citrobacter","Vibrio","Enterobacter","Salmonella","Klebsiella","Escherichia")
  nList <- readRDS(file = paste0(dir,"data/processed/global/nList.rds"))    
  
  remove <- c(1,6); species=species[-remove]; nList=nList[-remove]
  
  
  nSp <- length(species)
  nG <- 0
  dir <- "/home/projects/feariel/mishash/Development/ThreeFold/"
  # dir <- "/home/mishash@PHYSICS.weizmann.ac.il/Documents/Development/ThreeFold/" 
  
  dirname <- "addGammaSum"
  MLDanalytic <- readRDS( file = paste0(dir,"data/processed/",dirname,"/MLDanalytic.rds"))    
  MLDdirect <- readRDS(file = paste0(dir,"data/processed/",dirname,"/MLDdirect.rds"))   
  MLDghost <- readRDS(file = paste0(dir,"data/processed/",dirname,"/MLDghost.rds"))
  rhoM <- readRDS(file = paste0(dir,"data/processed/",dirname,"/rhoM.rds"))   
  ijk_from_set <- readRDS(file = paste0(dir,"data/processed/",dirname,"/ijk_from_set.rds"))  
  set_from_ijk <- readRDS(file = paste0(dir,"data/processed/",dirname,"/set_from_ijk.rds"))  
  MLD <- readRDS(file = paste0(dir,"data/processed/",dirname,"/MLD.rds"))
  gamma <- readRDS(file = paste0(dir,"data/processed/",dirname,"/gamma.rds"))
  L1 <- readRDS(file = paste0(dir,"data/processed/",dirname,"/L1.rds"))
  L2 <- readRDS(file = paste0(dir,"data/processed/",dirname,"/L2.rds"))
  
  
  library(igraph)
  pdf(paste0(dir,"data/processed/",dirname,"/network.pdf"),width=5,height=5)
  network <- graph_from_adjacency_matrix(rhoM/2, mode='plus', diag=F, weighted = TRUE)
  plot(network, layout=layout.circle, main="rho", edge.width = log10(1+10000*E(network)$weight/max(E(network)$weight)), edge.curved=0.25,edge.arrow.size = 0.3)
  plot(network, layout=layout.circle, main="rho", edge.width = (4*E(network)$weight/max(E(network)$weight)), edge.curved=0.25,edge.arrow.size = 0.3)
  w <- log10(1e-1+E(network)$weight)
  # col <- heat.colors(n = 100, rev = FALSE)[round(100*(w-min(w))/(max(w)-min(w)))]
  col <- brewer.pal(n = nSp, name = "YlOrRd")[1+round(8*(w-min(w))/(max(w)-min(w)))]
  plot(network, layout=layout.circle, main="rho", edge.width = log10(1+10000*E(network)$weight/max(E(network)$weight)), edge.curved=0,edge.color = col, edge.label = round(E(network)$weight, 2), vertex.label.cex=0.6, edge.label.cex=0.3)
  w <- E(network)$weight
  # col <- heat.colors(n = 100, rev = FALSE)[round(100*(w-min(w))/(max(w)-min(w)))]
  col <- brewer.pal(n = nSp+1, name = "YlOrRd")[1+round(8*(w-min(w))/(max(w)-min(w)))]
  plot(network, layout=layout.circle, main="rho", edge.width = w/max(w)*8, edge.curved=0,edge.color = col, edge.label = round(E(network)$weight, 2), vertex.label.cex=0.6, edge.label.cex=0.3)
  plot(network, layout=layout.circle, main="rho", edge.width = w/max(w)*8, edge.curved=0,edge.color = col, vertex.label.cex=0.6, edge.label.cex=0.3)

  rhoMh <- matrix(0,nSp+1,nSp+1)
  rhoMh[2:nrow(rhoMh),2:nrow(rhoMh)] <- rhoM
  rhoMh[1,2:nrow(rhoMh)] <- gamma/5
  rhoMh[2:nrow(rhoMh),1] <- gamma/5
  rownames(rhoMh) <- c("hub",species)
  colnames(rhoMh) <- c("hub",species)
  
  diag(rhoMh) <- 0
  for (i in 2:nrow(rhoMh))
  {
    rhoMh[i,i] <- rhoMh[i-1,i]*sum(rhoMh[i,-c(i,i-1)])/sum(rhoMh[i-1,-c(i,i-1)])
  }
  rhoMh[1,1] <- rhoMh[1+1,1]*sum(rhoMh[1,-c(1,2)])/sum(rhoMh[1+1,-c(1+1,1)])
  svd_res <- svd(rhoMh)
  u1 <- svd_res$u[, 1]
  v1 <- svd_res$v[, 1]
  s <- svd_res$d
  s1 <- svd_res$d[1]
  A_rank1 <- s1 * outer(u1, u1)
  s[1]^2 / sum(s^2)
  g <- -s1^0.5*u1
  names(g) <- rownames(rhoMh)

  # betterorder <- c( "Serratia", "Salmonella" ,"Proteus"    ,"Escherichia"  ,"Raoultella"  ,  "Enterobacter"   ,  "Citrobacter" ,  "Klebsiella" ,"hub" )
  # rhoMh <- rhoMh[betterorder,betterorder]
  
  network <- graph_from_adjacency_matrix(rhoMh/2, mode='plus', diag=F, weighted = TRUE)
  w <- log10(1e-1+E(network)$weight)
  w <- E(network)$weight
  col <- brewer.pal(n = 9, name = "YlOrRd")[3+round(6*(w-min(w))/(max(w)-min(w)))]
  coords <- layout_in_circle(network, order = c(1,5,6,9,7,4,3,2,8))
  rotate_layout <- function(layout, theta_degrees) {
    theta <- theta_degrees * pi / 180  # convert to radians
    rot_matrix <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2)
    t(rot_matrix %*% t(layout))
  }
  coords <- rotate_layout(coords, theta_degrees = 270)
  plot(network, layout=coords, main="rho", edge.width = w/max(w)*20, edge.curved=0,edge.color = col, edge.label = round(E(network)$weight, 2), vertex.label.cex=0.6, edge.label.cex=0.3)
  plot(network, layout=coords, main="rho", edge.width = w/max(w)*20, edge.curved=0,edge.color = col, vertex.label.cex=0.6, edge.label.cex=0.3)
  col <- brewer.pal(n = 9, name = "Blues")[3+round(6*(w-min(w))/(max(w)-min(w)))]
  plot(network, layout=coords, main="rho", edge.width = w/max(w)*20, edge.curved=0,edge.color = col, edge.label = round(E(network)$weight, 2), vertex.label.cex=0.6, edge.label.cex=0.3)
  plot(network, layout=coords, main="rho", edge.width = w/max(w)*20, edge.curved=0,edge.color = col, vertex.label.cex=0.6, edge.label.cex=0.3)
  
  barplot(A, horiz = TRUE, xlab = "X-axis",
          ylab = "Y-axis", main ="Horizontal Bar Chart"
  )
  
  nodes <- c("ghost", species)
  edges <- cbind(rep("ghost", nSp), species)  # Central node 1 to nodes 2-11
  weights <- gamma  # The rates represent the edge weights
  g <- graph_from_edgelist(edges, directed = TRUE)
  E(g)$weight <- weights
  w <- E(network)$weight
  # col <- heat.colors(n = 100, rev = FALSE)[round(100*(w-min(w))/(max(w)-min(w)))]
  col <- "blue"
  plot(g, 
       layout = layout.circle,  # Force-directed layout
       vertex.size = 10,
       vertex.label = nodes,
       vertex.label.cex = 0.5,
       edge.width = E(g)$weight * 0.01,
       main = "Network Visualization with Force-directed Layout",
       edge.color = col)
  
  p <- Heatmap(rhoM, cluster_rows = FALSE, cluster_columns = FALSE)
  plot(p)
  rhoNA <- rhoM
  rhoNA[nSp,] <- 1e-5
  rhoNA[,nSp] <- 1e-5
  dend <- hclust(as.dist(1/(rhoNA+1)))
  rhoNA <- rhoM
  rhoNA[rhoNA==0] <- NA
  col_fun = colorRamp2(c(0, 3), c( "blue", "red"))
  # col_fun(seq(-3, 3))
  p <- Heatmap(log10(rhoNA), cluster_columns = dend, cluster_rows = dend, col = col_fun)
  plot(p)
  col_fun = colorRamp2(c(-2, 4), c( "blue", "red"))
  # col_fun(seq(-3, 3))
  p <- Heatmap(log10(rhoNA), cluster_columns = dend, cluster_rows = dend, col = col_fun)
  plot(p)
  
  dev.off()
  
  diag(rhoM) <- 0
  
  
  
  alphaHill <- data.frame()
  # 2-fold with plots
  
  legends <- c()
  m2 <- rV*0
  colors <-  c("black", "red", "blue", "green", "orange", "purple", "cyan")
  markers <- 1:length(colors)
  
  fold <- 2
  pdf(paste0(dir,"data/processed/",dirname,"/MLD_",fold,".pdf"),height=6*3,width=5*3)
  mar <- 2
  par(mfrow=c(6,5),mar = c(mar, mar, mar, mar))
  for (i in 1:(nSp-nG-1))
  {
    for (j in (i+1):(nSp-nG))
    {
      Ind <- paste0(i,"_",j)
      s <- set_from_ijk[[fold]][Ind]
      ijk <- as.numeric(str_split(Ind,"_")[[1]])
      if (sum(MLD[[fold]][set_from_ijk[[fold]][Ind],])>0)
      {
        plot(log10(rV),log10(MLD[[fold]][s,]),pch=2,cex=0.5,col="blue",ylab="",xlim=c(3,log10(rV[which.max(which(MLD[[fold]][s,]>0))])), ann=FALSE)
        lines(log10(rV),log10(MLDanalytic[[fold]][s,]),col="blue")
        xlim <- par("usr")[1:2];ylim <- par("usr")[3:4]-1;text(x = xlim[2], y = ylim[2], labels = paste0(species[i],",",species[j]), pos = 2, xpd = TRUE)
        # lines(log10(rV),log10(MLDdirect[[fold]][s,]),col="grey",lty=2)
        # lines(log10(rV),log10(MLDghost[[fold]][s,]),col="grey")
        # lines(log10(rV),-log10(rV^(fold+1))+mean(log10(rV[rV<10^4.0 & MLD[[fold]][s,]>0]^(fold+1)))+mean(log10(MLD[[fold]][s,][rV<10^4.0 & MLD[[fold]][s,]>0])),col="black",lty=2)
      }
    }
  }
  dev.off()
  pdf(paste0(dir,"data/processed/",dirname,"/MLD_",fold,"_one.pdf"),height=10,width=5)
  plot(1, type = "n", xlim = c(1e3,1e5), ylim = c(1e-35,1e-1), log = "xy", axes = FALSE)
  for (i in 1:(nSp-nG-1))
  {
    for (j in (i+1):(nSp-nG))
    {
      Ind <- paste0(i,"_",j)
      s <- set_from_ijk[[fold]][Ind]
      ijk <- as.numeric(str_split(Ind,"_")[[1]])
      if (sum(MLD[[fold]][set_from_ijk[[fold]][Ind],])>0)
      {
        points((rV),MLD[[fold]][s,]/MLDanalytic[[fold]][s,1]/10^s,cex=0.5,col = colors[s %% length(colors) + 1], pch = markers[s %% length(markers) + 1])
        lines((rV),MLDanalytic[[fold]][s,]/MLDanalytic[[fold]][s,1]/10^s,col = colors[s %% length(colors) + 1])
      }
    }
  }
  xticks <- 10^(3:6)
  axis(1, at = xticks, labels = parse(text = paste0("10^", 3:6)))
  yticks <- 10^seq(-50,50,5)
  axis(2, at = yticks, labels = parse(text = paste0("10^", seq(-50,50,5))))
  box()
  dev.off()
  pdf(paste0(dir,"data/processed/",dirname,"/seq_freq_dist_",fold,".pdf"),height=6*3,width=5*3)
  mar <- 2
  par(mfrow=c(6,5),mar = c(mar, mar, mar, mar))
  for (i in 1:(nSp-nG-1))
  {
    for (j in (i+1):(nSp-nG))
    {
      Ind <- paste0(i,"_",j)
      s <- set_from_ijk[[fold]][Ind]
      ijk <- as.numeric(str_split(Ind,"_")[[1]])
      if (sum(MLD[[fold]][set_from_ijk[[fold]][Ind],])>0)
      {
        Comparison <- paste0(species[j],"_vs_",species[i])
        dat <- table(read.table(paste0("/home/projects/feariel/mishash/Development/ThreeFold/data/processed/Comparisons/",Comparison,"/mmseq2/cluster_cluster.tsv"), header = FALSE,fill=NA,sep="\t", comment.char="#",na.strings=".", stringsAsFactors=FALSE, quote="")$V1); 
        dat <- dat/nList[[i]]/nList[[j]]
        fB <- 10^seq(log10(min(dat))-0.2,log10(max(dat))+0.2,0.2)
        fV <- (fB[-1]*fB[-length(fB)])^0.5
        p <- hist(dat,breaks=fB,plot=FALSE)
        x <- log10(fV); y <- (p$counts/diff(fB)); 
        Ind <- which(y>0); y <- log10(y[Ind]); x <- x[Ind]
        plot(x,y)
        lines(x,mean(y)-1.5*(x-mean(x)))
        xlim <- par("usr")[1:2];ylim <- par("usr")[3:4]-1;text(x = xlim[2], y = ylim[2], labels = paste0(species[i],",",species[j]), pos = 2, xpd = TRUE)
      }
    }
  }
  dev.off()
  if (nSp-nG>=3)
  {
    # 3-fold
    m3 <- rV*0
    fold <- 3
    pdf(paste0(dir,"data/processed/",dirname,"/MLD_",fold,".pdf"),height=8*2.5,width=7*2.5)
    par(mfrow=c(8,7),mar = 2*c(1, 1, 1, 1))
    for (i in 1:(nSp-nG-2))
    {
      for (j in (i+1):(nSp-nG-1))
      {
        for (k in (j+1):(nSp-nG))
        {
          Ind <- paste0(i,"_",j,"_",k)
          s <- set_from_ijk[[fold]][Ind]
          ijk <- as.numeric(str_split(Ind,"_")[[1]])
          if (sum(MLD[[fold]][s,])>0)
          {
            # alpha <- 1 + sum(diff(rB)*MLD[[fold]][s,])/sum(diff(rB)*MLD[[fold]][s,]*log(rV/rTh))
            # alpha <- alphaHillraw[[fold]][s]
            # alphaTheory <- 1 + sum(diff(rB)*MLDsim[[fold]][s,])/sum(diff(rB)*MLDsim[[fold]][s,]*log(rV/rTh))
            # alphaTheory <- alpha1[[fold]][s]
            # alphaHill <- rbind(alphaHill,data.frame(Ind=Ind,n=fold,alpha,alphaTheory=alphaTheory,pair=paste0(species[i],"_vs_",species[j],"_vs_",species[k])))
            plot(log10(rV),log10(MLD[[fold]][s,]),pch=2,cex=0.5,col="red",ylab="",xlim=c(3,log10(rV[which.max(which(MLD[[fold]][s,]>0))])), ann=FALSE)
            lines(log10(rV),log10(MLDanalytic[[fold]][s,]),col="red")
            xlim <- par("usr")[1:2];ylim <- par("usr")[3:4]-1;text(x = xlim[2], y = ylim[2], labels = paste0(species[i],",",species[j],",\n",species[k]), pos = 2, xpd = TRUE)
            
            # lines(log10(rV),log10(MLDdirect[[fold]][s,]),col="grey",lty=2)
            # lines(log10(rV),log10(MLDghost[[fold]][s,]),col="grey")
          }
        }
      }
    }
    dev.off()
    pdf(paste0(dir,"data/processed/",dirname,"/MLD_",fold,"_one.pdf"),height=10,width=5)
    plot(1, type = "n", xlim = c(1e3,0.8e5), ylim = c(1e-63,1e-1), log = "xy", axes = FALSE)
    for (i in 1:(nSp-nG-2))
    {
      for (j in (i+1):(nSp-nG-1))
      {
        for (k in (j+1):(nSp-nG))
        {
          Ind <- paste0(i,"_",j,"_",k)
          s <- set_from_ijk[[fold]][Ind]
          ijk <- as.numeric(str_split(Ind,"_")[[1]])
          if (sum(MLD[[fold]][set_from_ijk[[fold]][Ind],])>0)
          {
            points((rV),MLD[[fold]][s,]/MLDanalytic[[fold]][s,1]/10^s,cex=0.5,col = colors[s %% length(colors) + 1], pch = markers[s %% length(markers) + 1])
            lines((rV),MLDanalytic[[fold]][s,]/MLDanalytic[[fold]][s,1]/10^s,col = colors[s %% length(colors) + 1])
            # lines(log10(rV),log10(MLDdirect[[fold]][s,]),col="grey",lty=2)
            # lines(log10(rV),log10(MLDghost[[fold]][s,]),col="grey")
            # lines(log10(rV),-log10(rV^(fold+1))+mean(log10(rV[rV<10^4.0 & MLD[[fold]][s,]>0]^(fold+1)))+mean(log10(MLD[[fold]][s,][rV<10^4.0 & MLD[[fold]][s,]>0])),col="black",lty=2)
          }
        }
      }
    }
    xticks <- 10^(3:6)
    axis(1, at = xticks, labels = parse(text = paste0("10^", 3:6)))
    yticks <- 10^seq(-500,500,10)
    axis(2, at = yticks, labels = parse(text = paste0("10^", seq(-500,500,10))))
    box()
    dev.off()
  }
  if (nSp-nG>=4)
  {
    # 4-fold
    fold <- 4
    pdf(paste0(dir,"data/processed/",dirname,"/MLD_",fold,".pdf"),height=10*2.5,width=7*2.5)
    par(mfrow=c(10,7),mar = 2*c(1, 1, 1, 1))
    for (i in 1:(nSp-nG-3))
    {
      for (j in (i+1):(nSp-nG-2))
      {
        for (k in (j+1):(nSp-nG-1))
        {
          for (l in (k+1):(nSp-nG))
          {
            Ind <- paste0(i,"_",j,"_",k,"_",l)
            s <- set_from_ijk[[fold]][Ind]
            ijk <- as.numeric(str_split(Ind,"_")[[1]])
            if (sum(MLD[[fold]][s,])>0)
            {
              plot(log10(rV),log10(MLD[[fold]][s,]),pch=2,cex=0.5,col="purple",ylab="",xlim=c(3,log10(rV[which.max(which(MLD[[fold]][s,]>0))])), ann=FALSE)
              lines(log10(rV),log10(MLDanalytic[[fold]][s,]),col="purple")
              xlim <- par("usr")[1:2];ylim <- par("usr")[3:4]-1;text(x = xlim[2], y = ylim[2], labels = paste0(species[i],",",species[j],",\n",species[k],",",species[l]), pos = 2, xpd = TRUE)
              
              # lines(log10(rV),log10(MLDdirect[[fold]][s,]),col="black",lty=2)
              # lines(log10(rV),log10(MLDghost[[fold]][s,]),col="grey")
            }
          }
        }
      }
    }
    dev.off()
    pdf(paste0(dir,"data/processed/",dirname,"/MLD_",fold,"_one.pdf"),height=10,width=5)
    plot(1, type = "n", xlim = c(1e3,0.8e5), ylim = c(1e-80,1e-1), log = "xy", axes = FALSE)
    for (i in 1:(nSp-nG-3))
    {
      for (j in (i+1):(nSp-nG-2))
      {
        for (k in (j+1):(nSp-nG-1))
        {
          for (l in (k+1):(nSp-nG))
          {
            Ind <- paste0(i,"_",j,"_",k,"_",l)
            s <- set_from_ijk[[fold]][Ind]
            ijk <- as.numeric(str_split(Ind,"_")[[1]])
            if (sum(MLD[[fold]][set_from_ijk[[fold]][Ind],])>0)
            {
              points((rV),MLD[[fold]][s,]/MLDanalytic[[fold]][s,1]/10^s,cex=0.5,col = colors[s %% length(colors) + 1], pch = markers[s %% length(markers) + 1])
              lines((rV),MLDanalytic[[fold]][s,]/MLDanalytic[[fold]][s,1]/10^s,col = colors[s %% length(colors) + 1])
            }
          }
        }
      }
    }
    xticks <- 10^(3:6)
    axis(1, at = xticks, labels = parse(text = paste0("10^", 3:6)))
    yticks <- 10^seq(-500,500,10)
    axis(2, at = yticks, labels = parse(text = paste0("10^", seq(-500,500,10))))
    box()
    dev.off()
  }
  if (nSp-nG>=5)
  {
    # 5-fold
    fold <- 5
    pdf(paste0(dir,"data/processed/",dirname,"/MLD_",fold,".pdf"),height=8*2.5,width=7*2.5)
    par(mfrow=c(8,7),mar = 2*c(1, 1, 1, 1))
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
              Ind <- paste0(i,"_",j,"_",k,"_",l,"_",mm)
              s <- set_from_ijk[[fold]][Ind]
              ijk <- as.numeric(str_split(Ind,"_")[[1]])
              if (sum(MLD[[fold]][s,])>0)
              {
                plot(log10(rV),log10(MLD[[fold]][s,]),pch=2,cex=0.5,col="orange",ylab="",xlim=c(3,log10(rV[which.max(which(MLD[[fold]][s,]>0))])), ann=FALSE)
                lines(log10(rV),log10(MLDanalytic[[fold]][s,]),col="orange")
                xlim <- par("usr")[1:2];ylim <- par("usr")[3:4]-1.5;text(x = xlim[2], y = ylim[2], labels = paste0(species[i],",",species[j],",\n",species[k],",",species[l],",\n",species[mm]), pos = 2, xpd = TRUE)
                
                # lines(log10(rV),log10(MLDdirect[[fold]][s,]),col="black",lty=2)
                # lines(log10(rV),log10(MLDghost[[fold]][s,]),col="grey")
              }
            }
          }
        }
      }
    }
    dev.off()
    pdf(paste0(dir,"data/processed/",dirname,"/MLD_",fold,"_one.pdf"),height=10,width=5)
    plot(1, type = "n", xlim = c(1e3,0.5e5), ylim = c(1e-66,1e-1), log = "xy", axes = FALSE)
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
              Ind <- paste0(i,"_",j,"_",k,"_",l,"_",mm)
              s <- set_from_ijk[[fold]][Ind]
              ijk <- as.numeric(str_split(Ind,"_")[[1]])
              if (sum(MLD[[fold]][set_from_ijk[[fold]][Ind],])>0)
              {
                points((rV),MLD[[fold]][s,]/MLDanalytic[[fold]][s,1]/10^s,cex=0.5,col = colors[s %% length(colors) + 1], pch = markers[s %% length(markers) + 1])
                lines((rV),MLDanalytic[[fold]][s,]/MLDanalytic[[fold]][s,1]/10^s,col = colors[s %% length(colors) + 1])
              }
            }
          }
        }
      }
    }
    xticks <- 10^(3:6)
    axis(1, at = xticks, labels = parse(text = paste0("10^", 3:6)))
    yticks <- 10^seq(-500,500,10)
    axis(2, at = yticks, labels = parse(text = paste0("10^", seq(-500,500,10))))
    box()
    dev.off()
  }
  if (nSp-nG>=6)
  {
    # 6-fold
    fold <- 6
    pdf(paste0(dir,"data/processed/",dirname,"/MLD_",fold,".pdf"),height=6*3,width=5*3)
    mar <- 2
    par(mfrow=c(6,5),mar = c(mar, mar, mar, mar))
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
                Ind <- paste0(i,"_",j,"_",k,"_",l,"_",mm,"_",pp)
                s <- set_from_ijk[[fold]][Ind]
                ijk <- as.numeric(str_split(Ind,"_")[[1]])
                if (sum(MLD[[fold]][s,])>0)
                {
                  # alpha <- 1 + sum(diff(rB)*MLD[[fold]][s,])/sum(diff(rB)*MLD[[fold]][s,]*log(rV/rTh))
                  # alpha <- alphaHillraw[[fold]][s]
                  # alphaTheory <- 1 + sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],])/sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],]*log(rV/rTh))
                  # alphaHill <- rbind(alphaHill,data.frame(Ind=Ind,n=fold,alpha,alphaTheory=alphaTheory,pair=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l],"_vs_",species[mm],"_vs_",species[pp])))
                  plot(log10(rV),log10(MLD[[fold]][s,]),pch=2,cex=0.5,col="green",ylab="",xlim=c(3,log10(rV[which.max(which(MLD[[fold]][s,]>0))])), ann=FALSE)
                  # ,main=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l],"_vs_",species[mm],"_vs_",species[pp]," a ",round(alpha,2)," ",round(alphaTheory,2), " M ", round(M[[fold]][set_from_ijk[[fold]][Ind]],2), " ", round(Mtheory[[fold]][set_from_ijk[[fold]][Ind]],2), " Om ", round(Om[[fold]][set_from_ijk[[fold]][Ind]],2), " ", round(Omtheory[[fold]][set_from_ijk[[fold]][Ind]],2)), cex.main=0.3)
                  # lines(log10(rV),log10(MLDsim[[fold]][s,]),col="turquoise")
                  lines(log10(rV),log10(MLDanalytic[[fold]][s,]),col="green")
                  xlim <- par("usr")[1:2];ylim <- par("usr")[3:4]-1.5;text(x = xlim[2], y = ylim[2], labels = paste0(species[i],",",species[j],",\n",species[k],",",species[l],",\n",species[mm],",",species[pp]), pos = 2, xpd = TRUE)
                  
                  # lines(log10(rV),log10(MLDdirect[[fold]][s,]),col="grey",lty=2)
                  # lines(log10(rV),log10(MLDghost[[fold]][s,]),col="grey")
                }
              }
            }
          }
        }
      }
    }
    dev.off()
    pdf(paste0(dir,"data/processed/",dirname,"/MLD_",fold,"_one.pdf"),height=10,width=5)
    plot(1, type = "n", xlim = c(1e3,0.5e5), ylim = c(1e-42,1e-1), log = "xy", axes = FALSE)
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
                Ind <- paste0(i,"_",j,"_",k,"_",l,"_",mm,"_",pp)
                s <- set_from_ijk[[fold]][Ind]
                ijk <- as.numeric(str_split(Ind,"_")[[1]])
                if (sum(MLD[[fold]][set_from_ijk[[fold]][Ind],])>0)
                {
                  points((rV),MLD[[fold]][s,]/MLDanalytic[[fold]][s,1]/10^s,cex=0.5,col = colors[s %% length(colors) + 1], pch = markers[s %% length(markers) + 1])
                  lines((rV),MLDanalytic[[fold]][s,]/MLDanalytic[[fold]][s,1]/10^s,col = colors[s %% length(colors) + 1])
                }
              }
            }
          }
        }
      }
    }
    xticks <- 10^(3:6)
    axis(1, at = xticks, labels = parse(text = paste0("10^", 3:6)))
    yticks <- 10^seq(-500,500,5)
    axis(2, at = yticks, labels = parse(text = paste0("10^", seq(-500,500,5))))
    box()
    dev.off()
  }
  if (nSp-nG>=7)
  {  
    # 7-fold
    fold <- 7
    pdf(paste0(dir,"data/processed/",dirname,"/MLD_",fold,".pdf"),height=3*3,width=3*3)
    par(mfrow=c(3,3),mar = 2*c(1.0, 1.0, 1.0, 1.0))
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
                  Ind <- paste0(i,"_",j,"_",k,"_",l,"_",mm,"_",pp,"_",q)
                  s <- set_from_ijk[[fold]][Ind]
                  ijk <- as.numeric(str_split(Ind,"_")[[1]])
                  if (sum(MLD[[fold]][s,])>0)
                  {
                    # alpha <- 1 + sum(diff(rB)*MLD[[fold]][s,])/sum(diff(rB)*MLD[[fold]][s,]*log(rV/rTh))
                    # alpha <- alphaHillraw[[fold]][s]
                    # alphaTheory <- 1 + sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],])/sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],]*log(rV/rTh))
                    # alphaHill <- rbind(alphaHill,data.frame(Ind=Ind,n=fold,alpha,alphaTheory=alphaTheory,pair=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l],"_vs_",species[mm],"_vs_",species[pp])))
                    plot(log10(rV),log10(MLD[[fold]][s,]),pch=2,cex=0.5,col="purple",ylab="",xlim=c(3,log10(rV[which.max(which(MLD[[fold]][s,]>0))])), ann=FALSE)
                    # ,main=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l],"_vs_",species[mm],"_vs_",species[pp],"_vs_",species[q]]), cex.main=0.3)
                    lines(log10(rV),log10(MLDanalytic[[fold]][s,]),col="purple")
                    xlim <- par("usr")[1:2];ylim <- par("usr")[3:4]-2.5;text(x = xlim[2], y = ylim[2], labels = paste0(species[i],",",species[j],",\n",species[k],",",species[l],",\n",species[mm],",",species[pp],",\n",species[q]), pos = 2, xpd = TRUE)
                    
                    # lines(log10(rV),log10(MLDdirect[[fold]][s,]),col="grey",lty=2)
                    # lines(log10(rV),log10(MLDghost[[fold]][s,]),col="grey")
                  }
                }
              }
            }
          }
        }
      }
    }
    dev.off()
    pdf(paste0(dir,"data/processed/",dirname,"/MLD_",fold,"_one.pdf"),height=10,width=5)
    plot(1, type = "n", xlim = c(1e3,0.5e5), ylim = c(1e-22,1e-1), log = "xy", axes = FALSE)
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
                  Ind <- paste0(i,"_",j,"_",k,"_",l,"_",mm,"_",pp,"_",q)
                  s <- set_from_ijk[[fold]][Ind]
                  ijk <- as.numeric(str_split(Ind,"_")[[1]])
                  if (sum(MLD[[fold]][set_from_ijk[[fold]][Ind],])>0)
                  {
                    points((rV),MLD[[fold]][s,]/MLDanalytic[[fold]][s,1]/10^s,cex=0.5,col = colors[s %% length(colors) + 1], pch = markers[s %% length(markers) + 1])
                    lines((rV),MLDanalytic[[fold]][s,]/MLDanalytic[[fold]][s,1]/10^s,col = colors[s %% length(colors) + 1])
                  }
                }
              }
            }
          }
        }
      }
    }
    xticks <- 10^(3:6)
    axis(1, at = xticks, labels = parse(text = paste0("10^", 3:6)))
    yticks <- 10^seq(-500,500,2)
    axis(2, at = yticks, labels = parse(text = paste0("10^", seq(-500,500,2))))
    box()
    dev.off()
    pdf(paste0(dir,"data/processed/",dirname,"/seq_freq_dist_",fold,".pdf"),height=3*3,width=3*3)
    par(mfrow=c(3,3),mar = 2*c(1.0, 1.0, 1.0, 1.0))
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
                  Ind <- paste0(i,"_",j,"_",k,"_",l,"_",mm,"_",pp,"_",q)
                  s <- set_from_ijk[[fold]][Ind]
                  ijk <- as.numeric(str_split(Ind,"_")[[1]])
                  if (sum(MLD[[fold]][set_from_ijk[[fold]][Ind],])>0)
                  {
                    Comparison <- paste0(species[q],"_vs_",species[pp],"_vs_",species[mm],"_vs_",species[l],"_vs_",species[k],"_vs_",species[j],"_vs_",species[i])
                    dat <- table(read.table(paste0("/home/projects/feariel/mishash/Development/ThreeFold/data/processed/Comparisons/",Comparison,"/mmseq2/cluster_cluster.tsv"), header = FALSE,fill=NA,sep="\t", comment.char="#",na.strings=".", stringsAsFactors=FALSE, quote="")$V1); 
                    dat <- dat/nList[[i]]/nList[[j]]
                    fB <- 10^seq(log10(min(dat))-0.2,log10(max(dat))+0.2,0.2)
                    fV <- (fB[-1]*fB[-length(fB)])^0.5
                    p <- hist(dat,breaks=fB,plot=FALSE)
                    x <- log10(fV); y <- (p$counts/diff(fB)); 
                    Ind <- which(y>0); y <- log10(y[Ind]); x <- x[Ind]
                    plot(x,y)
                    lines(x,mean(y)-1.5*(x-mean(x)))
                    xlim <- par("usr")[1:2];ylim <- par("usr")[3:4]-2.5;text(x = xlim[2], y = ylim[2], labels = paste0(species[i],",",species[j],",\n",species[k],",",species[l],",\n",species[mm],",",species[pp],",\n",species[q]), pos = 2, xpd = TRUE)
                  }
                }
              }
            }
          }
        }
      }
    }
    dev.off()
  }
  if (nSp-nG>=8)
  {  
    # 8-fold
    fold <- 8
    pdf(paste0(dir,"data/processed/",dirname,"/MLD_",fold,".pdf"),height=5*3,width=9*3)
    mar <- 2
    par(mfrow=c(5,9),mar = c(mar, mar, mar, mar))
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
                  for (rr in (q+1):(nSp-nG))
                  {
                    Ind <- paste0(i,"_",j,"_",k,"_",l,"_",mm,"_",pp,"_",q,"_",rr)
                    s <- set_from_ijk[[fold]][Ind]
                    ijk <- as.numeric(str_split(Ind,"_")[[1]])
                    if (sum(MLD[[fold]][s,])>0)
                    {
                      # alpha <- 1 + sum(diff(rB)*MLD[[fold]][s,])/sum(diff(rB)*MLD[[fold]][s,]*log(rV/rTh))
                      # alpha <- alphaHillraw[[fold]][s]
                      # alphaTheory <- 1 + sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],])/sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],]*log(rV/rTh))
                      # alphaHill <- rbind(alphaHill,data.frame(Ind=Ind,n=fold,alpha,alphaTheory=alphaTheory,pair=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l],"_vs_",species[mm],"_vs_",species[pp])))
                      # 
                      plot(log10(rV),log10(MLD[[fold]][s,]),pch=2,cex=0.5,col="red",ylab="",xlim=c(3,log10(rV[which.max(which(MLD[[fold]][s,]>0))])), ann=FALSE)
                      # ,main=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l],"_vs_",species[mm],"_vs_",species[pp],"_vs_",species[q]],"_vs_",species[rr]]), cex.main=0.3)
                      lines(log10(rV),log10(MLDanalytic[[fold]][s,]),col="red")
                      # lines(log10(rV),log10(MLDdirect[[fold]][s,]),col="grey",lty=2)
                      # lines(log10(rV),log10(MLDghost[[fold]][s,]),col="grey")
                      xlim <- par("usr")[1:2];ylim <- par("usr")[3:4]-2.5;text(x = xlim[2], y = ylim[2], labels = paste0(species[i],",",species[j],",\n",species[k],",",species[l],",\n",species[mm],",",species[pp],",\n",species[q],",",species[rr]), pos = 2, xpd = TRUE)
                      
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    dev.off()
    pdf(paste0(dir,"data/processed/",dirname,"/MLD_",fold,"_one.pdf"),height=10,width=5)
    plot(1, type = "n", xlim = c(1e3,0.5e5), ylim = c(1e-17,1e-1), log = "xy", axes = FALSE)
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
                  for (rr in (q+1):(nSp-nG))
                  {
                    Ind <- paste0(i,"_",j,"_",k,"_",l,"_",mm,"_",pp,"_",q,"_",rr)
                    s <- set_from_ijk[[fold]][Ind]
                    ijk <- as.numeric(str_split(Ind,"_")[[1]])
                    if (sum(MLD[[fold]][set_from_ijk[[fold]][Ind],])>0)
                    {
                      points((rV),MLD[[fold]][s,]/MLDanalytic[[fold]][s,1]/10^s,cex=0.5,col = colors[s %% length(colors) + 1], pch = markers[s %% length(markers) + 1])
                      lines((rV),MLDanalytic[[fold]][s,]/MLDanalytic[[fold]][s,1]/10^s,col = colors[s %% length(colors) + 1])
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    xticks <- 10^(3:6)
    axis(1, at = xticks, labels = parse(text = paste0("10^", 3:6)))
    yticks <- 10^seq(-500,500,1)
    axis(2, at = yticks, labels = parse(text = paste0("10^", seq(-500,500,1))))
    box()
    dev.off()
    pdf(paste0(dir,"data/processed/",dirname,"/seq_freq_dist_",fold,".pdf"),height=5*3,width=9*3)
    par(mfrow=c(5,9),mar = c(mar, mar, mar, mar))
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
                  for (rr in (q+1):(nSp-nG))
                  {
                    Ind <- paste0(i,"_",j,"_",k,"_",l,"_",mm,"_",pp,"_",q,"_",rr)
                    s <- set_from_ijk[[fold]][Ind]
                    ijk <- as.numeric(str_split(Ind,"_")[[1]])
                    if (sum(MLD[[fold]][set_from_ijk[[fold]][Ind],])>0)
                    {
                      Comparison <- paste0(species[rr],"_vs_",species[q],"_vs_",species[pp],"_vs_",species[mm],"_vs_",species[l],"_vs_",species[k],"_vs_",species[j],"_vs_",species[i])
                      dat <- table(read.table(paste0("/home/projects/feariel/mishash/Development/ThreeFold/data/processed/Comparisons/",Comparison,"/mmseq2/cluster_cluster.tsv"), header = FALSE,fill=NA,sep="\t", comment.char="#",na.strings=".", stringsAsFactors=FALSE, quote="")$V1); 
                      dat <- dat/nList[[i]]/nList[[j]]
                      fB <- 10^seq(log10(min(dat))-0.2,log10(max(dat))+0.2,0.2)
                      fV <- (fB[-1]*fB[-length(fB)])^0.5
                      p <- hist(dat,breaks=fB,plot=FALSE)
                      x <- log10(fV); y <- (p$counts/diff(fB)); 
                      Ind <- which(y>0); y <- log10(y[Ind]); x <- x[Ind]
                      plot(x,y)
                      lines(x,mean(y)-1.5*(x-mean(x)))
                      xlim <- par("usr")[1:2];ylim <- par("usr")[3:4]-2.5;text(x = xlim[2], y = ylim[2], labels = paste0(species[i],",",species[j],",\n",species[k],",",species[l],",\n",species[mm],",",species[pp],",\n",species[q],",",species[rr]), pos = 2, xpd = TRUE)
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    dev.off()
  }
  if (nSp-nG>=9)
  {
    # 9-fold
    fold <- 9
    pdf(paste0(dir,"data/processed/",dirname,"/MLD_",fold,".pdf"),height=2*3,width=5*3)
    mar <- 2
    par(mfrow=c(2,5),mar = c(mar, mar, mar, mar))
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
                    for (ss in (rr+1):(nSp-nG))
                    {
                      Ind <- paste0(i,"_",j,"_",k,"_",l,"_",mm,"_",pp,"_",q,"_",rr,"_",ss)
                      s <- set_from_ijk[[fold]][Ind]
                      ijk <- as.numeric(str_split(Ind,"_")[[1]])
                      # if (sum(MLD[[fold]][s,])>0)
                      {                    
                        # alpha <- 1 + sum(diff(rB)*MLD[[fold]][s,])/sum(diff(rB)*MLD[[fold]][s,]*log(rV/rTh))
                        # alpha <- alphaHillraw[[fold]][s]
                        # alphaTheory <- 1 + sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],])/sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],]*log(rV/rTh))
                        # alphaHill <- rbind(alphaHill,data.frame(Ind=Ind,n=fold,alpha,alphaTheory=alphaTheory,pair=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l],"_vs_",species[mm],"_vs_",species[pp])))
                        # 
                        plot(log10(rV),log10(MLD[[fold]][s,]),pch=2,cex=0.5,col="blue",ylab="",xlim=c(3,log10(rV[which.max(which(MLD[[fold]][s,]>0))])), ann=FALSE)
                        # ,main=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l],"_vs_",species[mm],"_vs_",species[pp],"_vs_",species[q]],"_vs_",species[rr]],"_vs_",species[ss]), cex.main=0.3)
                        lines(log10(rV),log10(MLDanalytic[[fold]][s,]),col="black")
                        lines(log10(rV),log10(MLDdirect[[fold]][s,]),col="grey",lty=2)
                        lines(log10(rV),log10(MLDghost[[fold]][s,]),col="grey")
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
    dev.off()
  }
  if (nSp-nG>=10)
  {
    # 10-fold
    fold <- 10
    pdf(paste0(dir,"data/processed/",dirname,"/MLD_",fold,".pdf"),height=1*4,width=1*4)
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
                      for (t in (ss+1):(nSp-nG))
                      {
                        Ind <- paste0(i,"_",j,"_",k,"_",l,"_",mm,"_",pp,"_",q,"_",rr,"_",ss,"_",t)
                        s <- set_from_ijk[[fold]][Ind]
                        ijk <- as.numeric(str_split(Ind,"_")[[1]])
                        # if (sum(MLD[[fold]][s,])>0)
                        {
                          # alpha <- 1 + sum(diff(rB)*MLD[[fold]][s,])/sum(diff(rB)*MLD[[fold]][s,]*log(rV/rTh))
                          # alpha <- alphaHillraw[[fold]][s]
                          # alphaTheory <- 1 + sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],])/sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],]*log(rV/rTh))
                          # alphaHill <- rbind(alphaHill,data.frame(Ind=Ind,n=fold,alpha,alphaTheory=alphaTheory,pair=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l],"_vs_",species[mm],"_vs_",species[pp])))
                          # 
                          plot(log10(rV),log10(MLD[[fold]][s,]),pch=2,cex=0.5,col="blue",ylab="",xlim=c(3,log10(rV[which.max(which(MLD[[fold]][s,]>0))])), ann=FALSE)
                          # ,main=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l],"_vs_",species[mm],"_vs_",species[pp],"_vs_",species[q]],"_vs_",species[rr]],"_vs_",species[ss]), cex.main=0.3)
                          lines(log10(rV),log10(MLDanalytic[[fold]][s,]),col="black")
                          lines(log10(rV),log10(MLDdirect[[fold]][s,]),col="grey",lty=2)
                          lines(log10(rV),log10(MLDghost[[fold]][s,]),col="grey")
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
    dev.off()
  }
  
  # separate plots
  fold <- 2
  pdf(paste0(dir,"data/processed/",dirname,"/MLDs.pdf"),height=5,width=5)
  for (i in 1:(nSp-nG-1))
  {
    for (j in (i+1):(nSp-nG))
    {
      Ind <- paste0(i,"_",j)
      s <- set_from_ijk[[fold]][Ind]
      ijk <- as.numeric(str_split(Ind,"_")[[1]])
      if (sum(MLD[[fold]][set_from_ijk[[fold]][Ind],])>0)
      {
        alpha <- 1 + sum(diff(rB)*MLD[[fold]][s,])/sum(diff(rB)*MLD[[fold]][s,]*log(rV/rTh))
        Ind <- which(MLD[[fold]][s,]>0 & rV<3e3); alpha <- lm(y~x,data=data.frame(x=log(rV[Ind]),y=log(MLD[[fold]][s,Ind])))$coefficients[2]
        # alphaTheory <- 1 + sum(diff(rB)*MLDsim[[fold]][s,])/sum(diff(rB)*MLDsim[[fold]][s,]*log(rV/rTh))
        # Ind <- which(MLDsim[[fold]][s,]>0 & rV<3e3); alphaTheory <- lm(y~x,data=data.frame(x=log(rV[Ind]),y=log(MLDsim[[fold]][s,Ind])))$coefficients[2]
        
        # alphaHill <- rbind(alphaHill,data.frame(Ind=Ind,n=fold,alpha,alphaTheory=alphaTheory,pair=paste0(species[j],"_vs_",species[i])))
        plot(log10(rV),log10(MLD[[fold]][s,]),pch=2,cex=0.5,col="blue",ylab="",xlim=c(3,log10(rV[which.max(which(MLD[[fold]][s,]>0))])),main=paste0(species[i],"_vs_",species[j]), cex.main=0.3)
        lines(log10(rV),log10(MLDanalytic[[fold]][s,]),col="black")
        lines(log10(rV),log10(MLDdirect[[fold]][s,]),col="grey",lty=2)
        lines(log10(rV),log10(MLDghost[[fold]][s,]),col="grey")
      }
    }
  }
  if (nSp-nG>=3)
  {
    # 3-fold
    m3 <- rV*0
    fold <- 3
    for (i in 1:(nSp-nG-2))
    {
      for (j in (i+1):(nSp-nG-1))
      {
        for (k in (j+1):(nSp-nG))
        {
          Ind <- paste0(i,"_",j,"_",k)
          s <- set_from_ijk[[fold]][Ind]
          ijk <- as.numeric(str_split(Ind,"_")[[1]])
          if (sum(MLD[[fold]][s,])>0)
          {
            alpha <- 1 + sum(diff(rB)*MLD[[fold]][set_from_ijk[[fold]][Ind],])/sum(diff(rB)*MLD[[fold]][set_from_ijk[[fold]][Ind],]*log(rV/rTh))
            # alphaTheory <- 1 + sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],])/sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],]*log(rV/rTh))
            # alphaHill <- rbind(alphaHill,data.frame(Ind=Ind,n=fold,alpha,alphaTheory=alphaTheory,pair=paste0(species[i],"_vs_",species[j],"_vs_",species[k])))
            plot(log10(rV),log10(MLD[[fold]][s,]),pch=2,cex=0.5,col="red",ylab="",xlim=c(3,log10(rV[which.max(which(MLD[[fold]][s,]>0))])),main=paste0(species[i],"_vs_",species[j],"_vs_",species[k]), cex.main=0.3)
            # lines(log10(rV),log10(MLDsim[[fold]][s,]),col="red")
          }
        }
      }
    }
  }
  if (nSp-nG>=4)
  {
    # 4-fold
    fold <- 4
    for (i in 1:(nSp-nG-3))
    {
      for (j in (i+1):(nSp-nG-2))
      {
        for (k in (j+1):(nSp-nG-1))
        {
          for (l in (k+1):(nSp-nG))
          {
            Ind <- paste0(i,"_",j,"_",k,"_",l)
            s <- set_from_ijk[[fold]][Ind]
            ijk <- as.numeric(str_split(Ind,"_")[[1]])
            if (sum(MLD[[fold]][s,])>0)
            {
              alpha <- 1 + sum(diff(rB)*MLD[[fold]][set_from_ijk[[fold]][Ind],])/sum(diff(rB)*MLD[[fold]][set_from_ijk[[fold]][Ind],]*log(rV/rTh))
              # alphaTheory <- 1 + sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],])/sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],]*log(rV/rTh))
              # alphaHill <- rbind(alphaHill,data.frame(Ind=Ind,n=fold,alpha,alphaTheory=alphaTheory,pair=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l])))
              plot(log10(rV),log10(MLD[[fold]][s,]),pch=2,cex=0.5,col="purple",ylab="",xlim=c(3,log10(rV[which.max(which(MLD[[fold]][s,]>0))])),main=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l]), cex.main=0.3)
              # lines(log10(rV),log10(MLDsim[[fold]][s,]),col="purple")
            }
          }
        }
      }
    }
  }
  if (nSp-nG>=5)
  {
    # 5-fold
    fold <- 5
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
              Ind <- paste0(i,"_",j,"_",k,"_",l,"_",mm)
              s <- set_from_ijk[[fold]][Ind]
              ijk <- as.numeric(str_split(Ind,"_")[[1]])
              if (sum(MLD[[fold]][s,])>0)
              {
                alpha <- 1 + sum(diff(rB)*MLD[[fold]][set_from_ijk[[fold]][Ind],])/sum(diff(rB)*MLD[[fold]][set_from_ijk[[fold]][Ind],]*log(rV/rTh))
                # alphaTheory <- 1 + sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],])/sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],]*log(rV/rTh))
                # alphaHill <- rbind(alphaHill,data.frame(Ind=Ind,n=fold,alpha,alphaTheory=alphaTheory,pair=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l],"_vs_",species[mm])))
                plot(log10(rV),log10(MLD[[fold]][s,]),pch=2,cex=0.5,col="orange",ylab="",xlim=c(3,log10(rV[which.max(which(MLD[[fold]][s,]>0))])),main=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l],"_vs_",species[mm]), cex.main=0.3)
                # lines(log10(rV),log10(MLDsim[[fold]][s,]),col="orange")
              }
            }
          }
        }
      }
    }
  }
  if (nSp-nG>=6)
  {
    # 6-fold
    fold <- 6
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
                Ind <- paste0(i,"_",j,"_",k,"_",l,"_",mm,"_",pp)
                s <- set_from_ijk[[fold]][Ind]
                ijk <- as.numeric(str_split(Ind,"_")[[1]])
                if (sum(MLD[[fold]][s,])>0)
                {
                  alpha <- 1 + sum(diff(rB)*MLD[[fold]][set_from_ijk[[fold]][Ind],])/sum(diff(rB)*MLD[[fold]][set_from_ijk[[fold]][Ind],]*log(rV/rTh))
                  # alphaTheory <- 1 + sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],])/sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],]*log(rV/rTh))
                  # alphaHill <- rbind(alphaHill,data.frame(Ind=Ind,n=fold,alpha,alphaTheory=alphaTheory,pair=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l],"_vs_",species[mm],"_vs_",species[pp])))
                  plot(log10(rV),log10(MLD[[fold]][s,]),pch=2,cex=0.5,col="black",ylab="",xlim=c(3,log10(rV[which.max(which(MLD[[fold]][s,]>0))])),main=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l],"_vs_",species[mm],"_vs_",species[pp]), cex.main=0.3)
                  # lines(log10(rV),log10(MLDsim[[fold]][s,]),col="black")
                }
              }
            }
          }
        }
      }
    }
  }
  if (nSp-nG>=7)
  {  
    # 7-fold
    fold <- 7
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
                  Ind <- paste0(i,"_",j,"_",k,"_",l,"_",mm,"_",pp,"_",q)
                  s <- set_from_ijk[[fold]][Ind]
                  ijk <- as.numeric(str_split(Ind,"_")[[1]])
                  if (sum(MLD[[fold]][s,])>0)
                  {
                    alpha <- 1 + sum(diff(rB)*MLD[[fold]][set_from_ijk[[fold]][Ind],])/sum(diff(rB)*MLD[[fold]][set_from_ijk[[fold]][Ind],]*log(rV/rTh))
                    # alphaTheory <- 1 + sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],])/sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],]*log(rV/rTh))
                    # alphaHill <- rbind(alphaHill,data.frame(Ind=Ind,n=fold,alpha,alphaTheory=alphaTheory,pair=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l],"_vs_",species[mm],"_vs_",species[pp])))
                    plot(log10(rV),log10(MLD[[fold]][s,]),pch=2,cex=0.5,col="green",ylab="",xlim=c(3,log10(rV[which.max(which(MLD[[fold]][s,]>0))])),main=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l],"_vs_",species[mm],"_vs_",species[pp],"_vs_",species[q]), cex.main=0.3)
                    # lines(log10(rV),log10(MLDsim[[fold]][s,]),col="green")
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (nSp-nG>=8)
  {  
    # 8-fold
    fold <- 8
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
                  for (rr in (q+1):(nSp-nG))
                  {
                    Ind <- paste0(i,"_",j,"_",k,"_",l,"_",mm,"_",pp,"_",q,"_",rr)
                    s <- set_from_ijk[[fold]][Ind]
                    ijk <- as.numeric(str_split(Ind,"_")[[1]])
                    if (sum(MLD[[fold]][s,])>0)
                    {
                      alpha <- 1 + sum(diff(rB)*MLD[[fold]][set_from_ijk[[fold]][Ind],])/sum(diff(rB)*MLD[[fold]][set_from_ijk[[fold]][Ind],]*log(rV/rTh))
                      # alphaTheory <- 1 + sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],])/sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],]*log(rV/rTh))
                      # alphaHill <- rbind(alphaHill,data.frame(Ind=Ind,n=fold,alpha,alphaTheory=alphaTheory,pair=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l],"_vs_",species[mm],"_vs_",species[pp])))
                      plot(log10(rV),log10(MLD[[fold]][s,]),pch=2,cex=0.5,col="red",ylab="",xlim=c(3,log10(rV[which.max(which(MLD[[fold]][s,]>0))])),main=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l],"_vs_",species[mm],"_vs_",species[pp],"_vs_",species[q],"_vs_",species[rr]), cex.main=0.3)
                      # lines(log10(rV),log10(MLDsim[[fold]][s,]),col="red")
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
  if (nSp-nG>=9)
  {
    # 9-fold
    fold <- 9
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
                    for (ss in (rr+1):(nSp-nG))
                    {
                      Ind <- paste0(i,"_",j,"_",k,"_",l,"_",mm,"_",pp,"_",q,"_",rr,"_",ss)
                      s <- set_from_ijk[[fold]][Ind]
                      ijk <- as.numeric(str_split(Ind,"_")[[1]])
                      if (sum(MLD[[fold]][s,])>0)
                      {
                        alpha <- 1 + sum(diff(rB)*MLD[[fold]][set_from_ijk[[fold]][Ind],])/sum(diff(rB)*MLD[[fold]][set_from_ijk[[fold]][Ind],]*log(rV/rTh))
                        # alphaTheory <- 1 + sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],])/sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],]*log(rV/rTh))
                        # alphaHill <- rbind(alphaHill,data.frame(Ind=Ind,n=fold,alpha,alphaTheory=alphaTheory,pair=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l],"_vs_",species[mm],"_vs_",species[pp])))
                        plot(log10(rV),log10(MLD[[fold]][s,]),pch=2,cex=0.5,col="blue",ylab="",xlim=c(3,log10(rV[which.max(which(MLD[[fold]][s,]>0))])),main=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l],"_vs_",species[mm],"_vs_",species[pp],"_vs_",species[q],"_vs_",species[rr],"_vs_",species[ss]), cex.main=0.3)
                        # lines(log10(rV),log10(MLDsim[[fold]][s,]),col="blue")
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
  if (nSp-nG>=10)
  {
    # 10-fold
    fold <- 10
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
                      for (t in (ss+1):(nSp-nG))
                      {
                        Ind <- paste0(i,"_",j,"_",k,"_",l,"_",mm,"_",pp,"_",q,"_",rr,"_",ss,"_",t)
                        s <- set_from_ijk[[fold]][Ind]
                        ijk <- as.numeric(str_split(Ind,"_")[[1]])
                        if (sum(MLD[[fold]][s,])>0)
                        {
                          alpha <- 1 + sum(diff(rB)*MLD[[fold]][set_from_ijk[[fold]][Ind],])/sum(diff(rB)*MLD[[fold]][set_from_ijk[[fold]][Ind],]*log(rV/rTh))
                          alpha <- 1 + sum(diff(rB)*MLD[[fold]][set_from_ijk[[fold]][Ind],])/sum(diff(rB)*MLD[[fold]][set_from_ijk[[fold]][Ind],]*log(rV/rTh))
                          # alphaTheory <- 1 + sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],])/sum(diff(rB)*MLDsim[[fold]][set_from_ijk[[fold]][Ind],]*log(rV/rTh))
                          # alphaHill <- rbind(alphaHill,data.frame(Ind=Ind,n=fold,alpha,alphaTheory=alphaTheory,pair=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l],"_vs_",species[mm],"_vs_",species[pp])))
                          plot(log10(rV),log10(MLD[[fold]][s,]),pch=2,cex=0.5,col="blue",ylab="",xlim=c(3,log10(rV[which.max(which(MLD[[fold]][s,]>0))])),main=paste0(species[i],"_vs_",species[j],"_vs_",species[k],"_vs_",species[l],"_vs_",species[mm],"_vs_",species[pp],"_vs_",species[q],"_vs_",species[rr],"_vs_",species[ss],"_vs_",species[t]), cex.main=0.3)
                          # lines(log10(rV),log10(MLDsim[[fold]][s,]),col="blue")
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
  dev.off()
  
  write.table(alphaHill,paste0(dir,"data/processed/",dirname,"/alpha.tsv"),row.names=FALSE);
  pdf(paste0(dir,"data/processed/",dirname,"/alphaHist.pdf"),width=5,height=5)
  ggplot(alphaHill, aes(x = as.factor(n), y = alpha,fill=as.factor(n))) + 
    geom_violin() +
    geom_jitter(shape=16, position=position_jitter(0.2),cex=0.5) +
    scale_y_continuous(breaks=seq(0,12,1), expand = c(0, 0),limits=c(2,12)) +
    scale_x_discrete() + 
    theme(panel.grid.major = element_line(color = "red",linewidth = 0.5,linetype = 2))
  print(p)
  # ggplot(alphaHill, aes(x = as.factor(n), y = alphaTheory,fill=as.factor(n))) + 
  #   geom_violin() +
  #   geom_jitter(shape=16, position=position_jitter(0.2),cex=0.5) +
  #   scale_y_continuous(breaks=seq(0,12,1), expand = c(0, 0),limits=c(2,12)) +
  #   scale_x_discrete() + 
  #   theme(panel.grid.major = element_line(color = "red",linewidth = 0.5,linetype = 2))
  # print(p)
  # p <- ggplot(alphaHill, aes(x=alpha,y=alphaTheory,col=as.factor(n))) +
  #   geom_point() +
  #   geom_line(aes(x = alpha, y = alpha)) + 
  #   scale_x_continuous(breaks=seq(0,12,1), expand = c(0, 0),limits=c(2,12)) +
  #   scale_y_continuous(breaks=seq(0,12,1), expand = c(0, 0),limits=c(2,12)) +
  #   theme(panel.grid.major = element_line(color = "red",linewidth = 0.5,linetype = 2))
  # print(p)
  p <- ggplot(alphaHill, aes(x=as.factor(n), y=alpha,col=as.factor(n))) + 
    geom_violin() + geom_jitter(shape=16,cex=0.3, position=position_jitter(0.2)) +
    scale_y_continuous(breaks=seq(0,11,1), expand = c(0, 0),limits=c(2,11.5)) +
    theme(panel.grid.major = element_line(color = "grey",linewidth = 0.5,linetype = 2)) +
    geom_abline(intercept = 2, slope = 1, color="grey",  linetype="dashed", linewidth=0.5) +
    geom_abline(intercept = 3, slope = 1, color="grey",  linetype="dashed", linewidth=0.5)
  print(p)
  # p <- ggplot(alphaHill, aes(x=as.factor(n), y=alphaTheory,col=as.factor(n))) + 
  #   geom_violin() + geom_jitter(shape=16,cex=0.3, position=position_jitter(0.2)) +
  #   scale_y_continuous(breaks=seq(0,11,1), expand = c(0, 0),limits=c(2,11.5)) +
  #   theme(panel.grid.major = element_line(color = "grey",linewidth = 0.5,linetype = 2)) +
  #   geom_abline(intercept = 2, slope = 1, color="grey",  linetype="dashed", linewidth=0.5) +
  #   geom_abline(intercept = 3, slope = 1, color="grey",  linetype="dashed", linewidth=0.5)
  # print(p)
  dev.off()
}
# source("/home/projects/feariel/mishash/Development/ThreeFold/src/plotResults.R"); plotResults()
# 
# Analytical functions----
getMLD <- function(rhoM,N,IndR)
{
  # calculate R matrix
  R <- list()
  for (n in 2:N)
  {
    R[[n]] <- matrix(0,nrow=length(set_from_ijk[[n]]),ncol=length(set_from_ijk[[n]]))
    for (s in 1:length(set_from_ijk[[n]]))
    {
      ijk <- ijk_from_set[[n]][s,]
      for (i in ijk)
      {
        for (j in ijk)
        {
          R[[n]][s,s] <- R[[n]][s,s] + rhoM[i,j]
        }
      }
      not_ijk <- setdiff(1:nSp,ijk)
      for (i in not_ijk)
      {
        for (j in ijk)
        {
          R[[n]][s,s] <- R[[n]][s,s] + rhoM[i,j]
          i_insteadof_j <- ijk; i_insteadof_j[which(i_insteadof_j==j)] <- i; i_insteadof_j <- sort(i_insteadof_j); s_i_insteadof_j <- set_from_ijk[[n]][paste(i_insteadof_j,collapse='_')]
          R[[n]][s,s_i_insteadof_j] <- R[[n]][s,s_i_insteadof_j] - rhoM[i,j]
        }
      }
    } 
  }
  # calculate tau probs in Laplace rV space
  P <- list()
  Pp <- list()
  Pm <- list()
  MLD <- list()
  P[[1]] <- matrix(1,length(set_from_ijk[[1]]),length(rV))
  Pp[[1]] <- matrix(1,length(set_from_ijk[[1]]),length(rV))
  Pm[[1]] <- matrix(1,length(set_from_ijk[[1]]),length(rV))
  for (n in 2:N)
  {
    P[[n]] <- foreach (ir=1:length(rV), .inorder=TRUE, .combine=cbind) %do%
    {
      brackets <- rep(0,length(set_from_ijk[[n]]))
      for (s in 1:length(set_from_ijk[[n]]))
      {
        ijk <- ijk_from_set[[n]][s,]
        for (j in ijk)
        {
          s_not_j <- ijk[ijk != j]; s_not_j <- set_from_ijk[[n-1]][paste(s_not_j,collapse='_')]
          for (i in ijk)
          {
            brackets[s] <- brackets[s] + rhoM[i,j]*(P[[n-1]][s_not_j,ir]-1)
          }
        }
      }
      brackets <- brackets + rowSums(R[[n]])
      if (length(brackets)==1)
      {
        return(1/(n*rV[ir]*diag(nrow(R[[n]]))+R[[n]]) * brackets)
      }
      else
      {
        myMat <- n*rV[ir]*diag(nrow(R[[n]]))+R[[n]]
        # myMatsym <- (myMat + t(myMat))/2
        # diag(myMatsym) <- diag(myMat)
        # if (is.positive.definite(myMat)) {return(armaInv_sympd(myMat) %*% matrix(brackets,ncol=1))}
        # else {return(MASS::ginv((myMat)) %*% matrix(brackets,ncol=1))}
        return(MASS::ginv((myMat)) %*% matrix(brackets,ncol=1))
        # return(solve(n*rV[IndR[ir]]*diag(nrow(R[[n]]))+R[[n]],matrix(brackets,ncol=1)))
        # return(chol2inv(chol(n*rV[IndR[ir]]*diag(nrow(R[[n]]))+R[[n]])) %*% matrix(brackets,ncol=1))
      }
    }
    Pp[[n]] <- foreach (ir=1:length(rV), .inorder=TRUE, .combine=cbind) %do%
    {
      brackets <- rep(0,length(set_from_ijk[[n]]))
      for (s in 1:length(set_from_ijk[[n]]))
      {
        ijk <- ijk_from_set[[n]][s,]
        for (j in ijk)
        {
          s_not_j <- ijk[ijk != j]; s_not_j <- set_from_ijk[[n-1]][paste(s_not_j,collapse='_')]
          for (i in ijk)
          {
            brackets[s] <- brackets[s] + rhoM[i,j]*(Pp[[n-1]][s_not_j,ir]-1)
          }
        }
      }
      brackets <- brackets + rowSums(R[[n]])
      if (length(brackets)==1)
      {
        return(1/(n*(rV[ir]+dr[ir])*diag(nrow(R[[n]]))+R[[n]]) * brackets)
      }
      else
      {
        myMat <- n*(rV[ir]+dr[ir])*diag(nrow(R[[n]]))+R[[n]]
        # myMatsym <- (myMat + t(myMat))/2
        # diag(myMatsym) <- diag(myMat)
        # if (is.positive.definite(myMat)) {return(armaInv_sympd(myMat) %*% matrix(brackets,ncol=1))}
        # else {return(MASS::ginv((myMat)) %*% matrix(brackets,ncol=1))}
        return(MASS::ginv((myMat)) %*% matrix(brackets,ncol=1))
        # return(solve(n*(rV[IndR[ir]]+dr[IndR[ir]])*diag(nrow(R[[n]]))+R[[n]],matrix(brackets,ncol=1)))  
        # return(chol2inv(chol(n*(rV[IndR[ir]]+dr[IndR[ir]])*diag(nrow(R[[n]]))+R[[n]])) %*% matrix(brackets,ncol=1))
      }
    }
    Pm[[n]] <- foreach (ir=1:length(rV), .inorder=TRUE, .combine=cbind) %do%
    {
      brackets <- rep(0,length(set_from_ijk[[n]]))
      for (s in 1:length(set_from_ijk[[n]]))
      {
        ijk <- ijk_from_set[[n]][s,]
        for (j in ijk)
        {
          s_not_j <- ijk[ijk != j]; s_not_j <- set_from_ijk[[n-1]][paste(s_not_j,collapse='_')]
          for (i in ijk)
          {
            brackets[s] <- brackets[s] + rhoM[i,j]*(Pm[[n-1]][s_not_j,ir]-1)
          }
        }
      }
      brackets <- brackets + rowSums(R[[n]])
      if (length(brackets)==1)
      { 
        return(1/(n*(rV[ir]-dr[ir])*diag(nrow(R[[n]]))+R[[n]]) * brackets)
      }
      else
      {
        myMat <- n*(rV[ir]-dr[ir])*diag(nrow(R[[n]]))+R[[n]]
        # myMatsym <- (myMat + t(myMat))/2
        # diag(myMatsym) <- diag(myMat)
        # if (is.positive.definite(myMat)) {return(armaInv_sympd(myMat) %*% matrix(brackets,ncol=1))}
        # else {return(MASS::ginv((myMat)) %*% matrix(brackets,ncol=1))}
        return(MASS::ginv((myMat)) %*% matrix(brackets,ncol=1))
        # return(solve(n*(rV[IndR[ir]]-dr[IndR[ir]])*diag(nrow(R[[n]]))+R[[n]],matrix(brackets,ncol=1)))
        # return(chol2inv(chol(n*(rV[IndR[ir]]-dr[IndR[ir]])*diag(nrow(R[[n]]))+R[[n]])) %*% matrix(brackets,ncol=1))
      }
    }
    MLD[[n]] <- (Pp[[n]] + Pm[[n]] - 2*P[[n]]) %b*% matrix(drm2,nrow=1)
  }
  return(MLD)
}

