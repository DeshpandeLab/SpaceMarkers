#' @import matrixTests
#' @import rstatix
#' @import matrixStats
#import description end
0



## author: Atul Deshpande
## email: adeshpande@jhu.edu

do_dunnTest = function(in.kruskal){
  out = array(NA,dim = 3, dimnames = c("pval","interaction-pattern1","interaction-pattern2"))
  temp.summary = summary(in.kruskal)
  out["pval"] = in.summary[[1]]$`Pr(>F)`[1]
}
perform.dunn.test = function(in.data){
  test <- dunn_test(in.data,residuals~region)
  ofInterest <- which(test$group1=="Interacting"|test$group2=="Interacting")
  #ofInterest <- 1:2
  return(unlist(test[,c("statistic","p.adj")]))
}
fast.dunn.test = function(in.data,region){
  region = factor(region)
  pattern1 = levels(region)[2]
  pattern2 = levels(region)[3]
  in.ranks <- rowRanks(in.data,cols = !is.na(region),ties.method = "max")
  rsub <- region[!is.na(region)]
  
  N <- length(rsub)
  N1 <- sum(rsub==pattern1)
  N2 <-sum(rsub==pattern2)
  NI <-sum(rsub=="Interacting")
  SEI1 <- sqrt(N*(N+1)/12*(1/NI+1/N1))
  SEI2 <- sqrt(N*(N+1)/12*(1/NI+1/N2))
  SE12 <- sqrt(N*(N+1)/12*(1/N1+1/N2))
  
  mInt<- apply(in.ranks[,rsub=="Interacting"],1,mean)
  mP1<- apply(in.ranks[,rsub==pattern1],1,mean)
  mP2<- apply(in.ranks[,rsub==pattern2],1,mean)

  zP1_Int <- (mP1-mInt)/SEI1
  zP2_Int <- (mP2-mInt)/SEI2
  zP2_P1 <- (mP2-mP1)/SE12
  zVals <-cbind(zP1_Int,zP2_Int,zP2_P1)
  
  pval_1_Int <- ifelse(2*(pnorm(zP1_Int,lower.tail = TRUE))<=1,2*(pnorm(zP1_Int,lower.tail = TRUE)),1)
  pval_2_Int <- ifelse(2*(pnorm(zP2_Int,lower.tail = TRUE))<=1,2*(pnorm(zP2_Int,lower.tail = TRUE)),1)
  pval_2_1 <- ifelse(zP2_P1>0,2*(pnorm(zP2_P1,lower.tail = FALSE)),2*(pnorm(zP2_P1,lower.tail = TRUE)))
  
  pvals<-cbind(pval_1_Int,pval_2_Int,pval_2_1)
  p.adj <- ifelse(pvals*3>1,1,pvals*3)
  return(cbind(zVals,p.adj))
}

#===================
#' find_genes_of_interest_nonparametric_fast
#' Calculate ...
#'
#' This function calculates ...
#'
#'
#' @param testMat 	...
#' @param goodGenes ...
#' @param region	...
#'
#'
#'
#'
#' @return a list of ...


find_genes_of_interest_nonparametric_fast <- function(testMat, goodGenes = NULL, region)
{
  region = factor(region)
  patnames <- levels(region)[which(levels(region)!="Interacting")]
  pattern1 = patnames[1]
  pattern2 = patnames[2]
  if (!is.null(goodGenes)){
    subset_goodGenes <- intersect(rownames(testMat),goodGenes)
    testMat = testMat[subset_goodGenes,]
  }
  residuals.kruskal = matrixTests::row_kruskalwallis(x = as.matrix(testMat), g = !is.na(region))
  qq <- qvalue::qvalue(residuals.kruskal$pvalue,fdr.level = 0.05, pfdr = F, pi0 = 1)
  residuals.kruskal = cbind(residuals.kruskal,p.adj = qq$qvalues)
  ind <- rownames(residuals.kruskal[which(residuals.kruskal$p.adj<5e-2),])
  residuals.dunn.test = fast.dunn.test(as.matrix(testMat[ind,]),region)
  rownames(residuals.dunn.test) <- ind
  interacting_over_pattern1 <- residuals.dunn.test[,"pval_1_Int"]<0.05
  interacting_over_pattern2 <- residuals.dunn.test[,"pval_2_Int"]<0.05
  interacting_over_both_patterns <- interacting_over_pattern1 & interacting_over_pattern2
  not_pattern1_diff_pattern2 <- residuals.dunn.test[,"pval_2_1"]>=0.05
  exclusive_interacting_over_pattern1 <- interacting_over_pattern1 & not_pattern1_diff_pattern2
  exclusive_interacting_over_pattern2 <- interacting_over_pattern2 & not_pattern1_diff_pattern2
  names(interacting_over_both_patterns) <- names(exclusive_interacting_over_pattern1) <- names(exclusive_interacting_over_pattern2) <- ind
  genes_interacting = matrix(FALSE, c(dim(testMat)[1],1))
  rownames(genes_interacting) = rownames(testMat)
  genes_interacting[ind[which(exclusive_interacting_over_pattern1)],1]=paste0("vs",pattern1)
  genes_interacting[ind[which(exclusive_interacting_over_pattern2)],1]=paste0("vs",pattern2)
  genes_interacting[ind[which(interacting_over_both_patterns)],1]="vsBoth"
  colnames(genes_interacting) = paste0(pattern1,' x ', pattern2)
  genes_interacting <- cbind(genes_interacting,p.adj = residuals.kruskal[,"p.adj"])
  genes_interacting <- genes_interacting[genes_interacting[,1]!="FALSE",]
  return(list(genes_interacting))
}
