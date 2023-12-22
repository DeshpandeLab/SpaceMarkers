#' @importFrom matrixTests row_kruskalwallis
#' @importFrom rstatix filter
#' @importFrom matrixStats count
#' @importFrom stats lag dist pnorm sd
#' @importFrom qvalue qvalue
#' @import utils
#import description end



## author: Atul Deshpande
## email: adeshpande@jhu.edu

row.dunn.test <- function(in.data,region){
    region <- factor(region)
    pattern1 <- levels(region)[2]
    pattern2 <- levels(region)[3]
    in.ranks <- matrixStats::rowRanks(in.data,cols = !is.na(region),
                                        ties.method = "average")
    rsub <- region[!is.na(region)]
    
    N <- length(rsub)
    N1 <- sum(rsub==pattern1)
    N2 <- sum(rsub==pattern2)
    NI <- sum(rsub=="Interacting")

    tiesStat <- apply(in.ranks,1,function(rr) sum(table(rr)^3-table(rr)))
    tiesStat2 <- sqrt(1 - tiesStat/N/(N-1)/(N+1))

    SEI1 <- sqrt(N*(N+1)/12*(1/NI+1/N1))
    SEI2 <- sqrt(N*(N+1)/12*(1/NI+1/N2))
    SE12 <- sqrt(N*(N+1)/12*(1/N1+1/N2))

    mInt<- apply(in.ranks[,rsub=="Interacting"],1,mean)
    mP1<- apply(in.ranks[,rsub==pattern1],1,mean)
    mP2<- apply(in.ranks[,rsub==pattern2],1,mean)

    zP1_Int <- (mP1-mInt)/SEI1/tiesStat2
    zP2_Int <- (mP2-mInt)/SEI2/tiesStat2
    zP2_P1 <- (mP2-mP1)/SE12/tiesStat2
    zVals <-cbind(zP1_Int,zP2_Int,zP2_P1)

    pval_1_Int <- ifelse(2*(pnorm(zP1_Int,lower.tail = TRUE))
                        <=1,2*(pnorm(zP1_Int,lower.tail = TRUE)),1)
    pval_2_Int <- ifelse(2*(pnorm(zP2_Int,lower.tail = TRUE))
                            <=1,2*(pnorm(zP2_Int,lower.tail = TRUE)),1)
    pval_2_1 <- ifelse(zP2_P1>0,2*(pnorm(zP2_P1,lower.tail = FALSE)),
                        2*(pnorm(zP2_P1,lower.tail = TRUE)))
    pvals<-cbind(pval_1_Int,pval_2_Int,pval_2_1)
return(cbind(zVals,pvals))
}

#===================
#' find_genes_of_interest
#' Identify genes associated with pattern interaction.

#' This function identifies genes exhibiting significantly higher values of 
#' testMat in the Interaction region of the two 
#' patterns compared #' to regions with exclusive influence from either pattern.
#' It uses Kruskal-Wallis test followed by
#' posthoc analysis using Dunn's Test to identify the genes.
#'
#' @usage
#' find_genes_of_interest(
#'     testMat,
#'     goodGenes = NULL,
#'     region,
#'     fdr.level = 0.05,
#'     analysis = c("enrichment", "overlap")
#' )
#' @param    testMat A matrix of counts with cells as columns and genes as rows
#' @param    goodGenes A vector of user specified genes expected to interact 
#' a priori. The default for this is NULL as the function can find these genes 
#' itself
#' @param    region A data frame of the reference pattern regions that overlap 
#' with the other patterns
#' @param    fdr.level False Discovery Rate. The default value is 0.05.
#' @param    analysis a character string that specifies the type of analysis to 
#' carry out, whether overlap or enrichment.
#'
#' @return a list of genes exhibiting significantly higher values of testMat in 
#' the Interaction region of the two #' patterns compared to regions with 
#' exclusive influence from either pattern.



find_genes_of_interest<-function(
        testMat,goodGenes=NULL,region, fdr.level=0.05,
        analysis=c("enrichment","overlap")) {
    region <- factor(region)
    patnames <- levels(region)[which(levels(region)!="Interacting")]
    pattern1 <- patnames[1]
    pattern2 <- patnames[2]
    if (!is.null(goodGenes)){
        subset_goodGenes <- intersect(rownames(testMat),goodGenes)
        testMat <- testMat[subset_goodGenes,]
        }
    res_kruskal<- matrixTests::row_kruskalwallis(x=as.matrix(testMat),g=region)
    qq <- qvalue::qvalue(res_kruskal$pvalue,fdr.level = fdr.level,
                            pfdr = FALSE, pi0 = 1)
    res_kruskal <- cbind(res_kruskal,p.adj = qq$qvalues)
    res_dunn_test <- row.dunn.test(as.matrix(testMat),region)
    rownames(res_dunn_test) <- rownames(res_kruskal)
    ind <- rownames(res_kruskal[which(res_kruskal$p.adj<fdr.level),])
    qDunn <- qvalue::qvalue(res_dunn_test[,4:6],
                            fdr.level = fdr.level, pfdr = FALSE, pi0 = 1)
    qq<-qvalue::qvalue(
        res_dunn_test[ind,4:6],fdr.level=fdr.level,pfdr=FALSE,pi0 = 1)
    qDunn$qvalues[ind,] <- qq$qvalue
    res_dunn_test <- cbind(res_dunn_test,qDunn$qvalues)
    colnames(res_dunn_test)[7:9] <- paste0(colnames(res_dunn_test)[7:9],".adj")
    interact_patt1 <- res_dunn_test[ind,"pval_1_Int"]<fdr.level
    interact_patt2 <- res_dunn_test[ind,"pval_2_Int"]<fdr.level
    interacting_over_both_patterns <- interact_patt1 & interact_patt2
    not_pattern1_diff_pattern2 <- res_dunn_test[ind,"pval_2_1"]>=fdr.level
    exc_interact_patt1 <- interact_patt1 & not_pattern1_diff_pattern2
    exc_interact_patt2 <- interact_patt2 & not_pattern1_diff_pattern2
    names(exc_interact_patt2)<-ind
    names(exc_interact_patt1)<-names(exc_interact_patt2)
    names(interacting_over_both_patterns)<-names(exc_interact_patt1)
    interact_genes<-matrix(FALSE,nrow=nrow(res_dunn_test),ncol=2,dimnames=list(
        rownames(res_dunn_test),c("Gene",paste0(pattern1,' x ',pattern2))))
    interact_genes[,1] <- rownames(interact_genes)
    interact_genes[ind[which(exc_interact_patt1)],2]<-paste0("vs",pattern1)
    interact_genes[ind[which(exc_interact_patt2)],2]<-paste0("vs",pattern2)
    interact_genes[ind[which(interacting_over_both_patterns)],2]<-"vsBoth"
    colnames(res_kruskal) <- paste0("KW.",colnames(res_kruskal))
    colnames(res_dunn_test) <- paste0("Dunn.",colnames(res_dunn_test))
    interact_genes <- cbind(interact_genes,res_kruskal,res_dunn_test)
    if (analysis=="overlap"){
        interact_genes <-interact_genes[interact_genes[,2]!="FALSE",] 
    }
    return(list(interact_genes))
}
