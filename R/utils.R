
#' @title getOverlapScores
#' @description Calculate the overlap scores between patterns in hotspots
#' @param hotspots A data frame with columns x, y, barcode and pattern names
#' @param patternList A character vector of pattern names to calculate overlap 
#' scores for
#' @param method The method to calculate overlap scores. Options are
#' "Szymkiewicz-Simpson", "Jaccard", "Sorensen-Dice", "Ochiai" and "absolute"
#' @details The function calculates the overlap scores between patterns hotspots
#' using the specified method. The default method is "Szymkiewicz-Simpson"
#' overlap coefficient.
#' @return A data frame with columns pattern1, pattern2 and overlapScore
#' @export
#' @examples
#' hotspots <- data.frame(x = c(1,2,3,4,5),
#'                         y = c(1,2,3,4,5),
#'                         barcode = c("A","B","C","D","E"),
#'                         pattern1 = c(1,0,1,0,1),
#'                         pattern2 = c(1,1,0,0,1))
#' getOverlapScores(hotspots)   
#' getOverlapScores(hotspots, c("pattern1","pattern2"))
#' @importFrom ggplot2 ggplot geom_tile geom_text theme_minimal
#' @importFrom reshape2 melt
#' @importFrom stats complete.cases
getOverlapScores <- function(hotspots,
                             patternList = NULL, method = c("Szymkiewicz-Simpson",
                                                            "Jaccard", "Sorensen-Dice",
                                                            "Ochiai", "absolute") ) {
    
    #stop if more than one method is supplied, do not warn by default
    if(length(method) > 1){
        method <- method[1]
        message("Only one method can be used at a time. Using ", method)}

    if (is.null(patternList)) {
        patternList <- setdiff(colnames(hotspots),c("x","y","barcode"))
    } else if (!all(patternList %in% colnames(hotspots))) {
        stop("Pattern names not found in hotspots")
    }

    binarized <- (!is.na(hotspots[,patternList]))*1
    intersects <- t(binarized) %*% binarized
    nHotspots <- colSums(binarized)
    nHotsP1 <- t(t(nHotspots)) %*% array(1, length(patternList))
    nHotsP2 <- t(nHotsP1)
    overlapScore <- switch(method,
        "Szymkiewicz-Simpson" = intersects/pmin(nHotsP1,nHotsP2),
        "Jaccard" = intersects/(nHotsP1 + nHotsP2 - intersects),
        "Sorensen-Dice" = 2*intersects/(nHotsP1 + nHotsP2),
        "Ochiai" = intersects/sqrt(nHotsP1*nHotsP2),
        "absolute" = intersects,
        stop("Method not supported")
    )

    overlapScore[upper.tri(overlapScore,diag = TRUE)] <- NA
  
    # Melt normalized Jaccard for output
    dfOverlap <- reshape2::melt(overlapScore)
    dfOverlap <- dfOverlap[stats::complete.cases(dfOverlap),]
    # Due to melting in lower triangular orientation, the column names are flipped
    colnames(dfOverlap) <- c("pattern2", "pattern1", "overlapScore")
    dfOverlap <- dfOverlap[,c(2,1,3)]
    return(dfOverlap)
}

#' @title plotOverlapScores
#' @description Plot the overlap scores between patterns in hotspots
#' @param df A data frame with columns pattern1, pattern2 and overlapScore
#' @param title The title of the plot
#' @param fontsize The font size of the plot
#' @param out The output path for the plot
#' @return A ggplot object
#' @export
#' @examples
#' df <- data.frame(pattern1 = c("pattern1","pattern1","pattern2","pattern2"), 
#'                  pattern2 = c("pattern1","pattern2","pattern1","pattern2"),
#'                  overlapScore = c(0.5,0.7,0.3,0.9))
#' plotOverlapScores(df)
#' plotOverlapScores(df, "Overlap Scores", "overlapScores.png", 15)
#' @import ggplot2
#'
plotOverlapScores <- function(df, title = "Spatial Overlap Scores", out = NULL,fontsize = 15) {
    p <- ggplot2::ggplot(data = df, aes(pattern1, pattern2, fill = overlapScore)) +
        geom_tile(color = "black", size = 0.8) +
        geom_text(aes(label = round(overlapScore, 2)), size = 6) +  # Display values on the plot
        scale_fill_gradient2(low = "#FFF7EC", mid = "#FDBB84", high = "#D7301F", midpoint = 0.5) +
        scale_y_discrete(limits = rev, guide = guide_axis(angle = 45)) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, vjust = 1, size = fontsize, hjust = 1),
            axis.text.y = element_text(size = fontsize)
        ) +
        coord_fixed() +
        labs(x = NULL, y = NULL, fill = "overlapScore", title = title) +
        theme(
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank()
        )
    if(!is.null(out)){
        ggplot2::ggsave(filename = out, plot = p)
    }
    return(p)
}

#' @title getIMScores
#' @description Get the interaction scores for SpaceMarkers
#' @param SpaceMarkers A list of SpaceMarkers objects
#' @return A data frame with columns Gene and SpaceMarkersMetric
#' @export
#' @examples
#' example(getPairwiseInteractingGenes)
#' getIMScores(SpaceMarkers)
#' @importFrom stats setNames
#' 
getIMScores <- function(SpaceMarkers){
    smi <- SpaceMarkers[which(sapply(SpaceMarkers, function(x) length(x[['interacting_genes']]))>0)]
    fields <- c('Gene', 'SpaceMarkersMetric')

    imscores <- lapply(seq_along(smi), function(x) {
        df <- smi[[x]][['interacting_genes']][[1]][,fields]
        #rename to metric to its parent item name
        stats::setNames(df, c('Gene', names(smi)[x]))
    })

    imscores <- Reduce(function(x, y) {
                merge(x, y, by="Gene", all=TRUE)
            }, x=imscores, right=FALSE)
    
    if(is.null(imscores)) {
        imscores <- data.frame(Gene=character(0))
    }
    return(imscores)
}


#' @title plotIMScores
#' @description Plot the top SpaceMarkers IMScores
#' @param df A data frame with columns Gene and SpaceMarkersMetric
#' @param interaction The interaction to plot
#' @param cutOff The cut off value for the plot
#' @param nGenes The number of genes to plot
#' @param geneText The font size for the gene text
#' @param metricText The font size for the metric text
#' @param increments The increments for the y-axis
#' @param out The output path for the plot
#' @export
#' @examples 
#' example(getPairwiseInteractingGenes)
#' plotIMScores(getIMScores(SpaceMarkers), "Pattern_1_Pattern_3")
#' @import ggplot2
#' @importFrom stats reorder
#' @importFrom utils head
plotIMScores <- function(df, interaction, cutOff = 0, nGenes = 20,
    geneText = 12, metricText = 12, increments = 1, out = NULL) {
    df$genes <- df$Gene
    df <- df[order(df[[interaction]], decreasing = TRUE),]
    df <- utils::head(df,nGenes)
    df[[interaction]] <- round(df[[interaction]],2)
    
    p1 <- ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(Gene,!!sym(interaction)), y = !!sym(interaction))) +
                geom_bar(stat = "identity",
                show.legend = FALSE,
                fill = "blue",      # Background color
                color = "blue") + 
        theme(axis.text.x = element_text(size = metricText, angle=0), axis.text.y = element_text(angle = 0, size = geneText)) +
        xlab("Genes") +
        ylab("SpaceMarkers Metric") +
        ggtitle(interaction) + 
        scale_y_continuous(breaks= seq(0, max(df[[interaction]]),
            by = increments),
            limits = c(0, max(df[[interaction]]) + 0.2)) +
        coord_flip()
    if (!is.null(out)){
        ggplot2::ggsave(filename = out,plot = p1)
    }
    return(p1)
}

#' @title calculate_gene_set_score
#' @description Calculate the mean interaction score for a set of genes
#' @param IMscores A matrix of interaction scores
#' @param genes A list of gene sets, where each set is a vector of gene names
#' @return A matrix of mean interaction scores for genes in each gene set, with 
#' attributes for log p-value sums and number of genes for later fisher combination
#' @export
calculate_gene_set_score <- function(IMscores, gene_sets, weighted = TRUE, method = c("geometric_mean", "arithmetic_mean")) {
    method <- match.arg(method[1], choices = c("geometric_mean", "arithmetic_mean"))

    IMS <- split(IMscores, IMscores$cell_interaction)

    # Pre-calculate gene overlap counts
    # How many complexes does each gene appear in?
    all_genes <- unique(unlist(gene_sets))
    gene_complex_count <- sapply(all_genes, function(g) {
        sum(sapply(gene_sets, function(gs) g %in% gs))
    })

    scores <- matrix(NA, nrow=length(gene_sets), ncol=length(IMS))
    rownames(scores) <- names(gene_sets)
    colnames(scores) <- names(IMS)

    for (i in seq_along(IMS)) {
        temp <- lapply(gene_sets, function(geneSet) {
            if (length(geneSet) == 0) return(NA)
            ims <- IMS[[i]]
            geneSet <- intersect(geneSet, ims$gene)
            if (length(geneSet) == 0) return(NA)

            subset_data <- ims[ims$gene %in% geneSet,]
            
                        # Calculate uniqueness weight for each gene
            # Weight = 1 / (number of complexes containing this gene)
            if (!weighted) {
                gene_weights <- rep(1, nrow(subset_data))
            } else {
                gene_weights <- 1 / gene_complex_count[subset_data$gene]
                # Normalize weights to sum to number of genes
                # (preserves interpretation of the score)
                gene_weights <- gene_weights * length(gene_weights) / sum(gene_weights)
            }

            # Calculate weighted effect size
            weighted_effect <- switch(method,
                "geometric_mean" = exp(weighted.mean(log(subset_data$effect_size),gene_weights)),
                "arithmetic_mean" = weighted.mean(subset_data$effect_size, gene_weights),
                stop("Unknown method for gene set score calculation")
                )

            weighted_effect

        })
        scores[,i] <- unlist(temp)
    }
    
    scores <- .call(scores, 2, as.numeric)
    colnames(scores) <- names(IMS)
    rownames(scores) <- names(gene_sets)
        # Add gene overlap information as an attribute
    attr(scores, "gene_complex_count") <- gene_complex_count
    
    return(scores)
}

#' @title calculate_lr_scores
#' @description Calculate L-R pair scores using Fisher's method
#' @param ligand_results Output from getGeneSetScore for ligands
#' @param receptor_results Output from getGeneSetScore for receptors
#' @param lr_pairs Data frame with columns 'ligand' and 'receptor'
#' @return Data frame with L-R scores and p-values
#' @export
calculate_lr_scores <- function(ligand_scores, receptor_scores, lr_pairs,
                              ligand_test = c("greater", "two.sided"), method = c("geometric_mean", "arithmetic_mean"),
                              weighted = TRUE) {

    method <- match.arg(method[1], choices = c("geometric_mean", "arithmetic_mean"))
       # parameter checks
    ligand_test <- match.arg(ligand_test[1], choices = c("greater", "two.sided"))
    
    ## Scoring ligand overexpression near target cell type and receptor overexpression in target cell type
    if (!all(grepl("near",colnames(receptor_scores)))) {
        target_cells <- gsub("^.*_","",colnames(ligand_scores))
        if (!any(target_cells %in% colnames(receptor_scores))) {
            stop("Receptor scores do not have expected cell type names in column names.")
        }
        ligand_cols <- colnames(ligand_scores)[grepl(paste0("near_(", paste(target_cells, collapse="|"), ")"), colnames(ligand_scores))]
        names(ligand_cols) <- ligand_cols
        # Map ligand columns to receptor columns
        mapped_receptor_cols <- sapply(ligand_cols, function(lc) gsub("^.*_","",lc))
        names(mapped_receptor_cols) <- ligand_cols
        lr_scores <- matrix(0, nrow=nrow(lr_pairs), ncol=length(ligand_cols))
        if (weighted) {
            ligand_counts <- table(lr_pairs$ligand.symbol)
            receptor_counts <- table(lr_pairs$receptor.symbol)
            ligand_weights <- 1 / ligand_counts[lr_pairs$ligand.symbol]
            receptor_weights <- 1 / receptor_counts[lr_pairs$receptor.symbol]
            ligand_weights <- matrix(ligand_weights * length(ligand_weights) / sum(ligand_weights), nrow=nrow(lr_pairs), ncol=1)
            receptor_weights <- matrix(receptor_weights * length(receptor_weights) / sum(receptor_weights), nrow=nrow(lr_pairs), ncol=1)
            rownames(ligand_weights) <- rownames(receptor_weights) <- rownames(lr_pairs)
            weights <- cbind(ligand_weights, receptor_weights)
            colnames(weights) <- c("ligand_weights", "receptor_weights")
            rownames(weights) <- rownames(lr_pairs)
        } else {
            weights <- matrix(1, nrow=nrow(lr_pairs), ncol=2)
            colnames(weights) <- c("ligand_weights", "receptor_weights")
            rownames(weights) <- rownames(lr_pairs)
        }

        lr_scores <- matrix(NA, nrow=nrow(lr_pairs), ncol=length(ligand_cols))
        rownames(lr_scores) <- rownames(lr_pairs)
        mapped_lr_cols <- gsub("near_","to_",ligand_cols)
        names(mapped_lr_cols) <- ligand_cols
        colnames(lr_scores) <- mapped_lr_cols
            
        for (i in seq_along(ligand_cols)) {
            lc <- ligand_cols[i]
            rc <- mapped_receptor_cols[i]
            scores <- cbind(ligand_scores[,lc], receptor_scores[,rc])
            colnames(scores) <- c("ligand_score", "receptor_score")
            rownames(scores) <- rownames(lr_pairs)
            lr_scores[,i] <- switch(method,
                "geometric_mean" = exp(sapply(rownames(scores), function(r) {weighted.mean(log(scores[r, ]), w=weights[r,])})),
                "arithmetic_mean" = sapply(rownames(scores), function(r) {weighted.mean(scores[r, ], w=weights[r,])}),
                stop("Unknown method for L-R score calculation")
            )
        }
        colnames(lr_scores) <- mapped_lr_cols
        rownames(lr_scores) <- rownames(lr_pairs)
        return(lr_scores)
    }

    ## Scoring ligand overexpression near target cell type and receptor overexpression near source cell type

 
     n_interactions <- length(intersect(colnames(ligand_scores), colnames(receptor_scores)))
    if (n_interactions == 0) {
        stop("No matching cell interactions between ligand and receptor scores.")
    }
    # Map columns 
    ligand_cols <- intersect(colnames(ligand_scores), colnames(receptor_scores))
    mapped_receptor_cols <- sapply(ligand_cols, function(x) {parts <- strsplit(x, "_near_")[[1]]; return(paste0(parts[2], "_near_", parts[1]))})
    mapped_lr_cols <- sapply(ligand_cols, function(x) {parts <- strsplit(x, "_near_")[[1]]; return(paste0(parts[1], "_to_", parts[2]))})


    # Initialize result matrix
    lr_scores <- matrix(NA, nrow=nrow(lr_pairs), ncol=n_interactions)
    lr_pvalues <- matrix(NA, nrow=nrow(lr_pairs), ncol=n_interactions)

    colnames(lr_scores) <- colnames(lr_pvalues) <- mapped_lr_cols
    rownames(lr_scores) <- rownames(lr_pvalues) <- rownames(lr_pairs)
    
    l_name <- gsub(lr_pairs$ligand.symbol, pattern = ", ", replacement = "|")
    r_name <- gsub(lr_pairs$receptor.symbol, pattern = ", ", replacement = "|")

    # Get ligand data
    l_scores <- ligand_scores
    l_log_p_sum <- attr(ligand_scores, "log_p_sums")
    l_n <- attr(ligand_scores, "n_genes")

    # Get receptor data
    r_scores <- receptor_scores
    r_log_p_sum <- attr(receptor_scores, "log_p_sums")
    r_n <- attr(receptor_scores, "n_genes")

    names(mapped_receptor_cols) <- ligand_cols
    names(mapped_lr_cols) <- ligand_cols

    # Calculate L-R score using a switch statement
    lr_scores <- switch(method,
         "signed_geometric_mean" = sqrt(abs(l_scores[,ligand_cols]) * abs(r_scores[,mapped_receptor_cols])) * sign(r_scores[,mapped_receptor_cols]),
         "arithmetic_mean" = (l_scores[,ligand_cols] + r_scores[,mapped_receptor_cols]) / 2,
         stop("Unknown method for L-R score calculation")
        )
    
    colnames(lr_scores) <- mapped_lr_cols
    rownames(lr_scores) <- rownames(lr_pairs)

    
    # Fisher's method: -2 * sum(log(p)) ~ chi-squared with 2k df
    chi_stat_l <- -2 * l_log_p_sum[,ligand_cols]
    chi_stat_r <- -2 * r_log_p_sum
    l_pval <- pchisq(chi_stat_l, df = 2 * l_n[,ligand_cols], lower.tail = FALSE) # two-sided-test for ligands
    if (ligand_test=="greater") ifelse(l_scores>0,l_pval <- l_pval / 2, 1)
    r_pval <- pchisq(chi_stat_r[,mapped_receptor_cols], df = 2 * r_n[,mapped_receptor_cols], lower.tail = FALSE)
    
    # Only keep L-R p-value if both ligand and receptor are sufficiently expressed
    lr_pvalues <- pmax(l_pval, r_pval)
    colnames(lr_pvalues) <- mapped_lr_cols
    # Handle any NA p-values (e.g., no genes in set)
    lr_pvalues[is.na(lr_pvalues)] <- 1

    rownames(lr_pvalues) <- rownames(lr_scores) <- rownames(lr_pairs)
    # Combine results
    result <- data.frame(
        ligand = l_name,
        receptor = r_name,
        lr_pair = rownames(lr_pairs)
    )
    
    # Add scores and p-values for each interaction
    score_df <- as.data.frame(lr_scores)
    names(score_df) <- paste0("score_", names(score_df))
    
    pval_df <- as.data.frame(lr_pvalues)
    names(pval_df) <- paste0("pval_", names(pval_df))
    
    result <- cbind(result, score_df, pval_df)

     if (adjust_pvals) {
        if (adjustment_scope == "global") {
            # Adjust across all L-R pairs and conditions together
            all_pvals <- as.vector(lr_pvalues)
            all_padj <- p_adjust(all_pvals, method = adjustment_method)
            lr_padj <- matrix(all_padj, nrow = nrow(lr_pvalues), ncol = ncol(lr_pvalues))
            
        } else if (adjustment_scope == "per_interaction") {
            # Adjust within each condition separately
            lr_padj <- lr_pvalues
            for (i in 1:ncol(lr_pvalues)) {
                lr_padj[, i] <- p_adjust(lr_pvalues[, i], method = adjustment_method)
            }
        }
        
        # Add adjusted p-values to results
        padj_df <- as.data.frame(lr_padj)
        names(padj_df) <- paste0("p_adj_", colnames(lr_pvalues))
        result <- cbind(result, padj_df)
        }

    # Add specificity weights as an attribute
    attr(result, "specificity_weights") <- calculate_specificity_weights(lr_pairs)

    # Add concordance information as an attribute
    lr_concordance <- sign(l_scores[,ligand_cols]) == sign(r_scores[,mapped_receptor_cols])
    colnames(lr_concordance) <- mapped_lr_cols
    rownames(lr_concordance) <- rownames(lr_pairs)
    
    attr(result, "concordance") <- lr_concordance[rownames(lr_pairs), , drop = FALSE]
    attr(result, "ligand_concordance") <- attr(ligand_scores, "concordant")[, , drop = FALSE]
    attr(result, "receptor_concordance") <- attr(receptor_scores, "concordant")[, , drop = FALSE]
    
    # Store adjustment info as attributes
    attr(result, "p_adjustment") <- list(
        adjusted = adjust_pvals,
        method = adjustment_method,
        scope = adjustment_scope,
        n_tests = length(as.vector(lr_pvalues))
    )
    
    return(result)
}

#' Calculate Gene Set Specificity Scores
#' @title calculate_gene_set_specificity
#' @description
#' This function computes specificity scores for given gene sets across cell types or spatial patterns.
#' It uses fold-change scores and p-values to weight gene contributions, and supports both geometric and arithmetic means.
#'
#' @param data A numeric matrix or data frame of gene expression values (genes x samples).
#' @param spPatterns A data frame or matrix containing spatial pattern information, with columns for cell types and optionally "x", "y", "barcode".
#' @param gene_sets A named list of character vectors, where each vector contains gene names for a gene set.
#' @param weighted Logical; if TRUE, gene scores are weighted by their occurrence in multiple gene sets.
#' @param method Character; specifies the aggregation method for gene set scores. Options are "geometric_mean" or "arithmetic_mean".
#'
#' @return A numeric matrix of gene set specificity scores (gene sets x cell types).
#'
#' @details
#' - Genes not present in the data are excluded.
#' - Genes with all zero expression are removed.
#' - Fold-change scores and p-values are calculated using `calculate_all_fc_scores`.
#' - Scores are normalized and weighted by p-value significance.
#' - For each gene set, scores are aggregated using the specified method and gene weights.
#'
#' @examples
#' # Example usage:
#' # gene_set_scores <- calculate_gene_set_specificity(expr_matrix, spPatterns, gene_sets)
#'
#' @seealso \code{\link{calculate_all_fc_scores}}
#' @export
calculate_gene_set_specificity <- function(data, spPatterns, gene_sets, weighted = TRUE, method = c("geometric_mean", "arithmetic_mean")) {
    method <- match.arg(method[1], choices = c("geometric_mean", "arithmetic_mean"))
    genes <- unique(unlist(gene_sets))
    genes <- intersect(genes, rownames(data))
    if (length(genes) == 0) {
        stop("No genes from gene sets found in data")
    }
    expr <- data[genes, , drop = FALSE]
    expr[is.na(expr)] <- 0

    expr <- expr[apply(expr, 1, function(x) any(x > 0)), , drop = FALSE]
    
    genes <- rownames(expr)

    gene_complex_count <- sapply(genes, function(g) {
        sum(sapply(gene_sets, function(gs) g %in% gs))
    })
    if (weighted) {
        gene_weights <- 1 / gene_complex_count[genes]
            # (preserves interpretation of the score)
        gene_weights <- gene_weights * length(gene_weights) / sum(gene_weights)
    } else {
        gene_weights <- rep(1, length(genes))
    }
    names(gene_weights) <- genes
            # (preserves interpretation of the score)

    if (nrow(expr) == 0) {
        return(matrix(0, nrow=length(gene_sets), ncol=ncol(spPatterns)))
    }

    # Group by cell type
    patnames <- setdiff(colnames(spPatterns),c("x","y","barcode"))
    cell_types <- unique(patnames)

    # Initialize result matrix
    gene_set_scores <- matrix(0, nrow = length(gene_sets), ncol = length(cell_types))
    rownames(gene_set_scores) <- names(gene_sets)
    colnames(gene_set_scores) <- cell_types

    gene_scores <- calculate_all_fc_scores(expr, spPatterns, low_thr = 0.2, high_thr = 0.8)


    for (gene_set in names(gene_sets)) {
        genes <- gene_sets[[gene_set]]
        valid_genes <- intersect(genes, rownames(gene_scores))
        if (length(valid_genes) > 1) {
            
            # Calculate weighted scores
            gene_set_scores[gene_set, ] <- switch(method,
                "geometric_mean" = exp(apply(log(gene_scores[valid_genes, ]), 2, weighted.mean, w = gene_weights[valid_genes])),
                "arithmetic_mean" = apply(gene_scores[valid_genes, ], 2, weighted.mean, w = gene_weights[valid_genes]),
                stop("Unknown method for gene set score calculation")
            )
        } else if (length(valid_genes) == 1) {
        gene_set_scores[gene_set, ] <- gene_scores[valid_genes, ]
        } else {
        gene_set_scores[gene_set, ] <- 0
        }
    }

    return(gene_set_scores)
}

calculate_fc_score <- function(expr, spPatterns, gene, ct,
                              low_thr = 0.2, high_thr = 0.8) {

  # Get thresholds
  low_thr <- quantile(spPatterns[, ct], probs = low_thr, na.rm = TRUE)
  high_thr <- quantile(spPatterns[, ct], probs = high_thr, na.rm = TRUE)

  # Get high and low bins
  high_bins <- which(spPatterns[, ct] > high_thr)
  low_bins <- which(spPatterns[, ct] < low_thr)
  
  # Calculate means
  mean_high <- mean(expr[gene, high_bins])
  mean_low <- mean(expr[gene, low_bins])
  
  # Log fold change
  lfc <- mean_high - mean_low
  
  # Get p-value
  if (all(expr[gene, c(high_bins, low_bins)] == 0)) {
    score <- NA
    attr(score, "p_value") <- 1
    return(score)
  }
  w_test <- wilcox.test(as.matrix(expr[gene, high_bins]), as.matrix(expr[gene, low_bins]))

  score <- lfc
  attr(score, "p_value") <- w_test$p.value
  return(score)
}

# Wrapper for all genes and cell types
# importFrom BiocParallel bplapply
calculate_all_fc_scores <- function(expr, spPatterns, 
                                    low_thr = 0.2, high_thr = 0.9, ..., workers = NULL) {

  genes <- rownames(expr)
  cell_types <- setdiff(colnames(spPatterns), c("x", "y", "barcode"))
  
  # Initialize results
  score_matrix <- matrix(NA, nrow = length(genes), ncol = length(cell_types))
  rownames(score_matrix) <- genes
  colnames(score_matrix) <- cell_types
  
  gene_scores <- lfc_matrix <- p_values <- score_matrix
  use_biocparallel <- requireNamespace("BiocParallel", quietly = TRUE)
  # Calculate for each combination
  if (!use_biocparallel) {
    for (i in seq_along(genes)) {
        for (ct in cell_types) {
            result <- calculate_fc_score(expr, spPatterns, 
                                    genes[i], ct,
                                    low_thr, high_thr)
            lfc_matrix[i, ct] <- result
            p_values[i, ct] <- attr(result, "p_value")
        }
    }
  } else {
    if (is.null(workers) || workers < 1) {
        workers <- BiocParallel::bpworkers(BiocParallel::bpparam())
    }
    BiocParallel::register(BiocParallel::MulticoreParam(workers))
    lfc <- BiocParallel::bplapply(cell_types, function(ct) {
        lfc_col <- p_val_col <- numeric(length(genes))
        names(lfc_col) <- names(p_val_col) <- genes
        for (i in seq_along(genes)) {
            result <- calculate_fc_score(expr, spPatterns, 
                                genes[i], ct,
                                low_thr, high_thr)
            lfc_col[i] <- result
            p_val_col[i] <- attr(result, "p_value")
        }
        return(list(lfc = lfc_col, p_value = p_val_col))
    })
    for (i in seq_along(cell_types)) {
        ct <- cell_types[i]
        lfc_matrix[, ct] <- lfc[[i]]$lfc
        p_values[, ct] <- lfc[[i]]$p_value
    }
  }
    p_adj <- apply(p_values,2,p.adjust,method="BH")
    p_weights <- 1 / (1 + exp(2 * (log10(p_adj) - log10(0.05))))

    norm_scores <- 1/(1 + exp(-40 * (lfc_matrix-median(lfc_matrix[lfc_matrix>0], na.rm = TRUE))))
    gene_scores <- norm_scores * p_weights
    gene_scores[is.na(gene_scores)] <- 0
    attr(gene_scores,"p_values") <- p_values
    attr(gene_scores,"p_adj") <- p_adj
    attr(gene_scores,"p_weights") <- p_weights
    attr(gene_scores,"lfc") <- lfc_matrix
  return(gene_scores)
}

.call <- function(x, MARGIN = NULL, FUN, ...) {
    # Save attributes
    saved_attrs <- attributes(x)

    # Apply function
    if (is.null(MARGIN)) {
        result <- FUN(x, ...)
    } else {
        result <- apply(x, MARGIN, FUN, ...)
    }

    # Restore custom attributes
    for (attr_name in names(saved_attrs)) {
        attr(result, attr_name) <- saved_attrs[[attr_name]]
    }
    
    return(result)
}
