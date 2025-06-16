# Ensure required packages are installed (if not already)
# install.packages(c("circlize", "RColorBrewer"))

# --- Helper Function 1: Data Preparation ---
.prepare_circos_data_logic <- function(lr_df,
                                       all_molecules_def, # Expected: data.frame with cell_type, molecule_name, molecule_type
                                       split_segments_for_links,
                                       link_buffer_fraction,
                                       scale_link_width_by_score,
                                       score_transform_for_width,
                                       inter_molecule_segment_gap) {
    # This function now takes all_molecules_def as an argument.
    # lr_df is the already filtered data frame of interactions for the current plot.

    if (nrow(all_molecules_def) == 0) {
      return(list(bed_data = data.frame(), links_data = data.frame()))
    }

    bed_data_intermediate <- all_molecules_def

    if (split_segments_for_links) {
        if (scale_link_width_by_score) {
            lr_df$link_actual_width <- score_transform_for_width(lr_df$score)
            lr_df$link_actual_width[lr_df$link_actual_width <= 1e-6] <- 0.1 # Min width
        } else {
            lr_df$link_actual_width <- 1 # Each interaction gets 1 unit width
        }

        source_lengths <- stats::aggregate(link_actual_width ~ source_cell_type + ligand, data = lr_df, FUN = sum)
        if(nrow(source_lengths) > 0) {
          colnames(source_lengths) <- c("cell_type", "molecule_name", "total_segment_length")
          source_lengths$molecule_type <- "ligand"
        } else { # Handle empty aggregation
            source_lengths <- data.frame(cell_type=character(), molecule_name=character(), total_segment_length=numeric(), molecule_type=character())
        }


        target_lengths <- stats::aggregate(link_actual_width ~ target_cell_type + receptor, data = lr_df, FUN = sum)
        if(nrow(target_lengths) > 0){
          colnames(target_lengths) <- c("cell_type", "molecule_name", "total_segment_length")
          target_lengths$molecule_type <- "receptor"
        } else {
            target_lengths <- data.frame(cell_type=character(), molecule_name=character(), total_segment_length=numeric(), molecule_type=character())
        }

        molecule_total_lengths <- rbind(source_lengths, target_lengths)

        if(nrow(molecule_total_lengths) > 0){
            bed_data_intermediate <- merge(all_molecules_def, molecule_total_lengths,
                                       by = c("cell_type", "molecule_name", "molecule_type"), all.x = TRUE)
            bed_data_intermediate$total_segment_length[is.na(bed_data_intermediate$total_segment_length)] <- 1 # Default length
        } else { # if no interactions, molecule_total_lengths is empty
            bed_data_intermediate$total_segment_length <- 1
        }
        bed_data_intermediate$segment_length_val <- bed_data_intermediate$total_segment_length
    } else { # Not splitting segments, each molecule segment has base length 1
        bed_data_intermediate$segment_length_val <- 1
    }

    # Assign cumulative start/end for molecule segments on tracks, incorporating inter-segment gaps
    bed_data_list <- lapply(split(bed_data_intermediate, bed_data_intermediate$cell_type), function(df_cell) {
        df_cell <- df_cell[order(df_cell$molecule_type, df_cell$molecule_name), ]
        if (nrow(df_cell) > 0) {
            current_pos <- 0
            starts <- numeric(nrow(df_cell))
            ends <- numeric(nrow(df_cell))
            for (i in 1:nrow(df_cell)) {
                if (i > 1) current_pos <- current_pos + inter_molecule_segment_gap
                starts[i] <- current_pos
                ends[i] <- current_pos + df_cell$segment_length_val[i]
                current_pos <- ends[i]
            }
            df_cell$start <- starts
            df_cell$end <- ends
        } else {
            df_cell$start <- numeric(0); df_cell$end <- numeric(0);
        }
        df_cell
    })
    bed_data <- do.call(rbind, bed_data_list)
    if(!is.null(bed_data)) rownames(bed_data) <- NULL # bed_data can be NULL if no segments

    # Prepare links_data
    links_data <- data.frame() # Initialize empty
    if (nrow(lr_df) > 0 && !is.null(bed_data) && nrow(bed_data) > 0) {
        if (split_segments_for_links) {
            lr_df_temp_links <- lr_df # Already has interaction_id and link_actual_width
            lr_df_temp_links <- lr_df_temp_links[order(lr_df_temp_links$source_cell_type, lr_df_temp_links$ligand, lr_df_temp_links$interaction_id), ]
            lr_df_temp_links$source_sub_offset <- ave(lr_df_temp_links$link_actual_width, lr_df_temp_links$source_cell_type, lr_df_temp_links$ligand,
                                                FUN = function(x) c(0, cumsum(x)[-length(x)]))
            lr_df_temp_links <- lr_df_temp_links[order(lr_df_temp_links$target_cell_type, lr_df_temp_links$receptor, lr_df_temp_links$interaction_id), ]
            lr_df_temp_links$target_sub_offset <- ave(lr_df_temp_links$link_actual_width, lr_df_temp_links$target_cell_type, lr_df_temp_links$receptor,
                                                FUN = function(x) c(0, cumsum(x)[-length(x)]))
            lr_df_temp_links <- lr_df_temp_links[order(lr_df_temp_links$interaction_id), ]

            links_data_list_apply <- apply(lr_df_temp_links, 1, function(row) {
                lig_info <- bed_data[bed_data$cell_type == row["source_cell_type"] & bed_data$molecule_name == row["ligand"] & bed_data$molecule_type == "ligand", ]
                rec_info <- bed_data[bed_data$cell_type == row["target_cell_type"] & bed_data$molecule_name == row["receptor"] & bed_data$molecule_type == "receptor", ]
                if (nrow(lig_info) == 1 && nrow(rec_info) == 1) {
                    actual_w = as.numeric(row["link_actual_width"])
                    s_start = lig_info$start + as.numeric(row["source_sub_offset"]) + (actual_w * link_buffer_fraction)
                    s_end   = lig_info$start + as.numeric(row["source_sub_offset"]) + actual_w * (1 - link_buffer_fraction)
                    t_start = rec_info$start + as.numeric(row["target_sub_offset"]) + (actual_w * link_buffer_fraction)
                    t_end   = rec_info$start + as.numeric(row["target_sub_offset"]) + actual_w * (1 - link_buffer_fraction)
                    if (s_start >= s_end) { s_end <- s_start + 1e-6 }
                    if (t_start >= t_end) { t_end <- t_start + 1e-6 }
                    data.frame(source_cell_type = row["source_cell_type"], source_start = s_start, source_end = s_end,
                               target_cell_type = row["target_cell_type"], target_start = t_start, target_end = t_end,
                               score = as.numeric(row["score"]), stringsAsFactors = FALSE)
                } else { NULL }
            })
        } else { # Not splitting, links use full segment width from bed_data
            links_data_list_apply <- apply(lr_df, 1, function(row) {
                lig_info <- bed_data[bed_data$cell_type == row["source_cell_type"] & bed_data$molecule_name == row["ligand"] & bed_data$molecule_type == "ligand",]
                rec_info <- bed_data[bed_data$cell_type == row["target_cell_type"] & bed_data$molecule_name == row["receptor"] & bed_data$molecule_type == "receptor",]
                if (nrow(lig_info) == 1 && nrow(rec_info) == 1) {
                    data.frame(source_cell_type = row["source_cell_type"], source_start = lig_info$start, source_end = lig_info$end,
                               target_cell_type = row["target_cell_type"], target_start = rec_info$start, target_end = rec_info$end,
                               score = as.numeric(row["score"]), stringsAsFactors = FALSE)
                } else { NULL }
            })
        }
        # Filter out NULLs from apply if any rows failed to produce link data
        links_data <- do.call(rbind, Filter(Negate(is.null), links_data_list_apply))
        if(!is.null(links_data)) rownames(links_data) <- NULL
    }
    return(list(bed_data = bed_data, links_data = links_data))
}

# --- Helper Function 2: Color Generation ---
.generate_circos_colors <- function(molecules_in_plot_df, # Should be unique molecules: molecule_name, molecule_type
                                    cell_types_in_plot_vector,
                                    use_individual_molecule_colors,
                                    default_ligand_color,
                                    default_receptor_color,
                                    individual_ligand_palette_generator,
                                    individual_receptor_palette_generator,
                                    cell_type_palette_name = "Set1") {
    n_cell_types <- length(cell_types_in_plot_vector)
    pal_ct_max <- RColorBrewer::brewer.pal.info[cell_type_palette_name, "maxcolors"]
    pal_ct <- RColorBrewer::brewer.pal(min(n_cell_types, pal_ct_max), cell_type_palette_name)
    cell_type_colors_vec <- if (n_cell_types > length(pal_ct)) grDevices::colorRampPalette(pal_ct)(n_cell_types) else pal_ct[1:n_cell_types]
    names(cell_type_colors_vec) <- cell_types_in_plot_vector

    name_cols_ligand <- NULL; name_cols_receptor <- NULL
    if(use_individual_molecule_colors) {
        gen_mol_cols_helper <- function(mol_names, pal_gen, pal_default_name, type) {
            n_unique <- length(mol_names)
            if (n_unique == 0) return(stats::setNames(character(0), character(0)))
            colors_vec <- if (is.null(pal_gen)) {
                max_brewer <- RColorBrewer::brewer.pal.info[pal_default_name, "maxcolors"]
                # Ensure n for brewer.pal is at least 3, or adjust if n_unique is less
                n_for_brewer <- if (n_unique < 3 && n_unique > 0) 3 else n_unique
                base_pal <- RColorBrewer::brewer.pal(min(n_for_brewer, max_brewer), pal_default_name)
                if (n_unique <= length(base_pal)) base_pal[1:n_unique] else grDevices::colorRampPalette(base_pal)(n_unique)
            } else { pal_gen(n_unique) }
            # Ensure colors_vec length matches mol_names length
            if(length(colors_vec) != n_unique) {
                warning(paste("Palette generator for", type, "did not return correct number of colors. Using default."))
                # Fallback to default logic
                max_brewer <- RColorBrewer::brewer.pal.info[pal_default_name, "maxcolors"]
                n_for_brewer <- if (n_unique < 3 && n_unique > 0) 3 else n_unique
                base_pal <- RColorBrewer::brewer.pal(min(n_for_brewer, max_brewer), pal_default_name)
                colors_vec <- if (n_unique <= length(base_pal)) base_pal[1:n_unique] else grDevices::colorRampPalette(base_pal)(n_unique)
            }
            stats::setNames(colors_vec, mol_names)
        }
        lig_names <- sort(unique(molecules_in_plot_df$molecule_name[molecules_in_plot_df$molecule_type == "ligand"]))
        rec_names <- sort(unique(molecules_in_plot_df$molecule_name[molecules_in_plot_df$molecule_type == "receptor"]))
        
        name_cols_ligand <- gen_mol_cols_helper(lig_names, individual_ligand_palette_generator, "Greens", "ligands")
        name_cols_receptor <- gen_mol_cols_helper(rec_names, individual_receptor_palette_generator, "Blues", "receptors")
    }
    return(list(cell_type_colors = cell_type_colors_vec,
                name_colors_ligand = name_cols_ligand,
                name_colors_receptor = name_cols_receptor))
}

# @import circlize
# @import RColorBrewer
# --- Helper Function 3: Drawing Circlize Elements ---
.draw_circos_plot_elements <- function(bed_data,
                                       links_data,
                                       plot_cell_order,
                                       cell_type_colors,
                                       name_colors_ligand, # Used if use_individual_molecule_colors is TRUE
                                       name_colors_receptor, # Used if use_individual_molecule_colors is TRUE
                                       # Plotting parameters passed from top-level function
                                       gap_degree_after_sector, track_height_molecules,
                                       molecule_label_cex, cell_label_cex,
                                       link_transparency, link_connection_rou=0.7,
                                       link_arrowhead_type,
                                       link_arrowhead_width, link_arrowhead_length,
                                       score_color_palette_fun,
                                       molecule_segment_border_col,
                                       use_individual_molecule_colors,
                                       default_ligand_color,
                                       default_receptor_color,
                                       show_score_legend = TRUE) {
    circlize::circos.clear()
    circlize::circos.par(
        gap.after = stats::setNames(rep(gap_degree_after_sector, length(plot_cell_order)), plot_cell_order),
        start.degree = 90, 
        track.margin = c(0.01, 0.01), 
        cell.padding = c(0.02, 0, 0.02, 0),
        canvas.xlim = c(-1.1, 1.1), # Ensure full circle
        canvas.ylim = c(-1.1, 1.1),
        points.overflow.warning = FALSE
        )
    sector_xlim_df <- stats::aggregate(end ~ cell_type, data = bed_data, FUN = max, na.rm = TRUE)
    ordered_sector_xlim <- data.frame(cell_type = plot_cell_order, stringsAsFactors = FALSE)
    ordered_sector_xlim <- merge(ordered_sector_xlim, sector_xlim_df, by = "cell_type", all.x = TRUE, sort=FALSE) # Keep order
    ordered_sector_xlim$end[is.na(ordered_sector_xlim$end)] <- 1
    ordered_sector_xlim$start <- 0
    
    circlize::circos.initialize(factors = ordered_sector_xlim$cell_type, xlim = as.matrix(ordered_sector_xlim[, c("start", "end")]))

    circlize::circos.track(ylim = c(0, 1), bg.border = NA, track.height = track_height_molecules,
        panel.fun = function(x, y) {
            sector_name = circlize::get.cell.meta.data("sector.index")
            circlize::circos.text(x = circlize::CELL_META$xcenter, y = circlize::CELL_META$cell.ylim[2] + circlize::uy(6, "mm"),
                labels = sector_name, facing = "bending.inside", niceFacing = TRUE, adj = c(0.5, 0), cex = cell_label_cex, col = cell_type_colors[sector_name])
            
            sector_molecules <- bed_data[bed_data$cell_type == sector_name, ]
            if (nrow(sector_molecules) > 0) {
                for (i in 1:nrow(sector_molecules)) {
                    mol <- sector_molecules[i, ]
                    segment_col <- if(use_individual_molecule_colors) {
                        cols <- if(mol$molecule_type == "ligand") name_colors_ligand else name_colors_receptor
                        # Fallback if molecule name not in generated palette (should not happen with current logic)
                        col_val <- cols[mol$molecule_name]; if(is.na(col_val) || length(col_val)==0) ifelse(mol$molecule_type == "ligand", default_ligand_color, default_receptor_color) else col_val
                    } else { if(mol$molecule_type == "ligand") default_ligand_color else default_receptor_color }
                    if(is.na(segment_col) || length(segment_col) == 0 || segment_col == "") segment_col <- "grey" 

                    circlize::circos.rect(xleft = mol$start, ybottom = circlize::CELL_META$ylim[1], xright = mol$end, ytop = circlize::CELL_META$ylim[2],
                                          col = segment_col, border = molecule_segment_border_col,
                                          sector.index = sector_name, track.index = circlize::get.current.track.index())
                    circlize::circos.text(x = (mol$start + mol$end) / 2, y = mean(circlize::CELL_META$ylim), labels = mol$molecule_name,
                                          facing = "clockwise", niceFacing = TRUE, cex = molecule_label_cex, adj = c(0.5, 0.5), col = "black",
                                          sector.index = sector_name, track.index = circlize::get.current.track.index())
                }
            }
        })

    if (!is.null(links_data) && nrow(links_data) > 0) {
        link_cols_vec <- rep(grDevices::adjustcolor("grey50", alpha.f = link_transparency), nrow(links_data))
        if (!is.null(score_color_palette_fun) && "score" %in% colnames(links_data)) {
             if (is.function(score_color_palette_fun)) {
                link_scores <- links_data$score
                score_domain <- try(attr(score_color_palette_fun, "breaks"), silent = TRUE) 
                if (!inherits(score_domain, "try-error") && !is.null(score_domain) && length(score_domain) > 0) {
                     link_scores[!is.na(link_scores)] <- pmax(link_scores[!is.na(link_scores)], min(score_domain, na.rm=TRUE), na.rm=TRUE)
                     link_scores[!is.na(link_scores)] <- pmin(link_scores[!is.na(link_scores)], max(score_domain, na.rm=TRUE), na.rm=TRUE)
                }
                valid_scores_idx <- !is.na(link_scores)
                if(any(valid_scores_idx)){
                    raw_cols <- score_color_palette_fun(link_scores[valid_scores_idx])
                    link_cols_vec[valid_scores_idx] <- grDevices::adjustcolor(raw_cols, alpha.f = link_transparency)
                }
              } else { warning("`score_color_palette_fun` is not a function. Using default link colors.") }
        }
        
        #arr_width_val <- if(is.null(link_arrowhead_width)) circlize:::parse_circos_par_from_cherries("arrow.width") else link_arrowhead_width
        #arr_length_val <- if(is.null(link_arrowhead_length)) circlize:::parse_circos_par_from_cherries("arrow.length") else link_arrowhead_length
        for (i in 1:nrow(links_data)) {
            link <- links_data[i, ]
            if (link$source_start >= link$source_end || link$target_start >= link$target_end) next
            circlize::circos.link(
                sector.index1 = link$source_cell_type, 
                point1 = c(link$source_start, link$source_end),
                sector.index2 = link$target_cell_type, 
                point2 = c(link$target_start, link$target_end),
                col = link_cols_vec[i],
                rou = link_connection_rou,
                border = NA,
                directional = 1, 
                arr.type = "big.arrow",
                arr.col = link_cols_vec[i], 
                arr.width = link_arrowhead_width, 
                arr.length = link_arrowhead_length
            )
        }
    }
        # After the main plot is drawn but BEFORE circos.clear()
    if (show_score_legend && !is.null(score_color_palette_fun) && is.function(score_color_palette_fun)) {
        
        # This check is to avoid trying to make a legend if there's no color function
        if(!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
            warning("Package 'ComplexHeatmap' is needed to show the legend. Please install it.")
        } else {
            lgd_links = ComplexHeatmap::Legend(
                col_fun = score_color_palette_fun,
                title = "Interaction Score",
                direction = "vertical"
            )
            ComplexHeatmap::draw(
                lgd_links,
                x = grid::unit(1, "npc") - grid::unit(5, "mm"), # Adjusted position
                y = grid::unit(5, "mm"),
                just = c("right", "bottom")
            )
        }
    }
    circlize::circos.clear()
}


# --- Top-Level Function 1: General Cell Interaction Plot ---
#' Plot Ligand-Receptor Interactions between Cell Types
#'
#' Visualizes ligand-receptor interactions as a circular plot, allowing for
#' selection of cell types, customization of segment and link appearance,
#' and different modes of representing interactions.
#'
#' @param lr_interactions_df A data.frame with columns: `ligand`, `receptor`,
#'   `source_cell_type`, `target_cell_type`, `score`.
#' @param selected_cell_types Optional character vector. If provided, only these cell types
#'   and interactions *between* them will be shown.
#' @param cell_order Optional character vector for specific cell type ordering. If NULL,
#'   alphabetical order is used for the (selected) cell types.
#' @param gap_degree_after_sector Numeric, gap in degrees after each sector. Default is 5.
#' @param track_height_molecules Numeric, height of the track for molecule segments. Default is 0.1.
#' @param molecule_label_cex Numeric, cex for molecule labels. Default is 0.6.
#' @param cell_label_cex Numeric, cex for cell type labels. Default is 0.9.
#' @param link_transparency Numeric, alpha for links (0 to 1). Default is 0.5.
#' @param link_connection_rou Numeric (0-1) or vector of two for `rou` in `circos.link`. Default is 0.7.
#' @param link_arrowhead_type Character. Type of arrowhead (e.g., "triangle", "big.arrow"). Default is "big.arrow".
#' @param link_arrowhead_width Numeric. Width of the arrowhead. Default uses circlize default.
#' @param link_arrowhead_length Numeric. Length of the arrowhead. Default uses circlize default.
#' @param score_color_palette_fun A function (e.g., from `circlize::colorRamp2`) to map scores to colors.
#' @param split_segments_for_links Logical. If TRUE, molecule segments are sized by interaction count/score
#'   and links connect to unique sub-segments. Default is TRUE.
#' @param link_buffer_fraction Numeric (0 to <0.5). Buffer around each link within its sub-segment.
#'   Applies if `split_segments_for_links` is TRUE. Default is 0.1.
#' @param scale_link_width_by_score Logical. If TRUE (and `split_segments_for_links` is TRUE),
#'   link/sub-segment width is proportional to score. Default is FALSE.
#' @param score_transform_for_width Function to transform scores for width scaling. Default is `function(s) s`.
#' @param use_individual_molecule_colors Logical. If TRUE, each unique ligand/receptor name gets a distinct color.
#'   Default is TRUE.
#' @param default_ligand_color Character, color for all ligand segments if `use_individual_molecule_colors` is FALSE.
#'   Default is "lightgreen".
#' @param default_receptor_color Character, color for all receptor segments if `use_individual_molecule_colors` is FALSE.
#'   Default is "lightblue".
#' @param individual_ligand_palette_generator Function (takes n, returns n colors) for unique ligands.
#'   Default generates a green palette.
#' @param individual_receptor_palette_generator Function (takes n, returns n colors) for unique receptors.
#'   Default generates a blue palette.
#' @param molecule_segment_border_col Color for the border of L/R segments. Default `NA` (no border).
#' @param inter_molecule_segment_gap Numeric, gap between L/R segments on the same track. Default is 0.1.
#'
#' @return Invisibly returns NULL. The function is called for its side effect of creating a plot.
#' @export
#' @examples
#' \dontrun{
#' if (requireNamespace("circlize", quietly = TRUE) &&
#'     requireNamespace("RColorBrewer", quietly = TRUE)) {
#'   lr_data <- data.frame(
#'     ligand = c("LGF1", "LGF1", "LGF2", "LGF3", "LGF4", "LGF1", "LGF5", "LGF6"),
#'     receptor = c("REC1", "REC2A", "REC1", "REC3B", "REC1", "REC4", "REC4", "REC2A"),
#'     source_cell_type = c("CellA", "CellA", "CellB", "CellC", "CellA", "CellD", "CellD", "CellB"),
#'     target_cell_type = c("CellB", "CellC", "CellB", "CellA", "CellD", "CellA", "CellD", "CellC"),
#'     score = c(4, 0.9, 0.7, 0.95, 0.6, 0.1, 0.88, 2.5),
#'     stringsAsFactors = FALSE
#'   )
#'   # Basic plot with defaults
#'   # plotCellInteractionCircos(lr_data)
#'
#'   # Plot with selected cell types and custom order
#'   # plotCellInteractionCircos(lr_data,
#'   #                           selected_cell_types = c("CellA", "CellB", "CellD"),
#'   #                           cell_order = c("CellD", "CellA", "CellB"))
#'  }
#' }
#' @import circlize
#' @import RColorBrewer
#' @export
plotCellInteractionCircos <- function(lr_interactions_df,
                                     selected_cell_types = NULL,
                                     cell_order = NULL,
                                     gap_degree_after_sector = 5,
                                     track_height_molecules = 0.1,
                                     molecule_label_cex = 0.6,
                                     cell_label_cex = 0.9,
                                     link_transparency = 0.5,
                                     link_connection_rou = 0.7,
                                     link_arrowhead_type = "big.arrow",
                                     link_arrowhead_width = NULL,
                                     link_arrowhead_length = NULL,
                                     score_color_palette_fun = NULL,
                                     split_segments_for_links = TRUE,
                                     link_buffer_fraction = 0.1,
                                     scale_link_width_by_score = FALSE,
                                     score_transform_for_width = function(s) s,
                                     use_individual_molecule_colors = TRUE,
                                     default_ligand_color = "lightgreen",
                                     default_receptor_color = "lightblue",
                                     individual_ligand_palette_generator = NULL,
                                     individual_receptor_palette_generator = NULL,
                                     molecule_segment_border_col = NA,
                                     inter_molecule_segment_gap = 0.1) {

    if (!requireNamespace("circlize", quietly = TRUE)) stop("Package 'circlize' is needed.", call. = FALSE)
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) stop("Package 'RColorBrewer' is needed.", call. = FALSE)
    if (nrow(lr_interactions_df) == 0) { warning("Input 'lr_interactions_df' is empty."); return(invisible(NULL)) }
    if (link_buffer_fraction >= 0.5 || link_buffer_fraction < 0) stop("'link_buffer_fraction' must be between 0 and < 0.5.", call. = FALSE)

    effective_lr_df <- lr_interactions_df
    if (!is.null(selected_cell_types) && length(selected_cell_types) > 0) {
        effective_lr_df <- lr_interactions_df[
            lr_interactions_df$source_cell_type %in% selected_cell_types &
            lr_interactions_df$target_cell_type %in% selected_cell_types,
        ]
        if (nrow(effective_lr_df) == 0) {
            warning("No interactions remain after filtering for selected_cell_types.")
            return(invisible(NULL))
        }
        current_cell_order <- cell_order
        if(!is.null(current_cell_order) && !all(current_cell_order %in% selected_cell_types)){
            warning("'cell_order' contains types not in 'selected_cell_types' or not in data after filtering. It will be adjusted.")
            current_cell_order <- current_cell_order[current_cell_order %in% unique(c(effective_lr_df$source_cell_type, effective_lr_df$target_cell_type))]
        }
        cell_order <- current_cell_order # Update cell_order based on filtering
    }
    if (nrow(effective_lr_df) == 0) { warning("No interaction data to plot after initial processing."); return(invisible(NULL)) }
    effective_lr_df$interaction_id <- 1:nrow(effective_lr_df)


    source_mols <- unique(effective_lr_df[, c("source_cell_type", "ligand")])
    if(nrow(source_mols)>0) {colnames(source_mols) <- c("cell_type", "molecule_name"); source_mols$molecule_type <- "ligand"
    } else { source_mols <- data.frame(cell_type=character(), molecule_name=character(), molecule_type=character())}

    target_mols <- unique(effective_lr_df[, c("target_cell_type", "receptor")])
    if(nrow(target_mols)>0) {colnames(target_mols) <- c("cell_type", "molecule_name"); target_mols$molecule_type <- "receptor"
    } else {target_mols <- data.frame(cell_type=character(), molecule_name=character(), molecule_type=character())}
    
    all_molecules_definition <- unique(rbind(source_mols, target_mols))
    if(nrow(all_molecules_definition) == 0 && nrow(effective_lr_df) > 0) { # Should not happen if effective_lr_df has rows
        warning("Could not define any molecules for segments, though interactions exist.")
        return(invisible(NULL))
    }


    prepared_data <- .prepare_circos_data_logic(
        lr_df = effective_lr_df, all_molecules_def = all_molecules_definition,
        split_segments_for_links = split_segments_for_links, link_buffer_fraction = link_buffer_fraction,
        scale_link_width_by_score = scale_link_width_by_score, score_transform_for_width = score_transform_for_width,
        inter_molecule_segment_gap = inter_molecule_segment_gap
    )

    if (is.null(prepared_data$bed_data) || nrow(prepared_data$bed_data) == 0) {
        warning("No molecule segments to plot after data preparation.")
        return(invisible(NULL))
    }

    # Determine plot_cell_order based on bed_data and user's cell_order preference
    unique_cells_in_bed <- sort(unique(prepared_data$bed_data$cell_type))
    current_plot_cell_order <- if (!is.null(cell_order) && length(cell_order) > 0) {
        # Ensure provided cell_order is valid for the cells actually in bed_data
        ordered_cells <- cell_order[cell_order %in% unique_cells_in_bed]
        # Append any missing cells from bed_data alphabetically
        missing_cells <- unique_cells_in_bed[!unique_cells_in_bed %in% ordered_cells]
        if (length(missing_cells) > 0) {
             warning(paste("Some cells in data were not in 'cell_order':", paste(missing_cells, collapse=", "), ". Appending them."))
        }
        c(ordered_cells, missing_cells)
    } else {
        unique_cells_in_bed # Alphabetical
    }
    if (length(current_plot_cell_order) == 0) { warning("No cell types to plot."); return(invisible(NULL)) }

    plot_colors <- .generate_circos_colors(
        molecules_in_plot_df = unique(prepared_data$bed_data[, c("molecule_name", "molecule_type")]),
        cell_types_in_plot_vector = current_plot_cell_order,
        use_individual_molecule_colors = use_individual_molecule_colors, default_ligand_color = default_ligand_color,
        default_receptor_color = default_receptor_color, individual_ligand_palette_generator = individual_ligand_palette_generator,
        individual_receptor_palette_generator = individual_receptor_palette_generator
    )

    .draw_circos_plot_elements(
        bed_data = prepared_data$bed_data, links_data = prepared_data$links_data,
        plot_cell_order = current_plot_cell_order, cell_type_colors = plot_colors$cell_type_colors,
        name_colors_ligand = plot_colors$name_colors_ligand, name_colors_receptor = plot_colors$name_colors_receptor,
        gap_degree_after_sector=gap_degree_after_sector, track_height_molecules=track_height_molecules,
        molecule_label_cex=molecule_label_cex, cell_label_cex=cell_label_cex,
        link_transparency=link_transparency, link_connection_rou=link_connection_rou,
        link_arrowhead_type=link_arrowhead_type,
        link_arrowhead_width=link_arrowhead_width, link_arrowhead_length=link_arrowhead_length,
        score_color_palette_fun=score_color_palette_fun, molecule_segment_border_col=molecule_segment_border_col,
        use_individual_molecule_colors=use_individual_molecule_colors,
        default_ligand_color=default_ligand_color, default_receptor_color=default_receptor_color
    )
    return(invisible(NULL))
}


# --- Top-Level Function 2: Source-to-Target Plot ---
#' Plot Ligand-Receptor Interactions from a Single Source to Target Cell Types
#'
#' Visualizes ligand-receptor interactions focusing on a single source cell type
#' and its outgoing interactions to a specified set of target cell types. Only
#' ligands from the source and relevant receptors on the targets are shown.
#' @inheritParams plotCellInteractionCircos
#' @param lr_interactions_df A data.frame with columns: `ligand`, `receptor`,
#'   `source_cell_type`, `target_cell_type`, `score`.
#' @param source_cell_name Character, name of the source cell type.
#' @param target_cell_names Character vector, names of target cell types to show.
#' @param cell_order Optional character vector for specific cell type ordering.
#'   If NULL, source cell is first, then targets alphabetically.
#'
#' @return Invisibly returns NULL. The function is called for its side effect of creating a plot.
#' @export
#' @examples
#' \dontrun{
#' if (requireNamespace("circlize", quietly = TRUE) &&
#'     requireNamespace("RColorBrewer", quietly = TRUE)) {
#'   lr_data <- data.frame(
#'     ligand = c("LGF1", "LGF1", "LGF2", "LGF3", "LGF4", "LGF1", "LGF5", "LGF6", "LGF7"),
#'     receptor = c("REC1", "REC2A", "REC1", "REC3B", "REC1", "REC4", "REC4", "REC2A", "REC5"),
#'     source_cell_type = c("CellA", "CellA", "CellB", "CellC", "CellA", "CellD", "CellD", "CellB", "CellA"),
#'     target_cell_type = c("CellB", "CellC", "CellB", "CellA", "CellD", "CellA", "CellD", "CellC", "CellE"),
#'     score = c(4,0.9,0.7,0.95,0.6,0.1,0.88,2.5,3.0),
#'     stringsAsFactors = FALSE
#'   )
#'   # Plot interactions from CellA to CellB and CellC
#'   # plotSourceToTargetCircos(lr_data,
#'   #                          source_cell_name = "CellA",
#'   #                          target_cell_names = c("CellB", "CellC"))
#'
#'   # Custom order and score-based link widths
#'   # score_pal <- circlize::colorRamp2(c(0, 4), c("white", "blue"))
#'   # plotSourceToTargetCircos(lr_data,
#'   #                          source_cell_name = "CellD",
#'   #                          target_cell_names = c("CellA", "CellD"), # Include autocrine
#'   #                          cell_order = c("CellD", "CellA"),
#'   #                          scale_link_width_by_score = TRUE,
#'   #                          score_color_palette_fun = score_pal)
#'  }
#' }
#' @import circlize
#' @import RColorBrewer
#' 
#' @export
plotSourceToTargetCircos <- function(lr_interactions_df,
                                     source_cell_name,
                                     target_cell_names,
                                     cell_order = NULL,
                                     gap_degree_after_sector = 5,
                                     track_height_molecules = 0.1,
                                     molecule_label_cex = 0.6,
                                     cell_label_cex = 0.9,
                                     link_transparency = 0.5,
                                     link_connection_rou = 0.7,
                                     link_arrowhead_type = "triangle",
                                     link_arrowhead_width = NULL,
                                     link_arrowhead_length = NULL,
                                     score_color_palette_fun = NULL,
                                     split_segments_for_links = TRUE,
                                     link_buffer_fraction = 0.1,
                                     scale_link_width_by_score = FALSE,
                                     score_transform_for_width = function(s) s,
                                     use_individual_molecule_colors = TRUE,
                                     default_ligand_color = "lightgreen",
                                     default_receptor_color = "lightblue",
                                     individual_ligand_palette_generator = NULL,
                                     individual_receptor_palette_generator = NULL,
                                     molecule_segment_border_col = NA,
                                     inter_molecule_segment_gap = 0.1) {

    if (nrow(lr_interactions_df) == 0) { warning("Input 'lr_interactions_df' is empty."); return(invisible(NULL)) }
    if (missing(source_cell_name) || length(source_cell_name) != 1) stop("'source_cell_name' must be a single string.")
    if (missing(target_cell_names) || length(target_cell_names) == 0) stop("'target_cell_names' must be provided.")

    effective_lr_df <- lr_interactions_df[
        lr_interactions_df$source_cell_type == source_cell_name &
        lr_interactions_df$target_cell_type %in% target_cell_names,
    ]
    if (nrow(effective_lr_df) == 0) {
        warning("No interactions found for the specified source and target cells.")
        return(invisible(NULL))
    }
    effective_lr_df$interaction_id <- 1:nrow(effective_lr_df)


    source_ligands <- unique(effective_lr_df[effective_lr_df$source_cell_type == source_cell_name, c("source_cell_type", "ligand")])
    if(nrow(source_ligands)>0) {colnames(source_ligands) <- c("cell_type", "molecule_name"); source_ligands$molecule_type <- "ligand"
    } else {source_ligands <- data.frame(cell_type=character(), molecule_name=character(), molecule_type=character())}


    target_receptors <- unique(effective_lr_df[effective_lr_df$target_cell_type %in% target_cell_names, c("target_cell_type", "receptor")])
    if(nrow(target_receptors)>0) {colnames(target_receptors) <- c("cell_type", "molecule_name"); target_receptors$molecule_type <- "receptor"
    } else {target_receptors <- data.frame(cell_type=character(), molecule_name=character(), molecule_type=character())}

    all_molecules_definition <- unique(rbind(source_ligands, target_receptors))
    if(nrow(all_molecules_definition) == 0){
        warning("No molecules defined for source/target cells after filtering for source-target view.")
        return(invisible(NULL))
    }

    prepared_data <- .prepare_circos_data_logic(
        lr_df = effective_lr_df, all_molecules_def = all_molecules_definition,
        split_segments_for_links = split_segments_for_links, link_buffer_fraction = link_buffer_fraction,
        scale_link_width_by_score = scale_link_width_by_score, score_transform_for_width = score_transform_for_width,
        inter_molecule_segment_gap = inter_molecule_segment_gap
    )
    if (is.null(prepared_data$bed_data) || nrow(prepared_data$bed_data) == 0) {
        warning("No molecule segments to plot after data preparation for source-target view.")
        return(invisible(NULL))
    }

    # Determine plot_cell_order for this mode
    unique_cells_in_bed <- sort(unique(prepared_data$bed_data$cell_type))
    default_plot_order <- c(source_cell_name[source_cell_name %in% unique_cells_in_bed],
                            sort(target_cell_names[target_cell_names %in% unique_cells_in_bed & target_cell_names != source_cell_name]))

    current_plot_cell_order <- if (!is.null(cell_order) && length(cell_order) > 0) {
        # Ensure provided cell_order is valid for the cells actually in bed_data for this mode
        ordered_cells <- cell_order[cell_order %in% unique_cells_in_bed]
        missing_cells <- unique_cells_in_bed[!unique_cells_in_bed %in% ordered_cells]
        if (length(missing_cells) > 0) {
             warning(paste("Some relevant cells were not in 'cell_order' for source-target plot:", paste(missing_cells, collapse=", "), ". Appending them."))
        }
        c(ordered_cells, missing_cells)
    } else {
        default_plot_order
    }
    current_plot_cell_order <- unique(current_plot_cell_order) # Ensure unique in case source was in targets
    if (length(current_plot_cell_order) == 0) { warning("No cell types to plot in source-target view."); return(invisible(NULL)) }


    plot_colors <- .generate_circos_colors(
        molecules_in_plot_df = unique(prepared_data$bed_data[, c("molecule_name", "molecule_type")]),
        cell_types_in_plot_vector = current_plot_cell_order,
        use_individual_molecule_colors = use_individual_molecule_colors, default_ligand_color = default_ligand_color,
        default_receptor_color = default_receptor_color, individual_ligand_palette_generator = individual_ligand_palette_generator,
        individual_receptor_palette_generator = individual_receptor_palette_generator
    )

    .draw_circos_plot_elements(
        bed_data = prepared_data$bed_data, links_data = prepared_data$links_data,
        plot_cell_order = current_plot_cell_order, cell_type_colors = plot_colors$cell_type_colors,
        name_colors_ligand = plot_colors$name_colors_ligand, name_colors_receptor = plot_colors$name_colors_receptor,
        gap_degree_after_sector=gap_degree_after_sector, track_height_molecules=track_height_molecules,
        molecule_label_cex=molecule_label_cex, cell_label_cex=cell_label_cex,
        link_transparency=link_transparency, link_connection_rou=link_connection_rou,
        link_arrowhead_type=link_arrowhead_type,
        link_arrowhead_width=link_arrowhead_width, link_arrowhead_length=link_arrowhead_length,
        score_color_palette_fun=score_color_palette_fun, molecule_segment_border_col=molecule_segment_border_col,
        use_individual_molecule_colors=use_individual_molecule_colors,
        default_ligand_color=default_ligand_color, default_receptor_color=default_receptor_color,
    )
    return(invisible(NULL))
}

#' Plot Ligand-Receptor Interactions from Multiple Source to a Single Target Cell Type
#'
#' Visualizes ligand-receptor interactions focusing on a single target cell type
#' and its incoming interactions from a specified set of source cell types. Only
#' relevant ligands from the source cells and receptors on the target cell are shown.
#' @inheritParams plotSourceToTargetCircos
#' @param lr_interactions_df A data.frame with columns: `ligand`, `receptor`,
#'   `source_cell_type`, `target_cell_type`, `score`.
#' @param source_cell_names Character vector, names of source cell types to show.
#' @param target_cell_name Character, name of the single target cell type.
#' @param cell_order Optional character vector for specific cell type ordering.
#'   If NULL, target cell is first, then sources alphabetically.
#'
#' @return Invisibly returns NULL. The function is called for its side effect of creating a plot.
#' @export
#' @import circlize
#' @import RColorBrewer
#' @examples
#' \dontrun{
#' if (requireNamespace("circlize", quietly = TRUE) &&
#'     requireNamespace("RColorBrewer", quietly = TRUE)) {
#'   # Use the same test data from previous examples
#'   test_lr_data <- data.frame(
#'     ligand = c("LGF1","LGF1","LGF2","LGF3","LGF4","LGF1","LGF5","LGF6","LGF5"),
#'     receptor = c("REC1","REC2","REC1","REC1","REC5","REC4","REC4","REC6","REC2"),
#'     source_cell_type = c("CellA","CellA","CellB","CellC","CellA","CellD","CellD","CellE","CellD"),
#'     target_cell_type = c("CellB","CellC","CellA","CellD","CellD","CellA","CellD","CellF","CellA"),
#'     score = c(4.0,2.5,3.0,2.2,1.5,1.0,3.5,0.5,2.8),
#'     stringsAsFactors = FALSE
#'   )
#'   # Plot interactions from CellB and CellD converging on CellA
#'   # plotTargetFromSourcesCircos(test_lr_data,
#'   #                             source_cell_names = c("CellB", "CellD"),
#'   #                             target_cell_name = "CellA")
#'
#'   # Another example: Interactions from CellA and CellC targeting CellD
#'   # plotTargetFromSourcesCircos(test_lr_data,
#'   #                             source_cell_names = c("CellA", "CellC"),
#'   #                             target_cell_name = "CellD",
#'   #                             scale_link_width_by_score = TRUE)
#'  }
#' }
plotTargetFromSourcesCircos <- function(lr_interactions_df,
                                        source_cell_names,
                                        target_cell_name,
                                        cell_order = NULL,
                                        # Pass through all other parameters
                                        gap_degree_after_sector = 5,
                                        track_height_molecules = 0.1,
                                        molecule_label_cex = 0.6,
                                        cell_label_cex = 0.9,
                                        link_transparency = 0.5,
                                        link_connection_rou = 0.7,
                                        link_arrowhead_type = "triangle",
                                        link_arrowhead_width = 0.1, # Explicit default
                                        link_arrowhead_length = 0.1, # Explicit default
                                        score_color_palette_fun = NULL,
                                        split_segments_for_links = TRUE,
                                        link_buffer_fraction = 0.1,
                                        scale_link_width_by_score = FALSE,
                                        score_transform_for_width = function(s) s,
                                        use_individual_molecule_colors = TRUE,
                                        default_ligand_color = "lightgreen",
                                        default_receptor_color = "lightblue",
                                        individual_ligand_palette_generator = NULL,
                                        individual_receptor_palette_generator = NULL,
                                        molecule_segment_border_col = NA,
                                        inter_molecule_segment_gap = 0.1) {

    # --- Parameter Validation & Package Checks ---
    if (nrow(lr_interactions_df) == 0) { warning("Input 'lr_interactions_df' is empty."); return(invisible(NULL)) }
    if (missing(source_cell_names) || length(source_cell_names) == 0) stop("'source_cell_names' must be provided.")
    if (missing(target_cell_name) || length(target_cell_name) != 1) stop("'target_cell_name' must be a single string.")

    # --- 1. Filter interactions for sources -> single target ---
    effective_lr_df <- lr_interactions_df[
        lr_interactions_df$source_cell_type %in% source_cell_names &
        lr_interactions_df$target_cell_type == target_cell_name,
    ]
    if (nrow(effective_lr_df) == 0) {
        warning("No interactions found for the specified source and target cells.")
        return(invisible(NULL))
    }
    effective_lr_df$interaction_id <- 1:nrow(effective_lr_df)

    # --- 2. Define specialized all_molecules_definition ---
    # Source cells: only ligands involved in the filtered interactions
    source_ligands <- unique(effective_lr_df[effective_lr_df$source_cell_type %in% source_cell_names, c("source_cell_type", "ligand")])
    if(nrow(source_ligands)>0) {
      colnames(source_ligands) <- c("cell_type", "molecule_name"); source_ligands$molecule_type <- "ligand"
    } else {
      source_ligands <- data.frame(cell_type=character(), molecule_name=character(), molecule_type=character())
    }

    # Target cell: only receptors involved in the filtered interactions
    target_receptors <- unique(effective_lr_df[effective_lr_df$target_cell_type == target_cell_name, c("target_cell_type", "receptor")])
    if(nrow(target_receptors)>0) {
      colnames(target_receptors) <- c("cell_type", "molecule_name"); target_receptors$molecule_type <- "receptor"
    } else {
      target_receptors <- data.frame(cell_type=character(), molecule_name=character(), molecule_type=character())
    }
    
    all_molecules_definition <- unique(rbind(source_ligands, target_receptors))
    if(nrow(all_molecules_definition) == 0){
        warning("No molecules defined for source/target cells after filtering.")
        return(invisible(NULL))
    }

    # --- 3. Prepare Data using Helper ---
    prepared_data <- .prepare_circos_data_logic(
        lr_df = effective_lr_df,
        all_molecules_def = all_molecules_definition,
        split_segments_for_links = split_segments_for_links,
        link_buffer_fraction = link_buffer_fraction,
        scale_link_width_by_score = scale_link_width_by_score,
        score_transform_for_width = score_transform_for_width,
        inter_molecule_segment_gap = inter_molecule_segment_gap
    )
    if (is.null(prepared_data$bed_data) || nrow(prepared_data$bed_data) == 0) {
        warning("No molecule segments to plot after data preparation.")
        return(invisible(NULL))
    }

    # --- 4. Determine Cell Order for Plotting ---
    unique_cells_in_bed <- sort(unique(prepared_data$bed_data$cell_type))
    actual_source_cells_in_plot <- sort(source_cell_names[source_cell_names %in% unique_cells_in_bed])
    default_plot_order <- c(target_cell_name, actual_source_cells_in_plot)
    
    current_plot_cell_order <- if (!is.null(cell_order) && length(cell_order) > 0) {
        expected_cells <- default_plot_order
        if(!all(expected_cells %in% cell_order) || !all(cell_order %in% expected_cells)){
             warning("Custom 'cell_order' is inconsistent. Using default: target then sorted sources.")
             default_plot_order
        } else { cell_order }
    } else {
        default_plot_order
    }
    current_plot_cell_order <- current_plot_cell_order[current_plot_cell_order %in% unique_cells_in_bed]
    if (length(current_plot_cell_order) == 0) { warning("No cell types to plot."); return(invisible(NULL)) }

    # --- 5. Generate Colors using Helper ---
    plot_colors <- .generate_circos_colors(
        molecules_in_plot_df = unique(prepared_data$bed_data[, c("molecule_name", "molecule_type")]),
        cell_types_in_plot_vector = current_plot_cell_order,
        use_individual_molecule_colors = use_individual_molecule_colors,
        default_ligand_color = default_ligand_color,
        default_receptor_color = default_receptor_color,
        individual_ligand_palette_generator = individual_ligand_palette_generator,
        individual_receptor_palette_generator = individual_receptor_palette_generator
    )

    # --- 6. Draw Plot using Helper ---
    .draw_circos_plot_elements(
        bed_data = prepared_data$bed_data,
        links_data = prepared_data$links_data,
        plot_cell_order = current_plot_cell_order,
        cell_type_colors = plot_colors$cell_type_colors,
        name_colors_ligand = plot_colors$name_colors_ligand,
        name_colors_receptor = plot_colors$name_colors_receptor,
        # Pass through all relevant plotting params
        gap_degree_after_sector=gap_degree_after_sector, track_height_molecules=track_height_molecules,
        molecule_label_cex=molecule_label_cex, cell_label_cex=cell_label_cex,
        link_transparency=link_transparency, link_connection_rou=link_connection_rou,
        link_arrowhead_type=link_arrowhead_type,
        link_arrowhead_width=link_arrowhead_width, link_arrowhead_length=link_arrowhead_length,
        score_color_palette_fun=score_color_palette_fun,
        molecule_segment_border_col=molecule_segment_border_col,
        use_individual_molecule_colors = use_individual_molecule_colors,
        default_ligand_color = default_ligand_color,
        default_receptor_color = default_receptor_color
    )
    return(invisible(NULL))
}