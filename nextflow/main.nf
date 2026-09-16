process SPACEMARKERS {
  tag "$meta.id"
  label 'process_high_memory'
  container 'ghcr.io/deshpandelab/spacemarkers@sha256:9c06f8f9340bb5c51300dbf3bc4e803613a15e1bd349eae43d5a129462a13f4e'

  input:
    tuple val(meta), path(features), path(data)
  output:
    tuple val(meta), path("${prefix}/sme.rds"),                 val(source),   emit: sme
    tuple val(meta), path("${prefix}/spPatterns.rds"),         val(source),   emit: spPatterns
    tuple val(meta), path("${prefix}/optParams.rds"),          val(source),   emit: optParams
    tuple val(meta), path("${prefix}/spaceMarkersObject.rds"), val(source),   emit: spaceMarkers
    tuple val(meta), path("${prefix}/hotspots.rds"),           val(source),   emit: hotspots
    tuple val(meta), path("${prefix}/overlapScores.csv"),      val(source),   emit: overlapScores
    path  "versions.yml",                                                     emit: versions

  script:
    def args = task.ext.args ?: ''
    source = features.simpleName
    prefix = task.ext.prefix ?: "${meta.id}/${source}"
    """
    #!/usr/bin/env Rscript
    dir.create("${prefix}", showWarnings = FALSE, recursive = TRUE)
    library("SpaceMarkers")
    set.seed(${params.seed})

    # ---- Setup: build the SpaceMarkersExperiment in one call ----
    sme <- load10X(
      visiumDir  = "$data",
      features   = "$features",
      method     = "CSV",
      resolution = "lowres"
    )

    # ---- Undirected SpaceMarkers (one line) ----
    sme <- SpaceMarkers(sme, directed = FALSE, cpus = $task.cpus, minOverlap = 10)

    # ---- Directed SpaceMarkers (one line), only if a ligand-receptor
    # reference CSV is configured via params.lr_reference (empty string by
    # default -- see nextflow.config). Runs that never set it keep
    # producing exactly the same outputs as before this change.
    lr_reference_path <- "${params.lr_reference}"
    run_directed <- nzchar(lr_reference_path) && file.exists(lr_reference_path)

    if (run_directed) {
      message("Running directed SpaceMarkers using lr_reference: ", lr_reference_path)
      LR_df <- read.csv(lr_reference_path, row.names = 1)
      LR_df[["ligand.symbol"]]   <- LR_df[["ligand"]]
      LR_df[["receptor.symbol"]] <- LR_df[["receptor"]]
      sme <- SpaceMarkers(sme, directed = TRUE, lr_pairs = LR_df)
    } else if (nzchar(lr_reference_path)) {
      warning("params.lr_reference was set to '", lr_reference_path,
              "' but that file does not exist; skipping directed SpaceMarkers.")
    }

    # ---- New output: the full SpaceMarkersExperiment (both analyses, if run) ----
    saveRDS(sme, file = "${prefix}/sme.rds")

    # ---- Backward-compatible legacy outputs (same file names/shapes as before) ----
    spPatterns <- data.frame(
      barcode = colnames(sme),
      y       = SpatialExperiment::spatialCoords(sme)[, "y"],
      x       = SpatialExperiment::spatialCoords(sme)[, "x"],
      as.data.frame(spatial_patterns(sme)),
      row.names = NULL,
      check.names = FALSE
    )
    saveRDS(spPatterns, file = "${prefix}/spPatterns.rds")

    optParams <- spatial_params(sme)
    saveRDS(optParams, file = "${prefix}/optParams.rds")

    hotspots_undirected <- hotspots(sme, "undirected")
    saveRDS(hotspots_undirected, file = "${prefix}/hotspots.rds")

    overlaps <- overlap_scores(sme)
    write.csv(overlaps, file = "${prefix}/overlapScores.csv", row.names = FALSE)

    #find genes that are differentially expressed in spatial patterns
    spaceMarkers <- interactions(sme)
    saveRDS(spaceMarkers, file = "${prefix}/spaceMarkersObject.rds")

    #save SpaceMarkers Interaction Scores
    IMScores <- get_im_scores(spaceMarkers)
    rownames(IMScores) <- IMScores[,"Gene"]
    IMScores[,"Gene"] <- NULL
    write.csv(IMScores, file = "${prefix}/IMScores.rds", row.names = FALSE)

    # Get the versions of the packages
    spaceMarkersVersion <- packageVersion("SpaceMarkers")
    rVersion <- packageVersion("base")
    cat(sprintf('"%s":\n  SpaceMarkers: %s\n  R: %s\n',
            "${task.process}", spaceMarkersVersion, rVersion),
        file = "versions.yml")
    """
    stub:
    def args = task.ext.args ?: ''
    source = features.simpleName
    prefix = task.ext.prefix ?: "${meta.id}/${source}"
    """
    mkdir -p "${prefix}"
    touch "${prefix}/sme.rds"
    touch "${prefix}/spPatterns.rds"
    touch "${prefix}/optParams.rds"
    touch "${prefix}/spaceMarkersObject.rds"
    touch "${prefix}/hotspots.rds"
    touch "${prefix}/overlapScores.csv"
    cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          SpaceMarkers: \$(Rscript -e 'print(packageVersion("SpaceMarkers"))' | awk '{print \$2}')
          R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
    END_VERSIONS
    """
}

process SPACEMARKERS_PLOTS {
  tag "$meta.id"
  label 'process_low'
  container 'ghcr.io/deshpandelab/spacemarkers@sha256:9c06f8f9340bb5c51300dbf3bc4e803613a15e1bd349eae43d5a129462a13f4e'

  input:
  tuple val(meta), path(spaceMarkers), path(overlapScores), val(source)

  output:
  tuple val(meta), path("${prefix}/overlapScores.png"),         val(source),     emit: overlapScores_plot
  tuple val(meta), path("${prefix}/*_interacting_genes.png"),   val(source),     emit: interaction_plots,   optional:true
  path  "versions.yml",                                                          emit: versions

  script:
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${meta.id}/${source}/plots"
  """
  #!/usr/bin/env Rscript
  dir.create("${prefix}", showWarnings = FALSE, recursive = TRUE)
  library("SpaceMarkers")
  overlaps <- read.csv("$overlapScores")

  set.seed(${params.seed})

  #getOverlapScores needs factors to be ordered
  overlaps[["pattern1"]] <- factor(overlaps[["pattern1"]], 
                                    levels = unique(overlaps[["pattern1"]]))
  overlaps[["pattern2"]] <- factor(overlaps[["pattern2"]], 
                                    levels = unique(overlaps[["pattern2"]]))
  plot <- plot_overlap_scores(overlaps) + ggplot2::labs(subtitle="$meta.id")
  ggplot2::ggsave("${prefix}/overlapScores.png", plot)

  #plot interaction plots
  sm <- readRDS("$spaceMarkers")
  plot_names <- colnames(sm)
  for (plot_name in plot_names) {
    plot <- plotIMScores(sm, plot_name) + ggplot2::labs(subtitle="$meta.id")
    ggplot2::ggsave(paste0("${prefix}/", plot_name, "_interacting_genes.png"), plot)
  }

  # Get the versions of the packages
  spaceMarkersVersion <- packageVersion("SpaceMarkers")
  rVersion <- packageVersion("base")
  cat(sprintf('"%s":\n  SpaceMarkers: %s\n  R: %s\n', 
        "${task.process}", spaceMarkersVersion, rVersion), 
        file = "versions.yml")
  """
  stub: 
  def args = task.ext.args ?: ''
  source = overlapScores.simpleName
  prefix = task.ext.prefix ?: "${meta.id}/${source}"
  """
  mkdir -p "${prefix}"
  touch "${prefix}/overlapScores.png"

  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SpaceMarkers: \$(Rscript -e 'print(packageVersion("SpaceMarkers"))' | awk '{print \$2}')
        R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
  END_VERSIONS
  """
}


process SPACEMARKERS_MQC {
  tag "$meta.id"
  label 'process_low'
  container 'ghcr.io/deshpandelab/spacemarkers@sha256:9c06f8f9340bb5c51300dbf3bc4e803613a15e1bd349eae43d5a129462a13f4e'

  input:
    tuple val(meta), path(spaceMarkers), val(source)
  output:
    tuple val(meta), path("${prefix}/spacemarkers_mqc.json"), emit: spacemarkers_mqc
    path  "versions.yml",                                     emit: versions

  script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}/${source}"
    mqc_sample = prefix.replaceAll("/", "_")
    """
    #!/usr/bin/env Rscript
    dir.create("${prefix}", showWarnings = FALSE, recursive = TRUE)

    set.seed(${params.seed})

    #[['']] notation needed to allow nextflow var susbtitution

    #init all report variables
    n_pairs_total <- NA
    n_pairs_interact <- NA
    min_spacemarker_metric <- NA
    max_spacemarker_metric <- NA
    min_genes <- NA
    max_genes <- NA
    avg_hotspot_area <- NA

    sm <- readRDS("$spaceMarkers")
    smi <- sm[which(sapply(sm, function(x) length(x[['interacting_genes']]))>0)]

    #interacting patterns stats
    n_pairs_total <- length(sm)
    n_pairs_interact <- length(smi)

    if(n_pairs_interact >0 ) {
      #spacemarker metric
      max_spacemarker_metric <- max(sapply(smi, function(x) {
        max(x[['interacting_genes']][[1]][['SpaceMarkersMetric']])
      }))
      min_spacemarker_metric <- min(sapply(smi, function(x) {
        min(x[['interacting_genes']][[1]][['SpaceMarkersMetric']])
      }))

      #average number of genes in each pair
      min_genes <- min(sapply(smi, function(x) {
        nrow(x[['interacting_genes']][[1]])
      }))

      #average number of genes in each pair
      max_genes <- max(sapply(smi, function(x) {
        nrow(x[['interacting_genes']][[1]])
      }))

      #average percent overlap across interacting patterns
      avg_hotspot_area <- mean(sapply(smi, function(x) {
        sum(!is.na(x[['hotspots']]))/length(x[['hotspots']][,1])
      }))
    }

    #report
    report_data <- list(
      "$mqc_sample" = list(
        'Pairs Total' = n_pairs_total,
        'Pairs Interact' = n_pairs_interact,
        'SpaceMarker Metric' = sprintf('%0.1f - %0.1f', min_spacemarker_metric, max_spacemarker_metric),
        'Gene Count' = sprintf('%0.f - %0.f', min_genes, max_genes),
        'Mean Hotspot Area' = avg_hotspot_area
      )
    )

    report <- list(
        id = "spacemarkers_mqc",
        section_name = "SpaceMarkers",
        description = "Tool to identify genes associated with latent space interactions in spatial transcriptomics.",
        plot_type = "table",
        pconfig = list(
            id = "custom_data_table",
            title = "SpacemMarkers Stats"
            ),
        data = report_data
    )
    jsonlite::write_json(
              x=report, 
              path = "${prefix}/spacemarkers_mqc.json", 
              auto_unbox = TRUE, 
              pretty = TRUE)
    
    # Get the versions of the packages
    spaceMarkersVersion <- packageVersion("SpaceMarkers")
    rVersion <- packageVersion("base")
    cat(sprintf('"%s":\n  SpaceMarkers: %s\n  R: %s\n', 
            "${task.process}", spaceMarkersVersion, rVersion), 
        file = "versions.yml")
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}/${source}"
    """
    mkdir -p "${prefix}"
    touch "${prefix}/spacemarkers_mqc.json"
    cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          SpaceMarkers: \$(Rscript -e 'print(packageVersion("SpaceMarkers"))' | awk '{print \$2}')
          R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
    END_VERSIONS
    """
}


// Nextflow pipeline to run SpaceMarkers
workflow {
    ch_versions = Channel.empty()

    ch_sm_inputs = Channel.fromPath(params.input)
    .splitCsv(header:true, sep: ",")
    .map { row-> [meta:[id:row.sample], features:file(row.annotation_file), data:file(row.data_dir)] }

    //spacemarkers - main
    SPACEMARKERS( ch_sm_inputs )
    ch_versions = ch_versions.mix(SPACEMARKERS.out.versions)


    //collate versions
    ch_versions
      .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'versions.yml', sort: true, newLine: true)
}