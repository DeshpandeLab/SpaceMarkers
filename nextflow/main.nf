process SPACEMARKERS {
  tag "$meta.id"
  label 'process_high_memory'
  container 'ghcr.io/deshpandelab/spacemarkers:main'

  input:
    tuple val(meta), path(features), path(data)
  output:
    tuple val(meta), path("${prefix}/spPatterns.rds"),         val(source),   emit: spPatterns
    tuple val(meta), path("${prefix}/optParams.rds"),          val(source),   emit: optParams
    tuple val(meta), path("${prefix}/spaceMarkersObject.rds"), val(source),   emit: spaceMarkers
    tuple val(meta), path("${prefix}/spaceMarkers.csv"),       val(source),   emit: spaceMarkersScores
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
    
    #load spatial coords from tissue positions, deconvolved patterns, and expression
    coords <- load10XCoords("$data")
    features <- getSpatialFeatures("$features")
    dataMatrix <- load10XExpr("$data")

    #add spatial coordinates to deconvolved data, only use barcodes present in data
    spPatterns <- merge(coords, features, by.x = "barcode", by.y = "row.names")
    spPatterns <- spPatterns[which(spPatterns[,"barcode"] %in% colnames(dataMatrix)),]
    saveRDS(spPatterns, file = "${prefix}/spPatterns.rds")

    #remove genes with low expression, only barcodes present in spatial data
    keepGenes <- which(apply(dataMatrix, 1, sum) > 10)
    keepBarcodes <- which(colnames(dataMatrix) %in% spPatterns[,"barcode"])
    dataMatrix <- dataMatrix[keepGenes, keepBarcodes]

    #compute optimal parameters for spatial patterns
    optParams <- getSpatialParameters(spPatterns, visiumDir="$data");
    saveRDS(optParams, file = "${prefix}/optParams.rds")

    #find hotspots in spatial patterns
    hotspots <- findAllHotspots(spPatterns);
    saveRDS(hotspots, file = "${prefix}/hotspots.rds");

    #find regions of overlapping spatial patterns
    overlaps <- getOverlapScores(hotspots);
    write.csv(overlaps, file = "${prefix}/overlapScores.csv", row.names = FALSE);

    #find genes that are differentially expressed in spatial patterns
    spaceMarkers <- getPairwiseInteractingGenes(data = dataMatrix,
                                                  optParams = optParams,
                                                  spPatterns = spPatterns,
                                                  hotspots = hotspots,
                                                  mode = "DE",
                                                  analysis = "enrichment",
                                                  minOverlap = 10,
					                                        workers=$task.cpus)
    saveRDS(spaceMarkers, file = "${prefix}/spaceMarkersObject.rds")

    #save SpaceMarkers Interaction Scores
    IMScores <- getIMScores(spaceMarkers)
    write.csv(IMScores, file = "${prefix}/spaceMarkers.csv", row.names = FALSE)

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
    touch "${prefix}/spPatterns.rds"
    touch "${prefix}/optParams.rds"
    touch "${prefix}/spaceMarkers.csv"
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
  container 'ghcr.io/deshpandelab/spacemarkers:main'

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

  #getOverlapScores needs factors to be ordered
  overlaps[["pattern1"]] <- factor(overlaps[["pattern1"]], 
                                    levels = unique(overlaps[["pattern1"]]))
  overlaps[["pattern2"]] <- factor(overlaps[["pattern2"]], 
                                    levels = unique(overlaps[["pattern2"]]))
  plot <- plotOverlapScores(overlaps) + ggplot2::labs(subtitle="$meta.id")
  ggplot2::ggsave("${prefix}/overlapScores.png", plot)

  #plot interaction plots
  sm <- read.csv("$spaceMarkers")
  plot_names <- names(sm[,(tolower(names(sm))!="gene")])
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
  container 'ghcr.io/deshpandelab/spacemarkers:main'

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
    .map { row-> tuple(meta:[id:row.sample], features:file(row.annotation_file), data:file(row.data_dir)) }

    //spacemarkers - main
    SPACEMARKERS( ch_sm_inputs )
    ch_versions = ch_versions.mix(SPACEMARKERS.out.versions)


    //collate versions
    ch_versions
      .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'versions.yml', sort: true, newLine: true)
}