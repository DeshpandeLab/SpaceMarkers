process SPACEMARKERS {
  tag "$meta.id"
  label 'process_high_memory'
  container 'ghcr.io/deshpandelab/spacemarkers:main'

  input:
    tuple val(meta), path(features), path(data)
  output:
    tuple val(meta), path("${prefix}/spPatterns.rds"), val(source),         emit: spPatterns
    tuple val(meta), path("${prefix}/optParams.rds"), val(source),          emit: optParams
    tuple val(meta), path("${prefix}/spaceMarkersObject.rds"), val(source), emit: spaceMarkers
    tuple val(meta), path("${prefix}/hotspots.rds"), val(source),           emit: hotspots
    path  "versions.yml",                                                   emit: versions

  stub:
    def args = task.ext.args ?: ''
    source = features.simpleName
    prefix = task.ext.prefix ?: "${meta.id}/${source}"
    """
    mkdir -p "${prefix}"
    touch "${prefix}/spPatterns.rds"
    touch "${prefix}/optParams.rds"
    touch "${prefix}/spaceMarkersObject.rds"
    cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          SpaceMarkers: \$(Rscript -e 'print(packageVersion("SpaceMarkers"))' | awk '{print \$2}')
          R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
    END_VERSIONS
    """
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

    #find genes that are differentially expressed in spatial patterns
    spaceMarkers <- getPairwiseInteractingGenes(data = dataMatrix,
                                                  optParams = optParams,
                                                  spPatterns = spPatterns,
                                                  hotspots = hotspots,
                                                  mode = "DE",
                                                  analysis="enrichment",
					                                        workers=$task.cpus)

    saveRDS(spaceMarkers, file = "${prefix}/spaceMarkersObject.rds")

    # Get the versions of the packages
    spaceMarkersVersion <- packageVersion("SpaceMarkers")
    rVersion <- packageVersion("base")
    cat(sprintf('"%s":\n  SpaceMarkers: %s\n  R: %s\n', 
            "${task.process}", spaceMarkersVersion, rVersion), 
        file = "versions.yml")
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

process SPACEMARKERS_IMSCORES {
  tag "$meta.id"
  label 'process_low'
  container 'ghcr.io/deshpandelab/spacemarkers:main'

  input:
    tuple val(meta), path(spaceMarkers), val(source)
  output:
    tuple val(meta), path("${prefix}/spaceMarkers.csv"), emit: spacemarkers_imscores
    path  "versions.yml",                                   emit: versions

  script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}/${source}"
    """
    #!/usr/bin/env Rscript
    dir.create("${prefix}", showWarnings = FALSE, recursive = TRUE)

    sm <- readRDS("$spaceMarkers")
    smi <- sm[which(sapply(sm, function(x) length(x[['interacting_genes']]))>0)]
  
    fields <- c('Gene', 'SpaceMarkersMetric')

    imscores <- lapply(seq_along(smi), function(x) {
        df <- smi[[x]][['interacting_genes']][[1]][,fields]
        #rename to metric to its parent item name
        setNames(df, c('Gene', names(smi)[x]))
    })

    imscores <- Reduce(function(x, y) {
                merge(x, y, by="Gene", all=TRUE)
            }, x=imscores, right=FALSE)
    
    if(is.null(imscores)) {
        imscores <- data.frame(Gene=character(0))
    }
  
    write.csv(imscores, file = "${prefix}/spaceMarkers.csv", row.names = FALSE)

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
    touch "${prefix}/imscores.csv"
    cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          SpaceMarkers: \$(Rscript -e 'print(packageVersion("SpaceMarkers"))' | awk '{print \$2}')
          R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
    END_VERSIONS
    """
}
