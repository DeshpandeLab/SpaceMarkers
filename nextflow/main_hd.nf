process SPACEMARKERS_HD {
  tag "$meta.id"
  label 'process_medium'
  container 'ghcr.io/deshpandelab/spacemarkers:hdutils'

  input:
    tuple val(meta), path(features), path(data)
  output:
    tuple val(meta), path("${prefix}/IMscores.rds"),       val(source),   emit: IMscores
    tuple val(meta), path("${prefix}/LRscores.rds"),       val(source),   emit: LRscores
    path  "versions.yml",                                                 emit: versions

   script:
     def args = task.ext.args ?: ''
     source = features.simpleName
     prefix = task.ext.prefix ?: "${meta.id}/${source}"
     template 'hdPipeline.R'

  stub:
    def args = task.ext.args ?: ''
    source = features.simpleName
    prefix = task.ext.prefix ?: "${meta.id}/${source}"
    """
    mkdir -p "${prefix}"
    touch "${prefix}/IMscores.rds"
    touch "${prefix}/LRscores.rds"
    touch "${prefix}/top_25_interactions.csv"
    mkdir -p "${prefix}/figures"
    touch "${prefix}/figures/circos_plot.png"
    touch "${prefix}/figures/circos_top10.png"

    cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          SpaceMarkers: \$(Rscript -e 'print(packageVersion("SpaceMarkers"))' | awk '{print \$2}')
          R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
    END_VERSIONS
    """

}

process SPACEMARKERS_HD_PLOTS{
  tag "$meta.id"
  label 'process_low'
  container 'ghcr.io/deshpandelab/spacemarkers:hdutils'

  input:
    tuple val(meta), path(lrscores), val(source)
  output:
    tuple val(meta), path("${prefix}/top*.csv"),           val(source),   emit: top_interactions, optional: true
    tuple val(meta), path("${prefix}/**.png"),             val(source),   emit: figures,          optional: true
    path  "versions.yml",                                                 emit: versions,         optional: true

   script:
     def args = task.ext.args ?: ''
     prefix = task.ext.prefix ?: "${meta.id}/${source}"
     template 'hdPlotting.R'
   stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}/${source}"
    """
    mkdir -p "${prefix}/figures"
    touch "${prefix}/figures/circos_plot.png"
    touch "${prefix}/figures/circos_top10.png"
    """
}

// Nextflow pipeline to run SpaceMarkers
workflow {
    ch_versions = Channel.empty()

    ch_sm_inputs = Channel.fromPath(params.input)
    .splitCsv(header:true, sep: ",")
    .map { row-> [meta:[id:row.sample], features:file(row.annotation_file), data:file(row.data_dir)] }

    //spacemarkers - main
    SPACEMARKERS_HD( ch_sm_inputs )
    ch_versions = ch_versions.mix(SPACEMARKERS_HD.out.versions)

    SPACEMARKERS_HD_PLOTS( SPACEMARKERS_HD.out.LRscores )
    ch_versions = ch_versions.mix(SPACEMARKERS_HD_PLOTS.out.versions)

    //collate versions
    ch_versions
      .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'versions.yml', sort: true, newLine: true)
}