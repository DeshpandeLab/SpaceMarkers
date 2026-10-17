nextflow.enable.dsl = 2

process SPACEMARKERS {
  tag "$meta.id"
  label 'process_medium'
  container 'ghcr.io/deshpandelab/spacemarkers@sha256:e13854a27622a04293fd8c26e8829a0407ab08a91d06259ece02eb440eab9ae2'

  input:
    tuple val(meta), path(spatial), path(features)
    path( lr_reference )

  output:
    tuple val(meta), path("${prefix}/sme.rds"),      val(source), emit: sme
    tuple val(meta), path("${prefix}/IMscores.rds"), val(source), emit: IMscores
    tuple val(meta), path("${prefix}/LRscores.rds"), val(source), emit: LRscores, optional: true
    path  "versions.yml",                                         emit: versions

  script:
    def args = task.ext.args ?: ''
    source = spatial.simpleName
    prefix = task.ext.prefix ?: "${meta.id}/${source}"
    template 'spacemarkers.R'

  stub:
    def args = task.ext.args ?: ''
    source = spatial.simpleName
    prefix = task.ext.prefix ?: "${meta.id}/${source}"
    """
    mkdir -p "${prefix}"
    touch "${prefix}/sme.rds"
    touch "${prefix}/IMscores.rds"
    touch "${prefix}/LRscores.rds"

    cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          SpaceMarkers: \$(Rscript -e 'print(packageVersion("SpaceMarkers"))' | awk '{print \$2}')
          R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
    END_VERSIONS
    """
}

// unnamed worklow to run on a samplesheet with anndata files
// example usage
// nextflow run SpaceMarkers/nextflow/main.nf -with-docker -resume -params-file params.yaml
workflow {
  samplesheet_ch = channel.fromPath(params.input)
  samplesheet_ch
    .splitCsv(header: true)
    .map { row -> tuple([id: row.sample], file(row.spatial), row.features ? file(row.features) : []) }
    .set { spacemarkers_ch }

  // optional user-supplied ligand-receptor reference; [] stages no file
  lr_reference_ch = params.sm_lr_reference ?
    channel.value(file(params.sm_lr_reference, checkIfExists: true)) :
    channel.value([])

  SPACEMARKERS(spacemarkers_ch, lr_reference_ch)
}
