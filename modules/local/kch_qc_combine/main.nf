process COMBINE_QC {

    tag { "combine_qc" }
    label "process_single"
    container "docker.io/seglh/snappy_python3_ngstools:1.0.2@sha256:c762d2ca67e52068b7bf9ef4ba4ff4466a44f8b3eead2ba5b91dd29ca7193bd1"
    debug true

    input:
    tuple val(meta), path(QC_json), path(QC_json_truth), path(depth_file), path(depth_file_truth)

    output:
    tuple val(meta), path("${meta}_consolidated_kch_qc*")
    
    script:
    coverage_settings = task.ext.args ?: ''

    """
    test1.sh ${QC_json} ${depth_file} ${meta} ${meta}_consolidated_kch_qc
    
    test1.sh ${QC_json_truth} ${depth_file_truth} ${meta} ${meta}_consolidated_kch_qc_truth 
    
    #cat <<-END_VERSIONS > versions.yml
    #"${task.process}":
    #    bcftools: \$(bcftools 2>&1 | grep Version | sed 's/^.*Version: //g' |  sed 's/ /_/g')
    #END_VERSIONS

    """
}
