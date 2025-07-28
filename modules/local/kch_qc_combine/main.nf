process COMBINE_KCH_QC {

    tag { "$meta" }
    label "process_single"
    container "docker.io/seglh/alex_validationtools@sha256:bec3658f87699b978eb8b2009f2a42f08cf328913bd58f58fa09b1a080890881"
    debug true

    input:
    tuple val(meta), path(QC_json), path(QC_json_truth), path(depth_file), path(depth_file_truth)

    output:
    tuple val(meta), path("${meta}*_consolidated_kch_qc*")
    
    script:
    coverage_settings = task.ext.args ?: ''

    """
    test1.sh ${QC_json} ${depth_file} ${meta}_consolidated_kch_qc obs
    
    test1.sh ${QC_json_truth} ${depth_file_truth} ${meta}_consolidated_kch_qc_truth exp

    paste -d "\\t" ${meta}_consolidated_kch_qc ${meta}_consolidated_kch_qc_truth |  tee -a ${meta}_combined_consolidated_kch_qc > /dev/null
    
    #cat <<-END_VERSIONS > versions.yml
    #"${task.process}":
    #    bcftools: \$(bcftools 2>&1 | grep Version | sed 's/^.*Version: //g' |  sed 's/ /_/g')
    #END_VERSIONS

    """
}
