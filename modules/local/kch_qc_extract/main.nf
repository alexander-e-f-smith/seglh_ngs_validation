process EXTRACT_KCH_QC {

    tag { "$meta" }
    label "process_single"
    container "docker.io/seglh/bw_validation_tools@sha256:99fb766c096bd828c9cd0a9c88dba97eec3aab8f1f4c31a99c691096e98f8e08"
    debug true

    input:
    tuple val(meta), path(QC_json), path(QC_json_truth), path(depth_file), path(depth_file_truth)
    path(githead)

    output:
    tuple val(meta), path("${meta}_combined_consolidated_kch_qc")    ,  emit: sample_kch_qc
    path  "versions.yml"                                             ,  emit: versions
    
    script:
    coverage_settings = task.ext.args ?: ''

    """
    test1.sh ${QC_json} ${depth_file} ${meta}_consolidated_kch_qc obs
    
    test1.sh ${QC_json_truth} ${depth_file_truth} ${meta}_consolidated_kch_qc_truth exp

    paste -d "\\t" ${meta}_consolidated_kch_qc ${meta}_consolidated_kch_qc_truth |  tee -a ${meta}_combined_consolidated_kch_qc > /dev/null
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jq: \$(jq --version 2>&1 |  sed 's/^.*Version: //g')
        git: \$(head -1 ORIG_HEAD 2>&1  )
        docker: \$(echo \$( grep 'docker run' .command.run 2>&1) )
    END_VERSIONS

    """
}
