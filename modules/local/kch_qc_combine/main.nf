process COMBINE_KCH_QC {

    tag { "combine_kch_qc" }
    label "process_single"
    container "docker.io/seglh/bw_validation_tools@sha256:99fb766c096bd828c9cd0a9c88dba97eec3aab8f1f4c31a99c691096e98f8e08"
    debug true

    input:
    path(QC_json_files)

    output:
    path("batch_kch_qc")
    path("*ggplot.pdf")
    path  "versions.yml"                                             ,  emit: versions
    
    script:
    loess_smoothing_level = task.ext.args ?: ''
    coverage_level = task.ext.args2 ?: ''
    read_number_required = task.ext.args3 ?: ''
    test_data = task.ext.args4 ?: ''
    truth_data = task.ext.args5 ?: ''
    
     
    if (params.qc_standalone == "yes")

    """
    cat $QC_json_files > temp_cat_qc
    header=\$(head -1  temp_cat_qc |  awk 'BEGIN{FS="\\t"; OFS="\\t"} {print \$0}')
    
    tee batch_kch_qc && printf "\$header\\n" | tee -a batch_kch_qc && grep -v "Sample_name" temp_cat_qc | tee -a batch_kch_qc > /dev/null
    
    ggplot_depth_vs_reads_1samplesets_with_duplication_by_size_july2025.R batch_kch_qc expt $test_data $loess_smoothing_level $coverage_level $read_number_required    
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1 | head -c 15)
        docker: \$(echo \$( grep 'docker run' .command.run 2>&1) )
    END_VERSIONS

    """
    
    else
    """
    cat $QC_json_files > temp_cat_qc
    header=\$(head -1  temp_cat_qc |  awk 'BEGIN{FS="\\t"; OFS="\\t"} {print \$0}')

    tee batch_kch_qc && printf "\$header\\n" | tee -a batch_kch_qc && grep -v "Sample_name" temp_cat_qc | tee -a batch_kch_qc > /dev/null

    ggplot_depth_vs_reads_2samplesets_with_duplication_by_size_july2025.R batch_kch_qc expt truth "$test_data vs $truth_data" $loess_smoothing_level $coverage_level $read_number_required


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1 | head -c 15)
        docker: \$(echo \$( grep 'docker run' .command.run 2>&1) )
    END_VERSIONS

    """
}
