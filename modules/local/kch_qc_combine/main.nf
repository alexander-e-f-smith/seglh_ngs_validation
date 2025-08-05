process COMBINE_KCH_QC {

    tag { "combine_kch_qc" }
    label "process_single"
    container "docker.io/seglh/alex_validationtools:1.0.0@sha256:bec3658f87699b978eb8b2009f2a42f08cf328913bd58f58fa09b1a080890881"
    debug true

    input:
    path(QC_json_files)

    output:
    path("batch_kch_qc")
    path("*ggplot.pdf")
    
    script:
    loess_smoothing_level = task.ext.args ?: '1'
    coverage_level = task.ext.args ?: '400'
    read_number_required = task.ext.args ?: '12000000'
    

    """
    cat $QC_json_files > temp_cat_qc
    header=\$(head -1  temp_cat_qc |  awk 'BEGIN{FS="\\t"; OFS="\\t"} {print \$0}')
    
    tee batch_kch_qc && printf "\$header\\n" | tee -a batch_kch_qc && grep -v "Sample_name" temp_cat_qc | tee -a batch_kch_qc > /dev/null
    
    ggplot_depth_vs_reads_2samplesets_with_duplication_by_size_july2025.R batch_kch_qc expt truth Truth_VAF_linearity_ $loess_smoothing_level $coverage_level $read_number_required    
    

    #cat <<-END_VERSIONS > versions.yml
    #"${task.process}":
    #    bcftools: \$(bcftools 2>&1 | grep Version | sed 's/^.*Version: //g' |  sed 's/ /_/g')
    #END_VERSIONS

    """
}
