process BATCH_CONCORDANT_VARIANTS_PLOTTING_VARDICT {

    tag { "batch_plot_concordant" }
    label "plot_concordant"
    container "docker.io/seglh/alex_validationtools:1.0.0@sha256:bec3658f87699b978eb8b2009f2a42f08cf328913bd58f58fa09b1a080890881"

    input:
    path(batch_sample_concordant)
    path(batch_truth_concordant)

    output:
    path("batch_combined_concordant_variants_headed.tsv")
    path("*_correlation.pdf")

    script:
    def target = task.ext.args ?: ''
    def expected_data_source = task.ext.args2 ?: ''
    def observed_data_source = task.ext.args3 ?: ''
    

    """
    #extract_essential_qc_and_depth_from_snappy2.sh combined_kch_qc
    head="CHROM\\tPOS\\tend\\tREF\\tALT\\tAF\\tVD\\tDP\\tSAMPLE\\tFILTER"
    headtruth="exp_CHROM\\texp_POS\\texp_end\\texp_REF\\texp_ALT\\texp_AF\\texp_VD\\texp_DP\\texp_SAMPLE\\texp_FILTER"
    for file in $batch_sample_concordant
    do
        id=\$(basename \$file .vcf)
        bgzip \${file} && \\
        tabix \${file}.gz && \\
        bcftools query -f '%CHROM\\t%POS\\tend\\t%REF\\t%ALT\\t[%AF\\t%VD\\t%DP\\t%SAMPLE]\\t%FILTER\\n'  \${file}.gz | sort -n -k1 -k2 >  \${id}_batch_sample_concordant_variants.tsv
    done
    echo -e \$head | tee -a batch_sample_concordant_variants_headed.tsv > /dev/null && \\
    cat *_batch_sample_concordant_variants.tsv  | tee -a batch_sample_concordant_variants_headed.tsv > /dev/null

   for file in $batch_truth_concordant
   do
        id=\$(basename \$file .vcf)
        bgzip \${file} && \\
        tabix \${file}.gz && \\
        bcftools query -f '%CHROM\\t%POS\\tend\\t%REF\\t%ALT\\t[%AF\\t%VD\\t%DP\\t%SAMPLE]\\t%FILTER\\n'  \${file}.gz | sort -n -k1 -k2 > \${id}_batch_truth_concordant_variants.tsv
    done
    echo -e \$headtruth | tee -a batch_truth_concordant_variants_headed.tsv > /dev/null && \\
    cat *_batch_truth_concordant_variants.tsv | tee -a batch_truth_concordant_variants_headed.tsv > /dev/null && \\

    paste -d '\\t' batch_sample_concordant_variants_headed.tsv batch_truth_concordant_variants_headed.tsv > batch_combined_concordant_variants_headed.tsv && \\

    correlation.R batch_combined_concordant_variants_headed.tsv Vardict VAF AF exp_AF ${expected_data_source} ${observed_data_source} batch 1 1
    correlation.R batch_combined_concordant_variants_headed.tsv Vardict VD VD exp_VD ${expected_data_source} ${observed_data_source} batch 5000 5000
    #cat <<-END_VERSIONS > versions.yml
    #"${task.process}":
    #    bcftools: \$(bcftools 2>&1 | grep Version | sed 's/^.*Version: //g' |  sed 's/ /_/g')
    #END_VERSIONS

    """

}
