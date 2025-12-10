process BATCH_CONCORDANT_VARIANTS_PLOTTING_VARDICT {

    tag { "batch_plot_concordant" }
    label "plot_concordant"
    container "docker.io/seglh/bw_validation_tools@sha256:99fb766c096bd828c9cd0a9c88dba97eec3aab8f1f4c31a99c691096e98f8e08"

    input:
    path(batch_sample_concordant_A)
    path(batch_truth_concordant_A)
    path(batch_sample_concordant_B)
    path(batch_truth_concordant_B)
    path(vcf_indexes)
    path(combined_concordant_variants_A)
    path(combined_concordant_variants_B)    

    output:
    path("batch_combined_concordant_variants_{A,B}_headed.tsv")    , emit: batch_concordant_combined_tsv
    path("batch_sample_concordant_{A,B}_merge.vcf.gz")             , emit: batch_sample_concordant_merged_vcf
    path("batch_sample_concordant_variants_{A,B}.tsv")             , emit: batch_sample_concordant_merged_tsv
    path("batch_truth_concordant_{A,B}_merge.vcf.gz")              , emit: batch_truth_concordant_merged_vcf
    path("*_correlation.pdf")                                      , emit: batch_correlation_pdfs
    path  "versions.yml"                                           , emit: versions

    script:
    def target = task.ext.args ?: ''
    def expected_data_source = task.ext.args2 ?: ''
    def observed_data_source = task.ext.args3 ?: ''
    def plot_depth = task.ext.args4 ?: ''
    def normalize_depth = task.ext.args5 ?: ''
    

    """
    ## merge input vcfs by catagory (truth vs test) and isec analysis type
    bcftools merge -m none $batch_sample_concordant_A -O z > batch_sample_concordant_A_merge.vcf.gz && \\
    bcftools query -H -f '%CHROM\\t%POS\\tend\\t%REF\\t%ALT\\t[%AF\\t%AD{1}\\t%DP\\t]\\n'  batch_sample_concordant_A_merge.vcf.gz  >  batch_sample_concordant_variants_A.tsv 
    bcftools merge -m none $batch_truth_concordant_A -O z > batch_truth_concordant_A_merge.vcf.gz
    bcftools merge -m none $batch_sample_concordant_B -O z > batch_sample_concordant_B_merge.vcf.gz && \\
    bcftools query -H -f '%CHROM\\t%POS\\tend\\t%REF\\t%ALT\\t[%AF\\t%AD{1}\\t%DP\\t]\\n'  batch_sample_concordant_B_merge.vcf.gz  >  batch_sample_concordant_variants_B.tsv
    bcftools merge -m none $batch_truth_concordant_B -O z > batch_truth_concordant_B_merge.vcf.gz

    ##generate header column variables 
    head="CHROM\\tPOS\\tend\\tREF\\tALT\\tAF\\tVD\\tDP\\tSAMPLE\\tFILTER\\texp_CHROM\\texp_POS\\texp_end\\texp_REF\\texp_ALT\\texp_AF\\texp_VD\\texp_DP\\texp_SAMPLE\\texp_FILTER\\texpt_norm_depth"

    
    ##concatenate input variant csv files and plot
    echo -e \$head | tee -a batch_combined_concordant_variants_A_headed.tsv > /dev/null && \\
    cat $combined_concordant_variants_A  | tee -a batch_combined_concordant_variants_A_headed.tsv > /dev/null
    echo -e \$head | tee -a batch_combined_concordant_variants_B_headed.tsv > /dev/null && \\
    cat $combined_concordant_variants_B | tee -a batch_combined_concordant_variants_B_headed.tsv > /dev/null && \\

    ##plot expected vs truth vafs/VDs/depths with regression line and peason correlation calculations

    if [[  $plot_depth == "yes" && $normalize_depth == "no" ]]; then
      correlation.R batch_combined_concordant_variants_A_headed.tsv Vardict VAF AF exp_AF ${expected_data_source}_A ${observed_data_source}_A batch 1 1
      correlation.R batch_combined_concordant_variants_A_headed.tsv Vardict Depth_at_variant_loci DP exp_DP ${expected_data_source}_A ${observed_data_source}_A batch 10000 10000
      correlation.R batch_combined_concordant_variants_B_headed.tsv Vardict VAF AF exp_AF ${expected_data_source}_B ${observed_data_source}_B batch 1 1
      correlation.R batch_combined_concordant_variants_B_headed.tsv Vardict Depth_at_variant_loci DP exp_DP ${expected_data_source}_B ${observed_data_source}_B batch 10000 10000
    elif [[  $plot_depth == "yes" && $normalize_depth == "yes" ]]; then
      correlation.R batch_combined_concordant_variants_A_headed.tsv Vardict VAF AF exp_AF ${expected_data_source}_A ${observed_data_source}_A batch 1 1
      correlation.R batch_combined_concordant_variants_A_headed.tsv Vardict Normalized_Depth_at_variant_loci DP expt_norm_depth ${expected_data_source}_A ${observed_data_source}_A batch 10000 10000
      correlation.R batch_combined_concordant_variants_B_headed.tsv Vardict VAF AF exp_AF ${expected_data_source}_B ${observed_data_source}_B batch 1 1
      correlation.R batch_combined_concordant_variants_B_headed.tsv Vardict Normalized_Depth_at_variant_loci DP expt_norm_depth ${expected_data_source}_B ${observed_data_source}_B batch 10000 10000
    elif [[  $plot_depth == "no" ]]; then
      correlation.R batch_combined_concordant_variants_A_headed.tsv Vardict VAF AF exp_AF ${expected_data_source}_A ${observed_data_source}_A batch 1 1
      correlation.R batch_combined_concordant_variants_B_headed.tsv Vardict VAF AF exp_AF ${expected_data_source}_B ${observed_data_source}_B batch 1 1
    fi
   
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools 2>&1 | grep Version | sed 's/^.*Version: //g' |  sed 's/ /_/g')
        R: \$(R --version 2>&1 | head -c 15)
        docker: \$(echo \$( grep 'docker run' .command.run 2>&1) )
    END_VERSIONS

    """

}
