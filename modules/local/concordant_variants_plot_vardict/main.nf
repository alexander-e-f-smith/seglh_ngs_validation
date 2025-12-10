process CONCORDANT_VARIANTS_PLOTTING_VARDICT {

    tag { "$meta" }
    label "plot_concordant"
    container "docker.io/seglh/bw_validation_tools@sha256:99fb766c096bd828c9cd0a9c88dba97eec3aab8f1f4c31a99c691096e98f8e08"

    input:
    tuple val(meta), path(sample_concordant_A), path(truth_concordant_A), path(sample_concordant_B), path(truth_concordant_B), path(vcf_indexes), val(norm_index)

    output:
    tuple val(meta), path("${meta}_*_concordant_variants_headed_{A,B}.tsv")
    tuple val(meta), path("${meta}*_correlation.pdf")
    tuple val(meta), path("${meta}_combined_concordant_variants_A.tsv")     , emit: combined_variants_A
    tuple val(meta), path("${meta}_combined_concordant_variants_B.tsv")     , emit: combined_variants_B
    path  "versions.yml"                                                    , emit: versions

    script:
    def target = task.ext.args ?: ''
    def expected_data_source = task.ext.args2 ?: ''
    def observed_data_source = task.ext.args3 ?: ''
    def plot_depth = task.ext.args4 ?: ''
    def normalize_depth = task.ext.args5 ?: ''
    def norm_factor = norm_index

    """
    #extract_essential_qc_and_depth_from_snappy2.sh combined_kch_qc
    head="CHROM\\tPOS\\tend\\tREF\\tALT\\tAF\\tVD\\tDP\\tSAMPLE\\tFILTER"
    headtruth="exp_CHROM\\texp_POS\\texp_end\\texp_REF\\texp_ALT\\texp_AF\\texp_VD\\texp_DP\\texp_SAMPLE\\texp_FILTER\\texpt_norm_depth"
    bcftools query  -f '%CHROM\\t%POS\\tend\\t%REF\\t%ALT\\t[%AF\\t%AD{1}\\t%DP\\t%SAMPLE]\\t%FILTER\\n'  ${sample_concordant_A} | sort -n -k1 -k2 -k4 -k5 > ${meta}_sample_concordant_variants_A.tsv
    bcftools query  -f '%CHROM\\t%POS\\tend\\t%REF\\t%ALT\\t[%AF\\t%AD{1}\\t%DP\\t%SAMPLE]\\t%FILTER\\n'  ${truth_concordant_A} | sort -n -k1 -k2 -k4 -k5 | \\
    awk -v norm=${norm_factor} '{FS="\\t"; OFS=FS} {print \$0, \$8*norm}' - > ${meta}_truth_concordant_variants_A.tsv && \\
    paste -d '\\t' ${meta}_sample_concordant_variants_A.tsv ${meta}_truth_concordant_variants_A.tsv > ${meta}_combined_concordant_variants_A.tsv && \\
    echo -e \$head | tee -a ${meta}_sample_concordant_variants_headed_A.tsv > /dev/null && \\
    echo -e \$headtruth | tee -a ${meta}_truth_concordant_variants_headed_A.tsv > /dev/null && \\
    cat ${meta}_sample_concordant_variants_A.tsv | tee -a ${meta}_sample_concordant_variants_headed_A.tsv > /dev/null        
    cat ${meta}_truth_concordant_variants_A.tsv |  tee -a ${meta}_truth_concordant_variants_headed_A.tsv > /dev/null && \\
    paste -d '\\t' ${meta}_sample_concordant_variants_headed_A.tsv ${meta}_truth_concordant_variants_headed_A.tsv > ${meta}_combined_concordant_variants_headed_A.tsv

    bcftools query  -f '%CHROM\\t%POS\\tend\\t%REF\\t%ALT\\t[%AF\\t%AD{1}\\t%DP\\t%SAMPLE]\\t%FILTER\\n'  ${sample_concordant_B} | sort -n -k1 -k2 -k4 -k5 > ${meta}_sample_concordant_variants_B.tsv
    bcftools query  -f '%CHROM\\t%POS\\tend\\t%REF\\t%ALT\\t[%AF\\t%AD{1}\\t%DP\\t%SAMPLE]\\t%FILTER\\n'  ${truth_concordant_B} | sort -n -k1 -k2 -k4 -k5 | \\
    awk -v norm=${norm_factor} '{FS="\\t"; OFS=FS} {print \$0, \$8*norm}' - > ${meta}_truth_concordant_variants_B.tsv && \\
    paste -d '\\t' ${meta}_sample_concordant_variants_B.tsv ${meta}_truth_concordant_variants_B.tsv > ${meta}_combined_concordant_variants_B.tsv && \\
    echo -e \$head | tee -a ${meta}_sample_concordant_variants_headed_B.tsv > /dev/null && \\
    echo -e \$headtruth | tee -a ${meta}_truth_concordant_variants_headed_B.tsv > /dev/null && \\
    cat ${meta}_sample_concordant_variants_B.tsv | tee -a ${meta}_sample_concordant_variants_headed_B.tsv > /dev/null
    cat ${meta}_truth_concordant_variants_B.tsv |  tee -a ${meta}_truth_concordant_variants_headed_B.tsv > /dev/null && \\
    paste -d '\\t' ${meta}_sample_concordant_variants_headed_B.tsv ${meta}_truth_concordant_variants_headed_B.tsv > ${meta}_combined_concordant_variants_headed_B.tsv && \\

    if [[  $plot_depth == "yes" && $normalize_depth == "no" ]]; then
      correlation.R ${meta}_combined_concordant_variants_headed_A.tsv Vardict VAF AF exp_AF ${expected_data_source}_A ${observed_data_source}_A ${meta} 1 1
      correlation.R ${meta}_combined_concordant_variants_headed_A.tsv Vardict Depth_at_variant_loci DP exp_DP ${expected_data_source}_A ${observed_data_source}_A ${meta} 10000 10000
      correlation.R ${meta}_combined_concordant_variants_headed_B.tsv Vardict VAF AF exp_AF ${expected_data_source}_B ${observed_data_source}_B ${meta} 1 1
      correlation.R ${meta}_combined_concordant_variants_headed_B.tsv Vardict Depth_at_variant_loci DP exp_DP ${expected_data_source}_B ${observed_data_source}_B ${meta} 10000 10000
    elif [[  $plot_depth == "yes" && $normalize_depth == "yes" ]]; then
      correlation.R ${meta}_combined_concordant_variants_headed_A.tsv Vardict VAF AF exp_AF ${expected_data_source}_A ${observed_data_source}_A ${meta} 1 1
      correlation.R ${meta}_combined_concordant_variants_headed_A.tsv Vardict Norm_Depth_at_variant_loci DP expt_norm_depth ${expected_data_source}_A ${observed_data_source}_A ${meta} 10000 10000
      correlation.R ${meta}_combined_concordant_variants_headed_B.tsv Vardict VAF AF exp_AF ${expected_data_source}_B ${observed_data_source}_B ${meta} 1 1
      correlation.R ${meta}_combined_concordant_variants_headed_B.tsv Vardict Norm_Depth_at_variant_loci DP expt_norm_depth ${expected_data_source}_B ${observed_data_source}_B ${meta} 10000 10000

    elif [[  $plot_depth == "no" ]]; then
      correlation.R ${meta}_combined_concordant_variants_headed_A.tsv Vardict VAF AF exp_AF ${expected_data_source}_A ${observed_data_source}_A ${meta} 1 1
      correlation.R ${meta}_combined_concordant_variants_headed_B.tsv Vardict VAF AF exp_AF ${expected_data_source}_B ${observed_data_source}_B ${meta} 1 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools 2>&1 | grep Version | sed 's/^.*Version: //g' |  sed 's/ /_/g')
        R: \$(R --version 2>&1 | head -c 15)
        docker: \$(echo \$( grep 'docker run' .command.run 2>&1) )
    END_VERSIONS

    """

}
