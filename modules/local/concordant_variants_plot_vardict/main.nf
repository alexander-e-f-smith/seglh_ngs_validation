process CONCORDANT_VARIANTS_PLOTTING_VARDICT {

    tag { "$meta" }
    label "plot_concordant"
    container "docker.io/seglh/alex_validationtools:1.0.0@sha256:bec3658f87699b978eb8b2009f2a42f08cf328913bd58f58fa09b1a080890881"

    input:
    tuple val(meta), path(sample_concordant_A)
    tuple val(meta), path(truth_concordant_A)
    tuple val(meta), path(sample_concordant_B)
    tuple val(meta), path(truth_concordant_B)
    tuple val(meta), path(vcf_indexes)

    output:
    tuple val(meta), path("${meta}_*_concordant_variants_headed_{A,B}.tsv")
    tuple val(meta), path("${meta}*_correlation.pdf")

    script:
    def target = task.ext.args ?: ''
    def expected_data_source = task.ext.args2 ?: ''
    def observed_data_source = task.ext.args3 ?: ''
    

    """
    #extract_essential_qc_and_depth_from_snappy2.sh combined_kch_qc
    head="CHROM\\tPOS\\tend\\tREF\\tALT\\tAF\\tVD\\tDP\\tSAMPLE]\\tFILTER"
    headtruth="exp_CHROM\\texp_POS\\texp_end\\texp_REF\\texp_ALT\\texp_AF\\texp_VD\\texp_DP\\texp_SAMPLE\\texp_FILTER" 
    bcftools query  -f '%CHROM\\t%POS\\tend\\t%REF\\t%ALT\\t[%AF\\t%VD\\t%DP\\t%SAMPLE]\\t%FILTER\\n'  ${sample_concordant_A} --output ${meta}_sample_concordant_variants_A.tsv
    bcftools query  -f '%CHROM\\t%POS\\tend\\t%REF\\t%ALT\\t[%AF\\t%VD\\t%DP\\t%SAMPLE]\\t%FILTER\\n'  ${truth_concordant_A} --output ${meta}_truth_concordant_variants_A.tsv && \\
    echo -e \$head | tee -a ${meta}_sample_concordant_variants_headed_A.tsv > /dev/null && \\
    echo -e \$headtruth | tee -a ${meta}_truth_concordant_variants_headed_A.tsv > /dev/null && \\
    sort -n -k1 -k2 ${meta}_sample_concordant_variants_A.tsv | tee -a ${meta}_sample_concordant_variants_headed_A.tsv > /dev/null        
    sort -n -k1 -k2 ${meta}_truth_concordant_variants_A.tsv |  tee -a ${meta}_truth_concordant_variants_headed_A.tsv > /dev/null && \\
    paste -d '\\t' ${meta}_sample_concordant_variants_headed_A.tsv ${meta}_truth_concordant_variants_headed_A.tsv > ${meta}_combined_concordant_variants_headed_A.tsv \\

    bcftools query  -f '%CHROM\\t%POS\\tend\\t%REF\\t%ALT\\t[%AF\\t%VD\\t%DP\\t%SAMPLE]\\t%FILTER\\n'  ${sample_concordant_B} --output ${meta}_sample_concordant_variants_B.tsv
    bcftools query  -f '%CHROM\\t%POS\\tend\\t%REF\\t%ALT\\t[%AF\\t%VD\\t%DP\\t%SAMPLE]\\t%FILTER\\n'  ${truth_concordant_B} --output ${meta}_truth_concordant_variants_B.tsv && \\
    echo -e \$head | tee -a ${meta}_sample_concordant_variants_headed_B.tsv > /dev/null && \\
    echo -e \$headtruth | tee -a ${meta}_truth_concordant_variants_headed_B.tsv > /dev/null && \\
    sort -n -k1 -k2 ${meta}_sample_concordant_variants_B.tsv | tee -a ${meta}_sample_concordant_variants_headed_B.tsv > /dev/null
    sort -n -k1 -k2 ${meta}_truth_concordant_variants_B.tsv |  tee -a ${meta}_truth_concordant_variants_headed_B.tsv > /dev/null && \\
    paste -d '\\t' ${meta}_sample_concordant_variants_headed_B.tsv ${meta}_truth_concordant_variants_headed_B.tsv > ${meta}_combined_concordant_variants_headed_B.tsv && \\

    correlation.R ${meta}_combined_concordant_variants_headed_A.tsv Vardict VAF AF exp_AF ${expected_data_source}_A ${observed_data_source}_A ${meta} 1 1
    correlation.R ${meta}_combined_concordant_variants_headed_A.tsv Vardict VD VD exp_VD ${expected_data_source}_A ${observed_data_source}_A ${meta} 5000 5000
    correlation.R ${meta}_combined_concordant_variants_headed_B.tsv Vardict VAF AF exp_AF ${expected_data_source}_B ${observed_data_source}_B ${meta} 1 1
    correlation.R ${meta}_combined_concordant_variants_headed_B.tsv Vardict VD VD exp_VD ${expected_data_source}_B ${observed_data_source}_B ${meta} 5000 5000
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools 2>&1 | grep Version | sed 's/^.*Version: //g' |  sed 's/ /_/g')
    END_VERSIONS

    """

}
