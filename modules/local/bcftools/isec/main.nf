process ISEC_VALIDATION {

    tag { "$meta" }
    label 'process_single'
    container "docker.io/seglh/snappy_python3_ngstools:1.0.2@sha256:c762d2ca67e52068b7bf9ef4ba4ff4466a44f8b3eead2ba5b91dd29ca7193bd1"
    debug true

    input:
    tuple val(meta), path(sample_vcf), path(truth_vcf), path(variant_vcf)
    path(bed)

    output:
    tuple val(meta), path("isec_out_{A,B}_${meta}/${meta}*{0000,0002,0003}.vcf.gz.tbi")    ,  emit: isec_vcf_indexes
    tuple val(meta), path("isec_out_A_${meta}/${meta}_A_0000.vcf.gz")                      ,  emit: isec_unique_sample_variants
    tuple val(meta), path("isec_out_B_${meta}/${meta}_B_0000.vcf.gz")                      ,  emit: isec_unique_truth_variants
    tuple val(meta), path("isec_out_B_${meta}/${meta}_B_0003.vcf.gz")                      ,  emit: isec_concordant_sample
    tuple val(meta), path("isec_out_B_${meta}/${meta}_B_0002.vcf.gz")                      ,  emit: isec_concordant_truth
    tuple val(meta), path("${meta}_variant_comparison_stats")                              ,  emit: isec_variant_comparison_stats
    path  "versions.yml"                                                                   ,  emit: versions
    

    script:
    def filter1_vcf = task.ext.args ?: ''
    def filter2_vcf = task.ext.args2 ?: ''
    def output1 = "isec_out_A_${meta}"
    def output2 = "isec_out_B_${meta}"

    """
    tabix ${sample_vcf} 
    tabix ${truth_vcf} \\
        && bcftools isec ${filter1_vcf} ${filter2_vcf} -T ${bed} ${sample_vcf} ${truth_vcf} -p ${output1} \\
        && bcftools isec ${filter1_vcf} ${filter2_vcf} -T ${bed} ${truth_vcf} ${sample_vcf} -p ${output2} \\
        && cd ${output1} && mv 0000.vcf ${meta}_A_0000.vcf && mv 0001.vcf ${meta}_A_0001.vcf \\
        && bgzip ${meta}_A_0000.vcf && tabix ${meta}_A_0000.vcf.gz && bgzip ${meta}_A_0001.vcf && tabix ${meta}_A_0001.vcf.gz  \\
        && no_unique_sample_variants=\$(bcftools stats ${meta}_A_0000.vcf.gz  | grep -v "^#" | grep "number of records" | awk 'BEGIN{OFS="\\t"} {print \$6}') &&  cd ../ \\
        && cd ${output2} && mv 0000.vcf ${meta}_B_0000.vcf && mv 0001.vcf ${meta}_B_0001.vcf \\
        && mv 0003.vcf ${meta}_B_0003.vcf && mv 0002.vcf ${meta}_B_0002.vcf \\
        && bgzip ${meta}_B_0000.vcf && tabix ${meta}_B_0000.vcf.gz && bgzip ${meta}_B_0001.vcf && tabix ${meta}_B_0001.vcf.gz \\
        && bgzip ${meta}_B_0002.vcf && tabix ${meta}_B_0002.vcf.gz && bgzip ${meta}_B_0003.vcf && tabix ${meta}_B_0003.vcf.gz \\
        && no_unique_truth_variants=\$(bcftools stats ${meta}_B_0000.vcf.gz  | grep -v "^#" | grep "number of records" | awk  'BEGIN{OFS="\\t"} {print \$6}') \\
        && no_concordant_variants=\$(bcftools stats ${meta}_B_0003.vcf.gz  | grep -v "^#" | grep "number of records" | awk 'BEGIN{OFS="\\t"} {print \$6}') && cd ../ \\
        && printf "Sample\\tNumber_unique_sample_variants\\tNumber_unique_truth_variants\\tNumber_concordant_variants\\n" | tee -a ${meta}_variant_comparison_stats \\
        && printf  "${meta}\\t\$no_unique_sample_variants\\t\$no_unique_truth_variants\\t\$no_concordant_variants\\n" | tee -a  ${meta}_variant_comparison_stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools 2>&1 | grep Version | sed 's/^.*Version: //g' |  sed 's/ /_/g')
    END_VERSIONS
    """
}
