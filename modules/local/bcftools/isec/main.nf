process ISEC_VALIDATION {

    tag { "$meta" }
    label 'process_single'
    container "docker.io/seglh/snappy_python3_ngstools:1.0.2@sha256:c762d2ca67e52068b7bf9ef4ba4ff4466a44f8b3eead2ba5b91dd29ca7193bd1"
    debug true

    input:
    tuple val(meta), path(sample_vcf), path(truth_vcf), path(variant_vcf)
    path(bed)

    output:
    tuple val(meta), path("*.tbi")                                                  ,  emit: vcf_index
    tuple val(meta), path("isec_out_{A,B}_${meta}/{${meta}_0000,${meta}_0001}.vcf")     ,  emit: isec_unique_variants
    tuple val(meta), path("isec_out_B_${meta}/${meta}_0003.vcf")                    ,  emit: isec_concordant_sample
    tuple val(meta), path("isec_out_B_${meta}/${meta}_0002.vcf")                    ,  emit: isec_concordant_truth
    path  "versions.yml"                                                            ,  emit: versions
    

    script:
    def filter1_vcf = task.ext.args ?: ''
    def filter2_vcf = task.ext.args2 ?: ''
    def output1 = "isec_out_A_${meta}"
    def output2 = "isec_out_B_${meta}"

    """
    tabix ${sample_vcf} & tabix ${truth_vcf} \\
        && bcftools isec ${filter1_vcf} ${filter2_vcf} -T ${bed} ${sample_vcf} ${truth_vcf} -p ${output1} \\
        && bcftools isec ${filter1_vcf} ${filter2_vcf} -T ${bed} ${truth_vcf} ${sample_vcf} -p ${output2} \\
        && cd ${output1} && mv 0000.vcf ${meta}_0000.vcf && mv 0001.vcf ${meta}_0001.vcf && cd ../ \\
        && cd ${output2} && mv 0000.vcf ${meta}_0000.vcf && mv 0001.vcf ${meta}_0001.vcf \\
        && mv 0003.vcf ${meta}_0003.vcf && mv 0002.vcf ${meta}_0002.vcf && cd ../ \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools 2>&1 | grep Version | sed 's/^.*Version: //g' |  sed 's/ /_/g')
    END_VERSIONS
    """
}
