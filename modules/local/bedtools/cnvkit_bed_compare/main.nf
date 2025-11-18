process COMPARE_CNVKIT_BED {

    tag { "$meta" }
    label "cnv_bed"
    container "docker.io/seglh/bw_validation_tools@sha256:99fb766c096bd828c9cd0a9c88dba97eec3aab8f1f4c31a99c691096e98f8e08"

    input:
    tuple val(meta), path(sample_cnv_bed), path(truth_cnv_bed)
    path(regions_bed)

    output:
    tuple val(meta), path("${meta}_cnv_bed_compare_{A,B}.tsv")
    path  "versions.yml"                                             ,  emit: versions

    script:
    def target = task.ext.args ?: ''
    def expected_data_source = task.ext.args2 ?: ''
    

    """
    headerA="chr_test\\tstart_test\\tend_test\\tsample_test\\tcn_test\\tchr_truth\\tstart_truth\\tend_truth\\tsample_truth\\tcn_truth\\tbcftools_overlap\\tsize_test\\tsize_truth\\ttest_loci_size/overlap\\ttruth_loci_size/overlap"
    headerB="chr_truth\\tstart_truth\\tend_truth\\tsample_truth\\tcn_truth\\tchr_test\\tstart_test\\tend_test\\tsample_test\\tcn_test\\tbcftools_overlap\\tsize_truth\\tsize_test\\ttruth_loci_size/overlap\\ttest_loci_size/overlap"
    echo -e \$headerA | tee -a ${meta}_cnv_bed_compare_A.tsv > /dev/null &&
    bedtools intersect -wa -a $sample_cnv_bed -b $regions_bed | bedtools intersect -wao  -a - -b  $truth_cnv_bed | awk 'BEGIN{OFS="\\t"} {print \$0,\$3-\$2,\$8-\$7}' - | awk 'BEGIN{OFS="\\t"} {print \$0,(\$12/\$11),(\$13/\$11)}' \\
    | tee -a ${meta}_cnv_bed_compare_A.tsv > /dev/null
    echo -e \$headerB | tee -a ${meta}_cnv_bed_compare_B.tsv > /dev/null &&
    bedtools intersect -wa -a $truth_cnv_bed -b $regions_bed | bedtools intersect -wao  -a - -b  $sample_cnv_bed | awk 'BEGIN{OFS="\\t"} {print \$0,\$3-\$2,\$8-\$7}' - | awk 'BEGIN{OFS="\\t"} {print \$0,(\$12/\$11)/1,(\$13/\$11)}' \\
    | tee -a ${meta}_cnv_bed_compare_B.tsv > /dev/null
        

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        docker: \$(echo \$( grep 'docker run' .command.run 2>&1) )
    END_VERSIONS

    """

}
