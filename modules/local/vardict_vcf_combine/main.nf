process MERGE_VCFS_PLOTTING {

    tag { "combine_vcf" }
    label "process_combine_vcf"
    container "docker.io/seglh/snappy_python3_ngstools:1.0.2@sha256:c762d2ca67e52068b7bf9ef4ba4ff4466a44f8b3eead2ba5b91dd29ca7193bd1"
    debug true

    input:
    path(vcf_files)
    path(truth_vcf_files)
    path(bedfile)

    output:
    path("combined.out")

    script:
    coverage_settings = task.ext.args ?: ''
    

    """
    #extract_essential_qc_and_depth_from_snappy2.sh combined_kch_qc
    vcf_list=\$(ls $vcf_files)
    echo \$vcf_list > test2.out
    head="CHROM\\tPOS\\tend\\tREF\\tALT\\tAF\\tVD\\tDP\\tSAMPLE]\\tFILTER"
    for file in $vcf_files
    do
        id=\$(basename \$file .vcf.gz)
        echo \$id | tee -a combined.out
        tabix \${file}
        bcftools query  -f '%CHROM\\t%POS\\tend\\t%REF\\t%ALT\\t[%AF\\t%VD\\t%DP\\t%SAMPLE]\\t%FILTER\\n' -R $bedfile  \${file} --output \${id}_vcf_extracted_data.tsv
    done
    echo -e \$head | tee -a headed_file 
    cat *_vcf_extracted_data.tsv | tee -a headed_file        
    #cat <<-END_VERSIONS > versions.yml
    #"${task.process}":
    #    bcftools: \$(bcftools 2>&1 | grep Version | sed 's/^.*Version: //g' |  sed 's/ /_/g')
    #END_VERSIONS

    """
    // #bcftools query  -f '%CHROM\\t%POS\\tend\\t%REF\\t%ALT\\t[%AF\\t%VD\\t%SAMPLE]\\t%FILTER\\n' -R $bedfile  ${file} --output ${ID}_

}
