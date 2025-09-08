process BATCH_SENSITIVITY_SPECIFICITY {

    tag { "batch_sen_spec" }
    label "batch_sen_spec"
    container "docker.io/seglh/bw_validation_tools@sha256:99fb766c096bd828c9cd0a9c88dba97eec3aab8f1f4c31a99c691096e98f8e08"

    input:
    path(batch_sample_unique)
    path(batch_truth_unique)
    path(batch_variant_stats)
    path(vcf_indexes)

    output:
    path("batch_*_unique_variants.tsv")
    path("batch_*_unique_variants.vcf.gz*")
    path("batch_*_unique_variants_only_samples_with_variants.vcf.gz*") 
    path("batch_variant_compartison_stats")

    script:
    def target = task.ext.args ?: ''
    def expected_data_source = task.ext.args2 ?: ''
    def observed_data_source = task.ext.args3 ?: ''
    

    """
    bcftools merge -m none -i DP:join,VD:join,AF:join,MQ:join,QUAL:join,NM:join,MSI:join,MSILEN:join,SN:join,PMEAN:join $batch_sample_unique -O z > batch_exp_unique_variants.vcf.gz && \\
    bcftools merge -m none -i DP:join,VD:join,AF:join,MQ:join,QUAL:join,NM:join,MSI:join,MSILEN:join,SN:join,PMEAN:join $batch_truth_unique -O z > batch_truth_unique_variants.vcf.gz && \\
    tabix batch_exp_unique_variants.vcf.gz && tabix batch_truth_unique_variants.vcf.gz  && \\
    bcftools query -H  -f  '%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t[\\t%DP][\\t%AF][\\t%VD]\\t%INFO/MSI\\t%INFO/MSILEN\\t%INFO/QUAL\\t%INFO/NM\\t%INFO/MQ\\t%INFO/SN\\t%INFO/PMEAN\\t%INFO/LSEQ\\t%INFO/RSEQ\\n' batch_exp_unique_variants.vcf.gz > batch_exp_unique_variants.tsv && \\
    bcftools query -H  -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t[\\t%DP][\\t%AF][\\t%VD]\\t%INFO/MSI\\t%INFO/MSILEN\\t%INFO/QUAL\\t%INFO/NM\\t%INFO/MQ\\t%INFO/SN\\t%INFO/PMEAN\\t%INFO/LSEQ\\t%INFO/RSEQ\\n' batch_truth_unique_variants.vcf.gz > batch_truth_unique_variants.tsv
    cat $batch_variant_stats | sed '2,\${/^Sample/d}' - > batch_variant_compartison_stats
    
    touch expt_sample_with_unique_variants
    touch truth_sample_with_unique_variants
    for line in $batch_sample_unique
    do
      no_unique_sample_variants=\$(bcftools stats \${line}  | grep -v "^#" | grep "number of records" | awk  'BEGIN{OFS="\\t"} {print \$6}')
      if [[ "\$no_unique_sample_variants" -gt 0 ]]; then
        echo \$line | tee -a expt_sample_with_unique_variants
      fi
    done

    for line in $batch_truth_unique
    do
      no_unique_truth_variants=\$(bcftools stats \${line}  | grep -v "^#" | grep "number of records" | awk  'BEGIN{OFS="\\t"} {print \$6}')
      if [[ "\$no_unique_truth_variants" -gt 0 ]]; then
        echo \$line | tee -a truth_sample_with_unique_variants
      fi
    done
        
    bcftools merge -m none -i DP:join,VD:join,AF:join,MQ:join,QUAL:join,NM:join,MSI:join,MSILEN:join,SN:join,PMEAN:join -l expt_sample_with_unique_variants -O z > batch_expt_unique_variants_only_samples_with_variants.vcf.gz || true
    bcftools merge -m none -i DP:join,VD:join,AF:join,MQ:join,QUAL:join,NM:join,MSI:join,MSILEN:join,SN:join,PMEAN:join -l truth_sample_with_unique_variants -O z > batch_truth_unique_variants_only_samples_with_variants.vcf.gz  || true  && \\
    tabix batch_truth_unique_variants_only_samples_with_variants.vcf.gz || true && tabix batch_expt_unique_variants_only_samples_with_variants.vcf.gz || true
    
    

    #cat <<-END_VERSIONS > versions.yml
    #"${task.process}":
    #    bcftools: \$(bcftools 2>&1 | grep Version | sed 's/^.*Version: //g' |  sed 's/ /_/g')
    #END_VERSIONS

    """

}
