#/bin/bash
input_json=$1
input_cov=$2
name=$3
outname=$4

basename1=$(basename $input_json .json)
header1="Sample_name\tTotal_paired_reads_(fastqc)\t>Q30_base_quality_proportion\tTotal_processed_paired_reads\tDuplication_rate\tFailed_primers\tOn_target_molecules_proportion\tMEDIAN_INSERT_SIZE\tMEAN_INSERT_SIZE\t"
echo -e ${header1} | tee -a ${name}_qc
jq1=$(jq '.[]  | select( .type |contains("fastqc")) | .data."Basic Statistics".contents'  $input_json | jq -s '.[0]."Total Sequences"')
jq2=$(jq '.[]  | select( .type |contains("base_quality_metrics")) | .data.metrics.contents.Q30' $input_json)
jq3=$(jq '.[]  | select( .type |contains("alignment_summary_metrics")) | .data.metrics.contents | .[0].TOTAL_READS' $input_json)
jq4=$(jq '.[]  | select( .type |contains("DUPLICATION")) | .data.metrics.contents.DUPLICATION ' $input_json)
jq5=$(jq '.[]  | select( .type |contains("DUPLICATION")) | .data.metrics.contents.FAILED_PRIMERS ' $input_json)
jq6=$(jq '.[]  | select( .type |contains("DUPLICATION")) | .data.metrics.contents.ON_TARGET_MOLECULES_FRAC ' $input_json)
jq7=$(jq '.[]  | select( .type |contains("insert_size_metrics")) | .data.metrics.contents.MEDIAN_INSERT_SIZE' $input_json)
jq8=$(jq '.[]  | select( .type |contains("insert_size_metrics")) | .data.metrics.contents.MEAN_INSERT_SIZE'  $input_json)

echo -e ${basename1}'\t'${jq1}'\t'${jq2}'\t'${jq3}'\t'${jq4}'\t'${jq5}'\t'${jq6}'\t'${jq7}'\t'${jq8}'\t' | tee -a ${name}_qc


# extract coverage for snappy exon coverage file
basename2=$(basename $input_cov _exoncoverage_metrics)
headercov=$(grep 'EXON_LENGTH' $input_cov |  awk 'BEGIN{FS="\t"; OFS="\t"} {print $6,$7,$8,$9,$10,$11,"sample"}')
echo -e ${headercov} | tee -a ${name}_coverage
data=$(grep  'EXON_LENGTH' -A 1   $input_cov | grep -v 'EXON_LENGTH'| awk -v var="$basename2" 'BEGIN{FS="\t"; OFS="\t"} {print $6,$7,$8,$9,$10,$11,var}') 
echo -e $data | tee -a ${name}_coverage

#head -1 ${outname}_qc_compendium > head_cov
#grep -v ^.*EXON_LENGTH ${outname}_catall_out_coverage  cat head_cov - > ${outname}_out_coverage_formatted
#rm ${outname}_catall_out_coverage
#rm head_cov

#paste different metrics outputs together
paste -d "\t" ${name}_qc ${name}_coverage > ${outname}
