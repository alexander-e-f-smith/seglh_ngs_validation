#/bin/bash
input_json=$1
input_cov=$2
outname=$3
type=$4

basename1=$(basename $input_json .json)
header1="${type}_Sample_name\t${type}_Total_paired_reads_(fastqc)\t${type}_>Q30_base_quality_proportion\t${type}_Total_processed_paired_reads\t${type}_Duplication_rate\t${type}_Failed_primers\t${type}_On_target_molecules_proportion\t${type}_MEDIAN_INSERT_SIZE\t${type}_MEAN_INSERT_SIZE\t"
echo -e ${header1} | tee -a ${basename1}_qc > /dev/null
jq1=$(jq '.[]  | select( .type |contains("fastqc")) | .data."Basic Statistics".contents'  $input_json | jq -s '.[0]."Total Sequences"')
jq2=$(jq '.[]  | select( .type |contains("base_quality_metrics")) | .data.metrics.contents.Q30' $input_json)
jq3=$(jq '.[]  | select( .type |contains("alignment_summary_metrics")) | .data.metrics.contents | .[0].TOTAL_READS' $input_json)
jq4=$(jq '.[]  | select( .type |contains("DUPLICATION")) | .data.metrics.contents.DUPLICATION ' $input_json)
jq5=$(jq '.[]  | select( .type |contains("DUPLICATION")) | .data.metrics.contents.FAILED_PRIMERS ' $input_json)
jq6=$(jq '.[]  | select( .type |contains("DUPLICATION")) | .data.metrics.contents.ON_TARGET_MOLECULES_FRAC ' $input_json)
jq7=$(jq '.[]  | select( .type |contains("insert_size_metrics")) | .data.metrics.contents.MEDIAN_INSERT_SIZE' $input_json)
jq8=$(jq '.[]  | select( .type |contains("insert_size_metrics")) | .data.metrics.contents.MEAN_INSERT_SIZE'  $input_json)

echo -e ${basename1}'\t'${jq1}'\t'${jq2}'\t'${jq3}'\t'${jq4}'\t'${jq5}'\t'${jq6}'\t'${jq7}'\t'${jq8}'\t' | tee -a ${basename1}_qc > /dev/null


# extract coverage for snappy exon coverage file
basename2=$(basename $input_cov .exoncoverage)
headercov=$(grep 'EXON_LENGTH' $input_cov |  awk -v var="${type}_" 'BEGIN{FS="\t"; OFS="\t"} {print var$6,var$7,var$8,var$9,var$10,var$11,var"sample"}')
echo -e "${headercov}" | tee -a ${basename2}_coverage > /dev/null
data=$(grep  'EXON_LENGTH' -A 1   $input_cov | grep -v 'EXON_LENGTH'| awk -v var="$basename2" 'BEGIN{FS="\t"; OFS="\t"} {print $6,$7,$8,$9,$10,$11,var}')
echo -e "$data" | tee -a ${basename2}_coverage > /dev/null

#head -1 ${outname}_qc_compendium > head_cov
#grep -v ^.*EXON_LENGTH ${outname}_catall_out_coverage  cat head_cov - > ${outname}_out_coverage_formatted
#rm ${outname}_catall_out_coverage
#rm head_cov

#paste different metrics outputs together
paste -d "\t" ${basename1}_qc ${basename2}_coverage > ${outname}
