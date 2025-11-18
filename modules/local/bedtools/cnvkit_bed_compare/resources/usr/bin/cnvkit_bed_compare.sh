#/bin/bash

read -p "list of cnvkit bed files, space separated ( and PATH) with comparison id in third position/column: " list
read -p "bed regions file for interogation ( and PATH): " bed

read -p "experiment and output dir name prefix: " expt


mkdir ${expt}_dir

OLDIFS=$IFS; IFS=$'\n'; for bed_pair in $(cat $list)
do
   echo $bed_pair
   tagbase=$(echo $bed_pair | awk '{print $1"_vs"$2}' -)
   pair_tag=$( echo ${tagbase} | sed -e 's/.bed//g' -)
   echo $pair_tag
   sample1=$(echo $bed_pair | awk '{print $1}' -)
   sample2=$(echo $bed_pair | awk '{print $2}' -)
   comparison_no=$(echo $bed_pair | awk '{print $3}' -)
   bedtools intersect -wa -a $sample1 -b $bed | bedtools intersect -wao  -a - -b  $sample2 > ${comparison_no}_${pair_tag}_compareA
   bedtools intersect -wa -a $sample2 -b $bed | bedtools intersect -wao  -a - -b  $sample1  > ${comparison_no}_${pair_tag}_compareB
done
