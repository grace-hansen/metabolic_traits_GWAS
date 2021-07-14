#!/bin/bash

trait=$1
tmp=$2

#For each sex, download v2 version if present, non-v2 version otherwise

cd $tmp

for sex in both_sexes female male; do
	prefix=${trait}.${sex}
	wget=$(awk -F '\t' '$5=="'$prefix'.v2.tsv.bgz" {print $6}' /project2/nobrega/grace/metabolic_traits_GWAS/UKBB_GWAS_Imputed_v3_File_Manifest_Release_20180731.txt)
	if [[ $wget != "" ]]; then 
		echo "not empty"
		$wget
		mv $prefix.v2.tsv.bgz $trait.$sex.tsv.bgz
	else
		wget=$(awk -F '\t' '$5=="'$prefix'.tsv.bgz" {print $6}' /project2/nobrega/grace/metabolic_traits_GWAS/UKBB_GWAS_Imputed_v3_File_Manifest_Release_20180731.txt)
		$wget
	fi	
done
