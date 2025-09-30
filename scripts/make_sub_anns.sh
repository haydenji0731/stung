#!/usr/bin/env bash

mat_fn="/ccb/salz7-data/ftp.ccb/pub/data/hg002-q100/v0.5/hg002.v1.1.loff.v0.5.mat.gff"
pat_fn="/ccb/salz7-data/ftp.ccb/pub/data/hg002-q100/v0.5/hg002.v1.1.loff.v0.5.pat.gff"
out_dir="/ccb/salz4-3/hji20/stung/results/sub_anns"

for i in {1..22}; do
	echo $i
	chr="chr${i}_MATERNAL"
	awk -F'\t' -v chr="$chr" '$1 == chr { print $0 }' "$mat_fn" > "${out_dir}/chr${i}_mat.gff"
	chr="chr${i}_PATERNAL"
	awk -F'\t' -v chr="$chr" '$1 == chr { print $0 }' "$pat_fn" > "${out_dir}/chr${i}_pat.gff"
done