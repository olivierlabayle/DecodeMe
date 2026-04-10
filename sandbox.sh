dx download file-J1GbX10J1ZJ6PG0y65q642Y8
dx download 
dx download project-Ggv5Q00JFKgG9P4Z7Fv71bxz:/ukb-imputation/impute5/ref_xcf_phased/chrX/rp_chrX_1_27440787.bcf
dx download project-Ggv5Q00JFKgG9P4Z7Fv71bxz:/ukb-imputation/impute5/ref_xcf_phased/chrX/rp_chrX_1_27440787.bcf.csi
dx download project-Ggv5Q00JFKgG9P4Z7Fv71bxz:/ukb-imputation/impute5/ref_xcf_phased/chrX/rp_chrX_24026048_50656185.bcf
dx download project-Ggv5Q00JFKgG9P4Z7Fv71bxz:/ukb-imputation/impute5/ref_xcf_phased/chrX/rp_chrX_24026048_50656185.bcf.csi
dx download project-Ggv5Q00JFKgG9P4Z7Fv71bxz:/ukb-imputation/impute5/ref_xcf_phased/chrX/rp_chrX_47743520_88764754.bcf
dx download project-Ggv5Q00JFKgG9P4Z7Fv71bxz:/ukb-imputation/impute5/ref_xcf_phased/chrX/rp_chrX_47743520_88764754.bcf.csi
dx download project-Ggv5Q00JFKgG9P4Z7Fv71bxz:/ukb-imputation/impute5/ref_xcf_phased/chrX/rp_chrX_87221448_115258184.bcf
dx download project-Ggv5Q00JFKgG9P4Z7Fv71bxz:/ukb-imputation/impute5/ref_xcf_phased/chrX/rp_chrX_87221448_115258184.bcf.csi
dx download project-Ggv5Q00JFKgG9P4Z7Fv71bxz:/ukb-imputation/impute5/ref_xcf_phased/chrX/rp_chrX_112836207_139441710.bcf
dx download project-Ggv5Q00JFKgG9P4Z7Fv71bxz:/ukb-imputation/impute5/ref_xcf_phased/chrX/rp_chrX_112836207_139441710.bcf.csi
dx download project-Ggv5Q00JFKgG9P4Z7Fv71bxz:/ukb-imputation/impute5/ref_xcf_phased/chrX/rp_chrX_137290831_1000000000.bcf
dx download project-Ggv5Q00JFKgG9P4Z7Fv71bxz:/ukb-imputation/impute5/ref_xcf_phased/chrX/rp_chrX_137290831_1000000000.bcf.csi


docker run -it --rm -v $PWD:/mnt/data olivierlabayle/genomicc:0.3.0 /bin/bash


# Make reference statistics for QC

for bcf_file in rp_chrX_1_27440787.bcf rp_chrX_24026048_50656185.bcf rp_chrX_47743520_88764754.bcf rp_chrX_87221448_115258184.bcf rp_chrX_112836207_139441710.bcf rp_chrX_137290831_1000000000.bcf; do
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n' ${bcf_file} > ${bcf_file%.bcf}.tsv
    echo "${bcf_file%.bcf}.tsv"
done > stats_files.txt

julia --project=/opt/genomicc-workflows/ --startup-file=no --sysimage=/opt/genomicc-workflows/GenomiccWorkflows.so --threads=auto \


bed_prefix=decodeme_chrx_batch1_preimputation
ukb_prefix=ukb22418_cX_b0_v2
sample_file=ukb_control_sample_ids.txt

# QC DecodeMe file
julia --project=/opt/genomicc-workflows/ --startup-file=no --sysimage=/opt/genomicc-workflows/GenomiccWorkflows.so --threads=auto \
DecodeMe/scripts/format_bim_with_ref.jl ${bed_prefix}.bim ${ukb_prefix}.bim variants_to_drop.txt

plink \
--nonfounders \
--bfile ${bed_prefix} \
--hwe 1e-12 \
--exclude variants_to_drop.txt \
--freq \
--make-bed \
--out ${bed_prefix}.hwe


# QC UKB

julia --project=/opt/genomicc-workflows/ --startup-file=no --sysimage=/opt/genomicc-workflows/GenomiccWorkflows.so --threads=auto \
DecodeMe/scripts/format_variant_id.jl ${ukb_prefix}.bim

plink \
--bfile ${ukb_prefix} \
--hwe 1e-9 \
--freq \
--make-bed \
--keep ${sample_file} \
--out ${ukb_prefix}.hwe

# Merge

for f in decodeme_chrx_batch1_preimputation.hwe.bed; do
    echo "${f%.bed}"
done > merge_list.txt
echo "${ukb_prefix}.hwe" >> merge_list.txt

julia --project=/opt/genomicc-workflows/ --startup-file=no --sysimage=/opt/genomicc-workflows/GenomiccWorkflows.so --threads=auto \
DecodeMe/scripts/make_shared_variants_list.jl merge_list.txt shared_variants.txt

plink \
    --biallelic-only \
    --merge-list merge_list.txt \
    --extract shared_variants.txt \
    --make-bed \
    --out chrX

plink2 \
--bfile chrX \
--mind 0.03 \
--geno 0.02 \
--make-bed \
--out chrX.qced