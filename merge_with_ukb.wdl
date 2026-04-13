version 1.0

import "structs.wdl"
import "utils.wdl"

workflow merge_with_ukb {
    input {
        String docker_image = "olivierlabayle/decodeme:main"

        Array[File] reference_panel_files
        Array[PLINKFileset] decodeme_genotypes
        PLINKFileset ukb_genotypes
        File ukb_samples

        String julia_use_sysimage = "true"
        String julia_threads = "auto"
    }

    # Get generic Julia command
    call utils.get_julia_cmd {
        input:
            use_sysimage = julia_use_sysimage,
            threads = julia_threads
    }

    call make_ref_panel_stats {
        input:
            docker_image = docker_image,
            julia_cmd = get_julia_cmd.julia_cmd,
            bcf_csi_files = reference_panel_files
    }

    # QC and Liftover DecodeMe batches
    scatter (decodeme_batch in decodeme_genotypes) {
        call QCAndLiftOverArray as QCAndLiftOverDecodeMe {
            input:
                docker_image = docker_image,
                chr = decodeme_batch.chr,
                bed_file = decodeme_batch.bed,
                bim_file = decodeme_batch.bim,
                fam_file = decodeme_batch.fam,
                hwe_pval = "1e-12",
                use_sample_file = false
        }
    }

    # Merge DecodeMe batches
    scatter (fileset in QCAndLiftOverDecodeMe.plink_fileset) {
        File decodeme_batches_bed_files = fileset.bed
    }

    call MergeGenotypingArrays as MergeDMEGenotypingArrays {
        input:
            docker_image = docker_image,
            julia_cmd = get_julia_cmd.julia_cmd,
            bed_files = decodeme_batches_bed_files,
            plink_filesets = QCAndLiftOverDecodeMe.plink_fileset,
            output_prefix = "decodeme.chrX.hg38.merged_batches"
    }

    # QC DecodeMe with Ref
    call QCWithRef as QCDecodeMeWithRef {
        input:
            docker_image = docker_image,
            julia_cmd = get_julia_cmd.julia_cmd,
            chr = MergeDMEGenotypingArrays.plink_fileset.chr,
            bed_file = MergeDMEGenotypingArrays.plink_fileset.bed,
            bim_file = MergeDMEGenotypingArrays.plink_fileset.bim,
            fam_file = MergeDMEGenotypingArrays.plink_fileset.fam,
            ref_stats = make_ref_panel_stats.ref_panel_stats
    }

    # QC and Liftover UKB Array
    call QCAndLiftOverArray as QCAndLiftOverUKB {
        input:
            docker_image = docker_image,
            chr = ukb_genotypes.chr,
            bed_file = ukb_genotypes.bed,
            bim_file = ukb_genotypes.bim,
            fam_file = ukb_genotypes.fam,
            hwe_pval = "1e-9",
            sample_keep_file=ukb_samples,
            use_sample_file = true
    }

    # QC UKB with Ref
    call QCWithRef as QCUKBWithRef {
        input:
            docker_image = docker_image,
            julia_cmd = get_julia_cmd.julia_cmd,
            chr = QCAndLiftOverUKB.plink_fileset.chr,
            bed_file = QCAndLiftOverUKB.plink_fileset.bed,
            bim_file = QCAndLiftOverUKB.plink_fileset.bim,
            fam_file = QCAndLiftOverUKB.plink_fileset.fam,
            ref_stats = make_ref_panel_stats.ref_panel_stats
    }

    # Merge DecodeMe and UKB
    Array[PLINKFileset] decodeme_and_ukb_filesets = [QCDecodeMeWithRef.plink_fileset, QCUKBWithRef.plink_fileset]
    scatter (fileset in decodeme_and_ukb_filesets) {
        File decodeme_and_ukb_bed_files = fileset.bed
    }
    call MergeGenotypingArrays as MergeUKBDecodeMe{
        input:
            docker_image = docker_image,
            julia_cmd = get_julia_cmd.julia_cmd,
            bed_files = decodeme_and_ukb_bed_files,
            plink_filesets = decodeme_and_ukb_filesets,
            output_prefix = "chrX.hg38.merged.ukb_decodeme"
    }

    call ConvertToBCF {
        input:
            docker_image = docker_image,
            chr=MergeUKBDecodeMe.plink_fileset.chr,
            bed_file=MergeUKBDecodeMe.plink_fileset.bed,
            bim_file=MergeUKBDecodeMe.plink_fileset.bim,
            fam_file=MergeUKBDecodeMe.plink_fileset.fam
    }

    output {
        File ref_panel_stats = make_ref_panel_stats.ref_panel_stats

        Array[PLINKFileset] lifted_over_decodeme_batches = QCAndLiftOverDecodeMe.plink_fileset
        PLINKFileset merged_decodeme_array = MergeDMEGenotypingArrays.plink_fileset
        PLINKFileset decodeme_qced_with_ref = QCDecodeMeWithRef.plink_fileset
        File decodeme_ref_action_stats = QCDecodeMeWithRef.action_stats

        PLINKFileset lifted_over_ukb = QCAndLiftOverUKB.plink_fileset
        PLINKFileset ukb_qced_with_ref = QCUKBWithRef.plink_fileset
        File ukb_ref_action_stats = QCUKBWithRef.action_stats

        PLINKFileset merged_qced_fileset = MergeUKBDecodeMe.plink_fileset
        File merged_qced_bcf = ConvertToBCF.bcf_file
    }
}

task QCWithRef {
    input {
        String docker_image
        String julia_cmd
        String chr
        File bed_file
        File bim_file
        File fam_file
        File ref_stats
    }

    String input_prefix = basename(bed_file, ".bed")

    command <<<
        bed_prefix=$(dirname "~{bed_file}")/$(basename "~{bed_file}" .bed)
        
        plink2 \
            --bfile ${bed_prefix} \
            --nonfounders \
            --freq \
            --out ${bed_prefix}

        ~{julia_cmd} /opt/decodeme/scripts/qc_with_ref.jl ${bed_prefix} ~{ref_stats} ~{input_prefix}

        plink \
            --bfile ${bed_prefix} \
            --flip ~{input_prefix}.to_flip.txt \
            --exclude ~{input_prefix}.to_drop.txt \
            --make-bed \
            --out temp

        plink \
            --bfile temp \
            --update-name ~{input_prefix}.new_ids.tsv \
            --make-bed \
            --out ~{input_prefix}.ref_qced
    >>>

    output {
        File action_stats = "${input_prefix}.action_stats.tsv"
        PLINKFileset plink_fileset = object {
            chr: chr,
            bed: "${input_prefix}.ref_qced.bed",
            bim: "${input_prefix}.ref_qced.bim",
            fam: "${input_prefix}.ref_qced.fam"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x8"
    }
}

task make_ref_panel_stats {
    input {
        String docker_image
        String julia_cmd
        Array[File] bcf_csi_files
    }

    command <<<
        for bcf_file in ~{sep=" " bcf_csi_files}; do
            if [[ "${bcf_file}" == *.bcf ]]; then
                bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n' ${bcf_file} > ${bcf_file%.bcf}.tsv
                echo "${bcf_file%.bcf}.tsv"
            fi
        done > stats_files.txt

        ~{julia_cmd} /opt/decodeme/scripts/merge_ref_panel_stats.jl stats_files.txt
    >>>

    output {
        File ref_panel_stats = "ref_panel_stats.tsv"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x8"
    }

}

task QCAndLiftOverArray {
    input {
        String docker_image
        String chr
        File bed_file
        File bim_file
        File fam_file
        String hwe_pval
        File? sample_keep_file
        Boolean use_sample_file = false
    }

    String input_prefix = basename(bed_file, ".bed")
    String? keep_opt = if (use_sample_file) then "--keep " + sample_keep_file else ""

    command <<<
        bed_prefix=$(dirname "~{bed_file}")/$(basename "~{bed_file}" .bed)

        # Apply HWE QC and format chromosome for liftover
        plink \
            --nonfounders \
            --bfile ${bed_prefix} \
            --hwe ~{hwe_pval} \
            --geno 0.02 \
            --mind 0.03 \
            --output-chr chrMT \
            --freq \
            --make-bed \
            --out ~{input_prefix}.hwe ~{keep_opt}
        # Make bed file for liftover
        awk '{
            chr = $1;
            start = $4 - 1;
            end = $4;
            print chr, start, end, $2
        }' OFS='\t' ~{input_prefix}.hwe.bim > variants.hg19.bed
        # Liftover
        liftOver variants.hg19.bed /opt/hg19ToHg38.over.chain.gz variants.hg38.bed variants.unlifted.bed
        # Update chromosomes and positions and only keep mapped snps
        awk '{
            id=$4;
            pos=$3;
            print id, pos
        }' variants.hg38.bed > update_map_positions.txt
        awk '{
            id=$4;
            chr=$1;
            print id, chr
        }' variants.hg38.bed > update_map_chr.txt
        cut -f4 variants.hg38.bed > kept_snps.txt
        plink \
            --bfile ~{input_prefix}.hwe \
            --update-map update_map_positions.txt \
            --update-chr update_map_chr.txt \
            --extract kept_snps.txt \
            --output-chr chrMT \
            --make-bed \
            --out ~{input_prefix}.hwe.hg38
    >>>

    output {
        PLINKFileset plink_fileset = object {
            chr: chr,
            bed: "${input_prefix}.hwe.hg38.bed",
            bim: "${input_prefix}.hwe.hg38.bim",
            fam: "${input_prefix}.hwe.hg38.fam"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x8"
    }
}

task MergeGenotypingArrays {
    input {
        String docker_image
        String julia_cmd
        Array[File] bed_files
        Array[PLINKFileset] plink_filesets
        String output_prefix
    }
    
    command <<<
        for f in ~{sep=" " bed_files}; do
            echo "${f%.bed}"
        done > merge_list.txt

        ~{julia_cmd} /opt/decodeme/scripts/make_shared_variants_list.jl merge_list.txt shared_variants.txt

        plink \
            --biallelic-only \
            --merge-list merge_list.txt \
            --output-chr chrMT \
            --extract shared_variants.txt \
            --make-bed \
            --out ~{output_prefix}
    >>>

    output {
        PLINKFileset plink_fileset = object {
            chr: "X",
            bed: "${output_prefix}.bed",
            bim: "${output_prefix}.bim",
            fam: "${output_prefix}.fam"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x8"
    }
}

task ConvertToBCF {
    input {
        String docker_image
        String chr
        File bed_file
        File bim_file
        File fam_file
    }

    String input_prefix = basename(bed_file, ".bed")

    command <<<
        bed_prefix=$(dirname "~{bed_file}")/$(basename "~{bed_file}" .bed)

        plink2 \
            --bfile ${bed_prefix} \
            --export bcf id-paste=iid \
            --out ${bed_prefix}
    >>>

    output {
        File bcf_file = "${input_prefix}.bcf"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x8"
    }
}