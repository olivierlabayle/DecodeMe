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

    scatter (decodeme_batch in decodeme_genotypes) {
        call QCAndLiftOverArray as QCAndLiftOverDecodeMe {
            input:
                docker_image = docker_image,
                chr = decodeme_batch.chr,
                bed_file = decodeme_batch.bed,
                bim_file = decodeme_batch.bim,
                fam_file = decodeme_batch.fam,
                hwe_pval = "1e-12"
        }
    }

    call QCAndLiftOverArray as QCAndLiftOverUKB {
        input:
            docker_image = docker_image,
            chr = ukb_genotypes.chr,
            bed_file = ukb_genotypes.bed,
            bim_file = ukb_genotypes.bim,
            fam_file = ukb_genotypes.fam,
            hwe_pval = "1e-9"
    }

    scatter (fileset in QCAndLiftOverDecodeMe.plink_fileset) {
        File decodeme_batches_bed_files = fileset.bed
    }

    call MergeGenotypingArrays as MergeDMEGenotypingArrays {
        input:
            docker_image = docker_image,
            bed_files = decodeme_batches_bed_files,
            plink_filesets = QCAndLiftOverDecodeMe.plink_fileset,
            output_prefix = "decodeme.chrX.hg38.merged_batches"
    }

    output {
        Array[PLINKFileset] lifted_over_decodeme_batches = QCAndLiftOverDecodeMe.plink_fileset
        PLINKFileset lifted_over_ukb = QCAndLiftOverUKB.plink_fileset
        File ref_panel_stats = make_ref_panel_stats.ref_panel_stats
        PLINKFileset merged_decodeme_array = MergeDMEGenotypingArrays.plink_fileset
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
    }

    String input_prefix = basename(bed_file, ".bed")

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
            --out ~{input_prefix}.hwe
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
        Array[File] bed_files
        Array[PLINKFileset] plink_filesets
        String output_prefix
    }
    
    command <<<
        for f in ~{sep=" " bed_files}; do
            echo "${f%.bed}"
        done > merge_list.txt

        plink \
            --biallelic-only \
            --merge-list merge_list.txt \
            --output-chr chrMT \
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


task MergeDecodeMeAndUKB {
    input {
        String docker_image
        String ukb_chr
        File ukb_bed_file
        File ukb_bim_file
        File ukb_fam_file
        File ukb_sample_file

        Array[File] decodeme_bed_files
        Array[File] decodeme_bim_files
        Array[File] decodeme_fam_files
    }

    String ukb_bed_prefix = basename(ukb_bed_file, ".bed")

    command <<<
    for f in ~{sep=" " decodeme_bed_files}; do
        echo "${f%.bed}"
    done > merge_list.txt
    echo "~{ukb_bed_prefix}" >> merge_list.txt

    plink \
        --biallelic-only \
        --merge-list merge_list.txt \
        --extract range shared_variants.txt
        --geno 0.02 \
        --mind 0.03 \
        --export bcf \
        --out chrX
    >>>
}