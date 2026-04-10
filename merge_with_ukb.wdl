version 1.0

import "structs.wdl"

workflow merge_with_ukb {
    input {
        String docker_image = "olivierlabayle/genomicc:0.3.0"

        Array[Files] reference_panel_files
        Array[PLINKFileset] decodeme_genotypes
        PLINKFileset ukb_genotypes

        julia_use_sysimage = "true"
        julia_threads = "auto"
    }

    # Get generic Julia command
    call get_julia_cmd as get_julia_cmd {
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

    # scatter (fileset in decodeme_genotypes) {
    #     call QCDecodeMeArray {
    #         input:
    #             docker_image = docker_image,
    #             chr = fileset.chr,
    #             bed_file = fileset.bed,
    #             bim_file = fileset.bim,
    #             fam_file = fileset.fam,
    #     }
    # }

    # call QCUKBArray {
    #     input:
    #         docker_image = docker_image,
    #         chr = ukb_genotypes.chr,
    #         bed_file = ukb_genotypes.bed,
    #         bim_file = ukb_genotypes.bim,
    #         fam_file = ukb_genotypes.fam,
    #         sample_file = ukb_samples
    # }


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

        ${julia_cmd} /opt/decodeme/scripts/merge_ref_panel_stats.jl stats_files.txt
    >>>

    output {
        File ref_panel_stats = "ref_panel_stats.tsv"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x8"
    }

}

task QCDecodeMeArray {
    input {
        String docker_image
        String chr
        File bed_file
        File bim_file
        File fam_file
    }

    bed_prefix = basename(bed_file, ".bed")

    command <<<
        plink \
        --bfile ~{bed_prefix} \
        --hwe 1e-12 \
        --make-bed \
        --out ~{bed_prefix}.hwe
    >>>

    output {
        PLINKFileset ld_pruned_fileset = object {
            chr: chr,
            bed: "${bed_prefix}.hwe.bed",
            bim: "${bed_prefix}.hwe.bim",
            fam: "${bed_prefix}.hwe.fam"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x8"
    }
}

task QCUKBArray {
    input {
        String docker_image
        String chr
        File bed_file
        File bim_file
        File fam_file
        File sample_file
    }

    bed_prefix = basename(bed_file, ".bed")

    command <<<
        plink \
        --bfile ~{bed_prefix} \
        --hwe 1e-9 \
        --make-bed \
        --keep ~{sample_file} \
        --out ~{bed_prefix}.hwe
    >>>

    output {
        PLINKFileset ld_pruned_fileset = object {
            chr: chr,
            bed: "${bed_prefix}.hwe.bed",
            bim: "${bed_prefix}.hwe.bim",
            fam: "${bed_prefix}.hwe.fam"
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

    ukb_bed_prefix = basename(ukb_bed_file, ".bed")

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