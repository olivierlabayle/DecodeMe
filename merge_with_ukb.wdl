version 1.0

import "structs.wdl"

workflow merge_with_ukb {
    input {
        String docker_image = "olivierlabayle/genomicc:0.3.0"
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

    scatter (fileset in decodeme_genotypes) {
        call QCDecodeMeArray {
            input:
                docker_image = docker_image,
                chr = fileset.chr,
                bed_file = fileset.bed,
                bim_file = fileset.bim,
                fam_file = fileset.fam,
        }
    }

    call QCUKBArray {
        input:
            docker_image = docker_image,
            chr = ukb_genotypes.chr,
            bed_file = ukb_genotypes.bed,
            bim_file = ukb_genotypes.bim,
            fam_file = ukb_genotypes.fam,
            sample_file = ukb_samples
    }


}

task get_julia_cmd {
    input {
        String use_sysimage = "true"
        String threads = "auto"
    }
    command <<<
        julia_cmd_string="julia --project=/opt/PopGen --startup-file=no"
        if [[ "~{use_sysimage}" == "true" ]]; then
            julia_cmd_string+=" --sysimage=/opt/PopGen/sysimage.so"
        fi
        if [[ "~{threads}" == "auto" ]]; then
            julia_cmd_string+=" --threads=auto"
        fi
        julia_cmd_string+=" /opt/PopGen/bin/wdl-gwas.jl"
        echo "$julia_cmd_string"
    >>>

    output {
        String julia_cmd = read_string(stdout())
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