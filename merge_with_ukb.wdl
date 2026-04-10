version 1.0

import "structs.wdl"

workflow merge_with_ukb {
    input {
        String docker_image = "olivierlabayle/genomicc:0.3.0"
        Array[PLINKFileset] decodeme_genotypes
        PLINKFileset ukb_genotypes
    }

    scatter (fileset in decodeme_genotypes) {
        call HWQC as DecodeMEHWQC{
            input:
                docker_image = docker_image,
                chr = chr,
                bed_file = bed_file,
                bim_file = bim_file,
                fam_file = fam_file,
                pval = "1e-12"

        }
    }
}


task QCDecodeMeArray {
    input {
        String docker_image
        String chr
        File bed_file
        File bim_file
        File fam_file
        String pval
    }

    bed_prefix = basename(bed_file, ".bed")
    String is_sample_file_provided = select_first([sample_file, ""])

    command <<<
        sample_filter_opt=""
        if [[ "~{is_sample_file_provided}" != "" ]]; then
            sample_filter_opt="--keep ~{sample_file}"
        fi

        plink \
        --bfile ~{bed_prefix} \
        --geno 0.03 \
        --hwe 1e-12 \
        --make-bed \
        --out ~{bed_prefix}.hwe ${sample_filter_opt}
    >>>

    output {
        PLINKFileset ld_pruned_fileset = object {
            chr: "all",
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
        --geno 0.03 \
        --hwe 1e-9 \
        --make-bed \
        --keep ~{sample_file} \
        --out ~{bed_prefix}.hwe
    >>>

    output {
        PLINKFileset ld_pruned_fileset = object {
            chr: "all",
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