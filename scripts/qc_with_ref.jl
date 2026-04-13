using CSV
using DataFrames

const ALLELE_COMPLEMENT = Dict("A" => "T", "T" => "A", "C" => "G", "G" => "C")

function load_bim_with_freq(bed_prefix)
    bim = CSV.read(string(bed_prefix, ".bim"), DataFrame; header=[:CHROM, :ID, :CM, :POS, :A1, :A2])
    freqs = CSV.read(string(bed_prefix, ".afreq"), DataFrame; delim="\t", select=[:ID, :ALT, :ALT_FREQS])
    leftjoin!(bim, freqs, on=:ID)
    transform!(bim, [:A1, :A2, :ALT, :ALT_FREQS] => ByRow((a1, a2, alt, alt_freq) -> alt_freq < 0.5 ? alt : filter(!=(alt), [a1, a2])[1]) => :MINOR_ALLELE)
    return bim
end

function get_allele_refmap(ref_stats_file, bim_positions)
    ref_stats = CSV.read(ref_stats_file, DataFrame)
    refmap = Dict()
    map(eachrow(unique(ref_stats, :POS, keep=:noduplicates))) do row
        if row.POS in bim_positions
            refmap[string("chrX:", row.POS)] = (
                allele_map = Dict(row.REF => row.ALT, row.ALT => row.REF), 
                ref=row.REF, 
                alt=row.ALT, 
                af=row.AF, 
                minor_allele = row.AF < 0.5 ? row.ALT : row.REF
            )
        end
    end
    return refmap
end

is_palyndromic(alleles) = 
    alleles == Set(["A", "T"]) || alleles == Set(["C", "G"])

function get_action(row, refmap)
    key = string(row.CHROM, ":", row.POS)
    if !haskey(refmap, key)
        return "DROP (VARIANT-NOT-IN-REF)"
    end
    allele_map, ref_ref, ref_alt, ref_alt_af, ref_minor_allele = refmap[key]
    ref_alleles = Set([ref_ref, ref_alt])
    bim_alleles = Set([row.A1, row.A2])
    
    # Check if alleles match
    if ref_alleles == bim_alleles
        # Check if minor allele match
        if ref_minor_allele == row.MINOR_ALLELE
            return "KEEP (FULL-MATCH)"
        # Otherwise, is the SNP palyndromic ?
        elseif is_palyndromic(ref_alleles)
            # Can it be resolved ?
            if row.ALT_FREQS > 0.4
                return "DROP (PALYNDROMIC-HIGH-MAF)"
            else
                return "FLIP (PALYNDROMIC-LOW-MAF)"
            end
        elseif (abs(row.ALT_FREQS - 0.5) < 0.05) && (abs(ref_alt_af - 0.5) < 0.05)
            return "KEEP (OPPOSITE-MINOR-ALLELE-BUT-NEAR-0.5)"
        # This should not occur much
        else
            return "DROP (OPPOSITE-MINOR-ALLELE)"
        end
    else
        if "0" in bim_alleles || "-" in bim_alleles
            return "DROP (UNKOWN-ALLELE)"
        else
            if (row.A1 in keys(ALLELE_COMPLEMENT)) && (row.A2 in keys(ALLELE_COMPLEMENT))
                new_minor_allele = ALLELE_COMPLEMENT[row.MINOR_ALLELE]
                new_bim_alleles = Set([ALLELE_COMPLEMENT[row.A1], ALLELE_COMPLEMENT[row.A2]])
                if new_bim_alleles == ref_alleles
                    if ref_minor_allele == new_minor_allele
                        return "FLIP (COMPLEMENT-FULL-MATCH)"
                    elseif (abs(row.ALT_FREQS - 0.5) < 0.05) && (abs(ref_alt_af - 0.5) < 0.05)
                        return "KEEP (COMPLEMENT-OPPOSITE-MINOR-ALLELE-BUT-NEAR-0.5)"
                    # This should not occur much
                    else
                        return "DROP (COMPLEMENT-OPPOSITE-MINOR-ALLELE)"
                    end
                else
                    return "DROP (ALLELES-NOT-MATCHING)"
                end
            else
                return "DROP (ALLELES-NOT-MATCHING)"
            end
        end
    end
end

# function fix_alleles(a1, a2, minor_allele, allele_map)
#     a1 = row.A1
#     a2 = row.A2
#     # Try fixing missing a1 with reference
#     if a1 == "-" || a1 == "0"
#         if haskey(allele_map, a2)
#             a1 = allele_map[a2]
#         end
#     end
#     # Try fixing missing a2 with reference
#     if a2 == "-" || a2 == "0"
#         if haskey(allele_map, a1)
#             a2 = allele_map[a1]
#         end
#     end
    
#     return a1, a2
# end

function main()
    bed_prefix = ARGS[1]
    ref_stats_file = ARGS[2]
    output_prefix =  ARGS[3]

    # Load data
    bim = load_bim_with_freq(bed_prefix)

    # Load reference and build reference allele map
    refmap = get_allele_refmap(ref_stats_file, Set(bim.POS))

    # QC
    transform!(bim, 
        AsTable(:) => ByRow(row -> get_action(row, refmap)) => :ACTION
    )
    
    # Write number of variants per action taken
    action_stats = combine(groupby(bim, :ACTION), nrow => :N_VARIANTS)
    CSV.write(string(output_prefix, ".action_stats.tsv"), action_stats, delim="\t")

    # Write drop list
    drop_list = subset(bim, :ACTION => x -> startswith.(x, "DROP")).ID
    open(string(output_prefix, ".to_drop.txt"), "w") do io
        for vid in drop_list
            println(io, vid)
        end
    end
    # Write flip list
    flip_list = subset(bim, :ACTION => x -> startswith.(x, "FLIP")).ID
    open(string(output_prefix, ".to_flip.txt"), "w") do io
        for vid in flip_list
            println(io, vid)
        end
    end
    # Write new variant ids
    bim.NEW_ID = map(zip(bim.CHROM, bim.POS, bim.ID)) do (chr, pos, old_id)
        ref_key = string(chr, ":", pos)
        newid = if haskey(refmap, ref_key)
            _, ref_ref, ref_alt, _, _ = refmap[ref_key]
            string(ref_key, ":", ref_ref, ":", ref_alt)
        else
            old_id
        end
    end
    CSV.write(string(output_prefix, ".new_ids.tsv"), bim[!, [:ID, :NEW_ID]]; header=false, delim="\t")

    return 0
end

main()