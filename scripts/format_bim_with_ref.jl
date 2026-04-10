using CSV
using DataFrames

function get_allele_refmap(ref_bim)
    refmap = Dict()
    map(eachrow(ref_bim)) do row
        refmap[string(row.CHROM, ":", row.POS)] = Dict(row.A1 => row.A2, row.A2 => row.A1)
    end
    return refmap
end

function fix_missing_allele(row, refmap)
    
    chr = row.CHROM
    pos = row.POS
    a1 = row.A1
    a2 = row.A2
    refmap_key = string(chr, ":", pos)
    if !haskey(refmap, refmap_key)
        return a1, a2
    end
    pos_refmap = refmap[refmap_key]
    # Try fixing missing a1 with reference
    if a1 == "-" || a1 == "0"
        if haskey(pos_refmap, a2)
            a1 = pos_refmap[a2]
        end
    end
    # Try fixing missing a2 with reference
    if a2 == "-" || a2 == "0"
        if haskey(pos_refmap, a1)
            a2 = pos_refmap[a1]
        end
    end
    
    return a1, a2
end

function main()
    bim_file = ARGS[1]
    ref_bim_file = ARGS[2]
    variant_to_drop_file = ARGS[3]

    # Load reference and build reference allele map
    ref_bim = CSV.read(ref_bim_file, DataFrame; header=[:CHROM, :ID, :CM, :POS, :A1, :A2])
    refmap = get_allele_refmap(ref_bim)

    # Try fix missing A1/A2 from reference
    # If they cannot be fixed it means they are either not in the reference or alleles do not match
    # We will drop these variants
    bim = CSV.read(bim_file, DataFrame; header=[:CHROM, :ID, :CM, :POS, :A1, :A2])
    transform!(bim, 
        AsTable(:) => ByRow(row -> fix_missing_allele(row, refmap)) => [:A1, :A2]
    )
    # Set variant ID to be chr:po:a1:a2 where a1 and a2 are just alphabetically sorted
    transform!(
        bim,
        [:CHROM, :POS, :A1, :A2] => ByRow((c, p, a1, a2) -> string(c, ":", p, ":", join(sort([a1, a2]), ":"))) => :ID
    )
    # Identify mismatching variants
    variant_ids_to_drop = Set(String[])
    map(eachrow(bim)) do row
        if (row.A1 == "-" || row.A2 == "-" || row.A1 == "0" || row.A2 == "0")
            push!(variant_ids_to_drop, row.ID)
        end
    end
    open(variant_to_drop_file, "w") do io
        for id in variant_ids_to_drop
            println(io, id)
        end
    end
    # Overwrite bim
    CSV.write(bim_file, bim; header=false, delim="\t")
    return 0
end

main()