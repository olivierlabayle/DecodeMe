using CSV
using DataFrames

function get_bim_ids(prefix)
    bim = CSV.read(string(prefix, ".bim"), DataFrame; header=[:CHROM, :ID, :CM, :POS, :A1, :A2])
    return bim.ID
end

function main()
    # Find shared variants
    prefixes = readlines(ARGS[1])
    first_prefix = pop!(prefixes)
    shared_ids = get_bim_ids(first_prefix)
    for prefix in prefixes
        bim_ids = get_bim_ids(prefix)
        intersect!(shared_ids, bim_ids)
    end
    # Write them to file
    open(ARGS[2], "w") do io
        for shared_id in shared_positions
            println(io, shared_id)
        end
    end
    return 0
end