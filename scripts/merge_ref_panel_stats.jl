using CSV
using DataFrames

function main()
    input_file = ARGS[1]
    stats_files = readlines(input_file)
    stats = mapreduce(
        file -> CSV.read(file, DataFrame; header=[:CHROM, :POS, :ID, :REF, :ALT, :AC, :AN]), 
        vcat, 
        stats_files
    )
    unique!(stats, :ID)
    stats.AF = stats.AC ./ stats.AN
    CSV.write("ref_panel_stats.tsv", stats; delim="\t")
    return 0
end

main()