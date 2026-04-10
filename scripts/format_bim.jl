using CSV
using DataFrames

function main()
    bim_file = ARGS[1]
    bim = CSV.read(bim_file, DataFrame; header=[:CHROM, :ID, :CM, :POS, :A1, :A2])
    transform!(
        bim,
        :A1 => (x -> replace.(x, "-" => "0")) => :A1,
        :A2 => (x -> replace.(x, "-" => "0")) => :A2,
        [:CHROM, :POS] => ByRow((c, p, a1, a2) -> string(c, ":", p)) => :ID
    )
    CSV.write(bim_file, bim; header=false, delim="\t")
    return 0
end

main()