# Generates golden reference output for the JMP package tests (JMP_PLUGIN.md §9.3).
#
# Usage (from packages/Jmp):
#   julia --project=@jmpreader tools/jmp2json.jl files/jmpreader files/golden
#
# One-time setup of the named environment, pinned to the commit the TS port follows:
#   julia -e 'using Pkg; Pkg.activate("jmpreader"; shared=true);
#             Pkg.add(url="https://github.com/jaakkor2/JMPReader.jl", rev="9fd68aa2651a8dde002ecdbd06af4d8d832b18f2");
#             Pkg.add(["JSON3", "DataFrames", "ColorTypes"])'
#
# Output: one <name>.json per .jmp file with column names, element types and values.
# Date-times are ISO strings to the millisecond, dates yyyy-mm-dd, durations seconds,
# row states "<marker> #RRGGBB", missing null.

using JMPReader, JSON3, Dates
using ColorTypes: RGB
hex(c) = join(string.(round.(Int, 255 .* (Float64(c.r), Float64(c.g), Float64(c.b))), base = 16, pad = 2))

encode(x::Missing) = nothing
encode(x::Number) = isnan(float(x)) ? nothing : x
encode(x::AbstractString) = String(x)
encode(x::DateTime) = Dates.format(x, dateformat"yyyy-mm-ddTHH:MM:SS.sss")
encode(x::Date) = Dates.format(x, dateformat"yyyy-mm-dd")
encode(x::Period) = Dates.value(convert(Millisecond, x)) / 1000
encode(x::JMPReader.Rowstate) = string(x.marker, " #", uppercase(hex(x.color)))
encode(x) = string(x)

function convert_file(src, dst)
    df = readjmp(src)
    cols = [(name = n, type = string(nonmissingtype(eltype(df[!, n]))), values = encode.(df[!, n]))
            for n in names(df)]
    open(dst, "w") do io
        JSON3.write(io, (nrows = size(df, 1), columns = cols))
    end
end

srcdir, dstdir = ARGS
mkpath(dstdir)
for f in sort(filter(endswith(".jmp"), readdir(srcdir)))
    try
        convert_file(joinpath(srcdir, f), joinpath(dstdir, replace(f, r"\.jmp$" => ".json")))
        println("ok    ", f)
    catch e
        println("FAIL  ", f, ": ", e)
    end
end
