#name: Julia Lines Count
#description: file lines count
#language: julia
#input: file file
#input: bool header
#input: string separator = ","
#input: string dec = "."
#condition: file.isFile && file.name.endsWith("csv")
#output: int num_lines

num_lines = 0
csv_reader = CSV.File(file, header=header, delim=separator)
for row in csv_reader
    num_lines += 1
end
