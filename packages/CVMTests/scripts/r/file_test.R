#name: R Lines Count
#description: file lines count
#language: r
#input: file file
#input: bool header
#input: string separator = ","
#input: string dec = "."
#condition: file.isFile && file.name.endsWith("csv")
#output: int num_lines

data = read.csv(file, header = header, sep = separator, dec = dec)
num_lines = nrow(data)
