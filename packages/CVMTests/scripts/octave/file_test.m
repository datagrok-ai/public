#name: Octave Lines Count
#description: file lines count
#language: octave
#input: file file
#input: bool header
#input: string separator = ","
#input: string dec = "."
#condition: file.isFile && file.name.endsWith("csv")
#output: int num_lines

num_lines = rows(csv2cell(file, separator));
