#name: Python Lines Count
#description: file lines count
#language: python
#input: file file
#condition: file.isFile && (file.name.endsWith("csv") || file.name.endsWith("yml") || file.name.endsWith("yaml"))
#output: int num_lines

num_lines = 0
with open(file, "r") as f:
    for line in f:
        if line != "\n":
            num_lines += 1
