//name: NodeJS Lines Count
//description: file lines count
//language: nodejs
//input: file file
//input: bool header
//input: string separator = ","
//input: string dec = "."
//condition: file.isFile && file.name.endsWith("csv")
//output: int num_lines

df = await DataFrame.fromCSV("file://" + path.resolve(file));
num_lines = df.dim()[0];
