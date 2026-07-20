//name: NodeJS Lines Count
//description: file lines count
//language: nodejs
//input: file file
//input: bool header
//input: string separator = ","
//input: string dec = "."
//condition: file.isFile && file.name.endsWith("csv")
//output: int num_lines

df = DG.DataFrame.fromCsv(await fs.readFile(file, 'utf-8'));
num_lines = df.rowCount;
