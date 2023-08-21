//name: JavaScript Lines Count
//description: file lines count
//language: javascript
//input: file file
//input: bool header
//input: string separator = ","
//input: string dec = "."
//condition: file.isFile && file.name.endsWith("csv")
//output: int num_lines

df = await DG.DataFrame.fromCsv(file);
num_lines = df.rowCount;
