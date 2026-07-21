//name: NodeJS Lines Count
//description: file lines count
//language: nodejs
//input: file file
//input: bool header
//input: string separator = ","
//input: string dec = "."
//condition: file.isFile && file.name.endsWith("csv")
//output: int num_lines

// TODO(GROK-20443): drop the legacy dataframe-js branch once the Node js-api
// runtime (reddata side) is on master and deployed.
if (typeof DG !== 'undefined') {
  df = DG.DataFrame.fromCsv(await fs.readFile(file, 'utf-8'));
  num_lines = df.rowCount;
}
else {
  df = await DataFrame.fromCSV("file://" + path.resolve(file));
  num_lines = df.dim()[0];
}
