//name: JavaScript Lines Count
//description: file lines count
//language: javascript
//input: file file
//input: bool header
//input: string separator = ","
//input: string dec = "."
//condition: file.isFile && file.name.endsWith("csv")
//output: int num_lines
let data = await new DG.FileSource().readAsText('System:AppData/' + file.path);
let d = await DG.DataFrame.fromCsv(data, {delimiter: separator, decimalSeparator: dec, headerRow: header});
let num_lines = d.rowCount;
