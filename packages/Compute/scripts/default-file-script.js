//name: Default import script
//description: desc
//language: javascript
//input: file uploadedFile
//output: dataframe df { viewer: Grid(); }
//editor: Compute:RichFunctionViewEditor

const wb = new ExcelJS.Workbook();
await wb.xlsx.load(await uploadedFile.arrayBuffer());

const df_ws = wb.getWorksheet("df");
let data = await df_ws.workbook.csv.writeBuffer();
df = grok.data.parseCsv(data.toString());
