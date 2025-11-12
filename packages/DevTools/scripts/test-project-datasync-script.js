//name: TestProjectDatasyncScript
//language: javascript
//output: dataframe df

const fileList = await grok.dapi.files.list('System:DemoFiles/chem', true, '');
const csvFiles = fileList.filter((fi) => fi.fileName.endsWith('.csv'));
if (csvFiles.length === 0)
  throw new Error('No CSV files found in System:DemoFiles/chem');

const file = csvFiles.find((fi) => fi.fileName === 'SPGI.csv') ?? csvFiles[0];
const csv = await grok.dapi.files.readAsText(file.fullPath);
df = DG.DataFrame.fromCsv(csv);
