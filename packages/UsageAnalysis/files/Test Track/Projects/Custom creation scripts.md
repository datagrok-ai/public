1. Run the script:
  ```
  //name: getLastCreatedFile  
  //language: javascript
  //output: dataframe df  

  const fileList = await grok.dapi.files.list('System:DemoFiles/chem', true, '');  
  const csvFiles = fileList.filter((fi) => fi.fileName.endsWith('.csv'));  
  if (csvFiles.length === 0)
  throw new Error('No CSV files found in System:DemoFiles/chem');  
    
  csvFiles.sort((a, b) => b.modified - a.modified);  
  const lastModifiedFile = csvFiles[0];  
  const csv = await grok.dapi.files.readAsText(lastModifiedFile.fullPath);  
  df = DG.DataFrame.fromCsv(csv);
  ```
2. Add some viewers and save the project with data sync enabled
1. Close all
1. (Optional) Update any CSV file in the 'System:DemoFiles/chem' folder
1. Open the saved project - verify, that:
  - No errors occur
  - The most recently created or modified file in the `chem` folder is loaded