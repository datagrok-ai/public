
1. Create a script 
    //name: test_Layout  
    //language: javascript
    //input: int idx=1
    //output: dataframe df  

    const fileList = await grok.dapi.files.list('System:DemoFiles/chem', true, '');  
    const csvFiles = fileList.filter((fi) => fi.fileName.endsWith('.csv'));  
    if (csvFiles.length === 0)
  throw new Error('No CSV files found in System:DemoFiles/chem');  
    
    csvFiles.sort((a, b) => b.createdOn - a.createdOn);  
    const lastModifiedFile = csvFiles[idx];
    const csv = await grok.dapi.files.readAsText(lastModifiedFile.fullPath);  
    df = DG.DataFrame.fromCsv(csv);
1. Go to the Layout tab
4. Click **Run script**
4. Add some viewers (with and without docking one over another)
1. Add coloring, change style, hide some columns
5. Save
6. Close All
7. Run the test_Layout script - the result should open with the new layout
7. Add some new viewers
7. Save the project
7. Close All
7. Open the project - check the layout
6. Go to the **Toolbox > File** and click **Refresh** - the layout shouldn't change
---
{
  "order": 12
}