1. Run the script:
  ```
 //name: getLastCreatedFile_by_numeric_suffix
//language: javascript
//output: dataframe result

const csvFiles = await grok.dapi.files.list('System:DemoFiles/chem', true, 'csv');

if (csvFiles.length === 0)
    throw new Error('No CSV files found in System:DemoFiles/chem');

// Function to extract a numeric suffix from the file name
function getNumberSuffix(name) {
    const match = name.match(/(\d+)(?=\.csv$)/);
    return match ? parseInt(match[1], 10) : -1; // -1 if no number is found
}

// Sort files by numeric suffix in ascending order
csvFiles.sort((a, b) => getNumberSuffix(a.fileName) - getNumberSuffix(b.fileName));

// Select the file with the highest numeric suffix
const lastFile = csvFiles[csvFiles.length - 1];

const csv = await grok.dapi.files.readAsText(lastFile.fullPath);
result = DG.DataFrame.fromCsv(csv);
  ```
2. Add some viewers and save the project with data sync enabled
1. Close all
1. Update any CSV file in the 'System:DemoFiles/chem' folder
1. Open the saved project - verify, that:
  - No errors occur
  - The most recently created or modified file in the `chem` folder is loaded