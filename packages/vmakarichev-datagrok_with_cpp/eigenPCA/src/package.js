/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//tags: init
export async function init() {
  await initNewModule();  
}

// A wrapper for PCA by Eigen
function pca(module, sourceArrays, numOfPrincipalComponents, principalComponents, approxData ){
 
  // sizes
  let height = sourceArrays.length;
  let width = sourceArrays[0].length;
  let sizeOfData = height * width;
  let sizeOfPrincComp = numOfPrincipalComponents * width; 
  let numOfBytes = sourceArrays[0].BYTES_PER_ELEMENT;

  // heap
  let shift = 2;
  let heap = module.HEAPF32;

  // allocate buffer for source data
  let bufForData = module._malloc(sizeOfData * numOfBytes);

  // assign allocated buffer with source data
  for(let i = 0; i < height; i++)
      heap.set(sourceArrays[i], (bufForData + i * width * numOfBytes) >> shift);
  
  // allocate buffers
  let bufForPrincComp = module._malloc(sizeOfPrincComp * numOfBytes);
  let bufForApprox = module._malloc(sizeOfData * numOfBytes);

  // EXPORTED METHOD CALL
  let res = module.ccall('principalComponentAnalysis', 
  'number', ['number', 'number', 'number', 'number', 'number', 'number'], 
  [bufForData, height, width, 
      numOfPrincipalComponents, bufForPrincComp, bufForApprox]);

  // get results: principal components    
  for(let i = 0; i < numOfPrincipalComponents; i++)
      principalComponents.push(new Float32Array(heap.buffer, 
          bufForPrincComp + i * width * numOfBytes, width));

  // get results: approximation of source data
  for(let i = 0; i < height; i++)
      approxData.push(new Float32Array(heap.buffer, 
          bufForApprox + i * width * numOfBytes, width));
  
  // clear allocted memory
  module._free(bufForData);
  module._free(bufForPrincComp);
  module._free(bufForApprox);
}

//name: testOfPCA
//input: dataframe df
//input: int numOfComponents = 3
export function testOfPCA(df, numOfComponents) {
  
  // structures for column data
  let data = [];
  let approx = [];
  let princComp = [];
  let names = [];

  // create dataset
  for(let col of df.columns){
    switch(col.type){
      case 'double':
        data.push(col.getRawData());
        break;
      case 'int':
        data.push(new Float32Array(col.getRawData()));
        break;
      default:
        alert('Wrong type of data!');
        return;
    }
    names.push(col.name);
  }

  // PRINCIPAL COMPONENT ANALYSIS
  pca(NewModule, data, numOfComponents, princComp, approx);

  // Create dataframe with approximation
  let approxCols = [];  
  for(let i = 0; i < data.length; i++)
    approxCols.push(DG.Column.fromFloat32Array(names[i] + '(approx)', approx[i]));  
  let tableWithApprox = DG.DataFrame.fromColumns(approxCols);
  tableWithApprox.name = 'Approximation';
  grok.shell.addTableView(tableWithApprox);

  // Create dataframe with principal components
  let princCompCols = [];
  for(let i = 1; i <= numOfComponents; i++)
    princCompCols.push(DG.Column.fromFloat32Array('comp. ' + i.toString(), princComp[i - 1]));
  let tableWithPrincComp = DG.DataFrame.fromColumns(princCompCols);
  tableWithPrincComp.name = 'PrincipalComponents';
  grok.shell.addTableView(tableWithPrincComp);
}
