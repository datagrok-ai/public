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

  // BUFFERS ALLOCATION

  let buffers = [];
  let excpectedNumOfBuffers = 3; // 1st - for source data, 
                                 // 2nd - for principal components
                                 // 3rd - for approximated data

  let bufForData = module._malloc(sizeOfData * numOfBytes);
  if(bufForData != 0)
      buffers.push(bufForData);

  let bufForPrincComp = module._malloc(sizeOfPrincComp * numOfBytes);
  if(bufForPrincComp != 0)
      buffers.push(bufForPrincComp);

  let bufForApprox = module._malloc(sizeOfData * numOfBytes);
  if(bufForApprox != 0)
      buffers.push(bufForApprox);
  
  console.log("Allocated buffers:");
  console.log(buffers);

  // check memory allocation
  if(buffers.length == excpectedNumOfBuffers) {

    console.log('Success: the required memory is allocated!');

      // assign allocated buffer with source data
      for(let i = 0; i < height; i++)
          heap.set(sourceArrays[i], (bufForData + i * width * numOfBytes) >> shift);
    
      // EXPORTED METHOD CALL
      let res = module.ccall('principalComponentAnalysis', 
      'number', ['number', 'number', 'number', 'number', 'number', 'number'], 
      [bufForData, height, width, numOfPrincipalComponents, bufForPrincComp, bufForApprox]);

      // get results: principal components    
      for(let i = 0; i < numOfPrincipalComponents; i++)
          principalComponents.push(new Float32Array(heap.buffer, 
              bufForPrincComp + i * width * numOfBytes, width));

      // get results: approximation of source data
      for(let i = 0; i < height; i++)
          approxData.push(new Float32Array(heap.buffer, 
              bufForApprox + i * width * numOfBytes, width));
  
      // evaluate error - maximum absolute deviation  <-- OPTIONAL
      console.log('Maximum absolute deviation: ' + module._error(bufForData, bufForApprox, sizeOfData));

  } // if
  else { // i.e. memory allocation fail
      console.log('FAIL:the required memory is not allocated!');
  }
  
  // clear allocted memory
  for(let i = 0; i < buffers.length; i++)
      module._free(buffers[i]);  
}

/*// A wrapper for PCA by Eigen
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
  let bufForPrincComp = module._malloc(sizeOfPrincComp * numOfBytes);
  let bufForApprox = module._malloc(sizeOfData * numOfBytes);

  // assign allocated buffer with source data
  for(let i = 0; i < height; i++)
      heap.set(sourceArrays[i], (bufForData + i * width * numOfBytes) >> shift);
    
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
  
  // evaluate error - maximum absolute deviation  <-- OPTIONAL
  console.log("Maximum absolute deviation: " + module._error(bufForData, bufForApprox, sizeOfData));
  
  // clear allocted memory
  module._free(bufForData);
  module._free(bufForPrincComp);
  module._free(bufForApprox);
}*/

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
  console.log('PCA ...');

  let start = new Date().getTime(); 

  pca(NewModule, data, numOfComponents, princComp, approx);

  let finish = new Date().getTime();

  console.log('Principal components and approximation of source data have been computed!');
  console.log('Time for PCA execution: ' + (finish - start) + ' ms.');

  // Create dataframe with approximation
  console.log('Creating dataframes from the arrays computed ...');

  start = new Date().getTime();
  
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

  finish = new Date().getTime();

  console.log('Dataframes and their views have been created!');
  console.log('Time for creating dataframes and their views: ' + (finish - start) + ' ms.');
}
