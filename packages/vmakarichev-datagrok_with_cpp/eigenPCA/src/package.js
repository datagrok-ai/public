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
  await initEigenPCA();  
}

// PCA by Eigen: approximation is computed and compared to source data
function pcaWithVerification(module, columns, numOfPrincipalComponents, principalComponents, approxData){
  
  let resultCode;
    
  // sizes
  let height = columns.length;
  let width = columns[0].length;
  let sizeOfData = height * width;
  let sizeOfPrincComp = numOfPrincipalComponents * width; 
  let numOfBytes = Float32Array.BYTES_PER_ELEMENT; // since we operate Float32

  // heap: here, we take into account an operating Float32Arrays
  let shift = 2;
  let heap = module.HEAPF32;

  // BUFFERS ALLOCATION

  let buffers = [];
  let excpectedNumOfBuffers = 3; // 1st - for source data, 
                                 // 2nd - for principal components
                                 // 3rd - for approximated data: OPTIONAL

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

    let namesOfColumns = [];

    // assign allocated buffer with source data
    for(let i = 0; i < columns.length; i++){
      let array = null;
      let col = columns[i];

      // get raw data from column
      switch(col.type){
        case 'double':
          array = col.getRawData();
          break;
        case 'int':
          array = new Float32Array(col.getRawData()); // conversion is applied, 
                                                      // since PCA requires floats
          break;
        default:
          console.log('Column ' + col.name + ' has uncorrect type: ' + col.type);
          height--; // i.e. actual number of columns, which are passed to PCA, 
                    // should be reduced
          break;
      }

      // check data array
      if(array != null) {
        heap.set(array, (bufForData + i * width * numOfBytes) >> shift);
        namesOfColumns.push(col.name);
      }
    }

    // check size of data to be processed:
    // height == 0 means that there is no column data in bufForData
    if(height > 0){ 
      
      // EXPORTED С++-FUNCTION CALL
      resultCode = module.ccall('pcaWithApproximation', 'number', 
                     ['number', 'number', 'number', 'number', 'number', 'number'],
                     [bufForData, height, width, numOfPrincipalComponents, 
                      bufForPrincComp, bufForApprox]);
     
      // check result of the exported C++-function call:
      // resultCode != 0 means some error (see PCA/PCA.h)
      if(resultCode == 0) {
        // get results: principal components
        for(let i = 0; i < numOfPrincipalComponents; i++)
          principalComponents.push(new Float32Array(heap.buffer, 
             bufForPrincComp + i * width * numOfBytes, width));

        // get results: approximation of source data 
        for(let i = 0; i < height; i++)
          approxData.push(DG.Column.fromFloat32Array(
            'approx (' + namesOfColumns[i] + ')',
            new Float32Array(heap.buffer, 
            bufForApprox + i * width * numOfBytes, width)));
        
        // evaluate error - maximum absolute deviation  <-- OPTIONAL
        console.log('Maximum absolute deviation: ' 
                + module._error(bufForData, bufForApprox, sizeOfData)); 
      }
    }
    else
      console.log('FAIL: no columns with int or double data are passed!');
  } 
  else { // i.e. memory allocation fail
      console.log('FAIL: the required memory is not allocated!');
  }  
  // clear allocted memory
  for(let i = 0; i < buffers.length; i++)
      module._free(buffers[i]);
  
  return resultCode;
}


//name: verifyPCA
//input: dataframe df
//input: int numOfComponents = 3
export function verifyPCA(df, numOfComponents) {
  // check input
  if(numOfComponents < 1) {
    console.log('Uncorrect number of principal components: ' + numOfComponents);
    return;
  }

  // structures for column data
  let approx = [];
  let princComp = [];
 
  // PRINCIPAL COMPONENT ANALYSIS
  console.log('PCA ...');

  let start = new Date().getTime(); 

  let resultCode = pcaWithVerification(EigenPCA, df.columns.toList(), 
                                       numOfComponents, princComp, approx);

  let finish = new Date().getTime();

  // check result of PCA computation
  if(resultCode == 0) {
    
    console.log('SUCCESS: principal components and approximation of source data have been computed!');
    console.log('Time for PCA execution, including int-to-float conversion: '
      + (finish - start) + ' ms.');

    // Create dataframe with approximation
    console.log('Creating dataframes from the arrays computed ...');

    start = new Date().getTime();
  
    let tableWithApprox = DG.DataFrame.fromColumns(approx);
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

    console.log('Dataframes have been created!');
    console.log('Time for creating dataframes: ' + (finish - start) + ' ms.');
  }
  else
    console.log('PCA computation fail: result code is ' + resultCode);
}

// PCA by Eigen
function pca(module, columns, numOfPrincipalComponents, principalComponents){
  
  let resultCode;
    
  // sizes
  let height = columns.length;
  let width = columns[0].length;
  let sizeOfData = height * width;
  let sizeOfPrincComp = numOfPrincipalComponents * width; 
  let numOfBytes = Float32Array.BYTES_PER_ELEMENT; // since we operate Float32

  // heap: here, we take into account an operating Float32Arrays
  let shift = 2;
  let heap = module.HEAPF32;

  // BUFFERS ALLOCATION

  let buffers = [];
  let excpectedNumOfBuffers = 3; // 1st - for source data, 
                                 // 2nd - for principal components
                                 // 3rd - for approximated data: OPTIONAL

  let bufForData = module._malloc(sizeOfData * numOfBytes);
  if(bufForData != 0)
      buffers.push(bufForData);

  let bufForPrincComp = module._malloc(sizeOfPrincComp * numOfBytes);
  if(bufForPrincComp != 0)
      buffers.push(bufForPrincComp);

  let bufForApprox = module._malloc(sizeOfData * numOfBytes);
  if(bufForApprox != 0)
      buffers.push(bufForApprox);      
   
  // check memory allocation
  if(buffers.length == excpectedNumOfBuffers) { 

    // assign allocated buffer with source data
    for(let i = 0; i < columns.length; i++){
      let array = null;
      let col = columns[i];

      // get raw data from column
      switch(col.type){
        case 'double':
          array = col.getRawData();
          break;
        case 'int':
          array = new Float32Array(col.getRawData()); // conversion is applied, 
                                                      // since PCA requires floats
          break;
        default:
          console.log('Column ' + col.name + ' has uncorrect type: ' + col.type);
          height--; // i.e. actual number of columns, which are passed to PCA, 
                    // should be reduced
          break;
      }

      // check data array
      if(array != null) {
        heap.set(array, (bufForData + i * width * numOfBytes) >> shift);
      }
    }

    // check size of data to be processed:
    // height == 0 means that there is no column data in bufForData
    if(height > 0){ 
      
      // EXPORTED С++-FUNCTION CALL
      resultCode = module.ccall('pcaWithApproximation', 'number', 
                     ['number', 'number', 'number', 'number', 'number', 'number'],
                     [bufForData, height, width, numOfPrincipalComponents, 
                      bufForPrincComp, bufForApprox]);
     
      // check result of the exported C++-function call:
      // resultCode != 0 means some error (see PCA/PCA.h)
      if(resultCode == 0) {
        // get results: principal components
        for(let i = 0; i < numOfPrincipalComponents; i++)
          principalComponents.push(new Float32Array(heap.buffer, 
             bufForPrincComp + i * width * numOfBytes, width));        
      }
    }
    else
      console.log('FAIL: no columns with int or double data are passed!');
  } 
  else { // i.e. memory allocation fail
      console.log('FAIL: the required memory is not allocated!');
  }  
  // clear allocted memory
  for(let i = 0; i < buffers.length; i++)
      module._free(buffers[i]);
  
  return resultCode;
}

//name: testPCA
//input: dataframe df
//input: int numOfComponents = 3
export function testPCA(df, numOfComponents) {
  // check input
  if(numOfComponents < 1) {
    console.log('Uncorrect number of principal components: ' + numOfComponents);
    return;
  }

  // structure for principal components  
  let princComp = [];
 
  // PRINCIPAL COMPONENT ANALYSIS
  let resultCode = pca(EigenPCA, df.columns.toList(), numOfComponents, princComp);

  // check result of PCA computation
  if(resultCode == 0) {     
  
    // Create dataframe with principal components
    let princCompCols = [];
  
    for(let i = 1; i <= numOfComponents; i++)
      princCompCols.push(DG.Column.fromFloat32Array('comp. ' + i.toString(), princComp[i - 1]));
   
    let tableWithPrincComp = DG.DataFrame.fromColumns(princCompCols);
  
    tableWithPrincComp.name = 'PrincipalComponents';
    grok.shell.addTableView(tableWithPrincComp);

  }
  else
    console.log('PCA computation fail: result code is ' + resultCode);
}
