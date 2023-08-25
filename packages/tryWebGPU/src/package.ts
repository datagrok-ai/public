/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {checkNavigator, accessGPU, buffers, 
  matrixMultiplicationExample, gpuMatrixProduct, cpuMatrixProduct} from './web-gpu-sandbox';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: First steps
export async function firstSteps() {
  //checkNavigator();
  //await accessGPU();
  //await buffers();

  await matrixMultiplicationExample();
}

//name: Test matrix product: GPU & CPU
export async function testMatrixProducts() {
  const firstMatrix = new Float32Array([
    3 /* rows */,
    4 /* columns */,
    1, 2, 3, 4,
    5, 6, 7, 8,
    1, 2, 1, 0
  ]);

  const secondMatrix = new Float32Array([
    4 /* rows */,
    3 /* columns */,
    1, 2, 3,
    3, 4, 5,
    5, 6, 7,
    7, 8, 8
  ]);  
  
  console.log('GPU product:');
  console.log(await gpuMatrixProduct(firstMatrix, secondMatrix));

  console.log('CPU product:');
  console.log(cpuMatrixProduct(firstMatrix, secondMatrix));
}

//name: GPU vs CPU
//description: Compare matrix product provided by GPU & CPU: two n x n matrices are multipled.
//input: int n = 512 {caption: rows & columns count}
export async function matrixProductPerformance(n: number) {
  const A = grok.data.demo.randomWalk(n, n);
  A.name = 'A';
  grok.shell.addTable(A);

  const B = grok.data.demo.randomWalk(n, n);
  B.name = 'B';
  grok.shell.addTable(B);

  /* Product is applied to matrices given in the form of Float32Array:
      1-st element: rows count,
      2-nd element: columns count,
      other elements: matrix data. 
  */

  // First multiplier
  const firstMatrix = new Float32Array(n * n + 2);
  firstMatrix[0] = n;
  firstMatrix[1] = n;  

  let idx = 2;
  
  for (const col of A.columns) {
    const arr = col.getRawData();

    for (let j = 0; j < n; ++j) {
      firstMatrix[idx] = arr[j];
      ++idx;
    }
  }

  idx = 2;

  // Second multiplier
  const secondMatrix = new Float32Array(n * n + 2);
  secondMatrix[0] = n;
  secondMatrix[1] = n;

  for (const col of B.columns) {
    const arr = col.getRawData();

    for (let j = 0; j < n; ++j) {
      secondMatrix[idx] = arr[j];
      ++idx;
    }
  }

  // GPU matrix product
  let start = new Date().getTime();
  await gpuMatrixProduct(firstMatrix, secondMatrix);
  let finish = new Date().getTime();
  console.log(`GPU: ${finish - start} ms.`);  

  // CPU matrix product
  start = new Date().getTime();
  cpuMatrixProduct(firstMatrix, secondMatrix);
  finish = new Date().getTime();
  console.log(`CPU: ${finish - start} ms.`);
} // productPerfromanceGPU
