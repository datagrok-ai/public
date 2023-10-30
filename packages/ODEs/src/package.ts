/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {initMatrOperApi, getInverseMatrix} from '../wasm/matrix-operations-api';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//tags: init
export async function init() {
  await initMatrOperApi();
}

//name: check
//input: int n
export function check(n: number) {
  const A = (n === 1) ? new Float64Array([10]) : 
  ((n === 2) ? new Float64Array([1,2,3,4]) : new Float64Array([1,1,1,0,1,1,0,0,1]));
  console.log('A:');
  console.log(A);
  console.log('inv A:');
  console.log(getInverseMatrix(A, n));
}

//name: performance
//input: int size = 15
//input: int times = 10000
export function performance(size: number, times: number) {
  const A = new Float64Array(size * size);

  for (let i = 0; i < size; ++i)
    for (let j = i; j < size; ++j)
      A[j + i * size] = Math.random();

  let sum = 0;

  for (let k = 0; k < times; ++k) {
    const start = new Date().getTime();
    getInverseMatrix(A, size);
    const finish = new Date().getTime();
    sum += finish - start;
  }

  console.log(`Total time: ${sum} ms.`);
  console.log(`Average time: ${sum/times} ms.`);
}
