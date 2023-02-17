/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: fib
//input: int num = 10
export async function fib(num) {
  var worker = new Worker(new URL('../wasm/workerGal.js', import.meta.url));
  
  worker.postMessage(num);

  worker.onmessage = function(e) {  
    alert(e.data);
  } 
}

//name: fillArray
export async function fillArray() {
  let len = 10;
  let array = new Float32Array(10);

  console.log('Array before:');
  console.log(array);

  var worker = new Worker(new URL('../wasm/workerFillArray.js', import.meta.url));
  
  worker.postMessage({'array': array});

  worker.onmessage = function(e) {  
    alert(e.data);
  } 
}

//name: testODE
export async function testODE() {
var worker = new Worker(new URL('../wasm/workerODE.js', import.meta.url));
  
  worker.postMessage(_package.webRoot);

  worker.onmessage = function(e) {  
    alert(e.data);
  } 
}

//name: testFoo
//input: int num = 10
export async function testFoo(num) {
  var worker = new Worker(new URL('../wasm/workerFoo.js', import.meta.url));
    
    worker.postMessage(num);
  
    worker.onmessage = function(e) {  
      alert(e.data);
    } 
}

//name: testFooO3
//input: int num = 10
export async function testFoo3(num) {
  var worker = new Worker(new URL('../wasm/workerFooO3.js', import.meta.url));
    
    worker.postMessage(num);
  
    worker.onmessage = function(e) {  
      alert(e.data);
    } 
}
