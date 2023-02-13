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
 // await initFib();
}

//name: fib1
//input: int num = 44
export async function fib1(num) {
  (async () => {
    const fetchPromise = fetch(new URL('../fibs/fib.wasm', import.meta.url));
    const module = await WebAssembly.compileStreaming(fetchPromise);
    const instance = await WebAssembly.instantiate(module);
    
    const func = instance.exports.fib; // this is the function required
    
    alert(func(num));
  })();  
}

//name: fib2
//input: int num = 44
export async function fib2(num) {
 // alert(Fib._fib(num));
}

//name: fibWW
//input: int num = 44
export async function fibWW(num) {  

  const worker = new Worker(new URL('../fibs/worker.js', import.meta.url));  

  worker.postMessage(num);

  worker.onmessage = function(e) {  
    alert(e.data);
  } 
}

//name: sum
export async function sum() {  

  const worker = new Worker(new URL('../sum/worker.js', import.meta.url));
  
  let arr = new Int32Array([1,2,3,4,5,6]);

  worker.postMessage(arr);

  worker.onmessage = function(e) {  
    alert(e.data);
  } 
}

//name: doubleArr
export async function doubleArr() {  

  const worker = new Worker(new URL('../double/worker.js', import.meta.url));
  
  let arr = new Int32Array([1,2,3,4,5,6]);

  console.log(arr);

  worker.postMessage(arr);

  worker.onmessage = function(e) {  
    console.log(e.data);
  } 
}

