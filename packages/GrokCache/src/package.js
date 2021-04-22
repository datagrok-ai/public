/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import createGrokCache from './grokcache';
import Worker from './grokcache.worker.js';
//import Worker from "worker-loader!./grokcache.worker.js";

export let _package = new DG.Package();
export let grokCache = null; // WASM module to be initialized later

// initialization promise to protect package from the init races
let initCachePromise = null;


// WASM Module init: obtain the entry point
async function grokCacheModuleInit() {
  grokCache = await createGrokCache();
  return grokCache;
}


//name: test
export function test() {
  grok.shell.info(_package.webRoot);
}

//name: init
//tags: autostart
export async function init() {
  if (!initCachePromise) {
    initCachePromise = grokCacheModuleInit(); // returns Promise
    initCachePromise.then((gC) => { // define a promise resolution functor
      console.log(gC.wasm_version());
    });
  } else {
    console.log('grokCache: init promise already created. waiting...');
  }
  await initCachePromise;
  if (grokCache == null) {
    console.error("grokCache initialization not fullfilled!");
  }
  return grokCache;
}


//name: diagnostics
export async function Diagnostics() {
  await initCachePromise;
  // draw a test view here
  // console is ok

  let greeting = grokCache.wasm_version();
  console.log(greeting);
}
