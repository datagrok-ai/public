/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import createGrokCache from './grokcache';

export let _package = new DG.Package();
export let grokCache = null; // WASM module to be initialized later


// WASM Module init: obtain the entry point
async function grokCacheModuleInit() {
  grokCache = await createGrokCache();
}


//name: test
export function test() {
  grok.shell.info(_package.webRoot);
}

//name: init
//tags: autostart
export async function init() {
  if (grokCache == null) {
    await grokCacheModuleInit();
    let greeting = grokCache.wasm_version();
    console.log("grokCache:" + greeting);
  }
}


//name: diagnostics
export async function grokCacheDiagnostics() {
  // draw a test view here
  // console is ok
  if (grokCache != null) {
    let greeting = grokCache.wasm_version();
    console.log(greeting);
  }
}
